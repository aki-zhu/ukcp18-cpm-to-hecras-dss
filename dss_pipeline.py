"""Metadata-driven rainfall preparation. This module is shared by both repositories.

Python 3 prepares portable payloads; HEC-DSSVue's Jython writes/reads DSS.
Source files are always opened read-only. All spatial arrays use increasing x/y.
"""
from __future__ import annotations

import argparse
import calendar
from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS
from rasterio import Affine, band
from rasterio.features import geometry_mask
from rasterio.io import MemoryFile
from rasterio.warp import calculate_default_transform, reproject, Resampling, transform_bounds


REGIONS = [
    ("anglian", "Region_Anglian"), ("clyde", "Region_Clyde"),
    ("dee", "Region_Dee"), ("forth", "Region_Forth"),
    ("humber", "Region_humber"), ("northeastscotland", "Region_NEScotland"),
    ("northhighland", "Region_NHighland"), ("northwestengland", "Region_NWEngland"),
    ("northumbria", "Region_Northumbria"), ("severn", "Region_Severn"),
    ("solway", "Region_Solway"), ("southeastengland", "Region_SEEngland"),
    ("southwestengland", "Region_SWEngland"), ("tay", "Region_Tay"),
    ("thames", "Region_Thames"), ("tweed", "Region_Tweed"),
    ("westernwales", "Region_WestWales"),
]
DEFAULT_WINDOWS = {
    "cpm": [("01", "2019-11-01", "2039-12-31"),
            ("02", "2039-11-01", "2059-12-31"),
            ("03", "2059-11-01", "2080-11-30")],
    "haduk": [("01", "1979-11-01", "1999-12-31"),
              ("02", "1999-11-01", "2019-12-31")],
}


def safe_name(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", str(value)):
        raise ValueError("Names/window IDs must contain letters, numbers, _, . or -: " + str(value))
    return str(value)


@dataclass(frozen=True)
class Window:
    name: str
    start: pd.Timestamp
    stop: pd.Timestamp  # exclusive

    @classmethod
    def parse(cls, name, start, end):
        a, b = pd.Timestamp(start), pd.Timestamp(end)
        if a.tzinfo is not None or b.tzinfo is not None:
            raise ValueError("Use timezone-free ISO dates/times (source model time).")
        # Date-only end includes that entire day; datetime end is exclusive.
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(end)):
            b += pd.Timedelta(days=1)
        if b <= a or a.second or a.microsecond or b.second or b.microsecond:
            raise ValueError("Window must have increasing, whole-minute endpoints.")
        return cls(safe_name(name), a, b)

    @property
    def tag(self):
        return f"{self.name}_{self.start:%Y%m%dT%H%M}_{self.stop:%Y%m%dT%H%M}"


def select_windows(args, kind):
    if args.window and (args.start or args.end):
        raise ValueError("Use either --start/--end or repeated --window ID START END.")
    if bool(args.start) != bool(args.end):
        raise ValueError("Supply both --start and --end.")
    entries = args.window or ([("custom", args.start, args.end)] if args.start else DEFAULT_WINDOWS[kind])
    windows = [Window.parse(*entry) for entry in entries]
    if len({w.name for w in windows}) != len(windows):
        raise ValueError("Window IDs must be unique.")
    if args.windows:
        unknown = set(args.windows) - {w.name for w in windows}
        if unknown:
            raise ValueError(f"Unknown window IDs: {sorted(unknown)}")
        windows = [w for w in windows if w.name in args.windows]
    return windows


def date_key(value):
    if isinstance(value, np.datetime64):
        value = pd.Timestamp(value)
    if getattr(value, "second", 0) or getattr(value, "microsecond", 0):
        raise ValueError("Sub-minute source timestamps are not supported by DSS grid times.")
    return tuple(int(getattr(value, k)) for k in ("year", "month", "day", "hour", "minute"))


def seconds_between(a, b):
    if isinstance(a, np.datetime64):
        return float((b - a) / np.timedelta64(1, "s"))
    return (b - a).total_seconds()


def shifted(value, seconds):
    if isinstance(value, np.datetime64):
        return value + np.timedelta64(int(seconds), "s")
    return value + timedelta(seconds=float(seconds))


def unit_scale(units, interval_seconds):
    """Convert water-equivalent amounts or rates to mm per source interval."""
    u = str(units).lower().strip().replace("**", "^").replace(" ", "").replace("^", "")
    amounts = {"mm": 1., "kgm-2": 1., "kg/m2": 1., "m": 1000.}
    if u in amounts:
        return amounts[u]
    for suffixes, duration in [(("/day", "d-1", "day-1"), 86400.),
                               (("/hour", "/hr", "/h", "h-1", "hour-1"), 3600.),
                               (("/s", "s-1"), 1.)]:
        for suffix in suffixes:
            if u.endswith(suffix) and u[:-len(suffix)] in amounts:
                return amounts[u[:-len(suffix)]] * interval_seconds / duration
    raise ValueError(f"Unsupported precipitation units {units!r}; supply a verified --units override.")


def axis_name(ds, da, axis):
    standards = {"x": {"projection_x_coordinate", "grid_longitude", "longitude"},
                 "y": {"projection_y_coordinate", "grid_latitude", "latitude"}}[axis]
    names = {"x": ["projection_x_coordinate", "grid_longitude", "x", "longitude", "lon"],
             "y": ["projection_y_coordinate", "grid_latitude", "y", "latitude", "lat"]}[axis]
    candidates = [name for name, coord in ds.coords.items() if coord.ndim == 1
                  and coord.dims[0] in da.dims and
                  (coord.attrs.get("standard_name") in standards or coord.attrs.get("axis", "").lower() == axis)]
    for name in names:
        if name in ds and ds[name].ndim == 1 and ds[name].dims[0] in da.dims:
            candidates.append(name)
    candidates = list(dict.fromkeys(candidates))
    if not candidates:
        raise ValueError("Expected two 1-D spatial axes. Regional averages/curvilinear meshes "
                         "cannot be interpreted as rainfall grids.")
    return candidates[0]


def spacing(values):
    delta = np.diff(values)
    if len(values) < 2 or not np.isfinite(values).all() or not (np.all(delta > 0) or np.all(delta < 0)):
        raise ValueError("Source axes must have at least two finite, strictly monotonic centres.")
    step = float(abs(values[-1]-values[0])/(len(values)-1))
    return step


def regular_axis(values):
    if len(values) < 2:
        return True
    step = abs(values[-1]-values[0])/(len(values)-1)
    # Native CPM has both float32 quantisation AND a genuinely stretched halo.
    return (np.allclose(np.abs(np.diff(values)), step, rtol=0, atol=step*1e-3)
            and np.allclose(values, np.linspace(values[0], values[-1], len(values)), rtol=0, atol=step*1e-3))


@dataclass
class Grid:
    x: np.ndarray
    y: np.ndarray
    dx: float
    dy: float
    crs: CRS

    @property
    def transform(self):
        # South-up, matching increasing y and DSS's lower-left-first array order.
        if not self.regular:
            raise ValueError("A stretched grid has no single affine transform; reproject using geolocation arrays.")
        return Affine(self.dx, 0, self.x[0] - self.dx / 2,
                      0, self.dy, self.y[0] - self.dy / 2)

    @property
    def bounds(self):
        left = (self.x[1]-self.x[0]) if len(self.x)>1 else self.dx
        right = (self.x[-1]-self.x[-2]) if len(self.x)>1 else self.dx
        bottom = (self.y[1]-self.y[0]) if len(self.y)>1 else self.dy
        top = (self.y[-1]-self.y[-2]) if len(self.y)>1 else self.dy
        return (self.x[0]-left/2, self.y[0]-bottom/2, self.x[-1]+right/2, self.y[-1]+top/2)

    @property
    def regular(self):
        return regular_axis(self.x) and regular_axis(self.y)

    def crop(self, rows, cols):
        x, y = self.x[cols], self.y[rows]
        return Grid(x, y, spacing(x) if len(x)>1 else self.dx,
                    spacing(y) if len(y)>1 else self.dy, self.crs)

    def same_as(self, other):
        return (self.crs.equals(other.crs) and self.x.shape == other.x.shape
                and self.y.shape == other.y.shape and np.allclose(self.x, other.x, rtol=0, atol=self.dx * 1e-5)
                and np.allclose(self.y, other.y, rtol=0, atol=self.dy * 1e-5))

    def dss_metadata(self):
        if not self.crs.is_projected or any(abs(a.unit_conversion_factor - 1) > 1e-8 for a in self.crs.axis_info[:2]):
            raise ValueError("DSS output requires a projected CRS in metres; set --target-crs and --target-cell-size.")
        if not self.regular:
            raise ValueError("DSS needs regular square cells; set --target-cell-size to resample the stretched grid.")
        if not np.isclose(self.dx, self.dy, rtol=1e-5):
            raise ValueError("DSS needs square cells; set --target-cell-size to resample rectangular cells.")
        # Zero offsets with the true southwest corner also support fractional origins.
        return dict(nx=len(self.x), ny=len(self.y), cell_size=self.dx,
                    lower_left_cell_x=0, lower_left_cell_y=0,
                    cell_zero_x=float(self.bounds[0]), cell_zero_y=float(self.bounds[1]),
                    srs_name=self.crs.name, srs_wkt=self.crs.to_wkt(version="WKT1_ESRI"))


def inspect_variable(ds, args):
    if args.variable not in ds:
        raise ValueError(f"Variable {args.variable!r} is missing.")
    da = ds[args.variable]
    selections = {}
    for item in args.select:
        dim, value = item.split("=", 1)
        if dim not in da.dims:
            raise ValueError(f"Selection dimension {dim} is absent.")
        coords = ds[dim].values
        matches = np.flatnonzero(np.asarray([str(v) == value for v in coords]))
        if not len(matches):
            try:
                matches = np.flatnonzero(np.asarray(coords, dtype=float) == float(value))
            except (ValueError, TypeError):
                pass
        if len(matches) != 1:
            raise ValueError(f"Selection {item} does not identify exactly one coordinate value.")
        da = da.isel({dim: int(matches[0])}, drop=True)
        selections[dim] = str(coords[matches[0]])
    xn, yn = axis_name(ds, da, "x"), axis_name(ds, da, "y")
    for dim in list(da.dims):
        if dim not in ("time", ds[xn].dims[0], ds[yn].dims[0]):
            if da.sizes[dim] != 1:
                raise ValueError(f"Dimension {dim} has {da.sizes[dim]} values; choose one with --select {dim}=VALUE.")
            selections[dim] = str(ds[dim].values[0]) if dim in ds else "singleton"
            da = da.isel({dim: 0}, drop=True)
    if "time" not in da.dims:
        raise ValueError("Rainfall must have a time dimension.")
    x, y = np.asarray(ds[xn].values, dtype=float), np.asarray(ds[yn].values, dtype=float)
    dx, dy = spacing(x), spacing(y)
    mapping_name = ds[args.variable].attrs.get("grid_mapping", "").split(":")[0]
    if args.source_crs:
        crs = CRS.from_user_input(args.source_crs)
    elif mapping_name in ds:
        attrs = ds[mapping_name].attrs
        wkt = attrs.get("spatial_ref", attrs.get("crs_wkt"))
        crs = CRS.from_wkt(wkt) if wkt else CRS.from_cf(attrs)
    else:
        raise ValueError("CF grid_mapping/CRS is missing. Supply the verified --source-crs (e.g. EPSG:27700).")
    if crs.is_geographic:
        # UKCP rotated longitudes use 353..370 degrees. Use the equivalent
        # -7..10 range so inverse transforms and source cropping agree.
        x = x - 360*round(float(np.median(x))/360)
    for name, values in ((xn, x), (yn, y)):
        units = ds[name].attrs.get("units", "").lower()
        if crs.is_projected and units in {"km", "kilometre", "kilometer", "kilometres", "kilometers"}:
            values *= 1000
            if name == xn:
                dx *= 1000
            else:
                dy *= 1000
        elif crs.is_projected and units and units not in {"m", "metre", "meter", "metres", "meters"}:
            raise ValueError(f"Unsupported projected axis units: {units}")
    if x[0] > x[-1]:
        x = x[::-1]
        da = da.isel({ds[xn].dims[0]: slice(None, None, -1)})
    if y[0] > y[-1]:
        y = y[::-1]
        da = da.isel({ds[yn].dims[0]: slice(None, None, -1)})
    return da.transpose("time", ds[yn].dims[0], ds[xn].dims[0]), Grid(x, y, dx, dy, crs), selections


def source_times(ds, args, kind):
    times = ds.time.values
    if len(times) == 0:
        raise ValueError("Empty source time coordinate.")
    cal = str(ds.time.encoding.get("calendar", "standard")).lower()
    if cal not in {"standard", "gregorian", "proleptic_gregorian", "360_day"}:
        raise ValueError(f"Unsupported source calendar {cal}; no silent calendar conversion is performed.")
    bounds_name = ds.time.attrs.get("bounds")
    bounds = ds[bounds_name].values if bounds_name in ds else None
    if bounds is not None:
        if bounds.shape != (len(times), 2):
            raise ValueError("Time bounds must have dimensions (time, 2).")
        durations = np.asarray([seconds_between(a, b) for a, b in bounds])
        step = float(durations[0])
        starts = list(bounds[:, 0])
        if not np.allclose(durations, step, rtol=0, atol=.001):
            raise ValueError("Variable length intervals (e.g. monthly/seasonal means) cannot be used as hourly/daily rainfall.")
        if any(not 0 <= seconds_between(a, t) <= step for t, (a, b) in zip(times, bounds)):
            raise ValueError("Time coordinates fall outside their CF bounds.")
    else:
        if args.time_label == "auto":
            raise ValueError("No CF time bounds: specify --time-label start|centre|end after checking the source.")
        step = float(args.interval_hours * 3600) if args.interval_hours else (
            seconds_between(times[0], times[1]) if len(times) > 1 else 0)
        offsets = {"start": 0, "centre": -.5, "end": -1}
        starts = [shifted(t, step * offsets[args.time_label]) for t in times]
    if step <= 0 or step > 86400 or step % 60 or 86400 % step:
        raise ValueError("Intervals must be whole minutes, at most one day, and divide 24 hours exactly.")
    if args.interval_hours and not np.isclose(step, args.interval_hours * 3600):
        raise ValueError("--interval-hours disagrees with CF time bounds.")
    policy = "cf_interval_bounds" if bounds is not None else "explicit_" + args.time_label
    if kind == "haduk" and step == 86400:
        # Preserve the original workflow's nominal daily labels from the downloaded files.
        starts = [pd.Timestamp(t).normalize().to_datetime64() for t in times]
        policy = "haduk_original_nominal_day_00_to_24"
        if cal == "360_day":
            raise ValueError("HadUK daily observations require a Gregorian calendar.")
    return [date_key(t) for t in starts], int(step), cal, policy


class Catalog:
    """Index only metadata; load rainfall in small time/spatial blocks on demand."""
    def __init__(self, paths, args, kind):
        self.args, self.kind, self.paths = args, kind, list(paths)
        self.records, self.file_info = {}, []
        self.grid = None
        self.invalid_cells = 0
        self.cache = None
        for path in paths:
            with xr.open_dataset(path, engine=args.engine, decode_times=True) as ds:
                da, grid, selections = inspect_variable(ds, args)
                keys, step, cal, policy = source_times(ds, args, kind)
                units = args.units or da.attrs.get("units", "")
                scale = unit_scale(units, step)
                signature = (step, cal, selections, str(units), policy)
                if self.grid is None:
                    self.grid, self.signature = grid, signature
                    self.step, self.calendar, self.units, self.policy = step, cal, str(units), policy
                elif not self.grid.same_as(grid) or self.signature != signature:
                    raise ValueError(f"Grid, units, calendar, interval or ensemble changed in {path.name}. "
                                     "Process different products/ensembles in separate runs.")
                if keys != sorted(keys) or len(keys) != len(set(keys)):
                    raise ValueError(f"Non-monotonic/duplicate time coordinates in {path.name}")
                if kind == "haduk" and step == 86400:
                    match = re.search(r"_(\d{8})-(\d{8})\.nc$", path.name)
                    if match:
                        expected = [date_key(t) for t in pd.date_range(*match.groups(), freq="D")]
                        if keys != expected:
                            raise ValueError(f"HadUK filename dates disagree with time coordinates: {path.name}")
                for index, key in enumerate(keys):
                    if key in self.records:
                        raise ValueError(f"Overlapping source records at {key}; narrow --input-pattern or --select.")
                    self.records[key] = (path, index, scale)
                self.file_info.append(dict(path=str(path.resolve()), size=path.stat().st_size,
                                           mtime_ns=path.stat().st_mtime_ns, records=len(keys)))
        if not self.records:
            raise ValueError("No source records.")
        phases = {(key[3] * 3600 + key[4] * 60) % self.step for key in self.records}
        if len(phases) != 1:
            raise ValueError("Source intervals are not aligned to a consistent time step.")
        self.phase = phases.pop()

    def close(self):
        if self.cache:
            self.cache[1].close()
            self.cache = None

    def read(self, keys, plan):
        missing = [key for key in keys if key not in self.records]
        if missing:
            raise ValueError(f"Missing source intervals, first: {missing[0]}. "
                             "360_day conversion requires complete source months, including edge months.")
        output = np.empty((len(keys), len(plan.source.y), len(plan.source.x)), dtype=np.float32)
        groups = defaultdict(list)
        for pos, key in enumerate(keys):
            path, index, scale = self.records[key]
            groups[(path, scale)].append((pos, index))
        for (path, scale), items in groups.items():
            if not self.cache or self.cache[0] != path:
                self.close()
                ds = xr.open_dataset(path, engine=self.args.engine, decode_times=True)
                da, _, _ = inspect_variable(ds, self.args)
                self.cache = path, ds, da
            da = self.cache[2]
            values = np.asarray(da.isel({"time": [i for _, i in items],
                                        da.dims[1]: plan.source_rows, da.dims[2]: plan.source_cols}).values,
                                dtype=np.float32) * scale
            invalid = ~np.isfinite(values)
            self.invalid_cells += int(invalid.sum())
            if self.args.missing_data == "error" and invalid.any():
                raise ValueError(f"Missing/nonfinite rainfall in {path.name}")
            values[invalid] = 0
            if (values < -1e-6).any():
                raise ValueError(f"Negative precipitation in {path.name}")
            values[values < 0] = 0
            output[[pos for pos, _ in items]] = values
        return plan.apply(output)


def dilate_mask(mask, cells):
    if cells < 0:
        raise ValueError("--buffer-cells cannot be negative.")
    result = mask.copy()
    for _ in range(cells):
        padded = np.pad(result, 1)
        result = np.logical_or.reduce([padded[r:r+mask.shape[0], c:c+mask.shape[1]]
                                       for r in range(3) for c in range(3)])
    return result


def load_polygons(path, crs, boundary_crs=None):
    import geopandas as gpd
    gdf = gpd.read_file(path)
    gdf = gdf[gdf.geometry.notna() & gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    if gdf.empty or not gdf.geometry.is_valid.all():
        raise ValueError(f"Boundary must contain valid polygons: {path}")
    if gdf.crs is None:
        if not boundary_crs:
            raise ValueError("Boundary CRS missing; supply --boundary-crs instead of guessing.")
        gdf = gdf.set_crs(boundary_crs)
    return gdf.to_crs(crs)


def find_boundary(base, folder):
    directory = Path(base) / folder
    overrides = {"Region_humber": "UKCP18_uk_river_hires_humber.shp",
                 "Region_Northumbria": "Northumbria/Northumbria_boundary.shp"}
    if folder in overrides and (directory / overrides[folder]).exists():
        return directory / overrides[folder]
    import geopandas as gpd
    candidates = []
    for path in directory.rglob("*.shp"):
        frame = gpd.read_file(path)
        if frame.geometry.geom_type.isin(["Polygon", "MultiPolygon"]).any():
            candidates.append(path)
    preferred = [p for p in candidates if "boundary" in str(p).lower() or "ukcp18_uk_river_hires" in p.name.lower()]
    candidates = preferred or candidates
    canonical = [p for p in candidates if not p.stem.lower().endswith(("_1", "_copy", "_backup"))]
    candidates = canonical or candidates
    if len(candidates) != 1:
        raise ValueError(f"Expected one polygon boundary under {directory}, found {candidates}; use --boundary explicitly.")
    return candidates[0]


class SpatialPlan:
    def __init__(self, source, args, boundary=None):
        target_crs = CRS.from_user_input(args.target_crs) if args.target_crs else source.crs
        self.resampling = Resampling[args.resampling]
        if args.target_cell_size is not None or not source.crs.equals(target_crs):
            if not args.target_cell_size or args.target_cell_size <= 0:
                raise ValueError("Set a positive --target-cell-size in metres when reprojecting (e.g. 2200).")
            if not target_crs.is_projected:
                raise ValueError("--target-crs must be a projected CRS in metres.")
            affine, nx, ny = calculate_default_transform(
                source.crs, target_crs, len(source.x), len(source.y),
                *source.bounds, resolution=args.target_cell_size)
            dx, dy = affine.a, abs(affine.e)
            target = Grid(affine.c + (np.arange(nx) + .5) * dx,
                          affine.f - (np.arange(ny)[::-1] + .5) * dy, dx, dy, target_crs)
        else:
            target = source
        target.dss_metadata()
        self.mask = None
        if boundary:
            frame = load_polygons(boundary, target.crs, args.boundary_crs)
            mask = geometry_mask(frame.geometry, (len(target.y), len(target.x)),
                                 target.transform, invert=True, all_touched=True)
            mask = dilate_mask(mask, args.buffer_cells)
            rows, cols = np.where(mask)
            if not len(rows):
                raise ValueError(f"Boundary {boundary} does not overlap the rainfall grid.")
            rs, cs = slice(rows.min(), rows.max()+1), slice(cols.min(), cols.max()+1)
            self.mask = mask[rs, cs]
            target = target.crop(rs, cs)
        self.target = target
        # Only read source cells near the output region, including interpolation neighbours.
        bounds = transform_bounds(target.crs, source.crs, *target.bounds, densify_pts=41)
        cols = np.flatnonzero((source.x >= bounds[0]-2*source.dx) & (source.x <= bounds[2]+2*source.dx))
        rows = np.flatnonzero((source.y >= bounds[1]-2*source.dy) & (source.y <= bounds[3]+2*source.dy))
        if not len(cols) or not len(rows):
            raise ValueError("No source cells overlap the requested target grid.")
        self.source_rows, self.source_cols = slice(rows[0], rows[-1]+1), slice(cols[0], cols[-1]+1)
        self.source = source.crop(self.source_rows, self.source_cols)
        self.direct = self.source.same_as(target)
        self.geoloc_file = None

    def close(self):
        if self.geoloc_file is not None:
            self.geoloc_file.close()
            self.geoloc_file = None

    def geolocation(self):
        if self.geoloc_file is None:
            self.geoloc_file = MemoryFile()
            with self.geoloc_file.open(driver="GTiff", width=len(self.source.x), height=len(self.source.y),
                                       count=2, dtype="float64", crs=self.source.crs) as dataset:
                xx, yy = np.meshgrid(self.source.x, self.source.y[::-1])
                dataset.write(xx, 1)
                dataset.write(yy, 2)
                # GDAL otherwise treats geolocation coordinates as pixel corners,
                # which would introduce a half-cell shift in the stretched halo.
                dataset.update_tags(GEOREFERENCING_CONVENTION="PIXEL_CENTER")
        return self.geoloc_file.name

    def apply(self, values):
        if self.direct:
            output = values
        else:
            output = np.zeros((len(values), len(self.target.y), len(self.target.x)), dtype=np.float32)
            for i, array in enumerate(values):
                # GDAL works north-up; reverse rows on both sides of this boundary.
                src_transform = Affine(self.source.dx, 0, self.source.bounds[0],
                                       0, -self.source.dy, self.source.bounds[3])
                dst_transform = Affine(self.target.dx, 0, self.target.bounds[0],
                                       0, -self.target.dy, self.target.bounds[3])
                north_up = np.zeros(output.shape[1:], dtype=np.float32)
                destination = dict(dst_transform=dst_transform, dst_crs=self.target.crs,
                                   resampling=self.resampling, dst_nodata=0, init_dest_nodata=True)
                if self.source.regular:
                    reproject(array[::-1].copy(), north_up, src_transform=src_transform,
                              src_crs=self.source.crs, **destination)
                else:
                    # Describe geolocation in the source dataset itself. This
                    # also works with GDAL releases predating SRC_GEOLOC_ARRAY.
                    geoloc = self.geolocation()
                    with MemoryFile() as memory:
                        with memory.open(driver="GTiff", width=len(self.source.x), height=len(self.source.y),
                                         count=1, dtype="float32", crs=self.source.crs) as dataset:
                            dataset.write(array[::-1], 1)
                            dataset.update_tags(ns="GEOLOCATION", SRS=self.source.crs.to_wkt(),
                                                X_DATASET=geoloc, Y_DATASET=geoloc, X_BAND="1", Y_BAND="2",
                                                PIXEL_OFFSET="0", LINE_OFFSET="0", PIXEL_STEP="1", LINE_STEP="1",
                                                GEOREFERENCING_CONVENTION="PIXEL_CENTER")
                            reproject(band(dataset, 1), north_up, SRC_METHOD="GEOLOC_ARRAY", **destination)
                output[i] = north_up[::-1]
        if self.mask is not None:
            output[:, ~self.mask] = 0
        if not np.isfinite(output).all() or (output < 0).any():
            raise ValueError("Invalid values after reprojection.")
        return output


def conservative_month(values, days):
    """Interpolate 30 days, then conserve each cell's monthly amount (one subdaily slot)."""
    if len(values) != 30:
        raise ValueError("A complete 30-day source month is required.")
    target = np.linspace(0, 29, days)
    lo, hi = np.floor(target).astype(int), np.ceil(target).astype(int)
    weight = (target-lo)[:, None, None]
    out = values[lo] * (1-weight) + values[hi] * weight
    total, interpolated = values.sum(axis=0, dtype=np.float64), out.sum(axis=0, dtype=np.float64)
    out *= np.divide(total, interpolated, out=np.zeros_like(total), where=interpolated > 0)
    return out.astype(np.float32)


def expected_dates(catalog, window):
    anchor = window.start.normalize() + pd.Timedelta(seconds=catalog.phase)
    dates = pd.date_range(anchor, window.stop, freq=pd.Timedelta(seconds=catalog.step), inclusive="left")
    dates = dates[dates >= window.start]
    # Reject clipping through an accumulation interval, which would invent rainfall.
    if len(dates) == 0 or dates[-1] + pd.Timedelta(seconds=catalog.step) > window.stop + pd.Timedelta(seconds=catalog.phase):
        raise ValueError("Window does not contain complete source intervals.")
    if (window.start - anchor).total_seconds() % catalog.step or (window.stop - anchor).total_seconds() % catalog.step:
        # Date-only daily requests include records with the source's non-midnight phase.
        if not (window.start == window.start.normalize() and window.stop == window.stop.normalize()):
            raise ValueError("Datetime endpoints must align with the source accumulation intervals.")
    return dates


def iter_values(catalog, window, plan, temp_dir):
    dates = expected_dates(catalog, window)
    chunk = catalog.args.chunk_records
    if catalog.calendar != "360_day":
        keys = [date_key(t) for t in dates]
        missing = [k for k in keys if k not in catalog.records]
        if missing:
            raise ValueError(f"Requested time range is not covered; first missing interval {missing[0]}")
        for i in range(0, len(dates), chunk):
            yield dates[i:i+chunk], catalog.read(keys[i:i+chunk], plan)
        return
    # A disk-backed month keeps hourly national grids out of RAM. Each hour/slot
    # is remapped independently across days, preserving the diurnal phase.
    for period in pd.period_range(dates[0].to_period("M"), dates[-1].to_period("M"), freq="M"):
        days = calendar.monthrange(period.year, period.month)[1]
        slots = 86400 // catalog.step
        first = period.start_time + pd.Timedelta(seconds=catalog.phase)
        month_dates = pd.date_range(first, periods=days*slots, freq=pd.Timedelta(seconds=catalog.step))
        with tempfile.TemporaryDirectory(prefix="calendar_", dir=temp_dir) as working:
            mapped = np.memmap(Path(working) / "month.bin", mode="w+", dtype="float32",
                               shape=(days, slots, len(plan.target.y), len(plan.target.x)))
            flat = None
            try:
                for slot in range(slots):
                    seconds = catalog.phase + slot*catalog.step
                    keys = [(period.year, period.month, day, seconds//3600, (seconds%3600)//60)
                            for day in range(1, 31)]
                    mapped[:, slot] = conservative_month(catalog.read(keys, plan), days)
                flat = mapped.reshape((-1, len(plan.target.y), len(plan.target.x)))
                indices = np.flatnonzero((month_dates >= window.start) & (month_dates < window.stop))
                for i in range(0, len(indices), chunk):
                    take = indices[i:i+chunk]
                    yield month_dates[take], np.array(flat[take], copy=True)
            finally:
                del flat
                mapped._mmap.close()


def fingerprint(catalog, plan, window, args, boundary):
    payload = dict(sources=catalog.file_info, selection=catalog.signature,
                   grid=plan.target.dss_metadata(), window=[str(window.start), str(window.stop)],
                   boundary=str(boundary.resolve()) if boundary else None,
                   boundary_hash=boundary_digest(boundary) if boundary else None,
                   buffer_cells=args.buffer_cells, resampling=args.resampling,
                   missing_data=args.missing_data, algorithm_version=1)
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def boundary_digest(path):
    paths = sorted(path.parent.glob(path.stem + ".*")) if path.suffix.lower() == ".shp" else [path]
    digest = hashlib.sha256()
    for part in paths:
        digest.update(part.name.encode())
        digest.update(part.read_bytes())
    return digest.hexdigest()


def prepare(catalog, plan, window, args, name, boundary=None):
    stem = f"pr_{safe_name(name)}_{window.tag}"
    out = Path(args.output_dir).resolve()
    payload = out / "payloads" / stem
    payload.mkdir(parents=True, exist_ok=True)
    destination = out / (stem + ".dss")
    meta_path = payload / "meta.json"
    if destination.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists: {destination}. Choose another output directory or --overwrite.")
    meta_path.unlink(missing_ok=True)  # An interrupted job must never leave valid-looking metadata.
    count = 0
    with (payload / "values.bin").open("wb") as binary, (payload / "dates.csv").open("w", newline="", encoding="ascii") as table:
        writer = csv.writer(table)
        writer.writerow(["date", "start", "end", "min", "max", "mean"])
        for dates, values in iter_values(catalog, window, plan, payload):
            values = np.asarray(values, dtype="<f4")
            if not np.isfinite(values).all():
                raise ValueError("Float32 overflow in prepared rainfall.")
            binary.write(values.tobytes(order="C"))
            for date, array in zip(dates, values):
                stop = date + pd.Timedelta(seconds=catalog.step)
                writer.writerow([date.strftime("%Y-%m-%d"), date.isoformat(), stop.isoformat(),
                                 float(array.min()), float(array.max()), float(array.mean(dtype=np.float64))])
                count += 1
    expected = len(expected_dates(catalog, window))
    grid_meta = plan.target.dss_metadata()
    if count != expected or (payload / "values.bin").stat().st_size != count*grid_meta["nx"]*grid_meta["ny"]*4:
        raise ValueError("Payload record/byte count mismatch.")
    meta = dict(**grid_meta, values_bin=(payload / "values.bin").as_posix(),
                dates_csv=(payload / "dates.csv").as_posix(), out_dss=destination.as_posix(),
                template_dss=Path(args.template_dss).resolve().as_posix() if args.template_dss else "",
                a_part="", b_part="", c_part="PRECIPITATION", f_part="",
                data_units="MM", record_count=count, interval_seconds=catalog.step,
                start=window.start.isoformat(), end_exclusive=window.stop.isoformat(),
                source_calendar=catalog.calendar, source_units=catalog.units,
                time_policy=catalog.policy, missing_data=args.missing_data,
                nonfinite_cell_reads=catalog.invalid_cells,
                calendar_conversion=("monthly_total_preserving_per_subdaily_slot" if catalog.calendar == "360_day" else "none"),
                resampling=args.resampling, source_crs=catalog.grid.crs.to_string(),
                overwrite_dss=args.overwrite, source_files=catalog.file_info,
                fingerprint=fingerprint(catalog, plan, window, args, boundary))
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=True), encoding="ascii")
    print(f"Prepared {name} / {window.name}: {count} intervals, {grid_meta['nx']} x {grid_meta['ny']}, "
          f"{grid_meta['cell_size']:g} m", flush=True)
    return meta_path


def run_writer(meta_paths, args):
    if args.prepare_only:
        return
    if not args.jython_exe or not args.template_dss:
        raise ValueError("DSS writing requires --jython-exe and --template-dss; use --prepare-only to inspect payloads.")
    for path in (args.jython_exe, args.template_dss):
        if not Path(path).is_file():
            raise FileNotFoundError(path)
    script = Path(__file__).with_name("hec_dss_grid_writer_jython.py")
    out = Path(args.output_dir).resolve()
    manifest = out / "dss_jobs_manifest.json"
    manifest.write_text(json.dumps(dict(meta_paths=[str(p) for p in meta_paths],
                                        validation_csv=(out / "dss_validation.csv").as_posix()), indent=2), encoding="ascii")
    with (out / "dss_write.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen([str(Path(args.jython_exe).resolve()), "-here", str(script), str(manifest)],
                                   cwd=str(Path(args.jython_exe).resolve().parent), stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, errors="replace")
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        if process.wait() != 0:
            raise RuntimeError(f"HEC-DSS writer/readback validation failed; see {out / 'dss_write.log'}")


def build_parser(kind):
    parser = argparse.ArgumentParser(description="Convert metadata-described CPM/HadUK rainfall to HEC-RAS DSS.")
    parser.add_argument("--config", help="JSON configuration; command-line options override its values")
    parser.add_argument("--input-dir")
    parser.add_argument("--input-pattern", default="pr*.nc" if kind == "cpm" else "rainfall*.nc")
    parser.add_argument("--output-dir")
    parser.add_argument("--variable", default="pr" if kind == "cpm" else "rainfall")
    parser.add_argument("--engine", default=None, help="Optional xarray backend, e.g. netcdf4 or h5netcdf")
    parser.add_argument("--start")
    parser.add_argument("--end")
    parser.add_argument("--window", action="append", nargs=3, metavar=("ID", "START", "END"))
    parser.add_argument("--windows", nargs="+", help="Select IDs from configured/default windows")
    parser.add_argument("--regions", nargs="+", help="Region slugs; defaults to all with --boundary-dir")
    parser.add_argument("--boundary-dir")
    parser.add_argument("--boundary", help="One polygon file instead of the 17-region boundary directory")
    parser.add_argument("--boundary-crs", help="Only needed for boundary files with missing CRS")
    parser.add_argument("--name", default=kind, help="Dataset/region name for a single ungrouped input")
    parser.add_argument("--buffer-cells", type=int, default=1)
    parser.add_argument("--source-crs", help="Explicit source CRS when metadata is absent")
    parser.add_argument("--target-crs", help="Projected target CRS, e.g. EPSG:27700")
    parser.add_argument("--target-cell-size", type=float, help="Output square cell size in metres")
    parser.add_argument("--resampling", choices=["nearest", "bilinear", "average"], default="average")
    parser.add_argument("--select", action="append", default=[], help="Select a coordinate, e.g. ensemble_member=1")
    parser.add_argument("--units", help="Verified precipitation units override")
    parser.add_argument("--time-label", choices=["auto", "start", "centre", "end"], default="auto")
    parser.add_argument("--interval-hours", type=float, help="Interval when bounds/sufficient times are absent")
    parser.add_argument("--chunk-records", type=int, default=24)
    parser.add_argument("--missing-data", choices=["zero", "error"], default="zero")
    parser.add_argument("--jython-exe")
    parser.add_argument("--template-dss")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--inspect", action="store_true", help="Read metadata and print coverage/grid; write nothing")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def parse_args(kind, argv=None):
    parser = build_parser(kind)
    preliminary, _ = parser.parse_known_args(argv)
    if preliminary.config:
        config_path = Path(preliminary.config).resolve()
        config = json.loads(config_path.read_text(encoding="utf-8"))
        known = {action.dest for action in parser._actions} - {"help", "config"}
        if not isinstance(config, dict) or set(config) - known:
            raise ValueError("Unknown configuration keys: " + str(set(config)-known))
        for key in ("input_dir", "output_dir", "boundary_dir", "boundary", "jython_exe", "template_dss"):
            if config.get(key):
                path = Path(config[key])
                config[key] = str(path if path.is_absolute() else config_path.parent / path)
        import sys
        explicit = list(argv) if argv is not None else sys.argv[1:]
        for key in ("window", "select"):
            if "--" + key in explicit:
                config.pop(key, None)
        parser.set_defaults(**config)
    args = parser.parse_args(argv)
    if not args.input_dir or (not args.output_dir and not args.inspect):
        parser.error("--input-dir and --output-dir are required (output is optional with --inspect).")
    if args.chunk_records < 1 or args.buffer_cells < 0:
        parser.error("chunk-records must be positive and buffer-cells must be non-negative.")
    if args.boundary and args.boundary_dir:
        parser.error("Use either --boundary or --boundary-dir.")
    if args.target_cell_size is not None and args.target_cell_size <= 0:
        parser.error("--target-cell-size must be positive.")
    return args


def input_groups(args):
    files = sorted(Path(args.input_dir).glob(args.input_pattern))
    if not files:
        raise FileNotFoundError(f"No files match {args.input_pattern} under {args.input_dir}")
    if args.boundary or args.boundary_dir:
        return {safe_name(args.name): files}
    grouped = defaultdict(list)
    for path in files:
        matches = [slug for slug, _ in REGIONS if path.stem.startswith("pr_" + slug + "_")]
        grouped[matches[0] if matches else safe_name(args.name)].append(path)
    if args.regions:
        missing = set(args.regions) - set(grouped)
        if missing:
            raise ValueError(f"Requested regions absent from input: {sorted(missing)}")
        grouped = {key: value for key, value in grouped.items() if key in args.regions}
    return grouped


def main(kind="cpm", argv=None, clip=False):
    args = parse_args(kind, argv)
    windows = select_windows(args, kind)
    if clip and not (args.boundary or args.boundary_dir):
        raise ValueError("Regional clipping requires --boundary or --boundary-dir.")
    if not args.inspect and not clip and not args.prepare_only:
        for label in ("jython_exe", "template_dss"):
            if not getattr(args, label) or not Path(getattr(args, label)).is_file():
                raise ValueError(f"Set --{label.replace('_', '-')} or use --prepare-only.")
    jobs, summaries = [], []
    for group, paths in input_groups(args).items():
        catalog = Catalog(paths, args, kind)
        try:
            keys = sorted(catalog.records)
            if clip and not (args.start or args.window or args.windows):
                first, last = keys[0], keys[-1]
                start = pd.Timestamp(*first)
                if catalog.calendar == "360_day":
                    start = pd.Timestamp(first[0], first[1], 1)
                    last_day = calendar.monthrange(last[0], last[1])[1]
                    stop = pd.Timestamp(last[0], last[1], last_day) + pd.Timedelta(days=1)
                else:
                    stop = pd.Timestamp(*last) + pd.Timedelta(seconds=catalog.step)
                windows = [Window("all", start, stop)]
            print(json.dumps(dict(dataset=group, files=len(paths), records=len(keys),
                                  first=keys[0], last=keys[-1], calendar=catalog.calendar,
                                  interval_seconds=catalog.step, units=catalog.units,
                                  nx=len(catalog.grid.x), ny=len(catalog.grid.y),
                                  dx=catalog.grid.dx, dy=catalog.grid.dy,
                                  regular_axes=catalog.grid.regular,
                                  crs=catalog.grid.crs.to_string(), time_policy=catalog.policy), indent=2), flush=True)
            if args.inspect:
                continue
            if args.boundary_dir:
                requested = set(args.regions or [slug for slug, _ in REGIONS])
                unknown = requested - {slug for slug, _ in REGIONS}
                if unknown:
                    raise ValueError(f"Unknown regions: {sorted(unknown)}")
                boundaries = [(slug, find_boundary(args.boundary_dir, folder))
                              for slug, folder in REGIONS if slug in requested]
            else:
                boundaries = [(group, Path(args.boundary) if args.boundary else None)]
            for name, boundary in boundaries:
                plan = SpatialPlan(catalog.grid, args, boundary)
                try:
                    for window in windows:
                        meta = prepare(catalog, plan, window, args, name, boundary)
                        jobs.append(meta)
                        if clip:
                            summaries.append(payload_to_netcdf(meta, args, name, window))
                finally:
                    plan.close()
        finally:
            catalog.close()
    if not args.inspect:
        if not clip:
            run_writer(jobs, args)
        Path(args.output_dir, "preparation_summary.json").write_text(
            json.dumps(dict(metadata=[str(p) for p in jobs], netcdf=summaries), indent=2), encoding="utf-8")
        print(f"Completed: {len(jobs)} prepared job(s).", flush=True)
    return jobs


def payload_to_netcdf(meta_path, args, name, window):
    """Export the SAME prepared Gregorian intervals as CF NetCDF, retaining CRS/bounds.

    Stage 1 converts 360_day once. Stage 2 detects Gregorian and does not remap again.
    Native rotated grids are reprojected using the explicitly requested output grid.
    """
    meta = json.loads(meta_path.read_text(encoding="ascii"))
    rows = list(csv.DictReader(Path(meta["dates_csv"]).open(encoding="ascii")))
    dates = pd.DatetimeIndex([row["start"] for row in rows])
    ends = pd.DatetimeIndex([row["end"] for row in rows])
    values = np.memmap(meta["values_bin"], dtype="<f4", mode="r", shape=(len(rows), meta["ny"], meta["nx"]))
    output = Path(args.output_dir) / f"pr_{name}_{window.tag}.nc"
    if output.exists() and not args.overwrite:
        values._mmap.close()
        raise FileExistsError(output)
    try:
        ds = xr.Dataset({
            args.variable: (("time", "projection_y_coordinate", "projection_x_coordinate"), values,
                            {"units": "mm", "standard_name": "lwe_thickness_of_precipitation_amount", "grid_mapping": "crs", "cell_methods": "time: sum"}),
            "crs": ((), 0, {"spatial_ref": meta["srs_wkt"], "crs_wkt": meta["srs_wkt"]}),
            "time_bnds": (("time", "bnds"), np.column_stack([dates.values, ends.values])),
        }, coords={"time": dates.values,
                   "projection_x_coordinate": meta["cell_zero_x"]+(np.arange(meta["nx"])+.5)*meta["cell_size"],
                   "projection_y_coordinate": meta["cell_zero_y"]+(np.arange(meta["ny"])+.5)*meta["cell_size"]},
            attrs={"Conventions": "CF-1.8", "source_calendar": meta["source_calendar"],
                   "calendar_conversion": meta["calendar_conversion"], "history": "Prepared by dss_pipeline.py"})
        ds.time.attrs["bounds"] = "time_bnds"
        for axis in ("x", "y"):
            ds[f"projection_{axis}_coordinate"].attrs.update(standard_name=f"projection_{axis}_coordinate", units="m", axis=axis.upper())
        encoding = {args.variable: {"zlib": True, "complevel": 4, "chunksizes": (1, meta["ny"], meta["nx"])}}
        for var in ("time", "time_bnds"):
            encoding[var] = {"units": "minutes since 1970-01-01", "calendar": "proleptic_gregorian"}
        ds.to_netcdf(output, engine=args.engine or "netcdf4", encoding=encoding)
        ds.close()
    finally:
        values._mmap.close()
    return str(output)
