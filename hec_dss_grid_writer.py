"""Write prepared precipitation grids directly to DSS-7 with pydsstools."""
from __future__ import annotations

import csv
from datetime import datetime, timedelta
import json
import os
from pathlib import Path
import uuid

import numpy as np
from pyproj import CRS
from pydsstools.heclib.dss import HecDss
from pydsstools.core.gridinfo import DataType, GridType, SpecifiedInfo


def dss_time(value, *, end=False):
    """Use locale-independent dates and DSS's end-of-day convention."""
    midnight = end and value.time() == datetime.min.time()
    if midnight:
        value -= timedelta(days=1)
    months = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN",
              "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
    clock = "2400" if midnight else value.strftime("%H%M")
    return f"{value.day:02d}{months[value.month-1]}{value.year:04d}:{clock}"


def pathname(meta, start, end):
    return "/{}/{}/{}/{}/{}/{}/".format(
        meta.get("a_part", ""), meta.get("b_part", ""), "PRECIPITATION",
        dss_time(start), dss_time(end, end=True), meta.get("f_part", ""))


def read_intervals(meta):
    with Path(meta["dates_csv"]).open(encoding="ascii", newline="") as handle:
        rows = list(csv.DictReader(handle))
    intervals = [(datetime.fromisoformat(row["start"]), datetime.fromisoformat(row["end"]))
                 for row in rows]
    if not intervals or len(intervals) != int(meta["record_count"]):
        raise ValueError("Date/payload record count mismatch.")
    previous = None
    for start, end in intervals:
        if (end - start).total_seconds() != meta["interval_seconds"] or (
                previous is not None and start != previous):
            raise ValueError("Accumulation intervals must be contiguous and have a consistent duration.")
        previous = end
    if (intervals[0][0].isoformat() != meta["start"] or
            intervals[-1][1].isoformat() != meta["end_exclusive"]):
        raise ValueError("Accumulation intervals disagree with the requested window.")
    count = len(intervals) * int(meta["nx"]) * int(meta["ny"])
    if Path(meta["values_bin"]).stat().st_size != count * 4:
        raise ValueError("Binary rainfall payload byte count mismatch.")
    return intervals


def validate_dss(path, meta, intervals, values):
    expected_paths = [pathname(meta, start, end) for start, end in intervals]
    expected_crs, seen_crs = CRS.from_wkt(meta["srs_wkt"]), set()
    with HecDss.Open(str(path), mode="r") as dss:
        actual = dss.search_path()
        if set(actual) != set(expected_paths) or len(actual) != len(intervals):
            raise ValueError("DSS catalog disagrees with the precipitation intervals.")
        for index, record in enumerate(expected_paths):
            read_values = index in {0, len(intervals)-1}
            grid = dss.read_grid(record, metadata_only=not read_values)
            info = grid.gridinfo
            if (info.grid_type != GridType.specified_time or
                    info.data_type != DataType.per_cum or
                    info.data_units.strip().upper() != "MM" or
                    tuple(info.shape) != (meta["ny"], meta["nx"]) or
                    tuple(info.lower_left_cell) != (meta["lower_left_cell_x"], meta["lower_left_cell_y"])):
                raise ValueError("DSS grid dimensions, offsets, units or type disagree: " + record)
            if not np.allclose([info.cell_size, *info.coords_cell0],
                               [meta["cell_size"], meta["cell_zero_x"], meta["cell_zero_y"]],
                               rtol=1e-7, atol=1e-3):
                raise ValueError("DSS grid spacing/origin disagree: " + record)
            if info.crs not in seen_crs:
                if not CRS.from_user_input(info.crs).equals(expected_crs):
                    raise ValueError("DSS coordinate reference system disagrees: " + record)
                seen_crs.add(info.crs)
            if read_values:
                # pydsstools returns north-up rasters; payload rows increase northwards.
                actual_values = grid.read()
                np.testing.assert_allclose(actual_values[::-1], values[index], rtol=1e-6, atol=1e-6,
                                           err_msg="DSS rainfall array readback disagrees with the payload.")


def write_one(meta_path):
    meta = json.loads(Path(meta_path).read_text(encoding="ascii"))
    intervals = read_intervals(meta)
    destination = Path(meta["out_dss"])
    if destination.exists() and not meta.get("overwrite_dss", False):
        raise FileExistsError(f"DSS already exists; use --overwrite: {destination}")
    temporary = destination.with_name(destination.stem + ".building-" + uuid.uuid4().hex + ".dss")
    info = SpecifiedInfo(
        grid_type=GridType.specified_time, data_type=DataType.per_cum, data_units="MM",
        shape=(meta["ny"], meta["nx"]), cell_size=meta["cell_size"],
        lower_left_cell=(meta["lower_left_cell_x"], meta["lower_left_cell_y"]),
        coords_cell0=(meta["cell_zero_x"], meta["cell_zero_y"]),
        crs=meta["srs_wkt"], crs_name=meta["srs_name"], nodata=-9999.,
        is_interval=True, time_stamped=True)
    values = np.memmap(meta["values_bin"], dtype="<f4", mode="r",
                       shape=(len(intervals), meta["ny"], meta["nx"]))
    try:
        with HecDss.Open(str(temporary), version=7) as dss:
            for index, (start, end) in enumerate(intervals):
                # The Cython writer needs a writable buffer; keep the memmap read-only.
                data = np.array(values[index], dtype=np.float32, copy=True, order="C")
                if not np.isfinite(data).all() or (data < 0).any():
                    raise ValueError("Invalid precipitation in the prepared payload.")
                # Both payload and DSS storage start at the southwest corner.
                dss.put_grid(data, pathname(meta, start, end), info,
                             flipud=False, normalize=False)
                number = index + 1
                if number == 1 or number % 500 == 0 or number == len(intervals):
                    print(f"Wrote {number}/{len(intervals)} {destination}", flush=True)
        validate_dss(temporary, meta, intervals, values)
        os.replace(temporary, destination)
    finally:
        values._mmap.close()
        temporary.unlink(missing_ok=True)
    return destination
