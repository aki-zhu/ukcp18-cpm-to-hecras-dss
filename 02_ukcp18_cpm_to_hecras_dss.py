# -*- coding: utf-8 -*-
"""将 UKCP18 CPM 360_day 降雨 NetCDF 转为 HEC-RAS DSS。
Convert UKCP18 CPM 360_day rainfall NetCDF files to HEC-RAS DSS grids.
"""

from __future__ import annotations
from pathlib import Path
import argparse, calendar, csv, json, subprocess, time, tempfile
import h5py
import numpy as np
import pandas as pd
import xarray as xr

INPUT_DIR = Path("PATH_TO_CLIPPED_REGION_NETCDF_DIRECTORY")
OUT_DIR = Path("PATH_TO_DSS_OUTPUT_DIRECTORY")
DSSVUE_HOME = Path("PATH_TO_HEC_DSSVUE_DIRECTORY")
JYTHON_EXE = DSSVUE_HOME / "Jython.exe"
TEMPLATE_DSS = Path("PATH_TO_TEMPLATE_DSS_WITH_A_PRECIPITATION_GRID")
RAS_GEOM_HDF = Path("PATH_TO_HECRAS_GEOMETRY_HDF_WITH_PROJECTION")
WRITER_SCRIPT = OUT_DIR / "write_region_dss_grid_origin00_mm.py"
LOG_FILE = OUT_DIR / "batch_region_nc_to_dss_20y.log"

WINDOWS = [
    ("01", "2019-11-01", "2039-12-31"),
    ("02", "2039-11-01", "2059-12-31"),
    ("03", "2059-11-01", "2080-11-30"),
]

JYTHON_WRITER = r'''
from hec.heclib.dss import HecDss
from java.nio import ByteBuffer, ByteOrder
from java.io import File
from jarray import zeros
import csv, json, sys
from datetime import datetime, timedelta

def dss_time(d):
    return d.strftime('%d%b%Y').upper() + ':0000'

def read_float_array(fh, n):
    b = fh.read(n * 4)
    if len(b) != n * 4:
        raise Exception('Expected %d bytes, got %d' % (n * 4, len(b)))
    bb = ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN)
    arr = zeros(n, 'f')
    for i in xrange(n):
        arr[i] = bb.getFloat()
    return arr

def find_template_grid(template_dss):
    src = HecDss.open(template_dss)
    try:
        for p in src.getPathnameList():
            s = str(p)
            if '/PRECIPITATION/' in s.upper():
                try:
                    obj = src.get(s)
                    if obj.__class__.__name__.endswith('GridContainer'):
                        print 'Template grid:', s
                        return obj
                except Exception, e:
                    pass
        raise Exception('No readable GridContainer template found in ' + template_dss)
    finally:
        src.close()

def make_path(meta, start_s, end_s):
    return '/%s/%s/%s/%s/%s/%s/' % (
        meta.get('a_part', ''),
        meta.get('b_part', ''),
        meta.get('c_part', 'PRECIPITATION'),
        start_s,
        end_s,
        meta.get('f_part', '')
    )

def main(meta_path):
    with open(meta_path, 'r') as f:
        meta = json.load(f)
    rows = list(csv.DictReader(open(meta['dates_csv'], 'r')))
    n = int(meta['nx']) * int(meta['ny'])

    out_file = File(meta['out_dss'])
    if out_file.exists():
        out_file.delete()

    obj = find_template_grid(meta['template_dss'])
    gi = obj.getGridInfo()
    gi.setGridType(430)
    gi.setCellInfo(int(meta['lower_left_cell_x']), int(meta['lower_left_cell_y']), int(meta['nx']), int(meta['ny']), float(meta['cell_size']))
    gi.setSpatialReference(meta['srs_name'], meta['srs_wkt'], float(meta['cell_zero_x']), float(meta['cell_zero_y']))
    gi.setCoordOfGridCellZero(float(meta['cell_zero_x']), float(meta['cell_zero_y']))
    gi.setDataType(1)  # PER-CUM
    gi.setDataUnits(meta.get('data_units', 'MM'))

    dst = HecDss.open(meta['out_dss'])
    written = 0
    try:
        fh = open(meta['values_bin'], 'rb')
        try:
            for row in rows:
                start = datetime.strptime(row['date'], '%Y-%m-%d')
                end = start + timedelta(days=1)
                start_s = dss_time(start)
                end_s = dss_time(end)

                gd = obj.getGridData()
                gd.setData(read_float_array(fh, n))
                gd.initGridInfo(gi)
                try:
                    gd.updateStatistics()
                except Exception, e:
                    pass

                gi.setGridTimes(start_s, end_s)
                gi.setDataInfo(float(row['min']), float(row['max']), float(row['mean']))

                full = make_path(meta, start_s, end_s)
                obj.fullName = full
                obj.setFullName(full)
                dst.put(obj)
                written += 1
                if written <= 1 or written % 1000 == 0 or written == len(rows):
                    print 'Wrote %d/%d %s' % (written, len(rows), full)
        finally:
            fh.close()
    finally:
        dst.close()

    chk = HecDss.open(meta['out_dss'])
    try:
        paths = chk.getPathnameList()
        print 'Catalog records:', len(paths)
        if len(paths) > 0:
            first = str(paths[0])
            g = chk.get(first).getGridInfo()
            print 'First path:', first
            print 'Grid:', g.getNumberOfCellsX(), 'x', g.getNumberOfCellsY(), 'cell', g.getCellSize(), 'units', g.getDataUnits(), 'type', g.getDataType()
            print 'Grid time:', g.getStartTime(), '|', g.getEndTime(), '|', g.getTimePeriod()
            print 'Cell zero:', g.getXCoordOfGridCellZero(), g.getYCoordOfGridCellZero()
            print 'Lower-left cell:', g.getLowerLeftCellX(), g.getLowerLeftCellY()
    finally:
        chk.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception('Usage: write_region_dss_grid_origin00_mm.py meta.json')
    main(sys.argv[1])
'''

def log(message):
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = "[%s] %s" % (stamp, message)
    print(line, flush=True)
    with LOG_FILE.open("a", encoding="utf-8") as f:
        f.write(line + "\n")

def as_posix(path):
    return str(path).replace("\\", "/")

def load_model_projection():
    with h5py.File(RAS_GEOM_HDF, "r") as hdf:
        proj = hdf.attrs["Projection"]
    return proj.decode("utf-8") if isinstance(proj, bytes) else str(proj)

def region_slug_from_nc(nc_path):
    name = nc_path.stem
    if name.startswith("pr_") and name.endswith("_2010_2080"):
        return name[3:-10]
    return name

def grid_origin_from_edges(xx, yy):
    dx = float(np.median(np.diff(xx)))
    dy = float(np.median(np.diff(yy)))
    if abs(dx - dy) > 1e-6:
        raise ValueError("DSS specified grids need square cells; got dx=%s, dy=%s" % (dx, dy))
    x_edge = float(xx[0] - dx / 2.0)
    y_edge = float(yy[0] - dy / 2.0)
    llx = int(round(x_edge / dx))
    lly = int(round(y_edge / dy))
    cell_zero_x = float(x_edge - llx * dx)
    cell_zero_y = float(y_edge - lly * dy)
    return dx, x_edge, y_edge, llx, lly, cell_zero_x, cell_zero_y

def build_month_index(time_values):
    months = {}
    for i, t in enumerate(time_values):
        months.setdefault((int(t.year), int(t.month)), []).append(i)
    return {key: np.asarray(value, dtype=np.int64) for key, value in months.items()}

def interp_360_month_to_gregorian(month_values, n_days):
    if month_values.shape[0] != 30:
        raise ValueError("Expected 30 source days in 360_day month, got %s" % month_values.shape[0])
    src = np.nan_to_num(month_values.astype("float32", copy=False), nan=0.0)
    target = np.linspace(0.0, 29.0, n_days, dtype=np.float32)
    lo = np.floor(target).astype(np.int16)
    hi = np.ceil(target).astype(np.int16)
    w = (target - lo).astype("float32")[:, None, None]
    out = src[lo] * (1.0 - w) + src[hi] * w
    monthly_total = src.sum(axis=0, dtype="float64")
    out_sum = out.sum(axis=0, dtype="float64")
    scale = np.divide(monthly_total, out_sum, out=np.zeros_like(monthly_total), where=out_sum > 0)
    return (out * scale.astype("float32")).astype("float32", copy=False)

def gregorian_window_values(da, month_index, start, end):
    start_ts = pd.Timestamp(start)
    end_ts = pd.Timestamp(end)
    months = pd.period_range(start_ts.to_period("M"), end_ts.to_period("M"), freq="M")
    date_chunks = []
    value_chunks = []
    for period in months:
        key = (int(period.year), int(period.month))
        if key not in month_index:
            raise KeyError("Source 360_day month not found: %s" % (key,))
        idx = month_index[key]
        month_raw = da.isel(time=idx).values
        n_days = calendar.monthrange(int(period.year), int(period.month))[1]
        month_values = interp_360_month_to_gregorian(month_raw, n_days)
        month_dates = pd.date_range(period.start_time.normalize(), period.end_time.normalize(), freq="D")
        keep = (month_dates >= start_ts) & (month_dates <= end_ts)
        date_chunks.append(month_dates[keep])
        value_chunks.append(month_values[keep])
    dates = date_chunks[0].append(date_chunks[1:]) if len(date_chunks) > 1 else date_chunks[0]
    values = np.concatenate(value_chunks, axis=0)
    return dates, values

def build_dss_filename(region_slug, window_id, start, end):
    return OUT_DIR / ("pr_%s_20y_%s_%s_%s.dss" % (region_slug, window_id, start.replace("-", ""), end.replace("-", "")))

def make_payload(nc_path, region_slug, da, month_index, xx, yy, window_id, start, end, srs_wkt):
    payload_dir = OUT_DIR / "payloads" / ("%s_%s_%s_%s" % (region_slug, window_id, start.replace("-", ""), end.replace("-", "")))
    payload_dir.mkdir(parents=True, exist_ok=True)
    out_dss = build_dss_filename(region_slug, window_id, start, end)
    if out_dss.exists():
        log("Skip existing %s" % out_dss)
        return payload_dir / "meta.json", out_dss, True

    dates, values = gregorian_window_values(da, month_index, start, end)
    values = np.nan_to_num(values.astype("float32", copy=False), nan=0.0)
    values_bin = payload_dir / "values_f32_le.bin"
    values.astype("<f4", copy=False).tofile(values_bin)

    dates_csv = payload_dir / "dates.csv"
    with dates_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["date", "min", "max", "mean"])
        for i, date in enumerate(dates):
            arr = values[i]
            writer.writerow([date.strftime("%Y-%m-%d"), float(arr.min()), float(arr.max()), float(arr.mean())])

    cell_size, x_edge, y_edge, llx, lly, cell_zero_x, cell_zero_y = grid_origin_from_edges(xx, yy)
    meta = {
        "region_slug": region_slug,
        "nx": int(len(xx)),
        "ny": int(len(yy)),
        "cell_size": float(cell_size),
        "x_edge_min": x_edge,
        "y_edge_min": y_edge,
        "x_edge_max": float(xx[-1] + cell_size / 2.0),
        "y_edge_max": float(yy[-1] + cell_size / 2.0),
        "lower_left_cell_x": llx,
        "lower_left_cell_y": lly,
        "cell_zero_x": cell_zero_x,
        "cell_zero_y": cell_zero_y,
        "values_bin": as_posix(values_bin),
        "dates_csv": as_posix(dates_csv),
        "out_dss": as_posix(out_dss),
        "template_dss": as_posix(TEMPLATE_DSS),
        "a_part": "",
        "b_part": "",
        "c_part": "PRECIPITATION",
        "f_part": "",
        "data_units": "MM",
        "srs_name": "OSGB36_British_National_Grid",
        "srs_wkt": srs_wkt,
        "source_nc": as_posix(nc_path),
        "source_calendar": "360_day",
        "calendar_conversion_method": "monthly_total_preserving_360_day_to_gregorian",
        "variable": "pr",
        "source_units": "mm/day",
        "stored_units": "MM",
        "value_scale": 1.0,
        "start": start,
        "end": end,
        "record_count": int(len(dates)),
        "first_date": dates[0].strftime("%Y-%m-%d"),
        "last_date": dates[-1].strftime("%Y-%m-%d"),
    }
    meta_path = payload_dir / "meta.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log("%s window %s: %s records, grid %sx%s, lower-left=(%s,%s), cell-zero=(%.1f,%.1f)" % (region_slug, window_id, meta["record_count"], meta["nx"], meta["ny"], llx, lly, cell_zero_x, cell_zero_y))
    return meta_path, out_dss, False

def run_writer(meta_path, out_dss, skipped):
    if skipped:
        return
    cmd = [str(JYTHON_EXE), "-here", str(WRITER_SCRIPT), str(meta_path)]
    result = subprocess.run(cmd, cwd=str(DSSVUE_HOME), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    with LOG_FILE.open("a", encoding="utf-8") as f:
        f.write(result.stdout + "\n")
    if result.returncode != 0:
        log(result.stdout[-4000:])
        raise RuntimeError("Jython writer failed with exit code %s for %s" % (result.returncode, out_dss))
    for line in result.stdout.splitlines():
        if any(key in line for key in ["Catalog records:", "First path:", "Grid:", "Grid time:", "Cell zero:", "Lower-left cell:"]):
            log(line)

def convert_region(nc_path, srs_wkt, only_windows=None):
    region_slug = region_slug_from_nc(nc_path)
    log("Opening %s: %s" % (region_slug, nc_path))
    ds = xr.open_dataset(nc_path, decode_times=True)
    try:
        da = ds["pr"]
        if "ensemble_member" in da.dims:
            da = da.isel(ensemble_member=0, drop=True)
        da = da.transpose("time", "projection_y_coordinate", "projection_x_coordinate").astype("float32")
        xx = ds["projection_x_coordinate"].values.astype(float)
        yy = ds["projection_y_coordinate"].values.astype(float)
        if xx[1] < xx[0]:
            da = da.isel(projection_x_coordinate=slice(None, None, -1)); xx = xx[::-1]
        if yy[1] < yy[0]:
            da = da.isel(projection_y_coordinate=slice(None, None, -1)); yy = yy[::-1]
        month_index = build_month_index(da["time"].values)
        outputs = []
        for window_id, start, end in WINDOWS:
            if only_windows is not None and window_id not in only_windows:
                continue
            meta_path, out_dss, skipped = make_payload(nc_path, region_slug, da, month_index, xx, yy, window_id, start, end, srs_wkt)
            run_writer(meta_path, out_dss, skipped)
            outputs.append(out_dss)
            log("Finished %s window %s: %s" % (region_slug, window_id, out_dss))
        return outputs
    finally:
        ds.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--regions", nargs="*")
    parser.add_argument("--windows", nargs="*", choices=[w[0] for w in WINDOWS])
    return parser.parse_args()

def main():
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not LOG_FILE.exists():
        LOG_FILE.write_text("", encoding="utf-8")
    WRITER_SCRIPT.write_text(JYTHON_WRITER, encoding="utf-8")
    for path in [INPUT_DIR, JYTHON_EXE, TEMPLATE_DSS, RAS_GEOM_HDF]:
        if not path.exists():
            raise FileNotFoundError(path)
    region_filter = {r.lower() for r in args.regions} if args.regions else None
    window_filter = set(args.windows) if args.windows else None
    nc_files = sorted(INPUT_DIR.glob("pr_*_2010_2080.nc"))
    if region_filter is not None:
        nc_files = [p for p in nc_files if region_slug_from_nc(p).lower() in region_filter]
    if not nc_files:
        raise RuntimeError("No matching region NetCDF files found.")
    log("Input NetCDF files: %s" % len(nc_files))
    log("Output folder: %s" % OUT_DIR)
    log("Naming: pr_<region>_20y_<window>_<start>_<end>.dss")
    srs_wkt = load_model_projection()
    outputs = []
    for nc_path in nc_files:
        outputs.extend(convert_region(nc_path, srs_wkt, only_windows=window_filter))
    log("All requested DSS outputs complete: %s" % len(outputs))
    for path in outputs:
        log("  " + str(path))

if __name__ == "__main__":
    main()
