# UKCP18 CPM to HEC-RAS DSS

Configurable date ranges, multiple time windows, daily/hourly rainfall and
metadata-driven grid resolution. Sources are read-only. DSS grids use MM,
PER-CUM and explicit accumulation bounds.

中文使用说明：[README_zh.md](README_zh.md)。

## Quick start

```bash
python 02_ukcp18_cpm_to_hecras_dss.py --input-dir "D:/climate/cpm" --inspect
python 02_ukcp18_cpm_to_hecras_dss.py --input-dir "D:/climate/cpm" --output-dir "D:/output/cpm" --start 2020-02-03 --end 2020-02-07 --prepare-only
python 02_ukcp18_cpm_to_hecras_dss.py --config config.example.json
```

For actual DSS output remove `--prepare-only` and supply `--jython-exe` and
`--template-dss`. The example configuration starts in preparation-only mode.

## Native 2.2 km CPM and optional regional NetCDF

```bash
python 02_ukcp18_cpm_to_hecras_dss.py --input-dir "D:/climate/cpm_2.2km" --output-dir "D:/output/cpm" --boundary "D:/boundaries/forth.shp" --name forth --target-crs EPSG:27700 --target-cell-size 2200 --start 2020-02-01 --end 2020-02-29 --prepare-only
python 02_ukcp18_cpm_to_hecras_dss.py --input-dir "D:/climate/cpm" --output-dir "D:/output/windows" --window event 2020-02-03 2020-02-10 --window future 2040-05-01 2040-05-20 --prepare-only
```

Gregorian inputs retain their dates. `360_day` data are converted month by
month: interpolate 30 days to the Gregorian month length, then scale each
output cell to preserve its monthly rainfall total. For hourly data, remap
each hour-of-day independently, preserving its monthly amount and diurnal
phase. This changes event chronology; it is not a physical date-equivalence
mapping. Complete source months are required even for a partial output month.

The original default DSS windows are 2019-11-01..2039-12-31,
2039-11-01..2059-12-31, and 2059-11-01..2080-11-30.

The optional `01_clip_ukcp18_cpm_regions.py` now uses the same preparation
engine and exports **Gregorian interval accumulation in mm** as CF NetCDF,
including CRS and time bounds. Without date options it exports one full-coverage
file per region. This replaces the old hardcoded `pr_<region>_2010_2080.nc`
clipping stage. Conversion from `360_day` occurs here once; stage 2 detects
the Gregorian output and does not convert it again. Use `--boundary` or
`--boundary-dir`; rotated inputs also need the explicit target grid.
When deliberately exporting overlapping NC windows, select one window's files
per DSS run to avoid duplicate source records.

## Setup

Python 3.11 or newer:

```bash
python -m pip install -r requirements.txt
```

Writing DSS additionally requires HEC-DSSVue 3.4.x `Jython.exe` and a template
DSS with a precipitation GridContainer. A complete bundled example can serve
as the template. Python preparation and tests do not require HEC-DSSVue.

## Time selection

Supply `--start YYYY-MM-DD --end YYYY-MM-DD` for any inclusive date range.
The output record on the last day ends the following day. For hourly/subdaily
data, datetime bounds are **start-inclusive, end-exclusive**, for example
`--start 2020-01-02T03:00 --end 2020-01-02T09:00` writes six hourly records.
Datetime endpoints must align with the source accumulation intervals.

Use repeated `--window ID START END` to write independent DSS files, including
overlapping windows. IDs must be unique. With no dates supplied, the original
conversion windows remain the defaults; `--windows 01 02` filters their IDs.
Dates and output names have no fixed year range or `20y` assumption.

## Grids, units and validation

- Detect spatial axes, CF `grid_mapping`, units, time bounds, calendar and native
  interval from metadata. Process different products/ensembles in separate runs.
- Native regular projected square grids in metres retain their resolution
  automatically, including 1, 2.2, 5, 12 or 25 km. Descending axes are reversed
  together with the data. Fractional grid origins are preserved.
- For geographic/rotated or stretched grids, explicitly specify a projected
  `--target-crs` and `--target-cell-size` in metres. Reprojection transforms
  values as well as coordinates. Native CPM stretched halo cells use GDAL
  geolocation metadata with pixel-centre registration, rather than a false affine.
- Resampling: `average` (default), `nearest`, or `bilinear`. Area averaging
  is a practical depth resampling method; it is not an exact spherical
  conservative-remapping guarantee. Spatial resampling changes local values.
- Amounts in `mm`, `m` or `kg m-2` and rates such as `mm/day`, `mm h-1` and
  `kg m-2 s-1` are converted to **MM per native accumulation interval**.
  Hourly rainfall is kept hourly; no automatic daily disaggregation is invented.
- Multi-member dimensions require `--select ensemble_member=VALUE`; singleton
  dimensions are removed and their values are checked across source files.
- Missing/nonfinite rainfall becomes zero by default, preserving the previous
  workflows. Set `--missing-data error` for strict rejection. Negative rainfall
  is rejected. Out-of-boundary cells are zero; `--buffer-cells` controls padding.
- Regional-mean products without two spatial axes cannot be used as spatial
  rainfall grids. Curvilinear/unstructured meshes are not supported.
- Missing source coverage, overlaps, changing grids/units/calendars, unknown
  CRS, and monthly/seasonal mean inputs fail clearly. No zero-filled time gaps.
  Old clipped NC files without CRS require a verified `--source-crs EPSG:27700`.
- If time bounds are absent, explicitly provide `--time-label start|centre|end`;
  a single-record file also needs `--interval-hours`.

## Regions and large inputs

Use `--boundary path/to/region.shp --name forth` for one polygon, or
`--boundary-dir path/to/region_folders --regions forth tay` for the original
17-region layout. Without either option, the input grid is used directly;
recognized regional `pr_<slug>_*.nc` files are grouped separately, independent
of their date suffix. A missing boundary CRS requires `--boundary-crs`.

Files are indexed by their time coordinates, sorted chronologically, and opened
read-only. Only one rainfall source file is cached. Preparation reads spatial
crops and `--chunk-records` (default 24); calendar conversion works on one
30-day subdaily slot at a time and uses a disk-backed output month. Reserve
temporary space for a month plus the requested binary payloads.

`--inspect` prints metadata without writing files. `--prepare-only` builds the
binary payloads, interval CSVs and metadata JSON without HEC-DSSVue. Payloads
are retained for audit/retry. Existing DSS files require explicit `--overwrite`.
The writer builds a separate temporary DSS, validates every record's catalog,
dimensions, offsets, origin, cell size, WKT, units and accumulation bounds,
and first/last arrays, then replaces the destination. A failed write preserves
the previous DSS. The template may round rainfall by approximately 0.01 mm.

## JSON and Notebook

Copy `config.example.json`, edit the input/output paths and windows, then run
the entry script with `--config your_config.json`. Relative paths in JSON are
resolved relative to that JSON file. CLI settings override JSON settings.

The notebooks call the same `dss_pipeline.py` as the command-line scripts.
Edit their CONFIG cell, inspect metadata, and set `RUN = True` to execute.
`PREPARE_ONLY = True` is the initial setting; turn it off and provide the
HEC-DSSVue executable/template to produce DSS.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests create small synthetic CF files and check dates, hourly intervals,
360-day monthly amounts, leap days, units, metadata, grid direction, regional
grouping, polygon edges, stretched/rotated grids and NetCDF round trips.
CI runs the same suite on Windows and Linux. HEC-DSSVue integration needs
a local installation and is recorded separately in `VALIDATION.md`.

## Complete DSS output examples

`examples/output_dss/` contains only complete, previously generated **Forth**
DSS files, copied byte-for-byte from the supplied local output directories.
These are original workflow outputs, not regenerated with the new code.
All original time windows are included. See `examples/README.md` for sizes,
dates, record counts and SHA-256 digests. No large NC inputs, model HDFs,
boundary data, logs or temporary payloads are bundled.
