# Validation record — 2026-10-27

## Automated regression suite

`python -m unittest discover -s tests -v`: **23 tests passed** locally.
The same shared module and suite are included in both standalone repositories.
Tests cover arbitrary and overlapping windows, inclusive date-only endpoints,
exclusive datetime endpoints, leap days, hourly accumulation bounds, native
unit conversion, monthly total conservation for 28/29/30/31 days, subdaily
slot conservation, partial output months, missing/duplicate coverage,
ensemble selection, axis reversal, 1/2.2/5/12/25 km grids, fractional origins,
polygon clipping, boundary-edge buffering, rotated grids, stretched-grid
pixel-centre registration, metadata-only inspection and CF NetCDF round trips.

## Real local data

Source files were opened read-only. Small spatial/temporal subsets were used
for functional checks; this was not a rerun of every region for all decades.

| Input | Selected range | Result |
|---|---|---|
| Native UKCP18 2.2 km, rotated coordinates with stretched halo, 360_day, mm/day | 2020-02-03 through 2020-02-07 | 5 daily records; 8 x 8 output cells at 2200 m after EPSG:27700 reprojection |
| UKCP18 5 km projected grid, 360_day, mm/day | 2020-02-03 through 2020-02-07 | 5 daily records; 4 x 4 output cells at 5000 m |
| HadUK-Grid 5 km daily amounts | 1969-11-12 through 1969-11-18 | 7 unchanged nominal daily labels; 4 x 4 output cells at 5000 m |

The native 2.2 km and HadUK cases were actually written with **HEC-DSSVue
3.4.15** and read back. Catalog paths, all grid headers, WKT, units, start/end
times, and first/last rainfall arrays passed the writer's validation.

A separate synthetic **360_day hourly** CF file in `kg m-2 s-1` was passed
through the real HEC-DSSVue writer. The requested Gregorian interval
2020-01-02 03:00 to 06:00 produced exactly three one-hour MM/PER-CUM records;
header times and first/last arrays passed readback. No real hourly source
file was available, so real-product hourly metadata remain an integration
case for the user to check with `--inspect` before a long run.

## Full original DSS examples

Five complete Forth files were copied from the supplied output directories.
They are historical outputs, not results of the revised implementation.

| Dataset/window | Inclusive daily range | Catalog records |
|---|---|---:|
| CPM 01 | 2019-11-01 to 2039-12-31 | 7366 |
| CPM 02 | 2039-11-01 to 2059-12-31 | 7366 |
| CPM 03 | 2059-11-01 to 2080-11-30 | 7701 |
| HadUK 01 | 1979-11-01 to 1999-12-31 | 7366 |
| HadUK 02 | 1999-11-01 to 2019-12-31 | 7366 |

All use 35 x 17 cells, 5000 m, MM and PER-CUM. Native HEC-DSSVue catalog
inspection confirmed the first/last dates. SHA-256 checksums and byte sizes
are recorded in `examples/checksums.json`; copied bytes match the originals.
Only the DSS files are placed under `examples/output_dss/`.

## Scope and scientific conventions

Monthly total preservation is checked in prepared floating-point rainfall.
DSS compression can round individual values (approximately 0.01 mm).
Calendar conversion changes event chronology. Average spatial resampling is
not an exact spherical conservative-remapping guarantee. HadUK's original
nominal-day labelling is deliberately retained, although its physical
observation bounds are 09:00 to 09:00.

## Implementation references

- [HEC grid storage order and square-cell metadata](https://www.hec.usace.army.mil/confluence/dssdocs/dssjavaprogrammer/appendix-c-grid-classes-in-hec-java-programs/extent-functions)
- [HEC specified-grid origins and offsets](https://www.hec.usace.army.mil/confluence/metdoc/metum/3.4.0/general-information-and-tips/dss-specifiedgridinfo-in-hec-metvue)
- [GDAL geolocation conventions](https://gdal.org/en/stable/api/gdal_alg.html)
- [Rasterio reprojection](https://rasterio.readthedocs.io/en/stable/api/rasterio.warp.html)
