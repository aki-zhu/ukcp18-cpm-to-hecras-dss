# UKCP18 CPM to HEC-RAS DSS

## Overview

This repository converts UKCP18 CPM 5 km daily gridded rainfall NetCDF files into HEC-DSS specified-grid precipitation inputs for HEC-RAS. It includes regional NetCDF clipping for 17 UK regions and a monthly-total-preserving `360_day`-to-Gregorian calendar conversion.

The original UKCP18 NetCDF files are read-only and are never modified. Data, boundary shapefiles, the template DSS, and the HEC-RAS geometry HDF are not included.

## Two-stage pipeline

### 1. Regional NetCDF

1. Read the 5 km British National Grid from the first UKCP18 CPM file.
2. Select one polygon boundary per region and transform it to EPSG:27700.
3. Rasterize and buffer each boundary by one grid cell.
4. Apply the mask and minimum rectangular crop to every source file.
5. Concatenate along time and write `pr_<region>_2010_2080.nc`.

### 2. NetCDF to DSS

1. Check `pr`, spatial coordinates, grid direction, and the `360_day` time axis.
2. Interpolate each 30-day model month to the actual Gregorian month length.
3. Scale each cell so the interpolated monthly accumulation exactly matches the source month.
4. Store NoData and out-of-region cells as `0 mm`.
5. Build daily little-endian float32 payloads.
6. Use HEC-DSSVue Jython to write `MM`, `PER-CUM` specified grids.
7. Log and verify the DSS catalog, grid dimensions, projection, units, and times.

## Files

| File | Description |
|---|---|
| `01_clip_ukcp18_cpm_regions.ipynb` | 17-region clipping notebook |
| `01_clip_ukcp18_cpm_regions.py` | Regional clipping Python script |
| `02_ukcp18_cpm_to_hecras_dss.ipynb` | Calendar conversion and DSS-writing notebook |
| `02_ukcp18_cpm_to_hecras_dss.py` | CLI script with region and window filters |
| `hec_dss_grid_writer_jython.py` | HEC-DSSVue Jython grid writer |

## Input assumptions

Source CPM files:

- Variable: `pr`
- Units: `mm/day`
- Calendar: `360_day`
- Coordinates: `projection_x_coordinate`, `projection_y_coordinate`
- Grid: 5000 m British National Grid
- Example pattern: `pr_rcp85_land-cpm_uk_5km_01_day_*.nc`

Regional filename pattern:

```text
pr_<region>_2010_2080.nc
```

Regions:

```text
Anglian, Clyde, Dee, Forth, Humber, Northeast Scotland,
North Highland, Northwest England, Northumbria, Severn, Solway,
Southeast England, Southwest England, Tay, Thames, Tweed, Western Wales
```

## Calendar and units

Source `pr` is already in `mm/day`. The workflow does not divide values by 24 or 1000. Each daily DSS record uses units `MM` and data type `PER-CUM`.

The `360_day` conversion is performed month by month: 30 source days are interpolated to 28, 29, 30, or 31 Gregorian days, then corrected cell-by-cell to conserve the monthly accumulation. This preserves monthly water volume while changing the shape of the daily sequence, so the method should be documented in model metadata.

Default output windows:

- `2019-11-01` to `2039-12-31`
- `2039-11-01` to `2059-12-31`
- `2059-11-01` to `2080-11-30`

The windows intentionally overlap, but no individual DSS file contains duplicate dates.

## Environment

Python 3 dependencies:

```text
xarray
h5netcdf
h5py
cftime
numpy
pandas
geopandas
rasterio
matplotlib
```

DSS writing also requires `Jython.exe` from HEC-DSSVue 3.4.x, a template DSS containing a precipitation `GridContainer`, and the target HEC-RAS geometry HDF. The geometry HDF is used only to read the target model projection WKT.

## Configuration

Replace the following placeholders in the notebooks or scripts:

```python
PATH_TO_UKCP18_CPM_DAILY_NETCDF_DIRECTORY
PATH_TO_19_REGION_BOUNDARY_DIRECTORY
PATH_TO_CLIPPED_REGION_NETCDF_OUTPUT_DIRECTORY
PATH_TO_CLIPPED_REGION_NETCDF_DIRECTORY
PATH_TO_DSS_OUTPUT_DIRECTORY
PATH_TO_HEC_DSSVUE_DIRECTORY
PATH_TO_TEMPLATE_DSS_WITH_A_PRECIPITATION_GRID
PATH_TO_HECRAS_GEOMETRY_HDF_WITH_PROJECTION
```

## Run

Create regional NetCDF files:

```bash
python 01_clip_ukcp18_cpm_regions.py
```

Create all DSS outputs:

```bash
python 02_ukcp18_cpm_to_hecras_dss.py
```

Filter regions or windows:

```bash
python 02_ukcp18_cpm_to_hecras_dss.py --regions forth tay --windows 01 03
```

## DSS output

- Data type: `PER-CUM`
- Units: `MM`
- Time: daily 00:00 to next-day 00:00
- Projection: target HEC-RAS geometry projection, normally OSGB36 / British National Grid
- NoData and out-of-boundary cells: `0 mm`
- Naming: `pr_<region>_20y_<window>_<start>_<end>.dss`

Splitting outputs into approximately 20-year files reduces the load when HEC-RAS imports many daily grids. If an import appears successful but disappears after reopening, check the geometry projection, DSS pathnames, saved event conditions, and precipitation linkage in the model HDF.
