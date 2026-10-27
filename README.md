# UKCP18 CPM to HEC-RAS DSS

Convert UKCP18 Convection-Permitting Model (CPM) precipitation NetCDF files
into gridded HEC-DSS rainfall inputs for HEC-RAS rain-on-grid models.

The workflow selects a region and time period, handles grid coordinates,
rainfall units and the `360_day` calendar, and writes precipitation in
millimetres per accumulation interval (`MM`, `PER-CUM`). Source files are read-only.

## Requirements

- Python 3.11 or newer and the packages in `requirements.txt`.
- UKCP18 CPM precipitation NetCDF files containing `pr` and spatial coordinates.
- A polygon boundary for regional output.
- HEC-DSSVue 3.4.x with `Jython.exe` and a precipitation grid DSS template.
  A [Forth example](examples/output_dss/pr_forth_20y_01_20191101_20391231.dss)
  can be used as the template.

## Quick start

```bash
python -m pip install -r requirements.txt
```

Edit `config.example.json`: set the input folder, output folder, boundary
and `jython_exe` paths, then choose the time windows. Use one ensemble member
per run; `input_pattern` selects the input files. Relative paths are resolved
from the configuration file's location.

Inspect the input metadata, then create the DSS files:

```bash
python 02_ukcp18_cpm_to_hecras_dss.py --config config.example.json --inspect
python 02_ukcp18_cpm_to_hecras_dss.py --config config.example.json
```

The example configuration writes DSS files directly. Set `prepare_only` to
`true` to generate intermediate rainfall data without running HEC-DSSVue.
Use a separate output folder for each dataset or ensemble member.

## Time and grid settings

- Edit the `window` list to create separate files for multiple periods, or
  use `--start 2020-02-03 --end 2020-02-07` with a configuration that omits
  `window`. Date-only end values include the whole day; datetime end values
  are exclusive and must align with the input interval.
- Regular projected square grids in metres retain their native resolution.
  For native rotated or geographic grids, specify a target projection and
  cell size, for example:

  ```bash
  python 02_ukcp18_cpm_to_hecras_dss.py --config config.example.json --target-crs EPSG:27700 --target-cell-size 2200
  ```

- Daily inputs remain daily and hourly inputs remain hourly. `360_day`
  conversion preserves monthly rainfall totals but changes event timing and
  intensities. Complete source months are required for this conversion.
- Missing rainfall is set to zero by default. Use `--missing-data error`
  for strict handling. Regional averages without a spatial grid are unsupported.

## Use the DSS rainfall in HEC-RAS

After conversion, the DSS file contains a time series of precipitation grids
that can be assigned to a HEC-RAS rain-on-grid model.

1. Open **Unsteady Flow Analysis**, select the intended **Plan**, and set its
   **Simulation Time Window** to the required start and end dates and times.
2. Open **Unsteady Flow Data** and select **Meteorological Data**. Set
   **Precipitation/Evapotranspiration** to **Enable** if needed. Under
   **Precipitation**, select **Gridded** as the **Mode** and **DSS** as the
   **Source**.
3. In **DSS Data**, use the open-file button to select the converted `.dss`
   file. In the **DSS Viewer**, double-click the first timestep of the intended
   **PRECIPITATION** grid series.
4. Set **Time Window** to **Plan Time**. If an **Apply** button is shown,
   apply the time-window setting, then confirm the selection with **OK**.
   This displays the data for the active Plan's time window. The selected
   record identifies the rainfall grid series used by the model.
5. Save the unsteady flow data and ensure that the active Plan uses this
   unsteady flow file. During computation, HEC-RAS reads the rainfall grids
   corresponding to the Plan's simulation period.
6. In **RAS Mapper**, display precipitation under **Event Conditions** to
   inspect its location, coverage and variation over time. Incremental and
   accumulated precipitation maps are available.
7. Return to **Unsteady Flow Analysis**, save the Plan and click **Compute**
   to run the model.

**Time coverage:** The DSS accumulation intervals must continuously cover the
entire Plan simulation window. The first interval must start at or before the
Plan start, and the final interval must end at or after the Plan end. Choose a
DSS output period that includes the complete simulation; **Plan Time** does
not extend the rainfall data available in the file.

See the official [Meteorological Data tutorial](https://www.hec.usace.army.mil/confluence/rasdocs/hgt/latest/tutorials/2d-unsteady-flow/meteorological-data),
[gridded precipitation guidance](https://www.hec.usace.army.mil/confluence/rasdocs/r2dum/6.6/boundary-and-initial-conditions-for-2d-flow-areas/global-boundary-conditions)
and [DSS Viewer time-window documentation](https://www.hec.usace.army.mil/confluence/rasdocs/rasum/7.0/viewing-results/viewing-data-contained-in-an-hec-dss-file).

## Files

| File | Purpose |
|---|---|
| `02_ukcp18_cpm_to_hecras_dss.py` | Convert NetCDF rainfall directly to DSS, with optional regional clipping. |
| `01_clip_ukcp18_cpm_regions.py` | Optionally export regional Gregorian NetCDF rainfall for reuse. |
| `*.ipynb` | Notebook interfaces with editable configuration cells. |
| `config.example.json` | Paths and three configurable output periods. |
| `examples/output_dss/` | Complete Forth DSS examples for all three periods. |

NetCDF output from the optional regional export contains interval rainfall
amounts in `mm`, a Gregorian time axis and CRS information. It can be used as
input to the DSS conversion without repeating calendar conversion.

See [examples](examples/README.md) for the complete output files, or run
`python 02_ukcp18_cpm_to_hecras_dss.py --help` for all options.
