# UKCP18 CPM to HEC-RAS DSS

Convert UKCP18 Convection-Permitting Model (CPM) precipitation NetCDF files
into gridded HEC-DSS rainfall inputs for HEC-RAS rain-on-grid models.

The workflow selects a time period, optionally clips to a polygon boundary,
handles grid coordinates, rainfall units and the `360_day` calendar, and writes precipitation in
millimetres per accumulation interval (`MM`, `PER-CUM`). Source files are read-only.

## Requirements

- Python 3.11 or newer and the packages in `requirements.txt`.
- UKCP18 CPM precipitation NetCDF files containing `pr` and spatial coordinates.
- An optional polygon shapefile for regional clipping.
- For DSS output: HEC-DSSVue 3.4.x with `Jython.exe` and a precipitation grid DSS template.
  A [Forth example](examples/output_dss/pr_forth_20y_01_20191101_20391231.dss)
  can be used as the template.

## Quick start

1. Open a terminal in the repository folder and install the dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

2. Save the following as `config.json` in that folder and replace the data,
   output, boundary and HEC-DSSVue paths. This example selects daily 5 km
   files from ensemble member `01`; adjust `input_pattern` for your files.

   ```json
   {
     "input_dir": "D:/climate/ukcp18-cpm",
     "input_pattern": "pr_rcp85_land-cpm_uk_5km_01_day_*.nc",
     "output_dir": "D:/hec-ras/precipitation",
     "boundary": "D:/boundaries/forth.shp",
     "name": "forth",
     "jython_exe": "C:/HEC/HEC-DSSVue/Jython.exe",
     "template_dss": "examples/output_dss/pr_forth_20y_01_20191101_20391231.dss",
     "window": [["event", "2020-02-03", "2020-02-07"]]
   }
   ```

   Use one dataset and ensemble member per run. Relative paths are resolved
   from the configuration file's location. For native rotated grids, add the
   target settings in the resolution table below.

3. Inspect the input, then run the conversion:

   ```bash
   python ukcp18_cpm_to_hecras_dss.py --config config.json --inspect
   python ukcp18_cpm_to_hecras_dss.py --config config.json
   ```

   The first command prints metadata. The second writes one `.dss` file per
   time window to `output_dir`. Use a separate output folder for each dataset
   or ensemble member.

To run the same conversion from Python:

```python
from ukcp18_cpm_to_hecras_dss import main

main(["--config", "config.json"])
```

In the supplied Notebook, set `CONFIG_PATH = Path("config.json").resolve()`
and `RUN = True`. Keep `INSPECT_ONLY = True` to inspect, then set it to
`False` to convert.

DSS is the default output. Use `--output-format netcdf` for reusable Gregorian
NetCDF rainfall without HEC-DSSVue, or `--prepare-only` for intermediate rainfall
data without writing DSS.

## Boundary selection (optional)

Set `boundary` to any polygon shapefile and `name` to the output label.
The example uses the broad Forth region from the UK river-basin regions
dataset. Replace it with a catchment or model-domain boundary as needed:

```bash
python ukcp18_cpm_to_hecras_dss.py --config config.json --boundary "D:/boundaries/my_catchment.shp" --name my_catchment
```

To browse for a shapefile, use `--choose-boundary`, or set `choose_boundary`
to `true` in the JSON configuration. This requires
a local desktop with Tkinter. Keep the `.shp`, `.shx`, `.dbf` and `.prj`
files together; boundary coordinates are transformed to the output grid CRS.

For full-grid output, set `boundary` to `null`, keep `choose_boundary`
`false` and leave `boundary_dir` unset.

## Time windows

Each `window` entry is `[ID, START, END]`. Add entries for separate output files,
or use `--windows ID` to select one. To replace the configured periods for one run:

```bash
python ukcp18_cpm_to_hecras_dss.py --config config.json --window event 2020-02-03 2020-02-07
```

Date-only ends include the whole day; datetime ends are exclusive and must align
with the input interval. Daily and hourly inputs retain their interval.
`360_day` conversion needs complete source months and preserves monthly totals
while changing event timing and intensities. Monthly and seasonal means are unsupported.

## Spatial resolution and grid handling

Grid handling follows the spatial coordinates, coordinate units and CF CRS
metadata. The target projection and cell size control any reprojection.

| Input grid | Handling | Target settings |
|---|---|---|
| Regular projected 5 km grid, such as British National Grid | Retains its CRS and 5000 m cells by default. | No override needed. |
| Regular projected 2.2 km grid | Retains its CRS and 2200 m cells by default. | No override needed. |
| Native 2.2 km rotated grid | Reprojects to a regular projected grid. Regular axes use an affine transform; stretched axes use their actual cell-centre coordinates. | Set `target_crs` and `target_cell_size`, e.g. `EPSG:27700` and `2200`. Use `5000` for a 5 km output grid. |
| Other regular projected square grids in metres | Infers cell spacing from the coordinates and retains the grid by default. | Override both settings to change the output grid. |
| Geographic or other rotated grids | Require a projected output grid with square cells in metres. | Set `target_crs` and `target_cell_size` explicitly. |
| Regional averages or curvilinear meshes without two 1-D spatial axes | Unsupported as gridded rainfall inputs. | Use a supported spatial grid. |

For example, to reproject native CPM data to 2200 m British National Grid cells:

```bash
python ukcp18_cpm_to_hecras_dss.py --config config.json --target-crs EPSG:27700 --target-cell-size 2200
```

Resampling defaults to `average` at every resolution; `--resampling nearest` or
`--resampling bilinear` selects another method. An output cell size sets the
target grid spacing; it does not increase the detail in the source rainfall.
If CRS metadata is absent, supply a verified `--source-crs`.
Missing rainfall is zero by default; use `--missing-data error` for strict handling.

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
| `ukcp18_cpm_to_hecras_dss.py` | Shared converter and command-line entry point for DSS and NetCDF output. |
| `ukcp18_cpm_to_hecras_dss.ipynb` | Notebook interface using the same converter and JSON configuration. |
| `hec_dss_grid_writer_jython.py` | DSS reader/writer run by HEC-DSSVue's separate Jython runtime. |
| `config.example.json` | Paths, output format and three configurable periods. |
| `examples/output_dss/` | Complete Forth DSS examples for all three periods. |

NetCDF exports contain interval rainfall
amounts in `mm`, a Gregorian time axis and CRS information. It can be used as
input to the DSS conversion without repeating calendar conversion.

See [examples](examples/README.md) for the complete output files, or run
`python ukcp18_cpm_to_hecras_dss.py --help` for all options.
