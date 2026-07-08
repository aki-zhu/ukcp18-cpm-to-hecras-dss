# UKCP18 CPM to HEC-RAS DSS / UKCP18 CPM 转 HEC-RAS DSS

## 项目简介 / Overview

本仓库将 UKCP18 CPM 5 km 逐日格点降雨 NetCDF 转换为 HEC-RAS 可读取的 HEC-DSS specified grid precipitation 数据。流程包括 17 个区域的 NetCDF 裁剪，以及保持月累计量的 `360_day` 至 Gregorian 日历转换。

This repository converts UKCP18 CPM 5 km daily gridded rainfall NetCDF files into HEC-DSS specified-grid precipitation inputs for HEC-RAS. It includes regional NetCDF clipping for 17 UK regions and a monthly-total-preserving `360_day`-to-Gregorian calendar conversion.

原始 UKCP18 NetCDF 只读，不会被修改。数据、boundary shapefile、模板 DSS 和 HEC-RAS geometry HDF 不包含在仓库中。

The original UKCP18 NetCDF files are read-only and are never modified. Data, boundary shapefiles, the template DSS, and the HEC-RAS geometry HDF are not included.

## 两阶段流程 / Two-stage pipeline

### 1. 区域 NetCDF / Regional NetCDF

1. 从首个 UKCP18 CPM 文件读取 5 km British National Grid。
2. 为每个区域选择唯一 polygon boundary，并转换到 EPSG:27700。
3. 将 boundary 栅格化并向外缓冲一个格点。
4. 对每个输入文件应用掩膜和最小矩形裁剪。
5. 沿 time 维拼接并输出 `pr_<region>_2010_2080.nc`。

1. Read the 5 km British National Grid from the first UKCP18 CPM file.
2. Select one polygon boundary per region and transform it to EPSG:27700.
3. Rasterize and buffer each boundary by one grid cell.
4. Apply the mask and minimum rectangular crop to every source file.
5. Concatenate along time and write `pr_<region>_2010_2080.nc`.

### 2. NetCDF 转 DSS / NetCDF to DSS

1. 检查 `pr`、空间坐标、格网方向和 `360_day` 时间轴。
2. 将每个 30 天模式月插值到真实 Gregorian 月长。
3. 对每个格点缩放插值结果，使该月累计降雨量与源月完全一致。
4. 将 NoData 和区域外格点写为 `0 mm`。
5. 生成逐日 little-endian float32 payload。
6. 使用 HEC-DSSVue Jython 写入 `MM`、`PER-CUM` specified grids。
7. 记录并核对 DSS catalog、格网尺寸、投影、单位和时间。

1. Check `pr`, spatial coordinates, grid direction, and the `360_day` time axis.
2. Interpolate each 30-day model month to the actual Gregorian month length.
3. Scale each cell so the interpolated monthly accumulation exactly matches the source month.
4. Store NoData and out-of-region cells as `0 mm`.
5. Build daily little-endian float32 payloads.
6. Use HEC-DSSVue Jython to write `MM`, `PER-CUM` specified grids.
7. Log and verify the DSS catalog, grid dimensions, projection, units, and times.

## 文件 / Files

| 文件 / File | 说明 / Description |
|---|---|
| `01_clip_ukcp18_cpm_regions.ipynb` | 17 区域裁剪 notebook / 17-region clipping notebook |
| `01_clip_ukcp18_cpm_regions.py` | 区域裁剪 Python 脚本 / Regional clipping Python script |
| `02_ukcp18_cpm_to_hecras_dss.ipynb` | 日历转换与 DSS 写入 notebook / Calendar conversion and DSS-writing notebook |
| `02_ukcp18_cpm_to_hecras_dss.py` | 支持区域和时间窗筛选的命令行脚本 / CLI script with region and window filters |
| `hec_dss_grid_writer_jython.py` | HEC-DSSVue Jython 格网写入器 / HEC-DSSVue Jython grid writer |

## 输入要求 / Input assumptions

原始 CPM 文件 / Source CPM files:

- 变量 / Variable: `pr`
- 单位 / Units: `mm/day`
- 日历 / Calendar: `360_day`
- 坐标 / Coordinates: `projection_x_coordinate`, `projection_y_coordinate`
- 网格 / Grid: 5000 m British National Grid
- 示例模式 / Example pattern: `pr_rcp85_land-cpm_uk_5km_01_day_*.nc`

区域文件名 / Regional filename pattern:

```text
pr_<region>_2010_2080.nc
```

17 regions:

```text
Anglian, Clyde, Dee, Forth, Humber, Northeast Scotland,
North Highland, Northwest England, Northumbria, Severn, Solway,
Southeast England, Southwest England, Tay, Thames, Tweed, Western Wales
```

## 日历与单位 / Calendar and units

源 `pr` 已经是 `mm/day`。转换不会除以 24，也不会除以 1000。DSS 中每日记录的单位是 `MM`，数据类型是 `PER-CUM`。

Source `pr` is already in `mm/day`. The workflow does not divide values by 24 or 1000. Each daily DSS record uses units `MM` and data type `PER-CUM`.

`360_day` 转换在每个月内部完成：先从 30 天插值到 28、29、30 或 31 天，再按格点校正，使月累计量守恒。该方法保留月尺度水量，但会改变逐日序列的形状，因此应在方法说明中明确记录。

The `360_day` conversion is performed month by month: 30 source days are interpolated to 28, 29, 30, or 31 Gregorian days, then corrected cell-by-cell to conserve the monthly accumulation. This preserves monthly water volume while changing the shape of the daily sequence, so the method should be documented in model metadata.

默认输出时间窗 / Default output windows:

- `2019-11-01` to `2039-12-31`
- `2039-11-01` to `2059-12-31`
- `2059-11-01` to `2080-11-30`

时间窗之间有意重叠，但每个 DSS 文件内部不包含重复日期。

The windows intentionally overlap, but no individual DSS file contains duplicate dates.

## 环境 / Environment

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

DSS 写入还需要 HEC-DSSVue 3.4.x 的 `Jython.exe`、一个含 precipitation `GridContainer` 的模板 DSS，以及目标 HEC-RAS geometry HDF。geometry HDF 仅用于读取目标模型投影 WKT。

DSS writing also requires `Jython.exe` from HEC-DSSVue 3.4.x, a template DSS containing a precipitation `GridContainer`, and the target HEC-RAS geometry HDF. The geometry HDF is used only to read the target model projection WKT.

## 配置 / Configuration

在 notebook 或脚本中替换以下占位符：

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

## 运行 / Run

先生成区域 NetCDF：

First create regional NetCDF files:

```bash
python 01_clip_ukcp18_cpm_regions.py
```

再生成全部 DSS：

Then create all DSS outputs:

```bash
python 02_ukcp18_cpm_to_hecras_dss.py
```

筛选区域或时间窗 / Filter regions or windows:

```bash
python 02_ukcp18_cpm_to_hecras_dss.py --regions forth tay --windows 01 03
```

## DSS 输出 / DSS output

- 数据类型 / Data type: `PER-CUM`
- 单位 / Units: `MM`
- 时间 / Time: daily 00:00 to next-day 00:00
- 投影 / Projection: target HEC-RAS geometry projection, normally OSGB36 / British National Grid
- NoData and out-of-boundary cells: `0 mm`
- 命名 / Naming: `pr_<region>_20y_<window>_<start>_<end>.dss`

按约 20 年分文件可降低 HEC-RAS 导入大量逐日 bands 时的负担。若 HEC-RAS 显示导入成功但重新打开后数据消失，应检查 geometry 投影、DSS pathname、event conditions 保存状态和模型 HDF 中的 precipitation linkage。

Splitting outputs into approximately 20-year files reduces the load when HEC-RAS imports many daily grids. If an import appears successful but disappears after reopening, check the geometry projection, DSS pathnames, saved event conditions, and precipitation linkage in the model HDF.
