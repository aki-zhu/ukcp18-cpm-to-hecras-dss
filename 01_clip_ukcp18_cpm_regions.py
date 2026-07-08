# -*- coding: utf-8 -*-
"""按区域裁剪 UKCP18 CPM 5 km 降雨 NetCDF。
Clip UKCP18 CPM 5 km rainfall NetCDF files to regional boundaries.
"""

# 路径配置：将 PATH_TO_... 替换为本机数据路径。
# Path configuration: replace PATH_TO_... with local data paths.
import xarray as xr
import geopandas as gpd
import rasterio
import rasterio.features
import numpy as np
import glob
from pathlib import Path

SOURCE_DIR = Path("PATH_TO_UKCP18_CPM_DAILY_NETCDF_DIRECTORY")
files = sorted(glob.glob(str(SOURCE_DIR / "pr_rcp85_land-cpm_uk_5km_01_day_*.nc")))

REGION_BASE = Path("PATH_TO_19_REGION_BOUNDARY_DIRECTORY")
OUTPUT_DIR = Path("PATH_TO_CLIPPED_REGION_NETCDF_OUTPUT_DIRECTORY")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

REGIONS = [
    {"number": 1, "name": "Anglian", "folder": "Region_Anglian", "slug": "anglian"},
    {"number": 3, "name": "Clyde", "folder": "Region_Clyde", "slug": "clyde"},
    {"number": 4, "name": "Dee", "folder": "Region_Dee", "slug": "dee"},
    {"number": 5, "name": "Forth", "folder": "Region_Forth", "slug": "forth"},
    {"number": 6, "name": "Humber", "folder": "Region_humber", "slug": "humber"},
    {"number": 8, "name": "Northeast Scotland", "folder": "Region_NEScotland", "slug": "northeastscotland"},
    {"number": 10, "name": "North Highland", "folder": "Region_NHighland", "slug": "northhighland"},
    {"number": 11, "name": "Northwest England", "folder": "Region_NWEngland", "slug": "northwestengland"},
    {"number": 13, "name": "Northumbria", "folder": "Region_Northumbria", "slug": "northumbria"},
    {"number": 15, "name": "Severn", "folder": "Region_Severn", "slug": "severn"},
    {"number": 16, "name": "Solway", "folder": "Region_Solway", "slug": "solway"},
    {"number": 17, "name": "Southeast England", "folder": "Region_SEEngland", "slug": "southeastengland"},
    {"number": 18, "name": "Southwest England", "folder": "Region_SWEngland", "slug": "southwestengland"},
    {"number": 19, "name": "Tay", "folder": "Region_Tay", "slug": "tay"},
    {"number": 20, "name": "Thames", "folder": "Region_Thames", "slug": "thames"},
    {"number": 21, "name": "Tweed", "folder": "Region_Tweed", "slug": "tweed"},
    {"number": 23, "name": "Western Wales", "folder": "Region_WestWales", "slug": "westernwales"},
]

BOUNDARY_OVERRIDES = {
    "Region_humber": "UKCP18_uk_river_hires_humber.shp",
    "Region_Northumbria": "Northumbria/Northumbria_boundary.shp",
}

print(f"Source files: {len(files)}")
print("Output folder:", OUTPUT_DIR)

# %% ------------------------------------------------------------

# === Step 2: Build source grid ONCE ===

# open one file just to get grid
ds0 = xr.open_dataset(files[0], decode_times=False)
pr0 = ds0["pr"]

# get coordinates
x = pr0["projection_x_coordinate"].values
y = pr0["projection_y_coordinate"].values

# resolution
res_x = x[1] - x[0]
res_y = y[1] - y[0]

# correct transform; kept the same as the original Forth workflow
transform = rasterio.transform.Affine(
    res_x, 0, x.min() - res_x / 2,
    0, res_y, y.min() - res_y / 2
)

ds0.close()

POLYGON_TYPES = {"Polygon", "MultiPolygon"}

def _contains_polygon(shp_path):
    try:
        gdf = gpd.read_file(shp_path)
    except Exception:
        return False
    geom_types = set(gdf.geometry.dropna().geom_type.unique())
    return bool(geom_types & POLYGON_TYPES)

def _choose_unique_polygon(paths, label):
    unique_paths = sorted(set(paths), key=lambda p: str(p).lower())
    polygon_paths = [p for p in unique_paths if _contains_polygon(p)]
    if len(polygon_paths) == 1:
        return polygon_paths[0]
    if len(polygon_paths) > 1:
        canonical = [p for p in polygon_paths if not p.stem.lower().endswith(("_1", "_copy", "_backup"))]
        if len(canonical) == 1:
            return canonical[0]
        lines = "\n".join(str(p) for p in polygon_paths)
        raise ValueError(f"More than one polygon boundary candidate found for {label}:\n{lines}")
    return None

def find_boundary_shapefile(region_dir):
    region_dir = Path(region_dir)
    label = region_dir.name

    override = BOUNDARY_OVERRIDES.get(label)
    if override is not None:
        shp = region_dir / override
        if not shp.exists():
            raise FileNotFoundError(shp)
        if not _contains_polygon(shp):
            raise ValueError(f"Override boundary is not polygon geometry: {shp}")
        return shp

    priority_dirs = sorted(
        [p for p in region_dir.rglob("*") if p.is_dir() and ("boundary" in p.name.lower() or "shapefile" in p.name.lower())],
        key=lambda p: (len(p.relative_to(region_dir).parts), str(p).lower())
    )

    # Priority 1: the single polygon shapefile inside a boundary/shapefile folder.
    priority_shps = []
    for folder in priority_dirs:
        priority_shps.extend(folder.rglob("*.shp"))
    chosen = _choose_unique_polygon(priority_shps, label + " priority folders")
    if chosen is not None:
        return chosen

    # Priority 2: the single polygon shapefile whose filename contains boundary.
    boundary_shps = [p for p in region_dir.rglob("*.shp") if "boundary" in p.name.lower()]
    chosen = _choose_unique_polygon(boundary_shps, label + " boundary filenames")
    if chosen is not None:
        return chosen

    # Final controlled fallback: regions such as Humber store the polygon region boundary as UKCP18_uk_river_hires_*.shp.
    ukcp18_shps = [p for p in region_dir.rglob("*.shp") if "ukcp18_uk_river_hires" in p.name.lower()]
    chosen = _choose_unique_polygon(ukcp18_shps, label + " UKCP18 boundary filenames")
    if chosen is not None:
        return chosen

    raise FileNotFoundError(f"No unique polygon boundary shapefile found under {region_dir}")

def load_boundary(region):
    region_dir = REGION_BASE / region["folder"]
    shp_path = find_boundary_shapefile(region_dir)
    gdf = gpd.read_file(shp_path)
    gdf = gdf[gdf.geometry.notna()].copy()
    gdf = gdf[gdf.geometry.geom_type.isin(POLYGON_TYPES)].copy()
    if gdf.empty:
        raise ValueError(f"No polygon geometry in {shp_path}")
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:27700")
    else:
        gdf = gdf.to_crs("EPSG:27700")
    return shp_path, gdf

# %% ------------------------------------------------------------

def build_buffered_mask(gdf):
    # create mask
    mask = rasterio.features.geometry_mask(
        gdf.geometry,
        transform=transform,
        invert=True,
        out_shape=(len(y), len(x))
    )

    # buffer mask by one grid cell, same as the original Forth workflow
    buffered_mask = mask.copy()
    for dy in [-1, 0, 1]:
        for dx in [-1, 0, 1]:
            shifted = np.roll(np.roll(mask, dy, axis=0), dx, axis=1)
            buffered_mask |= shifted

    # fix edges
    buffered_mask[0, :] = False
    buffered_mask[-1, :] = False
    buffered_mask[:, 0] = False
    buffered_mask[:, -1] = False

    return buffered_mask

# %% ------------------------------------------------------------

def get_crop_window(buffered_mask):
    rows, cols = np.where(buffered_mask)

    pad = 0  # already buffered above, no need extra

    if len(rows) == 0 or len(cols) == 0:
        raise ValueError("Mask does not overlap the UKCP18 grid")

    row_min = rows.min()
    row_max = rows.max()
    col_min = cols.min()
    col_max = cols.max()

    return row_min, row_max, col_min, col_max

# %% ------------------------------------------------------------

def clip_region(region):
    output_nc = OUTPUT_DIR / f"pr_{region['slug']}_2010_2080.nc"
    if output_nc.exists():
        print(f"\nRegion {region['number']} {region['name']} already exists, skipping:", output_nc)
        return {
            "number": region["number"],
            "name": region["name"],
            "slug": region["slug"],
            "output_nc": str(output_nc),
            "skipped": True,
        }

    shp_path, gdf = load_boundary(region)
    print(f"\nRegion {region['number']} {region['name']}")
    print("Boundary:", shp_path)

    buffered_mask = build_buffered_mask(gdf)
    row_min, row_max, col_min, col_max = get_crop_window(buffered_mask)
    print(f"Crop rows {row_min}:{row_max}, cols {col_min}:{col_max}")

    clipped_list = []

    for f in files:
        print(f"Processing: {f}")

        ds = xr.open_dataset(f, decode_times=False)
        pr = ds["pr"]

        # apply mask
        pr_clip = pr.where(buffered_mask)

        # crop
        pr_clip_small = pr_clip.isel(
            projection_y_coordinate=slice(row_min, row_max + 1),
            projection_x_coordinate=slice(col_min, col_max + 1)
        )

        clipped_list.append(pr_clip_small)

    pr_all = xr.concat(clipped_list, dim="time")

    encoding = {
        "pr": {
            "zlib": True,
            "complevel": 4
        }
    }

    pr_all.to_netcdf(output_nc, encoding=encoding)
    print("Written:", output_nc)
    return {
        "number": region["number"],
        "name": region["name"],
        "slug": region["slug"],
        "boundary": str(shp_path),
        "output_nc": str(output_nc),
        "rows": [int(row_min), int(row_max)],
        "cols": [int(col_min), int(col_max)],
        "time_steps": int(pr_all.sizes["time"]),
        "ny": int(pr_all.sizes["projection_y_coordinate"]),
        "nx": int(pr_all.sizes["projection_x_coordinate"]),
    }

# %% ------------------------------------------------------------

summaries = []

for region in REGIONS:
    summaries.append(clip_region(region))

summaries

# %% ------------------------------------------------------------

import matplotlib.pyplot as plt

# Final check plot for one output region
check_region = summaries[0]["slug"]
check_nc = OUTPUT_DIR / f"pr_{check_region}_2010_2080.nc"
check_ds = xr.open_dataset(check_nc, decode_times=False)
pr_map = check_ds["pr"].isel(time=5).squeeze()

fig, ax = plt.subplots(figsize=(6, 6))
mesh = ax.pcolormesh(
    check_ds["projection_x_coordinate"],
    check_ds["projection_y_coordinate"],
    pr_map,
    shading="auto"
)
plt.colorbar(mesh, ax=ax)
plt.title(f"Final check: {check_region}")
plt.show()
