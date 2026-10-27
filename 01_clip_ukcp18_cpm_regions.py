"""Prepare regional CF NetCDF with configurable dates and output grids."""
from dss_pipeline import main

if __name__ == "__main__":
    main("cpm", clip=True)
