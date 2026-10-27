"""Scientific invariants and end-to-end payload/CF round-trip regressions."""
import csv
from datetime import timedelta
import json
from pathlib import Path
import tempfile
import unittest

import cftime
import numpy as np
import pandas as pd
from pyproj import CRS
import xarray as xr

from dss_pipeline import (Catalog, Grid, SpatialPlan, Window, conservative_month,
                          dilate_mask, expected_dates, iter_values, main,
                          parse_args, payload_to_netcdf, prepare, unit_scale)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.source.mkdir()
        self.output = self.root / "output"

    def args(self, kind="cpm", extra=()):
        return parse_args(kind, ["--input-dir", str(self.source), "--output-dir", str(self.output),
                                 "--prepare-only", *extra])

    def make_source(self, name="pr_test.nc", calendar="standard", days=31, hours=24,
                    units="mm", size=5000, reverse=False, rotated=False,
                    variable="pr", start=None, values=None, nominal_hour=0):
        step = timedelta(hours=hours)
        if calendar == "360_day":
            first = start or cftime.Datetime360Day(2020, 1, 1)
        else:
            first = start or pd.Timestamp(2020, 1, 1).to_pydatetime()
        count = int(days*24/hours)
        starts = [first+i*step+timedelta(hours=nominal_hour) for i in range(count)]
        centres = [t+step/2 for t in starts]
        bounds = np.array([[t,t+step] for t in starts])
        x = np.arange(3)*size+size/2
        y = np.arange(4)*size+size/2
        if rotated:
            x, y = np.arange(3)*.02-1, np.arange(4)*.02-1
            crs_attrs = dict(grid_mapping_name="rotated_latitude_longitude", earth_radius=6371229.,
                             grid_north_pole_latitude=37.5, grid_north_pole_longitude=177.5,
                             north_pole_grid_longitude=0.)
            xn, yn = "grid_longitude", "grid_latitude"
        else:
            crs_attrs = CRS.from_epsg(27700).to_cf()
            xn, yn = "projection_x_coordinate", "projection_y_coordinate"
        if values is None:
            values = np.broadcast_to(np.arange(12).reshape(4,3)+1., (count,4,3)).copy()
        if reverse:
            x,y,values = x[::-1], y[::-1], values[:,::-1,::-1]
        ds = xr.Dataset({variable: (("time",yn,xn),values,{"units":units,"grid_mapping":"crs"}),
                         "crs": ((),0,crs_attrs), "time_bnds": (("time","bnds"),bounds)},
                        coords={"time":centres,xn:x,yn:y})
        ds.time.attrs["bounds"] = "time_bnds"
        for n in (xn,yn):
            ds[n].attrs.update(standard_name=n, units="degrees" if rotated else "m")
        encoding = {k:dict(units="minutes since 1900-01-01",calendar=calendar) for k in ("time","time_bnds")}
        path = self.source/name
        ds.to_netcdf(path, encoding=encoding)
        return path

    def catalog(self, paths, args=None, kind="cpm"):
        cat = Catalog(paths, args or self.args(), kind)
        self.addCleanup(cat.close)
        return cat

    def test_date_only_inclusive_and_datetime_exclusive(self):
        a = Window.parse("x","2020-02-02","2020-02-03")
        self.assertEqual(a.stop, pd.Timestamp("2020-02-04"))
        b = Window.parse("y","2020-02-02T03:00","2020-02-02T06:00")
        self.assertEqual((b.stop-b.start).total_seconds(),10800)
        with self.assertRaises(ValueError):
            Window.parse("../bad","2020-01-01","2020-01-02")

    def test_monthly_amounts_leap_and_dry_cells(self):
        rng=np.random.default_rng(9)
        values=rng.uniform(0,50,(30,4,3)).astype(np.float32)
        values[:,0,0]=0
        for days in (28,29,30,31):
            mapped=conservative_month(values,days)
            np.testing.assert_allclose(mapped.sum(axis=0,dtype=float),values.sum(axis=0,dtype=float),rtol=2e-7)
            self.assertTrue(np.isfinite(mapped).all())
            self.assertTrue((mapped[:,0,0]==0).all())

    def test_rate_units_daily_hourly(self):
        self.assertEqual(unit_scale("kg m-2 s-1",3600),3600)
        self.assertEqual(unit_scale("mm/day",3600),1/24)
        self.assertEqual(unit_scale("mm h^-1",86400),24)
        self.assertEqual(unit_scale("m",3600),1000)
        with self.assertRaises(ValueError):
            unit_scale("K",3600)

    def test_multiple_resolutions_descending_and_fractional_origin(self):
        for size in (1000,2200,5000,12000,25000):
            path=self.make_source(name=f"pr_{size}.nc",size=size,reverse=True,days=2)
            cat=self.catalog([path])
            plan=SpatialPlan(cat.grid,cat.args)
            values=cat.read([sorted(cat.records)[0]],plan)
            np.testing.assert_equal(values[0],np.arange(12).reshape(4,3)+1)
            self.assertEqual(plan.target.dss_metadata()["cell_size"],size)
        grid=Grid(np.array([1.25,2201.25]),np.array([2.5,2202.5]),2200,2200,CRS.from_epsg(27700))
        meta=grid.dss_metadata()
        self.assertEqual(meta["cell_zero_x"],-1098.75)

    def test_non_month_aligned_partial_range(self):
        path=self.make_source(calendar="360_day",days=30,units="mm/day")
        args=self.args()
        cat=self.catalog([path],args)
        plan=SpatialPlan(cat.grid,args)
        win=Window.parse("partial","2020-01-09","2020-01-13")
        chunks=list(iter_values(cat,win,plan,self.root))
        dates=chunks[0][0]
        self.assertEqual(len(dates),5)
        np.testing.assert_allclose(chunks[0][1][0],(np.arange(12).reshape(4,3)+1)*30/31,rtol=1e-6)

    def test_hourly_360_preserves_month_per_hour_and_times(self):
        raw=np.ones((720,4,3),dtype=float)
        raw[::24]*=6
        path=self.make_source(calendar="360_day",days=30,hours=1,units="mm h-1",values=raw)
        args=self.args(extra=["--chunk-records","17"])
        cat=self.catalog([path],args)
        plan=SpatialPlan(cat.grid,args)
        win=Window.parse("jan","2020-01-01","2020-01-31")
        chunks=list(iter_values(cat,win,plan,self.root))
        out=np.concatenate([v for d,v in chunks])
        self.assertEqual(len(out),744)
        self.assertTrue(all(len(v)<=17 for _,v in chunks))
        np.testing.assert_allclose(out.sum(axis=0,dtype=float),raw.sum(axis=0),rtol=2e-7)
        np.testing.assert_allclose(out[::24].sum(axis=0,dtype=float),raw[::24].sum(axis=0),rtol=2e-7)
        self.assertEqual(chunks[-1][0][-1],pd.Timestamp("2020-01-31T23:00"))

    def test_missing_360_edge_month_rejected_and_temp_released(self):
        path=self.make_source(calendar="360_day",days=29)
        cat=self.catalog([path])
        with self.assertRaisesRegex(ValueError,"Missing source"):
            list(iter_values(cat,Window.parse("x","2020-01-02","2020-01-03"),SpatialPlan(cat.grid,cat.args),self.root))
        self.assertFalse(list(self.root.glob("calendar_*")))

    def test_gregorian_leap_day_hourly_subset_and_missing(self):
        path=self.make_source(days=3,hours=1,start=pd.Timestamp("2020-02-28").to_pydatetime())
        cat=self.catalog([path])
        win=Window.parse("leap","2020-02-29T02:00","2020-02-29T06:00")
        chunks=list(iter_values(cat,win,SpatialPlan(cat.grid,cat.args),self.root))
        self.assertEqual(len(chunks[0][0]),4)
        with self.assertRaisesRegex(ValueError,"not covered"):
            list(iter_values(cat,Window.parse("bad","2020-03-02","2020-03-03"),SpatialPlan(cat.grid,cat.args),self.root))

    def test_haduk_original_filename_days_and_no_calendar_remapping(self):
        path=self.make_source(name="rainfall_hadukgrid_uk_5km_day_20200228-20200301.nc",
                              variable="rainfall",start=pd.Timestamp("2020-02-28").to_pydatetime(),
                              days=3,nominal_hour=9)
        args=self.args("haduk")
        cat=self.catalog([path],args,"haduk")
        self.assertIn((2020,2,29,0,0),cat.records)
        plan=SpatialPlan(cat.grid,args)
        chunks=list(iter_values(cat,Window.parse("x","2020-02-29","2020-03-01"),plan,self.root))
        np.testing.assert_equal(chunks[0][1][0],np.arange(12).reshape(4,3)+1)
        self.assertEqual(chunks[0][0][0].hour,0)

    def test_haduk_filename_metadata_mismatch_rejected(self):
        path=self.make_source(name="rainfall_uk_day_20200101-20200103.nc",days=2,variable="rainfall")
        with self.assertRaisesRegex(ValueError,"filename"):
            Catalog([path],self.args("haduk"),"haduk")

    def test_duplicate_and_changed_grid_rejected(self):
        a=self.make_source(name="pr_a.nc",days=2)
        b=self.make_source(name="pr_b.nc",days=2)
        with self.assertRaisesRegex(ValueError,"Overlapping"):
            Catalog([a,b],self.args(),"cpm")
        b=self.make_source(name="pr_b.nc",days=2,size=2200)
        with self.assertRaisesRegex(ValueError,"Grid"):
            Catalog([a,b],self.args(),"cpm")

    def test_nodata_zero_and_negative_rejected(self):
        raw=np.ones((2,4,3));raw[0,0,0]=np.nan
        path=self.make_source(days=2,values=raw)
        cat=self.catalog([path])
        plan=SpatialPlan(cat.grid,cat.args)
        self.assertEqual(cat.read([sorted(cat.records)[0]],plan)[0,0,0],0)
        strict=self.catalog([path],self.args(extra=["--missing-data","error"]))
        with self.assertRaisesRegex(ValueError,"Missing/nonfinite"):
            strict.read([sorted(strict.records)[0]],plan)
        cat.close();strict.close()
        raw[0,0,0]=-1
        path=self.make_source(days=2,values=raw)
        cat=self.catalog([path])
        with self.assertRaisesRegex(ValueError,"Negative"):
            cat.read([sorted(cat.records)[0]],SpatialPlan(cat.grid,cat.args))

    def test_buffer_keeps_edges_without_wraparound(self):
        mask=np.zeros((5,5),dtype=bool);mask[0,0]=True
        expanded=dilate_mask(mask,1)
        self.assertTrue(expanded[0,0])
        self.assertEqual(expanded.sum(),4)
        self.assertFalse(expanded[-1,-1])

    def test_rotated_grid_requires_target_and_reprojects(self):
        path=self.make_source(rotated=True,days=2,values=np.full((2,4,3),7.))
        cat=self.catalog([path])
        with self.assertRaisesRegex(ValueError,"projected CRS"):
            SpatialPlan(cat.grid,cat.args)
        args=self.args(extra=["--target-crs","EPSG:27700","--target-cell-size","2200"])
        plan=SpatialPlan(cat.grid,args)
        output=cat.read([sorted(cat.records)[0]],plan)
        self.assertTrue((output>0).any())
        np.testing.assert_allclose(output[output>0],7,rtol=1e-6)

    def test_two_windows_payload_bytes_stats_and_round_trip(self):
        self.make_source(days=5,hours=1)
        argv=["--input-dir",str(self.source),"--output-dir",str(self.output),"--prepare-only",
              "--window","one","2020-01-02T01:00","2020-01-02T04:00",
              "--window","two","2020-01-04","2020-01-04"]
        jobs=main("cpm",argv)
        self.assertEqual(len(jobs),2)
        meta=json.loads(jobs[0].read_text())
        data=np.fromfile(meta["values_bin"],dtype="<f4").reshape(3,4,3)
        np.testing.assert_equal(data[0],np.arange(12).reshape(4,3)+1)
        with open(meta["dates_csv"],encoding="ascii") as handle:
            rows=list(csv.DictReader(handle))
        self.assertEqual(rows[-1]["end"],"2020-01-02T04:00:00")
        self.assertEqual(float(rows[0]["mean"]),6.5)
        args=self.args()
        win=Window.parse("one","2020-01-02T01:00","2020-01-02T04:00")
        nc=payload_to_netcdf(jobs[0],args,"test",win)
        cat=self.catalog([Path(nc)])
        self.assertEqual(cat.calendar,"proleptic_gregorian")
        again=cat.read(sorted(cat.records),SpatialPlan(cat.grid,args))
        np.testing.assert_equal(again,data)

    def test_metadata_inspection_writes_nothing(self):
        self.make_source(days=2)
        main("cpm",["--input-dir",str(self.source),"--inspect"])
        self.assertFalse(self.output.exists())

    def test_json_config_paths_relative_and_cli_override(self):
        config=self.root/"settings.json"
        config.write_text(json.dumps(dict(input_dir="source",output_dir="output",chunk_records=2)))
        args=parse_args("cpm",["--config",str(config),"--chunk-records","5"])
        self.assertTrue(Path(args.input_dir).samefile(self.source))
        self.assertEqual(args.chunk_records,5)

    def test_stretched_grid_geolocation_has_correct_pixel_centres(self):
        source=Grid(np.array([500.,1500.,3500.]),np.array([500.,1500.,2500.,3500.]),
                    1500.,1000.,CRS.from_epsg(27700))
        args=self.args(extra=["--target-cell-size","1000","--resampling","nearest"])
        plan=SpatialPlan(source,args)
        self.addCleanup(plan.close)
        raw=np.broadcast_to(np.array([1.,2.,3.]),(1,4,3)).copy()
        result=plan.apply(raw)
        self.assertEqual(result[0,1,0],1)
        self.assertEqual(result[0,1,1],2)
        self.assertEqual(result[0,1,3],3)
        self.assertTrue(plan.target.regular)

    def test_no_bounds_requires_explicit_timestamp_convention(self):
        path=self.make_source(days=2,hours=1)
        with xr.open_dataset(path) as ds:
            copy=ds.drop_vars("time_bnds").load()
        copy.time.attrs.pop("bounds")
        copy.to_netcdf(path,mode="w")
        with self.assertRaisesRegex(ValueError,"No CF time bounds"):
            Catalog([path],self.args(),"cpm")
        cat=self.catalog([path],self.args(extra=["--time-label","centre"]))
        self.assertIn((2020,1,1,0,0),cat.records)

    def test_multi_member_selection_is_required(self):
        path=self.make_source(days=2)
        with xr.open_dataset(path) as ds:
            copy=ds.load()
        copy["pr"]=copy.pr.expand_dims(ensemble_member=[1,4])
        copy.to_netcdf(path,mode="w")
        with self.assertRaisesRegex(ValueError,"choose one"):
            Catalog([path],self.args(),"cpm")
        cat=self.catalog([path],self.args(extra=["--select","ensemble_member=4"]))
        self.assertEqual(cat.signature[2]["ensemble_member"],"4")

    def test_grouped_regional_files_keep_regions_separate(self):
        self.make_source(name="pr_forth_1990_2000.nc",days=3)
        self.make_source(name="pr_tay_1990_2000.nc",days=3,values=np.full((3,4,3),8.))
        jobs=main("cpm",["--input-dir",str(self.source),"--output-dir",str(self.output),
                          "--start","2020-01-01","--end","2020-01-02","--prepare-only"])
        self.assertEqual(len(jobs),2)
        self.assertIn("pr_forth_",str(jobs[0]))
        self.assertIn("pr_tay_",str(jobs[1]))

    def test_boundary_clipping_preserves_grid_direction(self):
        from shapely.geometry import box
        import geopandas as gpd
        boundary=self.root/"boundary.geojson"
        gpd.GeoDataFrame(geometry=[box(0,0,5000,5000)],crs=27700).to_file(boundary,driver="GeoJSON")
        path=self.make_source(days=2)
        args=self.args(extra=["--buffer-cells","0"])
        cat=self.catalog([path],args)
        plan=SpatialPlan(cat.grid,args,boundary)
        self.addCleanup(plan.close)
        values=cat.read([sorted(cat.records)[0]],plan)
        self.assertEqual(values[0,0,0],1)
        self.assertEqual(plan.target.bounds[0],0)
        self.assertEqual(plan.target.bounds[1],0)

    def test_clip_default_outputs_one_full_coverage_cf_file(self):
        from shapely.geometry import box
        import geopandas as gpd
        boundary=self.root/"boundary.geojson"
        gpd.GeoDataFrame(geometry=[box(0,0,5000,5000)],crs=27700).to_file(boundary,driver="GeoJSON")
        self.make_source(calendar="360_day",days=30)
        jobs=main("cpm",["--input-dir",str(self.source),"--output-dir",str(self.output),
                          "--boundary",str(boundary)],clip=True)
        self.assertEqual(len(jobs),1)
        paths=list(self.output.glob("*.nc"))
        self.assertEqual(len(paths),1)
        with xr.open_dataset(paths[0]) as ds:
            self.assertEqual(ds.sizes["time"],31)
            self.assertIn("crs",ds)
            self.assertIn("time_bnds",ds)


if __name__ == "__main__":
    unittest.main()
