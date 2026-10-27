# Forth DSS rainfall examples

Complete daily precipitation grids for the broad Forth region in the UK
river-basin regions dataset, using 5000 m cells on a 35 by 17 grid in
OSGB36 / British National Grid. Rainfall amounts are
stored as `MM`, `PER-CUM`, from 00:00 to 00:00 on the following day.

| DSS file | Inclusive daily dates | Size |
|---|---|---:|
| [Period 01](output_dss/pr_forth_20y_01_20191101_20391231.dss) | 2019-11-01 to 2039-12-31 | 22.3 MB |
| [Period 02](output_dss/pr_forth_20y_02_20391101_20591231.dss) | 2039-11-01 to 2059-12-31 | 22.3 MB |
| [Period 03](output_dss/pr_forth_20y_03_20591101_20801130.dss) | 2059-11-01 to 2080-11-30 | 23.3 MB |

The periods intentionally overlap. The `20y` filename label denotes the
period grouping; the exact coverage is given in the table.

These files illustrate the DSS output structure. For another catchment,
generate rainfall grids covering that catchment with its own input data
and boundary.
