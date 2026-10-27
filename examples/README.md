# Complete Forth DSS output examples

These are full original files, copied without truncation or recalculation from
the user-supplied local output directories. They illustrate the established
5 km MM / PER-CUM workflow. They are not outputs of the new implementation.
Only DSS files are stored in `output_dss/`; supporting notes stay here.

The windows intentionally overlap. Their precise inclusive dates are shown
below; "20y" is the original naming convention (the last CPM window is longer).
All five supplied Forth examples use 35 x 17 cells with 5000 m cell size.
Daily catalog counts: CPM 01=7366, 02=7366, 03=7701; HadUK 01=7366, 02=7366.
The first and last pathnames and headers were read with HEC-DSSVue 3.4.15.

| File | Inclusive daily dates | Bytes |
|---|---|---:|
| [pr_forth_20y_01_20191101_20391231.dss](output_dss/pr_forth_20y_01_20191101_20391231.dss) | 20191101–20391231 | 22,298,456 |
| [pr_forth_20y_02_20391101_20591231.dss](output_dss/pr_forth_20y_02_20391101_20591231.dss) | 20391101–20591231 | 22,304,512 |
| [pr_forth_20y_03_20591101_20801130.dss](output_dss/pr_forth_20y_03_20591101_20801130.dss) | 20591101–20801130 | 23,332,584 |

SHA-256 checksums are recorded in `checksums.json` and compared with the
local originals. Raw climate inputs, model files, logs and payloads are excluded.
One of these files can be used as a precipitation GridContainer template.
