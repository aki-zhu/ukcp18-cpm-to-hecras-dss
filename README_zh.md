# UKCP18 CPM to HEC-RAS DSS：使用说明

支持任意起止日期、多个时间段及从 metadata 读取网格分辨率。
脚本与 Notebook 共用转换模块，避免两套参数或逻辑不同步。

## 最常用的修改

1. 复制 `config.example.json`，修改 `input_dir`、`output_dir`。
2. 在 `window` 中修改或增删 `["编号", "开始", "结束"]`。
3. 日日期包含结束日；带时分的结束时间不包含该时刻，例如 03:00 到 09:00 为六条小时记录。
4. 原生旋转经纬度网格设置 `target_crs="EPSG:27700"`、`target_cell_size=2200`；规则投影网格可保留原分辨率。
5. `prepare_only=true` 只准备并检查数据；正式写 DSS 时设为 false，并填写 `jython_exe` 与 `template_dss`。

```bash
python 02_ukcp18_cpm_to_hecras_dss.py --config config.example.json
python 02_ukcp18_cpm_to_hecras_dss.py --input-dir "D:/climate/cpm" --inspect
```

可用 `--start 2020-02-03 --end 2020-02-07` 选择一个时间段；
也可重复 `--window 编号 开始日期 结束日期` 输出多个文件。
原本的默认时间段仍保留；自定义时间段会替换默认值。

## 网格与时间处理

规则平面网格自动识别 1、2.2、5、12、25 km 等分辨率，不固定文件名中的分辨率。
原生 CPM 的旋转经纬度、0–360° 经度和外围拉伸网格会按实际坐标处理。
目标 DSS 使用以米为单位的投影正方形网格；指定目标网格时会实际重投影数值。
小时输入保留小时累计时段，降雨率按原始时间步长换算为 mm。

CPM 的 `360_day` 日历按完整月转换，逐格保持月降雨总量；小时数据分别保持每个小时槽的月总量。选择月中日期也需要该月完整源数据。可选 NC 裁剪阶段现在输出含 CRS、时间边界和 Gregorian 日历的区间累计量，后续 DSS 阶段不会重复转换。

不同行政/河流区域平均值产品缺少二维网格时，不能当作空间降雨场。
缺时间、重复时间、混用成员/网格/单位、负降雨或缺少 CRS 会明确报错。
旧裁剪文件若丢失 CRS，可明确指定 `--source-crs EPSG:27700`。
缺测值默认按既有流程置零；严格检查用 `--missing-data error`。

## 原有完整输出示例

`examples/output_dss/` 仅存放原有 Forth 区域的完整 DSS 文件。
CPM 包含三个原有时段，HadUK 包含两个原有时段；未经截短或重算。
详细日期、文件大小、记录数和文件校验值见 `examples/README.md`。
这些示例属于之前的计算结果，新功能验证另见 `VALIDATION.md`。
