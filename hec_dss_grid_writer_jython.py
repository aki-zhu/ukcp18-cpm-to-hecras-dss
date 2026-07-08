# -*- coding: utf-8 -*-
"""HEC-DSSVue Jython grid writer.
"""

from hec.heclib.dss import HecDss
from java.nio import ByteBuffer, ByteOrder
from java.io import File
from jarray import zeros
import csv, json, sys
from datetime import datetime, timedelta

def dss_time(d):
    return d.strftime('%d%b%Y').upper() + ':0000'

def read_float_array(fh, n):
    b = fh.read(n * 4)
    if len(b) != n * 4:
        raise Exception('Expected %d bytes, got %d' % (n * 4, len(b)))
    bb = ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN)
    arr = zeros(n, 'f')
    for i in xrange(n):
        arr[i] = bb.getFloat()
    return arr

def find_template_grid(template_dss):
    src = HecDss.open(template_dss)
    try:
        for p in src.getPathnameList():
            s = str(p)
            if '/PRECIPITATION/' in s.upper():
                try:
                    obj = src.get(s)
                    if obj.__class__.__name__.endswith('GridContainer'):
                        print 'Template grid:', s
                        return obj
                except Exception, e:
                    pass
        raise Exception('No readable GridContainer template found in ' + template_dss)
    finally:
        src.close()

def make_path(meta, start_s, end_s):
    return '/%s/%s/%s/%s/%s/%s/' % (
        meta.get('a_part', ''),
        meta.get('b_part', ''),
        meta.get('c_part', 'PRECIPITATION'),
        start_s,
        end_s,
        meta.get('f_part', '')
    )

def main(meta_path):
    with open(meta_path, 'r') as f:
        meta = json.load(f)
    rows = list(csv.DictReader(open(meta['dates_csv'], 'r')))
    n = int(meta['nx']) * int(meta['ny'])

    out_file = File(meta['out_dss'])
    if out_file.exists():
        out_file.delete()

    obj = find_template_grid(meta['template_dss'])
    gi = obj.getGridInfo()
    gi.setGridType(430)
    gi.setCellInfo(int(meta['lower_left_cell_x']), int(meta['lower_left_cell_y']), int(meta['nx']), int(meta['ny']), float(meta['cell_size']))
    gi.setSpatialReference(meta['srs_name'], meta['srs_wkt'], float(meta['cell_zero_x']), float(meta['cell_zero_y']))
    gi.setCoordOfGridCellZero(float(meta['cell_zero_x']), float(meta['cell_zero_y']))
    gi.setDataType(1)  # PER-CUM
    gi.setDataUnits(meta.get('data_units', 'MM'))

    dst = HecDss.open(meta['out_dss'])
    written = 0
    try:
        fh = open(meta['values_bin'], 'rb')
        try:
            for row in rows:
                start = datetime.strptime(row['date'], '%Y-%m-%d')
                end = start + timedelta(days=1)
                start_s = dss_time(start)
                end_s = dss_time(end)

                gd = obj.getGridData()
                gd.setData(read_float_array(fh, n))
                gd.initGridInfo(gi)
                try:
                    gd.updateStatistics()
                except Exception, e:
                    pass

                gi.setGridTimes(start_s, end_s)
                gi.setDataInfo(float(row['min']), float(row['max']), float(row['mean']))

                full = make_path(meta, start_s, end_s)
                obj.fullName = full
                obj.setFullName(full)
                dst.put(obj)
                written += 1
                if written <= 1 or written % 1000 == 0 or written == len(rows):
                    print 'Wrote %d/%d %s' % (written, len(rows), full)
        finally:
            fh.close()
    finally:
        dst.close()

    chk = HecDss.open(meta['out_dss'])
    try:
        paths = chk.getPathnameList()
        print 'Catalog records:', len(paths)
        if len(paths) > 0:
            first = str(paths[0])
            g = chk.get(first).getGridInfo()
            print 'First path:', first
            print 'Grid:', g.getNumberOfCellsX(), 'x', g.getNumberOfCellsY(), 'cell', g.getCellSize(), 'units', g.getDataUnits(), 'type', g.getDataType()
            print 'Grid time:', g.getStartTime(), '|', g.getEndTime(), '|', g.getTimePeriod()
            print 'Cell zero:', g.getXCoordOfGridCellZero(), g.getYCoordOfGridCellZero()
            print 'Lower-left cell:', g.getLowerLeftCellX(), g.getLowerLeftCellY()
    finally:
        chk.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception('Usage: write_region_dss_grid_origin00_mm.py meta.json')
    main(sys.argv[1])
