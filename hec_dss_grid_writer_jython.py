# -*- coding: utf-8 -*-
"""HEC-DSSVue Jython 2.7 writer; timestamps are explicit accumulation bounds."""
from __future__ import print_function
import csv
from datetime import datetime, timedelta
import json
import os
import sys
import uuid

from hec.heclib.dss import HecDss
from hec.heclib.util import HecTime
from java.nio import ByteBuffer, ByteOrder
from java.nio.file import Files, Paths, StandardCopyOption
from jarray import zeros


def dss_time(value):
    # Avoid locale-dependent %b (DSS pathnames always use English month names).
    months = ('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC')
    return '%02d%s%04d:%02d%02d' % (value.day, months[value.month-1], value.year, value.hour, value.minute)


def interval(row):
    if row.get('start') and row.get('end'):
        return (datetime.strptime(row['start'], '%Y-%m-%dT%H:%M:%S'),
                datetime.strptime(row['end'], '%Y-%m-%dT%H:%M:%S'))
    start = datetime.strptime(row['date'], '%Y-%m-%d')
    return start, start + timedelta(days=1)


def make_path(meta, start, end):
    return '/%s/%s/%s/%s/%s/%s/' % (
        meta.get('a_part', ''), meta.get('b_part', ''),
        meta.get('c_part', 'PRECIPITATION'), dss_time(start), dss_time(end), meta.get('f_part', ''))


def read_float_array(handle, count):
    raw = handle.read(count*4)
    if len(raw) != count*4:
        raise ValueError('Truncated rainfall payload')
    buffer = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN)
    values = zeros(count, 'f')
    for i in range(count):
        values[i] = buffer.getFloat()
    return values


def find_template_grid(path):
    source = HecDss.open(path)
    try:
        for pathname in source.getPathnameList():
            if '/PRECIPITATION/' not in str(pathname).upper():
                continue
            container = source.get(str(pathname))
            if container.__class__.__name__.endswith('GridContainer'):
                return container
        raise ValueError('No readable precipitation GridContainer in ' + path)
    finally:
        source.close()


def read_rows(meta):
    handle = open(meta['dates_csv'], 'r')
    try:
        rows = list(csv.DictReader(handle))
    finally:
        handle.close()
    if not rows or len(rows) != int(meta['record_count']):
        raise ValueError('Date/payload record count mismatch')
    previous = None
    for row in rows:
        start, end = interval(row)
        if end <= start or (previous is not None and start != previous):
            raise ValueError('Non-contiguous or duplicate accumulation intervals')
        previous = end
    cells = int(meta['nx'])*int(meta['ny'])
    if os.path.getsize(meta['values_bin']) != len(rows)*cells*4:
        raise ValueError('Binary payload byte count mismatch')
    return rows


def validate_dss(path, meta, rows):
    dss = HecDss.open(path)
    try:
        actual = set(str(p) for p in dss.getPathnameList())
        expected = set(make_path(meta, *interval(row)) for row in rows)
        if actual != expected or len(expected) != len(rows):
            raise ValueError('DSS catalog dates/pathnames disagree with the payload')
        for row in rows:
            pathname = make_path(meta, *interval(row))
            container = dss.get(pathname)
            info = container.getGridInfo()
            checks = [
                (info.getNumberOfCellsX(), meta['nx']), (info.getNumberOfCellsY(), meta['ny']),
                (info.getLowerLeftCellX(), meta['lower_left_cell_x']),
                (info.getLowerLeftCellY(), meta['lower_left_cell_y']),
                (info.getDataType(), 1),
            ]
            if any(int(a) != int(b) for a,b in checks):
                raise ValueError('DSS grid dimensions/offsets/type disagree: ' + pathname)
            checks = [(info.getCellSize(), meta['cell_size']),
                      (info.getXCoordOfGridCellZero(), meta['cell_zero_x']),
                      (info.getYCoordOfGridCellZero(), meta['cell_zero_y'])]
            if any(abs(float(a)-float(b)) > max(.001, abs(float(b))*1e-7) for a,b in checks):
                raise ValueError('DSS cell size/origin disagree: ' + pathname)
            if str(info.getDataUnits()).strip().upper() != 'MM':
                raise ValueError('DSS units are not MM: ' + pathname)
            if ''.join(str(info.getSpatialReferenceSystem()).split()) != ''.join(meta['srs_wkt'].split()):
                raise ValueError('DSS projection WKT disagrees: ' + pathname)
            start, end = interval(row)
            for actual_time, expected_time in [(info.getStartTime(), start), (info.getEndTime(), end)]:
                actual_hec, expected_hec = HecTime(), HecTime()
                actual_hec.set(str(actual_time))
                expected_hec.set(dss_time(expected_time))
                if actual_hec.value() != expected_hec.value():
                    raise ValueError('DSS header accumulation times disagree: ' + pathname)
        # Read back actual first/last arrays, allowing the template's mm rounding.
        with open(meta['values_bin'], 'rb') as binary:
            cells = int(meta['nx'])*int(meta['ny'])
            for index in sorted(set([0, len(rows)-1])):
                binary.seek(index*cells*4)
                expected_values = read_float_array(binary, cells)
                actual_values = dss.get(make_path(meta, *interval(rows[index]))).getGridData().getData()
                if len(actual_values) != cells or any(
                        abs(float(a)-float(b)) > max(.011, abs(float(b))*1e-6)
                        for a,b in zip(actual_values,expected_values)):
                    raise ValueError('DSS rainfall array readback disagrees with payload')
        return len(actual)
    finally:
        dss.close()


def write_one(meta):
    rows = read_rows(meta)
    destination = meta['out_dss']
    if os.path.exists(destination) and not meta.get('overwrite_dss', False):
        raise ValueError('DSS already exists; use --overwrite: ' + destination)
    temporary = destination[:-4] + '.building-' + str(uuid.uuid4()) + '.dss'
    container = find_template_grid(meta['template_dss'])
    info = container.getGridInfo()
    info.setGridType(430)
    info.setCellInfo(int(meta['lower_left_cell_x']), int(meta['lower_left_cell_y']),
                     int(meta['nx']), int(meta['ny']), float(meta['cell_size']))
    info.setSpatialReference(meta['srs_name'], meta['srs_wkt'], float(meta['cell_zero_x']), float(meta['cell_zero_y']))
    info.setCoordOfGridCellZero(float(meta['cell_zero_x']), float(meta['cell_zero_y']))
    info.setDataType(1)
    info.setDataUnits('MM')
    cells = int(meta['nx'])*int(meta['ny'])
    dss = HecDss.open(temporary)
    try:
        handle = open(meta['values_bin'], 'rb')
        try:
            for number, row in enumerate(rows, 1):
                start, end = interval(row)
                info.setGridTimes(dss_time(start), dss_time(end))
                info.setDataInfo(float(row['min']), float(row['max']), float(row['mean']))
                data = container.getGridData()
                data.setData(read_float_array(handle, cells))
                data.initGridInfo(info)
                data.updateStatistics()
                pathname = make_path(meta, start, end)
                container.fullName = pathname
                container.setFullName(pathname)
                dss.put(container)
                if number == 1 or number % 500 == 0 or number == len(rows):
                    print('Wrote %d/%d %s' % (number, len(rows), destination))
        finally:
            handle.close()
    finally:
        dss.close()
    count = validate_dss(temporary, meta, rows)
    # Replace the destination only after every new record has passed readback.
    Files.move(Paths.get(temporary), Paths.get(destination), StandardCopyOption.REPLACE_EXISTING)
    return count


def main(path):
    with open(path) as handle:
        document = json.load(handle)
    paths = document.get('meta_paths', [path])
    results = []
    for meta_path in paths:
        with open(meta_path) as handle:
            meta = json.load(handle)
        try:
            count = write_one(meta)
            results.append(dict(file=meta['out_dss'], status='written', expected_records=meta['record_count'],
                                actual_records=count, nx=meta['nx'], ny=meta['ny'], units='MM', data_type=1, error=''))
        except Exception as error:
            print('FAILED:', str(error))
            results.append(dict(file=meta['out_dss'], status='failed', expected_records=meta['record_count'],
                                actual_records=0, nx=meta['nx'], ny=meta['ny'], units='', data_type='', error=str(error)))
    if document.get('validation_csv'):
        with open(document['validation_csv'], 'wb') as handle:
            writer = csv.DictWriter(handle, ['file','status','expected_records','actual_records',
                                            'nx','ny','units','data_type','error'])
            writer.writeheader()
            writer.writerows(results)
    if any(row['status'] == 'failed' for row in results):
        sys.exit(2)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Usage: Jython.exe -here hec_dss_grid_writer_jython.py manifest_or_meta.json')
    main(sys.argv[1])
