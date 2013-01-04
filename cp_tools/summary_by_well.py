import sys
import os
import sqlite3
import xlwt

debug = False

def get_grouping(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM Per_Image")
    # ignore result, we just need the description
    column_names = [col[0] for col in cursor.description]
    return [c for c in ['Image_Metadata_Plate', 'Image_Metadata_Well',
                        'Image_Metadata_Row', 'Image_Metadata_Column',
                        'Image_Metadata_Barcode']
            if c in column_names]

def get_cell_cols(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM Per_Object")
    # ignore result, we just need the description
    column_names = [col[0] for col in cursor.description]
    column_names.remove('ImageNumber')
    column_names.remove('ObjectNumber')
    return column_names

def get_image_count_cols(conn):
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM Per_Image')
    column_names = [col[0] for col in cursor.description]
    if debug:
        print "Counters", column_names
    column_names = [c for c in column_names if c.startswith('Image_Count_')]
    return column_names

def count_objects(conn, grouping_keys, count_columns):
    cursor = conn.cursor()
    prefix = ", ".join("Per_Image.%s" % k for k in grouping_keys)
    counter = ", ".join("SUM(Per_Image.%s)" % c for c in count_columns)
    grouper = ", ".join("Per_Image.%s" % k for k in grouping_keys)
    cursor.execute("SELECT %s, %s FROM Per_Image "
                   "GROUP BY %s" % (prefix, counter, grouper)) 
    counts = {}
    for v in cursor.fetchall():
        counts[v[:len(grouping_keys)]] = v[len(grouping_keys):]
    return counts

def summarize_data(conn, grouping_keys, summary_columns, count_columns, counts):
    cursor = conn.cursor()
    prefix = ", ".join("Per_Image.%s" % k for k in grouping_keys)
    all_groups = cursor.execute("SELECT %s from Per_Image" % prefix).fetchall()
    summarizer = ", ".join("AVG(Per_Object.%s)" % c for c in summary_columns)
    grouper = ", ".join("Per_Image.%s" % k for k in grouping_keys)
    if debug:
        print ("SELECT %s, %s FROM Per_Image, Per_Object "
               "WHERE Per_Image.ImageNumber == Per_Object.ImageNumber "
               "GROUP BY %s" % (prefix, summarizer, grouper))
    cursor.execute("SELECT %s, %s FROM Per_Image, Per_Object "
                   "WHERE Per_Image.ImageNumber == Per_Object.ImageNumber "
                   "GROUP BY %s" % (prefix, summarizer, grouper))
    missing_prefixes = set(all_groups)
    while True:
        v = cursor.fetchone()
        if not v:
            break
        missing_prefixes -= set([v[:len(grouping_keys)]])
        vcounts = counts.get(v[:len(grouping_keys)], [0] * len(count_columns))
        yield v + vcounts
    for m in missing_prefixes:
        print "No Cells:", " ".join(m)
        yield m + tuple([0] * len(summary_columns + count_columns))

if __name__ == '__main__':
    in_db = sys.argv[1]
    out_xls = sys.argv[2]
    if not os.path.exists(in_db):
        raise ValueError('No such file ' + in_db)
    if os.path.exists(out_xls):
        raise RuntimeError('Will not overwrite ' + out_xls)

    conn = sqlite3.connect(in_db)
    grouping_keys = get_grouping(conn)
    assert ('Image_Metadata_Well' in grouping_keys) or \
        (('Image_Metadata_Row' in grouping_keys) and
         ('Image_Metadata_Column' in grouping_keys)), \
         "Need either Well or Row & Column in Image Metadata."

    summary_columns = get_cell_cols(conn)
    count_columns = get_image_count_cols(conn)
    counts = count_objects(conn, grouping_keys, count_columns)

    outbook = xlwt.Workbook(encoding='utf-8')
    sheet = outbook.add_sheet('Summary')
    for colidx, header in enumerate(grouping_keys + ['Mean ' + c for c in summary_columns] + count_columns):
        sheet.write(0, colidx, header)

    for rowidx, rowvals in enumerate(summarize_data(conn, grouping_keys, summary_columns, count_columns, counts)):
        for colidx, v in enumerate(rowvals):
            sheet.write(1 + rowidx, colidx, v)

    outbook.save(out_xls)
