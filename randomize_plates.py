import numpy as np
import xlrd
import xlwt
import sys
import collections
import os.path

# Input:
# Header row, followed by columns of:
#  siRNA ID, Control?, Gene, Plate number, Row, Column, Plate Name
# Controls will not be moved

ID = 0
CTRL = 1
GENE = 2
BOUNDARY = 3  # everything before this is randomized
PNUM = 3
PROW = 4
PCOL = 5
PNAME = 6
NUM_COLS = 7  # keep this accurate

# for output
col_types = [unicode, float, unicode, float, unicode, float, unicode]  # excel stores ints as floats.  yay.

book = xlrd.open_workbook(sys.argv[1])
sheet = book.sheet_by_index(0)
rows = np.array([sheet.row_values(ridx, end_colx=NUM_COLS) for ridx in range(1, sheet.nrows)])

# sort by control status
rows = rows[np.argsort(rows[:, CTRL].astype(int).flatten()), :]

# only randomize the non-controls
num_non_control = sum(rows[:, CTRL].astype(int) == 0)

# this slice is a view on the original
nc_rows = rows[:num_non_control, :]
nc_rows[:, :BOUNDARY] = nc_rows[np.argsort(np.random.uniform(0, 1, num_non_control)), :BOUNDARY]

# Find any genes with more than one sirna on the border.  Randomly
# swap them with other locations until there aren't any genes with
# more than one such siRNA.
while True:
    border = ((nc_rows[:, PROW] == max(rows[:, PROW])) |
              (nc_rows[:, PROW] == min(rows[:, PROW])) |
              (nc_rows[:, PCOL] == max(rows[:, PCOL])) |
              (nc_rows[:, PCOL] == min(rows[:, PCOL])))
    genecounts = collections.Counter(nc_rows[border, GENE])
    totswaps = 0
    for g in genecounts:
        if genecounts[g] == 1:
            continue
        # swap all but the first non-border case of this gene
        idxs = np.nonzero(border & (nc_rows[:, GENE] == g))[0][1:]
        for i in idxs:
            ridx = int(np.random.uniform(0, nc_rows.shape[0]))
            nc_rows[i, :BOUNDARY], nc_rows[ridx, :BOUNDARY] = nc_rows[ridx, :BOUNDARY], nc_rows[i, :BOUNDARY]
            totswaps += 1
    print totswaps, "swaps"
    if totswaps == 0:
        break

# write out result
outbook = xlwt.Workbook(encoding='utf-8')
outsheet = outbook.add_sheet(sheet.name)
for ridx in range(sheet.nrows):
    for cidx in range(NUM_COLS):
        if ridx == 0:
            outsheet.write(ridx, cidx, sheet.cell_value(ridx, cidx))
        else:
            outsheet.write(ridx, cidx, col_types[cidx](rows[ridx - 1, cidx]))
assert not os.path.exists(sys.argv[2]), "Won't overwrite existing file " + str(sys.argv[2])
outbook.save(sys.argv[2])
