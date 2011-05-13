import xlrd
import xlutils.copy
import numpy as np
from welltools import extract_row, extract_col
from scipy.stats import nanmean, nanmedian
import os.path
import datetime
import getpass


DETECT = 'Detect'

TRANSFORM_NONE = 'None'
TRANSFORM_LOGARITHM = 'Logarithm (counts)'
TRANSFORM_LOGIT_FRACTION = 'Logit (fractions)'
TRANSFORM_LOGIT_PERCENT = 'Logit (percentages)'
TRANSFORMATIONS = [TRANSFORM_NONE,
                   TRANSFORM_LOGARITHM,
                   TRANSFORM_LOGIT_FRACTION,
                   TRANSFORM_LOGIT_PERCENT]

ALIGN_NEVER = "Never"
ALIGN_ONCE = "Once"
ALIGN_EACH = "Each iteration"
ALIGN_WHEN = [ALIGN_NEVER,
              ALIGN_ONCE,
              ALIGN_EACH]

ALIGN_NEGATIVE_CONTROLS = "Negative Controls"
ALIGN_POPULATION = "All But Controls"
ALIGN_EVERYTHING = "All Wells"
ALIGNMENT_METHODS = [ALIGN_NEGATIVE_CONTROLS,
                     ALIGN_POPULATION,
                     ALIGN_EVERYTHING]

CONTROL_POPULATION = "Tested Population"
CONTROL_NEGATIVE = "Negative Control"
CONTROL_POSITIVE = "Positive Control"
CONTROL_TYPES = [CONTROL_POPULATION,
                 CONTROL_NEGATIVE,
                 CONTROL_POSITIVE]

def safe_float(s):
    try:
        return float(s)
    except:
        return float('nan')

def logit(v):
    # make 0/1 representable, linear interpolate between 0.0001 and 0.9999
    v = 0.0001 + v * 0.9998
    lv = np.log(v) / np.log(1.0 - v)
    # deal with infinities
    lv[v <= 0.0] = - np.inf
    lv[v >= 1.0] = np.inf
    return lv

def fix_nans(shift_vals):
    # some rows or columns might be entirely controls, in which case
    # we take the normalization from the nearest non-nan rows/columns
    shift_vals = shift_vals.copy()
    bad = ~ np.isfinite(shift_vals)
    assert not all(bad)
    while any(bad):
        up = np.roll(shift_vals, 1)
        down = np.roll(shift_vals, -1)
        up[0] = np.nan
        down[-1] = np.nan
        shift_vals[bad] = nanmean(np.vstack((up, down)), 0)[bad]
        bad = ~ np.isfinite(shift_vals)
    return shift_vals

class Normalization(object):
    '''This object communicates the parameters (including input and output files) for a normalization'''
    def __init__(self):
        self.input_file = ''
        self.output_file = ''
        self.shape = DETECT
        self.detected_384 = False
        self.plate_column = ()
        self.well_column = ()
        self.wellrow_column = ()
        self.wellcol_column = ()
        self.combined_wellrowcol = True
        self.gene_column = ()
        self.num_replicates = 0
        self.replicate_features = {}
        self.transformation = TRANSFORMATIONS[0]
        self.align_when = ALIGN_WHEN[0]
        self.alignment_method = ALIGNMENT_METHODS[0]
        self.num_iterations = 10
        self.gene_to_control_type = {}
        self.combine_replicates = False
        self.iterations = 10
        self.book = None

        # used for fetching columns
        self.cached_book = None
        self.cached_values = {}

        self.need_renorm = True

        self.file_listeners = []
        self.parsing_listeners = []
        self.feature_selection_listeners = []
        self.parameter_change_listeners = []

    def set_input_file(self, val):
        self.input_file = val
        try:
            self.book = xlrd.open_workbook(self.input_file)
            self.update_file_listeners()
        except:
            # XXX - report error
            pass

    def set_output_file(self, val):
        self.output_file = val

    def update_file_listeners(self):
        for f in self.file_listeners:
            f()

    def parsing_finished(self):
        for f in self.parsing_listeners:
            f()

    def feature_selection_finished(self):
        self.need_renorm = True
        for f in self.feature_selection_listeners:
            f()

    def parameter_changed(self):
        self.need_renorm = True
        for f in self.parameter_change_listeners:
            f()

    def get_column_values(self, column_specifier):
        if self.book != self.cached_book:
            self.cached_values = {}
            self.cached_book = self.book
        if column_specifier not in self.cached_values:
            self.cached_values[column_specifier] = [cell.value for cell in self.book.sheet_by_index(column_specifier[0]).col(column_specifier[1])[1:]]
        return self.cached_values[column_specifier]

    def fetch_plates(self):
        return self.get_column_values(self.plate_column)

    def fetch_rows(self):
        if self.combined_wellrowcol:
            return [extract_row(v) for v in self.get_column_values(self.well_column)]
        else:
            return self.get_column_values(self.wellrow_column)

    def fetch_cols(self):
        if self.combined_wellrowcol:
            return [extract_col(v) for v in self.get_column_values(self.well_column)]
        else:
            return self.get_column_values(self.wellcol_column)

    def fetch_genes(self):
        return self.get_column_values(self.gene_column)

    def fetch_control_types(self):
        return [self.gene_to_control_type.get(g, CONTROL_POPULATION) for g in self.fetch_genes()]

    def set_control_type(self, gene, control_type):
        self.gene_to_control_type[gene] = control_type

    def ready(self):
        if self.num_replicates == 0:
            return False
        for i in range(self.num_replicates):
            if self.replicate_features.get(i, None) is None:
                return False
        return True

    def set_replicate_feature(self, index, val):
        self.replicate_features[index] = val
        if self.ready():
            self.feature_selection_finished()

    def set_iterations(self, val):
        self.iterations = val
        self.parameter_changed()

    def set_transformation(self, trans):
        assert trans in TRANSFORMATIONS
        self.transformation = trans
        self.parameter_changed()

    def set_alignment_when(self, w):
        assert w in ALIGN_WHEN
        self.align_when = w
        self.parameter_changed()

    def set_alignment_method(self, m):
        assert m in ALIGNMENT_METHODS
        self.alignment_method = m
        self.parameter_changed()

    def transform_data(self, vals):
        if self.transformation == TRANSFORM_NONE:
            return vals
        elif self.transformation == TRANSFORM_LOGARITHM:
            return np.log(vals)
        elif self.transformation == TRANSFORM_LOGIT_FRACTION:
            return logit(vals)
        elif self.transformation == TRANSFORM_LOGIT_PERCENT:
            return logit(vals / 100.0)

    def get_replicate_data(self, repindex, transformed=False):
        vals = np.array([safe_float(v) for v in self.get_column_values(self.replicate_features[repindex])])
        if transformed:
            vals = self.transform_data(vals)
        return vals

    def num_plates(self):
        return len(set(self.get_column_values(self.plate_column)))

    def plate_names(self):
        return sorted(set(self.get_column_values(self.plate_column)))

    def plate_shape(self):
        if self.shape == DETECT:
            return "384" if self.detected_384 else "96"
        return self.shape

    def plate_dims(self):
        return (16, 24) if self.plate_shape() == "384" else (8, 12)

    def plate_array(self, plate_name, repindex, transformed=False):
        if isinstance(plate_name, int):
            plate_name = self.plate_names()[plate_name]
        plate_mask = np.array([v == plate_name for v in self.get_column_values(self.plate_column)])
        indices = plate_mask.nonzero()[0]
        rows = np.array([ord(r) - ord('A') for r in self.fetch_rows()])
        cols = np.array([int(c) - 1 for c in self.fetch_cols()])
        vals = np.array(self.get_replicate_data(repindex, transformed))
        output = np.zeros(self.plate_dims(), np.float)
        output[rows[indices], cols[indices]] = vals[indices]
        return output

    # NORMALIZATION CODE
    # - self.normalization_* variables hold temporary results.
    # In particular, the steps of normalization are kept in
    # self.normalization_plate_values, a dictionary indexed by
    # (platename, repindex).

    def normalization_align_plates(self):
        # XXX - should not shift plates that are more than half filled by controls
        # compute an offset per-plate and per-replicate
        offsets = {}
        for (plate, repindex), values in self.normalization_plate_values.iteritems():
            control_map = self.normalization_control_maps[plate]
            if self.alignment_method == ALIGN_NEGATIVE_CONTROLS:
                align_values = values[control_map == CONTROL_NEGATIVE]
            elif self.alignment_method == ALIGN_POPULATION:
                align_values = values[control_map == CONTROL_POPULATION]
            elif self.alignment_method == ALIGN_EVERYTHING:
                align_values = values
            else:
                assert False, "Unknown normalization method: %s"%(self.alignment_method)
            assert len(align_values) > 0, "No valid wells on plate %s and replicate %d for alignment using %s"%(plate, repindex + 1, self.alignment_method)
            offsets[plate, repindex] = np.median(align_values)
        # shift offsets to zero-median to keep things identifiable
        if self.combine_replicates:
            global_shift = np.median(offsets.values())

            print "	", offsets
            print "	global: ", global_shift
            return dict(((plate, repindex), values - (offsets[plate, repindex] - global_shift))
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())
        else:
            replicate_indices = np.array([repindex for _, repindex in offsets])
            offset_vals = np.array(offsets.values())
            per_replicate_shifts = dict((repindex, np.median(offset_vals[replicate_indices == repindex])) for repindex in range(self.num_replicates))
            print "	", offsets
            print "	per-rep: ", per_replicate_shifts
            return dict(((plate, repindex), values - (offsets[plate, repindex] - per_replicate_shifts[repindex]))
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())

    def normalization_shift_rows_or_cols(self, axis, stacker):
        endshape = [-1, -1]
        endshape[axis] = 1
        if self.combine_replicates:
            all_plates = stacker(self.normalization_plate_values.keys())
            controls = (stacker([self.normalization_control_maps[pl] for pl, _ in self.normalization_plate_values.keys()]) != CONTROL_POPULATION)
            all_plates[controls] = np.nan
            offsets = fix_nans(nanmedian(all_plates, axis))
            # shift offsets to zero-median to keep things identifiable
            offsets -= np.median(offsets)
            return dict(((plate, repindex), values - offsets.reshape(endshape))
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())
        else:
            offsets = {}
            for repindex in range(self.num_replicates):
                rep_plates = stacker([v for (_, rep), v in self.normalization_plate_values.iteritems() if repindex == rep])
                controls = (stacker([self.normalization_control_maps[pl] for pl, rep in self.normalization_plate_values.keys() if repindex == rep]) != CONTROL_POPULATION)
                rep_plates[controls] = np.nan
                offsets[repindex] = fix_nans(nanmedian(rep_plates, axis))
                # shift offsets to zero-median to keep things identifiable
                offsets[repindex] -= np.median(offsets[repindex])
            return dict(((plate, repindex), values - offsets[repindex].reshape(endshape))
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())

    def normalization_shift_rows(self):
        return self.normalization_shift_rows_or_cols(1, np.hstack)

    def normalization_shift_columns(self):
        return self.normalization_shift_rows_or_cols(0, np.vstack)

    def run_normalization(self):
        if not self.ready():
            return
        if self.need_renorm == False:
            return
        self.normalization_plate_values = {}
        self.normalization_control_maps = {}
        rows = np.array([ord(r) - ord('A') for r in self.fetch_rows()])
        cols = np.array([int(c) - 1 for c in self.fetch_cols()])

        plate_names = np.array([v for v in self.get_column_values(self.plate_column)], dtype=object)
        for repindex in range(self.num_replicates):
            vals = np.array(self.get_replicate_data(repindex, True))
            for plate_name in set(plate_names):
                # fetch transformed values
                temp = np.zeros(self.plate_dims(), dtype=np.float)
                mask = (plate_names == plate_name)
                temp[rows[mask], cols[mask]] = vals[mask]
                self.normalization_plate_values[plate_name, repindex] = temp

        control_types = np.array(self.fetch_control_types(), dtype=object)
        for plate_name in set(plate_names):
            temp = np.zeros(self.plate_dims(), dtype=object)
            mask = (plate_names == plate_name)
            temp[rows[mask], cols[mask]] = control_types[mask]
            self.normalization_control_maps[plate_name] = temp

        for iteration in range(self.num_iterations):
            if (self.align_when == ALIGN_EACH) or (self.align_when == ALIGN_ONCE and iteration == 0):
                self.normalization_plate_values = self.normalization_align_plates()
            if iteration == 0:
                self.normalization_first_alignment = self.normalization_plate_values

            # XXX - shift control populations
            self.normalization_plate_values = self.normalization_shift_rows()
            self.normalization_plate_values = self.normalization_shift_columns()

        self.need_renorm = False

    def save(self):
        assert not os.path.exists(self.output_file), "Will not overwrite file %s"%(self.output_file)
        self.run_normalization()

        # first, copy all pages from the existing book
        outbook, results_sheet, provenance_sheet = duplicate_xlbook(self.book)

        # write out normalization values in order they appear in the plate/well columns
        plates = self.fetch_plates()
        rows = self.fetch_rows()
        cols = self.fetch_cols()
        genes = self.fetch_genes()

        def write_row(rowidx, pl, r, c, g, vals):
            results_sheet.write(rowidx, 0, pl)
            results_sheet.write(rowidx, 1, r)
            results_sheet.write(rowidx, 2, c)
            results_sheet.write(rowidx, 3, g)
            for vidx, v in enumerate(vals):
                results_sheet.write(rowidx, 4 + vidx, v)

        def get_normalized_value(pl, r, c, repidx):
            try:
                return self.normalization_plate_values[pl, repidx][ord(r) - ord('A'), int(c) - 1]
            except:
                return ""

        write_row(0, "Plate", "Row", "Column", "Gene", ["rep%d"%(r + 1) for r in range(self.num_replicates)])
        for rowidx, (pl, r, c, g) in enumerate(zip(plates, rows, cols, genes)):
            vals = [get_normalized_value(pl, r, c, ridx) for ridx in range(self.num_replicates)]
            write_row(rowidx + 1, pl, r, c, g, vals)

        nrows = len(provenance_sheet.rows)
        print "Provenance", nrows
        features = [self.replicate_features[repidx] for repidx in range(self.num_replicates)]
        feature_names = [(self.book.sheet_names()[f[0]], self.book.sheet_by_index(f[0]).row(0)[f[1]].value) for f in features]
        provenance = ([["Normalization", datetime.date.today().isoformat(), getpass.getuser()],
                       ["Original file", os.path.abspath(self.input_file)],
                       ["Features"]] +
                      [["", "Sheet:", f[0], "Column:", f[1]] for f in feature_names] +
                      [["Number of iterations:", self.num_iterations],
                       ["Transformation:", self.transformation],
                       ["When to align?", self.align_when],
                       ["How to align?", self.alignment_method],
                       ["Control(s):"]] +
                      [["", g, t] for g, t in self.gene_to_control_type.iteritems() if t != CONTROL_POPULATION])

        for idx, l in enumerate(provenance):
            for c, v in enumerate(l):
                provenance_sheet.write(nrows + 1 + idx, c, v)

        outbook.save(self.output_file)


def duplicate_xlbook(book):
    # add two sheets, one for normalization results, one for provenance tracking
    existing_sheets = book.sheet_names()
    count = 1
    while 'Normalization %d'%(count) in existing_sheets:
        count += 1
    normalization_sheet_name = 'Normalization %d'%(count)

    provenance_sheet_name = 'BFX Provenance Tracking'
    add_provenance_sheet = (provenance_sheet_name not in existing_sheets)

    # place the results sheet last, unless the BFX provenance sheet is
    # last, in which case, just before it.
    normalization_sheet_last = (existing_sheets[-1] != provenance_sheet_name)

    from xlutils.filter import process,XLRDReader,XLWTWriter
    import xlwt

    class WrapWT(XLWTWriter):
        def sheet(self, rdsheet, wtsheet_name):
            # if we're just adding a normalization sheet, add it before the provenance sheet
            if (not add_provenance_sheet) and (wtsheet_name == provenance_sheet_name):
                self.normalization_sheet = self.wtbook.add_sheet(normalization_sheet_name)
            # write the sheet as requested
            XLWTWriter.sheet(self, rdsheet, wtsheet_name)

        def finish(self):
            # add two new sheets at the end if we're adding normalization and results
            if add_provenance_sheet:
                self.normalization_sheet = self.wtbook.add_sheet(normalization_sheet_name)
                self.wtbook.add_sheet(provenance_sheet_name)
            XLWTWriter.finish(self)

    # hackity hack - need to force non-ascii encoding
    orig_Workbook_init = xlwt.Workbook.__init__
    def replacement(self, encoding='utf-8', style_compression=0):
        orig_Workbook_init(self, encoding, style_compression)
    xlwt.Workbook.__init__ = replacement
    w = WrapWT()
    process(XLRDReader(book, 'unknown.xls'), w)
    outbook = w.output[0][1]
    xlwt.Workbook.__init__ = orig_Workbook_init

    for provenance_sheet_idx in range(len(existing_sheets) + 2):
        if outbook.get_sheet(provenance_sheet_idx).name == provenance_sheet_name:
            break

    return w.output[0][1], w.normalization_sheet, outbook.get_sheet(provenance_sheet_idx)
