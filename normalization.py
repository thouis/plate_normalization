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
    lv = np.log(v) - np.log(1.0 - v)
    # deal with infinities
    lv[v <= 0.0] = - np.inf
    lv[v >= 1.0] = np.inf
    return lv

def inv_logit(lv):
    v = np.exp(lv) / (1 + np.exp(lv))
    v = (v - 0.0001) / 0.9998
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
        self.book = None

        # used for fetching columns
        self.cached_book = None
        self.cached_values = {}
        self.cached_replicate_data = {}

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
            self.cached_replicate_data = {}
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
        self.need_renorm = True
        if control_type != CONTROL_POPULATION:
            self.gene_to_control_type[gene] = control_type
        elif gene in self.gene_to_control_type:
            del self.gene_to_control_type[gene]
        self.parameter_changed()

    def ready(self):
        if self.num_replicates == 0:
            return False
        for i in range(self.num_replicates):
            if self.replicate_features.get(i, None) is None:
                return False
        return True

    def set_replicate_feature(self, index, val):
        self.need_renorm = True
        self.replicate_features[index] = val
        self.cached_replicate_data = {}
        self.parameter_changed()
        if self.ready():
            self.feature_selection_finished()

    def set_iterations(self, val):
        if self.num_iterations == val:
            return
        self.num_iterations = val
        self.parameter_changed()

    def set_transformation(self, trans):
        assert trans in TRANSFORMATIONS
        if self.transformation == trans:
            return
        self.transformation = trans
        self.parameter_changed()

    def set_alignment_when(self, w):
        assert w in ALIGN_WHEN
        if self.align_when == w:
            return
        self.align_when = w
        self.parameter_changed()

    def set_alignment_method(self, m):
        assert m in ALIGNMENT_METHODS
        if self.alignment_method == m:
            return
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
        if self.book != self.cached_book:
            self.cached_replicate_data = {}
            self.cached_values = {}
            self.cached_book = self.book
        if repindex not in self.cached_replicate_data:
            self.cached_replicate_data[repindex] = np.array([safe_float(v) for v in self.get_column_values(self.replicate_features[repindex])])
        vals = self.cached_replicate_data[repindex]
        if transformed:
            vals = self.transform_data(vals)
        return vals

    def get_transformed_values(self, repindex, cleaned=False):
        if not cleaned:
            vals = np.hstack([self.transformed_initial_plate_values[pl, rep]
                              for pl, rep in sorted(self.transformed_initial_plate_values.keys())
                              if rep == repindex]).flatten()
        else:
            vals = np.hstack([self.normalization_plate_values[pl, rep]
                              for pl, rep in sorted(self.normalization_plate_values.keys())
                              if rep == repindex]).flatten()
        return vals

    def get_transformed_control_groups(self):
        # needs to be in the same order as get_transformed_values()
        return (np.hstack([self.normalization_control_groups[pl] for pl in self.plate_names()]).flatten(),
                self.normalization_control_labels)

    def fetch_control_groups(self):
        control_keys = sorted(list(self.gene_to_control_type.keys())) # sorted so DMSO, GL2 before KIF11
        control_groups = [control_keys.index(g) + 1 if g in self.gene_to_control_type else 0
                          for g in self.fetch_genes()]
        return np.array(control_groups), ['tested'] + control_keys

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
        plate_mask = (np.array(self.get_column_values(self.plate_column)) == plate_name)
        indices = plate_mask.nonzero()[0]
        rows = np.array(self.fetch_rows()).view('uint32') - ord('A')
        cols = np.array(self.fetch_cols()).astype(int) - 1
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

            # XXX - should not shift plates that are more than half filled by controls
            # compute an offset per-plate and per-replicate
            if len(align_values) > 0:
                offsets[plate, repindex] = np.median(align_values)
            else:
                offsets[plate, repindex] = np.nan
        # shift offsets to zero-median to keep things identifiable
        if self.combine_replicates:
            # keep overall shift at 0
            global_shift = nanmedian(offsets.values())
            for plate, repindex in offsets:
                if np.isnan(offsets[plate, repindex]):
                    offsets[plate, repindex] = 0.0
                else:
                    offsets[plate, repindex] -= global_shift
                self.normalization_total_plate_shifts[plate, repindex] += offsets[plate, repindex]
            return dict(((plate, repindex), values - offsets[plate, repindex])
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())
        else:
            replicate_indices = np.array([repindex for _, repindex in offsets])
            offset_vals = np.array(offsets.values())
            per_replicate_shifts = dict((repindex, nanmedian(offset_vals[replicate_indices == repindex])) for repindex in range(self.num_replicates))
            for plate, repindex in offsets:
                if np.isnan(offsets[plate, repindex]):
                    offsets[plate, repindex] = 0.0
                else:
                    offsets[plate, repindex] -= per_replicate_shifts[repindex]
                self.normalization_total_plate_shifts[plate, repindex] += offsets[plate, repindex]
            return dict(((plate, repindex), values - offsets[plate, repindex])
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())

    def normalization_shift_rows_or_cols(self, axis, stacker, history):
        endshape = [-1, -1]
        endshape[axis] = 1
        if self.combine_replicates:
            all_plates = stacker(self.normalization_plate_values.keys())
            controls = (stacker([self.normalization_control_maps[pl] for pl, _ in self.normalization_plate_values.keys()]) != CONTROL_POPULATION)
            all_plates[controls] = np.nan
            # XXX - if a column is more than half controls, we should probably use adjacent columns to estimate its median
            offsets = fix_nans(nanmedian(all_plates, axis)).reshape(endshape)
            # shift offsets to zero-median to keep things identifiable
            offsets -= np.median(offsets)
            history += offsets
            return dict(((plate, repindex), values - offsets)
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())
        else:
            offsets = {}
            for repindex in range(self.num_replicates):
                rep_plates = stacker([v for (_, rep), v in self.normalization_plate_values.iteritems() if repindex == rep])
                controls = (stacker([self.normalization_control_maps[pl] for pl, rep in self.normalization_plate_values.keys() if repindex == rep]) != CONTROL_POPULATION)
                rep_plates[controls] = np.nan
                # XXX - if a column is more than half controls, we should probably use adjacent columns to estimate its median
                offsets[repindex] = fix_nans(nanmedian(rep_plates, axis)).reshape(endshape)
                # shift offsets to zero-median to keep things identifiable
                offsets[repindex] -= np.median(offsets[repindex])
                history[repindex] += offsets[repindex]
            return dict(((plate, repindex), values - offsets[repindex])
                        for ((plate, repindex), values) in self.normalization_plate_values.iteritems())

    def normalization_shift_rows(self):
        return self.normalization_shift_rows_or_cols(1, np.hstack, self.normalization_total_row_shifts)

    def normalization_shift_columns(self):
        return self.normalization_shift_rows_or_cols(0, np.vstack, self.normalization_total_col_shifts)

    def run_normalization(self):
        if not self.ready():
            return
        if self.need_renorm == False:
            return
        self.normalization_plate_values = {}
        self.normalization_control_maps = {}
        self.normalization_control_groups = {}
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

        # snapshot initial transformed state
        self.transformed_initial_plate_values = self.normalization_plate_values

        control_types = np.array(self.fetch_control_types(), dtype=object)
        for plate_name in set(plate_names):
            temp = np.zeros(self.plate_dims(), dtype=object)
            mask = (plate_names == plate_name)
            temp[rows[mask], cols[mask]] = control_types[mask]
            self.normalization_control_maps[plate_name] = temp

        control_groups, self.normalization_control_labels = self.fetch_control_groups()
        control_groups = np.array(control_groups)
        for plate_name in set(plate_names):
            temp = np.zeros(self.plate_dims(), dtype=object)
            mask = (plate_names == plate_name)
            temp[rows[mask], cols[mask]] = control_groups[mask]
            self.normalization_control_groups[plate_name] = temp


        # total shifts
        self.normalization_total_plate_shifts = dict(((pl, repindex), 0.0) for (pl, repindex) in self.normalization_plate_values)
        if self.combine_replicates:
            self.normalization_total_row_shifts = np.zeros((self.plate_dims[0], 1))
            self.normalization_total_col_shifts = np.zeros((1, self.plate_dims[1]))
        else:
            self.normalization_total_row_shifts = dict((idx, np.zeros((self.plate_dims()[0], 1))) for idx in range(self.num_replicates))
            self.normalization_total_col_shifts = dict((idx, np.zeros((1, self.plate_dims()[1]))) for idx in range(self.num_replicates))

        for iteration in range(self.num_iterations):
            if (self.align_when == ALIGN_EACH) or (self.align_when == ALIGN_ONCE and iteration == 0):
                self.normalization_plate_values = self.normalization_align_plates()
            if iteration == 0:
                self.normalization_first_alignment = self.normalization_plate_values

            # XXX - shift control populations
            self.normalization_plate_values = self.normalization_shift_rows()
            self.normalization_plate_values = self.normalization_shift_columns()

        # XXX - record shifts applied

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

        def write_row(rowidx, basecol, *args):
            for idx, v in enumerate(args):
                results_sheet.write(rowidx, basecol + idx, v)

        def get_normalized_value(pl, r, c, repidx):
            try:
                return self.normalization_plate_values[pl, repidx][ord(r) - ord('A'), int(c) - 1]
            except:
                return ""

        # XXX - should write Well = Row+Column
        write_row(0, 0, "Plate", "Row", "Column", "Gene", *["rep%d"%(r + 1) for r in range(self.num_replicates)])
        for rowidx, (pl, r, c, g) in enumerate(zip(plates, rows, cols, genes)):
            vals = [get_normalized_value(pl, r, c, ridx) for ridx in range(self.num_replicates)]
            vals = [v if np.isfinite(v) else "" for v in vals]
            write_row(rowidx + 1, 0, pl, r, c, g, *vals)

        # Write per-replicate median and MAD of all wells, population, negative controls, positive controls
        dest_column = len(vals) + 4 + 2 
        dest_row = [3]
        def write_summary(*vals):
            write_row(dest_row[0], dest_column, *vals)
            dest_row[0] += 1

        # fetch values
        allplates = list(set(plates))
        allvals = {}
        for ridx in range(self.num_replicates):
            allvals[ridx] = np.dstack([self.normalization_plate_values[pl, ridx] for pl in allplates])
        all_control_maps = np.dstack([self.normalization_control_maps[pl] for pl in allplates])

        # all wells
        write_summary("", "Median / MAD")
        write_summary("", *["R%d"%(ridx + 1) for ridx in range(self.num_replicates)])
        meds = []
        sigma_MADs = []
        for ridx in range(self.num_replicates):
            vals = allvals[ridx].flatten()
            meds.append(np.median(vals))
            sigma_MADs.append(1.4826 * np.median(np.abs(vals - meds[-1])))
        write_summary("Median All Wells", *meds)
        write_summary("sigma_MAD All Wells", *sigma_MADs)

        # population
        write_summary("")
        meds = []
        sigma_MADs = []
        for ridx in range(self.num_replicates):
            vals = allvals[ridx][all_control_maps == CONTROL_POPULATION]
            meds.append(np.median(vals))
            sigma_MADs.append(1.4826 * np.median(np.abs(vals - meds[-1])))
        write_summary("Median Population", *meds)
        write_summary("sigma_MAD Population", *sigma_MADs)

        # negative controls
        write_summary("")
        meds = []
        sigma_MADs = []
        for ridx in range(self.num_replicates):
            vals = allvals[ridx][all_control_maps == CONTROL_NEGATIVE]
            meds.append(np.median(vals))
            sigma_MADs.append(1.4826 * np.median(np.abs(vals - meds[-1])))
        write_summary("Median Neg. Controls", *meds)
        write_summary("sigma_MAD Neg. Controls", *sigma_MADs)

        # positive controls
        write_summary("")
        meds = []
        sigma_MADs = []
        for ridx in range(self.num_replicates):
            vals = allvals[ridx][all_control_maps == CONTROL_POSITIVE]
            meds.append(np.median(vals))
            sigma_MADs.append(1.4826 * np.median(np.abs(vals - meds[-1])))
        write_summary("Median Pos. Controls", *meds)
        write_summary("sigma_MAD Pos. Controls", *sigma_MADs)

        # Provenance sheet
        nrows = len(provenance_sheet.rows)
        features = [self.replicate_features[repidx] for repidx in range(self.num_replicates)]
        feature_names = [(self.book.sheet_names()[f[0]], self.book.sheet_by_index(f[0]).row(0)[f[1]].value) for f in features]
        provenance = ([["Normalization", results_sheet.name, datetime.date.today().isoformat(), getpass.getuser()],
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
