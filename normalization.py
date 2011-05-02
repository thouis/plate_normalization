import xlrd
import numpy as np
from welltools import extract_row, extract_col

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
        self.num_iterations = 1
        self.gene_to_control_type = {}
        self.combine_replicates = False

        self.file_listeners = []
        self.parsing_listeners = []
        self.feature_selection_listeners = []

    def set_input_file(self, val):
        self.input_file = val
        try:
            self.book = xlrd.open_workbook(self.input_file)
            self.update_file_listeners()
        except:
            # XXX - report error
            pass

    def update_file_listeners(self):
        for f in self.file_listeners:
            f()

    def parsing_finished(self):
        for f in self.parsing_listeners:
            f()

    def feature_selection_finished(self):
        for f in self.feature_selection_listeners:
            f()

    def get_column_values(self, column_specifier):
        # import pdb
        # pdb.set_trace()
        return [cell.value for cell in self.book.sheet_by_index(column_specifier[0]).col(column_specifier[1])[1:]]

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

    def ready(self):
        for i in range(self.num_replicates):
            if self.replicate_features.get(i, None) is None:
                return False
        return True

    def set_replicate_feature(self, index, val):
        self.replicate_features[index] = val
        if self.ready():
            self.feature_selection_finished()

    def set_transformation(self, trans):
        assert trans in TRANSFORMATIONS
        self.transformation = trans

    def set_alignment_when(self, w):
        assert w in ALIGN_WHEN
        self.align_when = w

    def set_alignment_method(self, m):
        assert m in ALIGNMENT_METHODS
        self.alignment_method = m

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
        print "trna", vals.min(), vals.max()
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
                assert False, "Unknown normalization method: %s"%(self.alignement_method)
            assert len(align_values) > 0, "No valid wells on plate %s and replicate %d for alignment %s"%(plate, repindex + 1, self.alignement_method)
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

    def normalization_shift_rows(self):
        return self.normalization_plate_values

    def normalization_shift_columns(self):
        return self.normalization_plate_values

    def run_normalization(self):
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
            print self.align_when, iteration
            if (self.align_when == ALIGN_EACH) or (self.align_when == ALIGN_ONCE and iteration == 0):
                print "prep align"
                self.normalization_plate_values = self.normalization_align_plates()
                if iteration == 0:
                    self.normalization_first_alignment = self.normalization_plate_values
            # XXX - shift control populations
            self.normalization_plate_values = self.normalization_shift_rows()
            self.normalization_plate_values = self.normalization_shift_columns()
