import wxversion
wxversion.select("2.8")
import wx, wx.html, wx.lib.scrolledpanel
import sys
import random
import os.path
import re
from normalization import Normalization, DETECT, TRANSFORMATIONS, ALIGN_WHEN, ALIGNMENT_METHODS, ALIGN_NEVER, CONTROL_POPULATION, CONTROL_CONTROL
from wrapfilename import wrap_filename
import wxplotpanel
import traceback
import numpy as np
from  matplotlib.patches import Rectangle
from errordialog import display_error_dialog
try:
    from collections import OrderedDict
except:
    from ordereddict import OrderedDict  # available on pypi

app_name = "Plate Normalizer"
aboutText = """<p>Plate normalizer v0.1.</p>"""

class HtmlWindow(wx.html.HtmlWindow):
    def __init__(self, parent, id, size=(600, 400)):
        wx.html.HtmlWindow.__init__(self, parent, id, size=size)

    def OnLinkClicked(self, link):
        wx.LaunchDefaultBrowser(link.GetHref())

class AboutBox(wx.Dialog):
    def __init__(self):
        wx.Dialog.__init__(self, None, -1, "About...",
                           style=wx.DEFAULT_DIALOG_STYLE | wx.THICK_FRAME)
        hwin = HtmlWindow(self, -1, size=(400, 200))
        hwin.SetPage(aboutText)
        self.SetClientSize(hwin.GetSize())
        self.CentreOnParent(wx.BOTH)
        self.SetFocus()


class TabPanel(wx.Panel):
    #----------------------------------------------------------------------
    def __init__(self, parent, *args):
        """ """
        wx.Panel.__init__(self, parent=parent)

        colors = ["red", "blue", "gray", "yellow", "green"]
        self.SetBackgroundColour(random.choice(colors))

        btn = wx.Button(self, label="Press Me")
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(btn, 0, wx.ALL, 10)
        self.SetSizer(sizer)

class ColumnSelector(wx.Panel):
    def __init__(self, parent, callback, substring_hints, normalization, callback_args=[]):
        """ """
        wx.Panel.__init__(self, parent=parent)
        self.callback = callback
        self.substring_hints = substring_hints
        self.normalization = normalization
        self.sheet_idx = 0
        self.prev_colidx = None
        self.callback_args = callback_args

        self.sheet_selector = wx.ComboBox(self, -1, choices=[], style=wx.CB_READONLY)
        self.column_selector = wx.ComboBox(self, -1, choices=[], style=wx.CB_READONLY)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.sheet_selector, 1, wx.EXPAND | wx.RIGHT, 5)
        sizer.Add(self.column_selector, 1, wx.EXPAND)
        self.SetSizer(sizer)

        self.sheet_selector.Bind(wx.EVT_COMBOBOX, self.select_sheet)
        self.column_selector.Bind(wx.EVT_COMBOBOX, self.select_column)

        normalization.file_listeners.append(self.input_file_updated)
        if normalization.input_file:
            self.input_file_updated()

    def input_file_updated(self):
        # find sheet and column names
        self.sheet_selector.Clear()
        names = self.normalization.book.sheet_names()
        self.sheet_selector.AppendItems(names)
        self.sheet_selector.Value = names[0]
        newevt = wx.PyCommandEvent(wx.wxEVT_COMMAND_COMBOBOX_SELECTED)
        newevt.SetInt(0)
        wx.PostEvent(self.sheet_selector, newevt)

    def select_sheet(self, evt):
        sheet_idx = evt.GetSelection()
        self.sheet_idx = sheet_idx
        self.column_selector.Clear()
        column_names = [c.value for c in self.normalization.book.sheet_by_index(sheet_idx).row(0)]
        self.column_selector.AppendItems(column_names)
        self.prev_colidx = None
        default = min([idx for idx, name in enumerate(column_names)
                       if any([hint.lower() in name.lower() for hint in self.substring_hints])]
                      or [-1])
        if default > -1:
            self.column_selector.Selection = default
            newevt = wx.PyCommandEvent(wx.wxEVT_COMMAND_COMBOBOX_SELECTED)
            newevt.SetInt(default)
            wx.PostEvent(self.column_selector, newevt)

    def set_default_sheet(self, idx):
        # set a default sheet and trigger a select event
        self.sheet_selector.Selection = idx
        newevt = wx.PyCommandEvent(wx.wxEVT_COMMAND_COMBOBOX_SELECTED)
        newevt.SetInt(idx)
        wx.PostEvent(self.sheet_selector, newevt)

    def select_column(self, evt):
        colidx = evt.GetSelection()
        if colidx != self.prev_colidx:
            self.prev_colidx = colidx
            self.callback((self.sheet_idx, colidx), *self.callback_args)

class PlateLayout(wx.Panel):
    def __init__(self, parent, normalization):
        wx.Panel.__init__(self, parent=parent)
        self.normalization = normalization

        input_box = wx.StaticBox(self, wx.ID_ANY, 'Input')
        input_sizer = wx.StaticBoxSizer(input_box, wx.VERTICAL)
        lower_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.input_text = wx.TextCtrl(self, -1, normalization.input_file, style=wx.TE_RIGHT)
        self.bfx_format = wx.CheckBox(self, label='BFX format?')
        input_browse = wx.Button(self, label="Browse")
        lower_sizer.Add(self.bfx_format)
        lower_sizer.AddStretchSpacer()
        lower_sizer.Add(input_browse)
        input_sizer.Add(self.input_text, flag=wx.EXPAND)
        input_sizer.Add(lower_sizer, flag=wx.EXPAND | wx.TOP, border=5)

        shape_box = wx.StaticBox(self, wx.ID_ANY, 'Plate shape')
        shape_sizer = wx.StaticBoxSizer(shape_box, wx.HORIZONTAL)
        shapeb1 = wx.RadioButton(self, -1, '96', style=wx.RB_GROUP)
        shapeb2 = wx.RadioButton(self, -1, '384')
        shapeb3 = wx.RadioButton(self, -1, DETECT)
        shape_sizer.Add((1, 1), 2)
        shape_sizer.Add(shapeb1, 0)
        shape_sizer.Add((1, 1), 1)
        shape_sizer.Add(shapeb2, 0)
        shape_sizer.Add((1, 1), 1)
        shape_sizer.Add(shapeb3, 0)
        shape_sizer.Add((1, 1), 2)
        shapeb3.Value = 1

        plate_column_box = wx.StaticBox(self, wx.ID_ANY, 'Plate column in spreadsheet')
        plate_column_sizer = wx.StaticBoxSizer(plate_column_box, wx.VERTICAL)
        plate_column_selector = ColumnSelector(self, self.set_plate_column, ['plate'], self.normalization)
        plate_column_sizer.Add(plate_column_selector, 0, wx.EXPAND)

        well_column_box = wx.StaticBox(self, wx.ID_ANY, 'Well column(s) in spreadsheet')
        self.well_column_sizer = wx.StaticBoxSizer(well_column_box, wx.VERTICAL)
        wells_combined = self.wells_combined = wx.RadioButton(self, -1, 'Wells in single spreadsheet column', style=wx.RB_GROUP)
        wells_separate = wx.RadioButton(self, -1, 'Rows && columns in separate spreadsheet columns')
        self.well_selector = ColumnSelector(self, self.set_well_column, ['well'], self.normalization)
        self.wellrow_selector = ColumnSelector(self, self.set_wellrow_column, ['row'], self.normalization)
        self.wellcol_selector = ColumnSelector(self, self.set_wellcol_column, ['col'], self.normalization)
        self.well_column_sizer.Add(wells_combined, 0)
        self.well_column_sizer.Add(wells_separate, 0, wx.TOP, 5)
        self.well_column_sizer.Add(self.well_selector, 0, wx.EXPAND | wx.TOP, 5)
        self.well_column_sizer.Add(self.wellrow_selector, 0, wx.EXPAND | wx.TOP, 5)
        self.well_column_sizer.Add(self.wellcol_selector, 0, wx.EXPAND | wx.TOP, 5)

        gene_column_box = wx.StaticBox(self, wx.ID_ANY, 'Gene/Chemical column in spreadsheet')
        gene_column_sizer = wx.StaticBoxSizer(gene_column_box, wx.VERTICAL)
        self.gene_column_selector = ColumnSelector(self, self.set_gene_column, ['gene', 'compound', 'treatment', 'content'],
                                                   self.normalization)
        gene_column_sizer.Add(self.gene_column_selector, 0, wx.EXPAND)

        status_box = wx.StaticBox(self, wx.ID_ANY, 'Status')
        status_sizer = wx.StaticBoxSizer(status_box, wx.VERTICAL)
        self.status_text = wx.StaticText(self, -1, "Choose settings above...")
        status_sizer.Add(self.status_text, 0, wx.EXPAND)

        # make things align
        self.status_text.SetFont(wx.Font(pointSize=-1, family=wx.FONTFAMILY_MODERN,
                                         style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL))

        wells_combined.Value = True
        self.well_column_sizer.Hide(self.wellrow_selector)
        self.well_column_sizer.Hide(self.wellcol_selector)

        input_browse.Bind(wx.EVT_BUTTON, self.browse_input)
        self.bfx_format.Bind(wx.EVT_CHECKBOX, self.set_bfx_format)
        shapeb1.Bind(wx.EVT_RADIOBUTTON, self.set_shape)
        shapeb2.Bind(wx.EVT_RADIOBUTTON, self.set_shape)
        shapeb3.Bind(wx.EVT_RADIOBUTTON, self.set_shape)

        wells_combined.Bind(wx.EVT_RADIOBUTTON, self.set_wells_combined)
        wells_separate.Bind(wx.EVT_RADIOBUTTON, self.set_wells_combined)

        self.topsizer = sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(input_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(shape_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(plate_column_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(self.well_column_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(gene_column_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(status_sizer, 0, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(sizer)

        self.normalization.file_listeners.append(self.set_file)

    # XXX - handle text field editing
    def browse_input(self, evt):
        dlg = wx.FileDialog(self, "Choose an input file (.XLS)", wildcard="*.xls", style=wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.normalization.set_input_file(dlg.GetPath())
        dlg.Destroy()

    def set_shape(self, evt):
        self.normalization.shape = evt.EventObject.Label
        self.preflight()

    def set_plate_column(self, val):
        self.normalization.plate_column = val
        # fill well and gene selectors, as well
        for sel in [self.well_selector, self.wellrow_selector, self.wellcol_selector, self.gene_column_selector]:
            sel.set_default_sheet(val[0])
        self.preflight()

    def set_well_column(self, val):
        self.normalization.well_column = val
        self.preflight()

    def set_wellrow_column(self, val):
        self.normalization.wellrow_column = val
        self.preflight()

    def set_wellcol_column(self, val):
        self.normalization.wellcol_column = val
        self.preflight()

    def set_wells_combined(self, evt):
        if evt.EventObject == self.wells_combined:
            self.normalization.combined_wellrowcol = True
            self.well_column_sizer.Show(self.well_selector)
            self.well_column_sizer.Hide(self.wellrow_selector)
            self.well_column_sizer.Hide(self.wellcol_selector)
        else:
            self.normalization.combined_wellrowcol = False
            self.well_column_sizer.Hide(self.well_selector)
            self.well_column_sizer.Show(self.wellrow_selector)
            self.well_column_sizer.Show(self.wellcol_selector)
        self.preflight()
        self.topsizer.Layout()

    def set_gene_column(self, val):
        self.normalization.gene_column = val
        self.preflight()

    def set_file(self):
        # should check existence, possibly pre-parse
        self.input_text.Value = self.normalization.input_file
        self.TopLevelParent.update_title()
        bfx_row = self.normalization.get_row_values((0, 1))
        pat = re.compile('[0-9]*')
        self.bfx_format.Value = all(pat.match(unicode(s)) for s in bfx_row)
        self.normalization.set_bfx_format(self.bfx_format.Value)
        self.Refresh()

    def set_bfx_format(self, evt):
        self.normalization.set_bfx_format(self.bfx_format.Value)

    def preflight(self):
        # Attempt to parse the data, report what we find
        self.valid = False
        try:
            # guess if we are BFX format

            # fetch (plate, row, col, gene) for every entry in the XLS
            # given the columns we have.
            current_data = zip(self.normalization.fetch_plates(),
                               self.normalization.fetch_rows(),
                               self.normalization.fetch_cols(),
                               self.normalization.fetch_genes())

            self.detected_384 = False
            num_rows = 0
            gene_counts = {}
            wells_per_plate = {}
            for rowidx, (plate, row, col, gene) in enumerate(current_data):
                # XXX - should we allow empty genes?
                if plate == '':
                    # +2 because rowidx is 0-based, plus one for the head.
                    assert row == col == gene == '', "Incomplete data at row %d (plate = '%s', row = '%s', col = '%s', gene = '%s')" % (rowidx + 2, plate, row, col, gene)
                    continue
                if self.normalization.shape == '96':
                    assert row in 'ABCDEFGH', "Bad row for 96 well plate - '%s' at spreadsheet row %d'" % (row, rowidx + 1)
                    assert 1 <= int(col) <= 12, "Bad col for 96 well plate - '%s' at spreadsheet row %d'" % (col, rowidx + 1)
                else:
                    assert 'A' <= row <= 'P', "Bad row for 384 well plate - '%s' at spreadsheet row %d'" % (row, rowidx + 1)
                    assert 1 <= int(col) <= 24, "Bad col for 384 well plate - '%s' at spreadsheet row %d'" % (col, rowidx + 1)
                    if (row > 'H') or (int(col) > 12):
                        self.detected_384 = True
                num_rows += 1
                wells_per_plate[plate] = wells_per_plate.get(plate, 0) + 1
                gene_counts[gene] = gene_counts.get(gene, 0) + 1

            # XXX - check multiple plates for matching well/gene information
            # XXX - check uniqueness of wells (plate, row, col)
            # 
            # count number of plates, wells, and wells per gene
            self.valid = True
            plate_shape = self.normalization.shape
            autodetected = ""
            if plate_shape == DETECT:
                plate_shape = "384" if self.detected_384 else "96"
                autodetected = " (autodetected)"
                self.normalization.detected_384 = self.detected_384

            if min(wells_per_plate.values()) == max(wells_per_plate.values()):
                plate_counts_text = "Wells per plate: %d" % (wells_per_plate.values()[0])
            else:
                plate_counts_text = ["Variable number of wells per plate:"]
                plate_counts_text += sorted(["    %s : %d" % (p, v) for p, v in wells_per_plate.iteritems()])
                plate_counts_text = "\n".join(plate_counts_text)

            # report top 10 genes by count
            countgenes = sorted([(c, g) for g, c in gene_counts.iteritems()])[-10:][::-1]
            gene_counts_text = "\n".join(["Number of wells per gene, top 10:"] +
                                         ["%12s : %d" % (g.replace('\n', ' - '), c) for (c, g) in countgenes])

            status = "\n".join(["Plate Type: %s%s" % (plate_shape, autodetected),
                                "Number of plates: %d" % (len(wells_per_plate)),
                                plate_counts_text,
                                gene_counts_text])

            self.normalization.gene_counts = gene_counts
            self.status_text.Label = status
            self.topsizer.Layout()

            self.normalization.parsing_finished()

        except AssertionError, e:
            self.status_text.Label = "Parsing error:\n" + e.message
            self.valid = False
        except Exception, e:
            self.status_text.Label = "Choose settings above...\n" + str(e)
            self.valid = False

class Controls(wx.Panel):
    def __init__(self, parent, normalization):
        wx.Panel.__init__(self, parent=parent)
        self.normalization = normalization
        self.row_controls = []

        box = wx.StaticBox(self, wx.ID_ANY, 'Choose control populations')
        box_sizer = wx.StaticBoxSizer(box, wx.VERTICAL)

        self.scroll_window = wx.lib.scrolledpanel.ScrolledPanel(self, -1)
        box_sizer.Add(self.scroll_window, 1, wx.EXPAND)

        self.row_sizer = wx.BoxSizer(wx.VERTICAL)
        self.scroll_window.SetSizer(self.row_sizer)

        self.scroll_window.SetupScrolling(False, True)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(box_sizer, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Layout()

        self.normalization.parsing_listeners.append(self.update)

        self.Bind(wx.EVT_RADIOBUTTON, self.set_gene_control_type)

    class GeneControlButton(wx.RadioButton):
        def __init__(self, population, gene, *args, **kwargs):
            wx.RadioButton.__init__(self, *args, **kwargs)
            self.population = population
            self.gene = gene

    def set_gene_control_type(self, evt):
        self.normalization.set_control_type(evt.EventObject.gene, evt.EventObject.population)

    def update(self):
        # populate with genes, counts, radiobuttons
        self.row_controls = []
        self.row_sizer.DeleteWindows()

        def make_row(panel, g, c, tested, control):
            self.row_controls.append([g, c, tested, control])
            g.Wrap(150)
            sizer = wx.BoxSizer(wx.HORIZONTAL)
            sizer.Add(g, 1, wx.ALIGN_CENTER)
            sizer.Add(c, 1, wx.ALIGN_CENTER)
            sizer.Add(tested, 1, wx.ALIGN_CENTER)
            sizer.Add(control, 1, wx.ALIGN_CENTER)
            panel.SetSizer(sizer)
            return panel

        # sort by count, then name
        countgenes = sorted([(-c, g) for g, c in self.normalization.gene_counts.iteritems()])
        countgenes = [(-c, g) for c, g in countgenes]

        # truncate to no more than 50 genes
        countgenes = countgenes[:50]

        self.scroll_window.Freeze()
        panel = wx.Panel(self.scroll_window, -1)
        panel.BackgroundColour = "light blue"
        # XXX - this should be outside the scrolled area.
        self.row_sizer.Add(make_row(panel,
                                    wx.StaticText(panel, -1, "Name"),
                                    wx.StaticText(panel, -1, "Count"),
                                    wx.StaticText(panel, -1, "Tested Population"),
                                    wx.StaticText(panel, -1, "Control")),
                           0, wx.EXPAND)

        for idx, (count, gene) in enumerate(countgenes):
            panel = wx.Panel(self.scroll_window, -1)
            # panel.BackgroundColour = "white" if (idx % 5) else "light grey"
            self.row_sizer.Add(make_row(panel,
                                        wx.StaticText(panel, -1, gene.strip().replace('\n', ' - ')),
                                        wx.StaticText(panel, -1, "%d" % (count)),
                                        Controls.GeneControlButton(CONTROL_POPULATION, gene, panel, -1, style=wx.RB_GROUP),
                                        Controls.GeneControlButton(CONTROL_CONTROL, gene, panel, -1)),
                               0, wx.EXPAND)
            if idx % 5 == 4:
                self.row_sizer.Add(wx.StaticLine(self.scroll_window), 0, wx.EXPAND | wx.ALL, 1)

        for g, c, tested, control in self.row_controls[1:]:
            tested.Value = True

        self.scroll_window.Thaw()
        self.row_sizer.Layout()
        self.scroll_window.VirtualSize = self.scroll_window.BestVirtualSize

class Feature(wx.Panel):
    def __init__(self, parent, normalization):
        wx.Panel.__init__(self, parent=parent)
        self.normalization = normalization
        self.num_replicates = 1
        normalization.num_replicates = self.num_replicates

        feature_column_box = wx.StaticBox(self, wx.ID_ANY, 'Feature columns for each replicate')
        self.feature_column_sizer = wx.StaticBoxSizer(feature_column_box, wx.VERTICAL)
        feature_column_selector = ColumnSelector(self, self.set_feature_column, [], self.normalization, callback_args=(0,))
        self.feature_column_sizer.Add(feature_column_selector, 0, wx.EXPAND)
        self.feature_column_sizer.Add((1, 10), 1)

        add_replicate_button = wx.Button(self, label="Add Replicate")
        remove_replicate_button = wx.Button(self, label="Remove Last Replicate")
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        button_sizer.Add(add_replicate_button, 0)
        button_sizer.Add((10, 1), 0)
        button_sizer.Add(remove_replicate_button, 0)
        button_sizer.Add((1, 1), 2)

        add_replicate_button.Bind(wx.EVT_BUTTON, self.add_replicate)
        remove_replicate_button.Bind(wx.EVT_BUTTON, self.remove_replicate)

        self.feature_column_sizer.Add(button_sizer, 0, wx.EXPAND)

        self.topsizer = sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.feature_column_sizer, 0, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(sizer)

    def set_feature_column(self, val, replicate_index):
        self.normalization.set_replicate_feature(replicate_index, val)

    def add_replicate(self, evt):
        default_column = self.feature_column_sizer.Children[(self.num_replicates - 1) * 2].GetWindow().column_selector.StringSelection
        if default_column.endswith('%d' % self.num_replicates):
            w = len('%d' % self.num_replicates)
            hints = [default_column[:-w] + ('%d' % (self.num_replicates + 1))]
        else:
            hints = []
        feature_column_selector = ColumnSelector(self, self.set_feature_column, hints, self.normalization, callback_args=(self.num_replicates,))
        self.feature_column_sizer.Insert(self.num_replicates * 2, feature_column_selector, 0, wx.EXPAND)
        self.feature_column_sizer.Insert(self.num_replicates * 2 + 1, (1, 10), 0)
        default_sheet = self.feature_column_sizer.Children[(self.num_replicates - 1) * 2].GetWindow().sheet_idx
        self.num_replicates += 1
        self.normalization.num_replicates = self.num_replicates
        feature_column_selector.set_default_sheet(default_sheet)
        self.Layout()

    def remove_replicate(self, evt):
        if self.num_replicates == 1:
            return
        self.num_replicates -= 1
        self.normalization.num_replicates = self.num_replicates
        win = self.feature_column_sizer.Children[self.num_replicates * 2].GetWindow()
        self.feature_column_sizer.Detach(self.num_replicates * 2)
        self.feature_column_sizer.Detach(self.num_replicates * 2)
        win.Destroy()
        self.Layout()

class Parameters(wx.Panel):
    def __init__(self, parent, normalization, preview_callback):
        wx.Panel.__init__(self, parent=parent)
        self.normalization = normalization
        self.preview_callback = preview_callback

        iterations_box = wx.StaticBox(self, wx.ID_ANY, '')
        iterations_sizer = wx.StaticBoxSizer(iterations_box, wx.HORIZONTAL)
        iterations_control = wx.SpinCtrl(self, -1, value=str(self.normalization.num_iterations), min=0, max=20, initial=self.normalization.num_iterations)
        iterations_sizer.Add((1, 1), 1)
        iterations_sizer.Add(wx.StaticText(self, -1, "Number of iterations:"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 2)
        iterations_sizer.Add(iterations_control, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 2)
        iterations_sizer.Add((1, 1), 1)

        transform_box = wx.StaticBox(self, wx.ID_ANY, 'Data transformation')
        transform_sizer = wx.StaticBoxSizer(transform_box, wx.HORIZONTAL)
        transform_buttons = ([wx.RadioButton(self, -1, TRANSFORMATIONS[0], style=wx.RB_GROUP)]
                            + [wx.RadioButton(self, -1, label) for label in TRANSFORMATIONS[1:]])
        transform_sizer.Add((1, 1), 1)
        for b in transform_buttons:
            transform_sizer.Add(b, 0)
            transform_sizer.Add((1, 1), 1)
        transform_buttons[0].Value = 1

        plate_alignment_box = wx.StaticBox(self, wx.ID_ANY, 'Align plates?')
        plate_alignment_sizer = wx.StaticBoxSizer(plate_alignment_box, wx.VERTICAL)

        align_when_sizer = wx.BoxSizer(wx.HORIZONTAL)
        align_when_buttons = ([wx.RadioButton(self, -1, ALIGN_WHEN[0], style=wx.RB_GROUP)]
                                + [wx.RadioButton(self, -1, label) for label in ALIGN_WHEN[1:]])
        align_when_sizer.Add(wx.StaticText(self, -1, "When to align:"), 0)
        align_when_sizer.Add((1, 1), 1)
        for b in align_when_buttons:
            align_when_sizer.Add(b, 0)
            align_when_sizer.Add((1, 1), 1)
        align_when_buttons[0].Value = 1

        align_how_sizer = wx.BoxSizer(wx.HORIZONTAL)
        align_how_buttons = ([wx.RadioButton(self, -1, ALIGNMENT_METHODS[0], style=wx.RB_GROUP)]
                            + [wx.RadioButton(self, -1, label) for label in ALIGNMENT_METHODS[1:]])
        align_how_sizer.Add(wx.StaticLine(self), 0, wx.EXPAND | wx.ALL, 1)
        align_how_sizer.Add(wx.StaticText(self, -1, "How to align:"), 0)
        align_how_sizer.Add((1, 1), 1)
        for b in align_how_buttons:
            align_how_sizer.Add(b, 0)
            align_how_sizer.Add((1, 1), 1)
        align_how_buttons[0].Value = 1

        plate_alignment_sizer.Add(align_when_sizer, 0, wx.EXPAND | wx.ALL, 2)
        plate_alignment_sizer.Add(align_how_sizer, 0, wx.EXPAND | wx.ALL, 2)

        preview_save_box = wx.StaticBox(self, wx.ID_ANY, '')
        preview_save_sizer = wx.StaticBoxSizer(preview_save_box, wx.HORIZONTAL)
        preview_button = wx.Button(self, -1, 'Preview normalization')
        save_button = wx.Button(self, -1, 'Save normalized values')
        save_figures_button = wx.Button(self, -1, 'Save figures')
        preview_save_sizer.Add((1, 1), 1)
        preview_save_sizer.Add(preview_button, 0)
        preview_save_sizer.Add((1, 1), 1)
        preview_save_sizer.Add(save_button, 0)
        preview_save_sizer.Add((1, 1), 1)
        preview_save_sizer.Add(save_figures_button, 0)
        preview_save_sizer.Add((1, 1), 1)

        self.topsizer = sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(iterations_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(transform_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(plate_alignment_sizer, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(preview_save_sizer, 0, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(sizer)

        iterations_control.Bind(wx.EVT_SPINCTRL, self.update_iterations)
        iterations_control.Bind(wx.EVT_TEXT, self.update_iterations)

        for b in transform_buttons:
            b.Bind(wx.EVT_RADIOBUTTON, self.update_transform)

        for b in align_when_buttons:
            b.Bind(wx.EVT_RADIOBUTTON, self.update_alignment_when)

        for b in align_how_buttons:
            b.Bind(wx.EVT_RADIOBUTTON, self.update_alignment_method)

        # XXX - preview/save buttons should be inactive until the normalization is ready
        preview_button.Bind(wx.EVT_BUTTON, self.preview)
        save_button.Bind(wx.EVT_BUTTON, self.save)
        save_figures_button.Bind(wx.EVT_BUTTON, self.save_figures)

        self.default_dir = '.'
        self.default_file = 'normalized.csv'
        self.normalization.file_listeners.append(self.set_default_output)

    def update_iterations(self, evt):
        try:
            int(evt.EventObject.Value)
        except:
            return
        self.normalization.set_iterations(int(evt.EventObject.Value))

    def update_transform(self, evt):
        self.normalization.set_transformation(evt.EventObject.Label)

    def update_alignment_when(self, evt):
        self.normalization.set_alignment_when(evt.EventObject.Label)

    def update_alignment_method(self, evt):
        self.normalization.set_alignment_method(evt.EventObject.Label)

    def preview(self, evt):
        self.preview_callback()

    def save(self, evt):
        dlg = wx.FileDialog(self, "Choose output file", wildcard="*.xls", style=wx.FD_SAVE,
                            defaultDir=self.default_dir, defaultFile=self.default_file)
        output_file = None
        if dlg.ShowModal() == wx.ID_OK:
            output_file = dlg.GetPath()
        dlg.Destroy()
        if not output_file:
            return
        if os.path.exists(output_file):
            dlg = wx.MessageDialog(self, "Will not overwrite existing file %s" % (output_file), "Warning", wx.OK|wx.ICON_WARNING)
            dlg.ShowModal()
            dlg.Destroy()
            return self.save(evt)
        else:
            self.normalization.set_output_file(output_file)
            self.normalization.save()

    def save_figures(self, evt):
        dlg = wx.FileDialog(self, "Choose figures file", wildcard="*.pdf", style=wx.FD_SAVE,
                            defaultDir=self.default_dir, defaultFile=self.default_file.replace('xls', 'pdf'))
        plots_output_file = None
        if dlg.ShowModal() == wx.ID_OK:
            plots_output_file = dlg.GetPath()
        dlg.Destroy()
        if not plots_output_file:
            return
        if os.path.exists(plots_output_file):
            dlg = wx.MessageDialog(self, "Will not overwrite existing file %s" % (plots_output_file), "Warning", 
                                   wx.OK | wx.ICON_WARNING)
            dlg.ShowModal()
            dlg.Destroy()
            return self.save_figures(evt)
        else:
            self.preview_callback()
            self.TopLevelParent.plots.save_plots(plots_output_file)

    def set_default_output(self):
        self.default_dir = os.path.dirname(self.normalization.input_file)
        self.default_file = wrap_filename(os.path.basename(self.normalization.input_file), 'normalized')

class Plot(wxplotpanel.PlotPanel):
    ''' shared superclass with common __init__ '''
    def __init__(self, parent, normalization, **kwargs):
        self.normalization = normalization
        self.parent = parent
        # initiate plot window
        wxplotpanel.PlotPanel.__init__(self, parent, **kwargs)

    def draw(self):
        self.figure.clear()
        if not self.normalization.ready():
            subplot = self.figure.add_subplot(111)
            subplot.annotate('waiting\nfor settings...', (0, 0),
                             horizontalalignment='center',
                             multialignment='center')
            subplot.axis([-1, 1, -1, 1])
        else:
            self.do_draw()

class PlatePlot(Plot):
    def do_draw(self):
        self.pre_draw()
        plotidx = 0
        bad_data = False
        self.grid = grid = self.image_grid(self.get_num_plates(), self.get_num_replicates())
        lo = np.inf
        hi = -np.inf
        for plate_index in range(self.get_num_plates()):
            for rep in range(self.get_num_replicates()):
                try:
                    im = self.get_plate(plate_index, rep)
                    im = im[np.isfinite(im)]
                    lo = min(np.min(im), lo)
                    hi = max(np.max(im), hi)
                except:
                    pass
        norm = self.get_norm(lo, hi)
        xticklocs = range(0, 24, 2) if (self.normalization.plate_shape() == '384') else range(0, 12)
        yticklocs = range(0, 16, 2) if (self.normalization.plate_shape() == '384') else range(0, 8)
        xticklabels = ["%02d" % (i+ 1) for i in xticklocs]
        yticklabels = ['ABCDEFGHIJKLMNOP'[i] for i in yticklocs]
        for plate_index in range(self.get_num_plates()):
            for rep in range(self.get_num_replicates()):
                try:
                    im = self.get_plate(plate_index, rep)
                    if np.any(~ np.isfinite(im)):
                        bad_data = True
                    mappable = grid[plotidx].imshow(im, interpolation='nearest', norm=norm, origin='upper')
                    grid[plotidx].set_xticks(xticklocs)
                    grid[plotidx].set_yticks(yticklocs)
                    grid[plotidx].set_xticklabels(xticklabels, size='xx-small')
                    grid[plotidx].set_yticklabels(yticklabels, size='x-small')
                    plotidx += 1
                except:
                    print "Exception", plate_index, rep
                    traceback.print_exc()
                    pass
        grid[0].set_xlim(-0.5, 23.5 if (self.normalization.plate_shape() == '384') else 11.5)
        grid[0].set_ylim(15.5 if (self.normalization.plate_shape() == '384') else 7.5, -0.5)
        grid[0].cax.colorbar(mappable)
        grid[0].cax.toggle_label(True)
        self.post_draw(bad_data)

    def pre_draw(self):
        pass
    def get_plate(self, plate_index, repindex):
        raise NotImplementedError('need get_plate()')
    def post_draw(self, bad_data):
        pass

    def get_num_replicates(self):
        # can be overriden for plots that want to show only one column
        return self.normalization.num_replicates

    def get_num_plates(self):
        # can be overriden for plots that want to show combined plates
        return self.normalization.num_plates()

class OriginalPlates(PlatePlot):
    def get_plate(self, plate_index, rep):
        return self.normalization.plate_array(plate_index, rep)

    def post_draw(self, bad_data):
        self.figure.suptitle('original%s' % (' (invalid values discarded)' if bad_data else ''))


class Histograms(Plot):
    def do_draw(self):
        self.pre_draw()
        bad_data = False

        # find control status, and unique groups (by treatment) in controls
        control_groups, control_names = self.normalization.get_transformed_control_groups()
        num_control_groups = max(control_groups) + 1

        for rep in range(self.normalization.num_replicates):
            subplot = self.figure.add_subplot(self.normalization.num_replicates, 1, rep + 1)
            if self.original:
                data = self.normalization.get_orig_values(rep)
            else:
                data = self.normalization.get_transformed_values(rep, cleaned=self.cleaned)
            good_mask = np.isfinite(data)
            bad_data = bad_data or np.any(~ good_mask)
            # black first = tested population
            colors = list('kgrbcmy')
            for control_idx in range(num_control_groups):
                ctrl_good_mask = (control_groups == control_idx) & good_mask
                if not np.any(ctrl_good_mask):
                    continue
                subplot.hist(data[ctrl_good_mask], 20,
                             color=(colors.pop(0) if colors else None),
                             alpha=0.5,
                             label=control_names[control_idx] or '---')
        self.align_subplots()

        # add legend above first plot
        legend = subplot.legend(loc='upper center', scatterpoints=1,
                                ncol=4, mode="expand", frameon=False,
                                bbox_to_anchor=(0, 0, 1, 1), bbox_transform=self.figure.transFigure,
                                title='%s%s' % (self.title,
                                       ' (invalid values discarded)' if bad_data else ''))

        # set legend font size
        for t in legend.get_texts():
            t.set_fontsize('small')

    def pre_draw(self):
        self.normalization.run_normalization()

class OriginalHistograms(Histograms):
    title = 'original'
    cleaned = False
    original = True

class TransformedHistograms(Histograms):
    title = 'transformed'
    original = False
    cleaned = False

class CleanedTransformedHistograms(Histograms):
    title = 'cleaned transformed'
    original = False
    cleaned = True

class TransformedPlates(PlatePlot):
    def get_plate(self, plate_index, rep):
        return self.normalization.plate_array(plate_index, rep, transformed=True)

    def post_draw(self, bad_data):
        self.figure.suptitle('transformed%s' % (' (invalid values discarded)' if bad_data else ''))

class Agreement(Plot):
    def do_draw(self):
        if self.normalization.num_replicates < 2:
            return
        self.pre_draw()
        nreps = self.normalization.num_replicates
        npairs = (nreps * (nreps - 1)) / 2
        bad_data = False
        pairidx = 1

        # find control status, and unique groups (by treatment) in controls
        control_groups, control_names = self.normalization.get_transformed_control_groups()
        num_control_groups = max(control_groups) + 1

        lo = np.inf
        hi = -np.inf

        subplots = []

        for rep_a in range(self.normalization.num_replicates):
            data_a = self.normalization.get_transformed_values(rep_a, cleaned=self.cleaned)
            for rep_b in range(rep_a):
                subplot = self.figure.add_subplot(npairs, 1, pairidx, aspect='equal')
                subplots.append(subplot)
                pairidx += 1
                data_b = self.normalization.get_transformed_values(rep_b, cleaned=self.cleaned)
                good_mask = np.isfinite(data_a + data_b)
                bad_data = np.any(~ good_mask)
                if any(good_mask):
                    lo = min(min(data_a[good_mask]), min(data_b[good_mask]), lo)
                    hi = max(max(data_a[good_mask]), max(data_b[good_mask]), hi)
                    # black first = tested population
                    colors = list('kgrbcmy')
                    for control_idx in range(num_control_groups):
                        ctrl_good_mask = (control_groups == control_idx) & good_mask
                        if not np.any(ctrl_good_mask):
                            continue
                        subplot.scatter(data_b[ctrl_good_mask],
                                        data_a[ctrl_good_mask],
                                        s=10,
                                        marker='o',
                                        color=(colors.pop(0) if colors else None),
                                        zorder=-control_idx,  # tested population on top
                                        label=control_names[control_idx] or '---',
                                        alpha=0.3)
                    subplot.set_xlabel('replicate %d' % (rep_b + 1))
                    subplot.set_ylabel('replicate %d' % (rep_a + 1))
                    if pairidx <= npairs:  # turn off xtick labels for all but last plot
                        subplot.set_xticklabels([])
                        subplot.get_xaxis().labelpad = 0

        # add legend above first plot
        legend = subplot.legend(loc='upper center', scatterpoints=1,
                                ncol=4, mode="expand", frameon=False,
                                bbox_to_anchor=(0, 0, 1, 1), bbox_transform=self.figure.transFigure,
                                title='%stransformed%s' % ('cleaned ' if self.cleaned else '',
                                                           ' (invalid values discarded)' if bad_data else ''))

        # set legend font size
        for t in legend.get_texts():
            t.set_fontsize('small')

        # give some space around data
        border = (hi - lo) / 20
        lo -= border
        hi += border
        # make every subplot have identical limits
        for subplot in subplots:
            subplot.axis([lo, hi, lo, hi])

    def pre_draw(self):
        self.normalization.run_normalization()

class TransformedAgreement(Agreement):
    cleaned = False

class CleanedAgreement(Agreement):
    cleaned = True

class AlignedPlates(PlatePlot):
    def do_draw(self):
        if self.normalization.align_when != ALIGN_NEVER:
            PlatePlot.do_draw(self)

    def pre_draw(self):
        self.normalization.run_normalization()

    def get_plate(self, plate_index, rep):
        plate_name = self.normalization.plate_names()[plate_index]
        return self.normalization.normalization_first_alignment[plate_name, rep]

    def post_draw(self, bad_data):
        self.figure.suptitle('transformed and aligned%s' % (' (invalid values discarded)' if bad_data else ''))

class ZScorePlates(PlatePlot):
    def __init__(self, *args, **kwargs):
        self.before_normalization = kwargs.pop('before_normalization', False)
        PlatePlot.__init__(self, *args, **kwargs)

    def do_draw(self):
        PlatePlot.do_draw(self)

    def pre_draw(self):
        self.normalization.run_normalization()
        # compute median of non-controls

    def get_plate(self, plate_index, rep):
        plate_name = self.normalization.plate_names()[plate_index]
        if self.before_normalization:
            return np.clip(((self.normalization.normalization_first_alignment[plate_name, rep] -
                             self.normalization.first_alignment_medians[rep]) /
                            self.normalization.first_alignment_MAD_sigmas[rep]),
                           -4, 4)
        else:
            return np.clip(((self.normalization.normalization_plate_values[plate_name, rep] -
                             self.normalization.normalization_medians[rep]) /
                            self.normalization.normalization_MAD_sigmas[rep]),
                           -4, 4)

    def post_draw(self, bad_data):
        if self.before_normalization:
            self.figure.suptitle('transformed & aligned robust Z-scores\n(clamped at -/+ 4)%s' % (' (invalid values discarded)' if bad_data else ''))
        else:
            self.figure.suptitle('cleaned robust Z-scores\n(clamped at -/+ 4)%s' % (' (invalid values discarded)' if bad_data else ''))

class MergedPlates(PlatePlot):
    def __init__(self, when, *args, **kwargs):
        assert when in ['before', 'after']
        self.when = when
        PlatePlot.__init__(self, *args, **kwargs)

    def do_draw(self):
        if self.normalization.num_replicates > 1:
            PlatePlot.do_draw(self)

    def pre_draw(self):
        self.normalization.run_normalization()

    def get_num_replicates(self):
        return 1

    def get_plate(self, plate_index, rep):
        plate_name = self.normalization.plate_names()[plate_index]
        if self.when == 'before':
            vals = [self.normalization.normalization_first_alignment[plate_name, r] for r in range(self.normalization.num_replicates)]
        else:
            vals = [self.normalization.normalization_plate_values[plate_name, r] for r in range(self.normalization.num_replicates)]
        return np.median(np.dstack(vals), 2)

    def post_draw(self, bad_data):
        self.figure.suptitle('Merged %s' % (self.when))

class MergedInReplicatesPlates(PlatePlot):
    def __init__(self, when, *args, **kwargs):
        assert when in ['before', 'after']
        self.when = when
        PlatePlot.__init__(self, *args, **kwargs)

    def pre_draw(self):
        self.normalization.run_normalization()

    def get_num_replicates(self):
        return 1

    def get_num_plates(self):
        return self.normalization.num_replicates

    def get_plate(self, plate_index, rep):
        # plate_index is actually rep index
        if self.when == 'before':
            vals = [v for (pl, r), v in self.normalization.normalization_first_alignment.iteritems() if r == plate_index]
        else:
            vals = [v for (pl, r), v in self.normalization.normalization_plate_values.iteritems() if r == plate_index]
        return np.median(np.dstack(vals), 2)

    def post_draw(self, bad_data):
        self.figure.suptitle('Merged %s' % (self.when))

class MedianZScorePlate(PlatePlot):
    def pre_draw(self):
        self.normalization.run_normalization()

    def get_num_replicates(self):
        return 1

    def get_num_plates(self):
        return self.normalization.num_replicates + 1

    def get_plate(self, plate_index, rep):
        if plate_index < self.normalization.num_replicates:
            # fetch vals
            vals = np.dstack([v for (pl, r), v in self.normalization.normalization_plate_values.iteritems() if r == plate_index])
            # convert to absolute-robust-Z-scores
            vals = vals - np.median(vals.flatten())
            vals = abs(vals / (1.4826 * np.median(abs(vals.flatten()))))
            return np.median(vals, 2)
        return np.median(np.dstack([self.get_plate(idx, 0) for idx in range(self.normalization.num_replicates)]), 2)

    def post_draw(self, bad_data):
        self.figure.suptitle('Median-abs-z-score (last plate is median of MAZS)')

class CleanedPlates(PlatePlot):
    def pre_draw(self):
        self.normalization.run_normalization()

    def get_plate(self, plate_index, rep):
        plate_name = self.normalization.plate_names()[plate_index]
        return self.normalization.normalization_plate_values[plate_name, rep]

    def post_draw(self, bad_data):
        self.figure.suptitle('cleaned transformed %s' % (' (invalid values discarded)' if bad_data else ''))

class ControlPlatemap(PlatePlot):
    def do_draw(self):
        PlatePlot.do_draw(self)

        # find control groups
        control_indices, control_names = self.normalization.fetch_control_groups()

        # get the colormap, and colors for each value
        cmap = self.grid[0].get_images()[0].get_cmap()
        rgb_values = cmap(np.linspace(0, 1, len(control_names)))
        # create proxies
        proxies = [Rectangle((0, 0), 1, 1, fc=rgb) for rgb in rgb_values]
        control_names = [c.replace('\n', '') for c in control_names]
        # add legend above first plot
        legend = self.grid[0].legend(proxies, control_names, 
                                     loc='upper center', scatterpoints=1,
                                     ncol=2, mode="expand", frameon=False,
                                     bbox_to_anchor=(0, 0, 1, 1), bbox_transform=self.figure.transFigure,
                                     title='Controls')

        # set legend font size
        for t in legend.get_texts():
            t.set_fontsize('small')

    def pre_draw(self):
        self.normalization.run_normalization()

    def get_num_replicates(self):
        return 1

    def get_plate(self, plate_index, rep):
        return self.normalization.normalization_control_groups[self.normalization.plate_names()[plate_index]]


class Plots(wx.Panel):
    def __init__(self, parent, normalization):
        wx.Panel.__init__(self, parent=parent)

        self.normalization = normalization

        # XXX - add control layout map

        # XXX - add old vs. new scatterplot

        self.scroll_window = wx.lib.scrolledpanel.ScrolledPanel(self, -1)
        self.subpanel = subpanel = wx.Panel(self.scroll_window, -1)
        self.panels = OrderedDict()
        self.panels['original data'] = OriginalHistograms(subpanel, normalization)
        self.panels['original platemaps'] = OriginalPlates(subpanel, normalization)
        self.panels['transformed data'] = TransformedHistograms(subpanel, normalization)
        self.panels['transformed platemaps'] = TransformedPlates(subpanel, normalization)
        self.panels['transformed agreement'] = TransformedAgreement(subpanel, normalization)
        self.panels['aligned platemaps'] = AlignedPlates(subpanel, normalization)
        self.panels['aligned zscore platemaps'] = ZScorePlates(subpanel, normalization, before_normalization=True)
        self.panels['merged before'] = MergedPlates('before', subpanel, normalization)
        self.panels['merged before inrep'] = MergedInReplicatesPlates('before', subpanel, normalization)
        self.panels['merged after inrep'] = MergedInReplicatesPlates('after', subpanel, normalization)
        self.panels['merged after'] = MergedPlates('after', subpanel, normalization)
        self.panels['medabsz'] = MedianZScorePlate( subpanel, normalization)
        self.panels['cleaned agreement'] = CleanedAgreement(subpanel, normalization)
        self.panels['cleaned plates'] = CleanedPlates(subpanel, normalization)
        self.panels['cleaned data'] = CleanedTransformedHistograms(subpanel, normalization)
        self.panels['cleaned zscore'] = ZScorePlates(subpanel, normalization, before_normalization=False)
        self.panels['controls'] = ControlPlatemap(subpanel, normalization)

        sizer = self.panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.panels['original data'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['original platemaps'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['transformed data'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['transformed platemaps'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['transformed agreement'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['aligned platemaps'], 1, wx.ALL | wx.EXPAND, 1, )
        sizer.Add(self.panels['aligned zscore platemaps'], 1, wx.ALL | wx.EXPAND, 1, )
        sizer.Add(self.panels['merged before'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['merged before inrep'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['merged after inrep'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['merged after'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['medabsz'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['cleaned agreement'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['cleaned plates'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['cleaned data'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['cleaned zscore'], 1, wx.ALL | wx.EXPAND, 1)
        sizer.Add(self.panels['controls'], 1, wx.ALL | wx.EXPAND, 1)
        subpanel.SetSizer(sizer)

        self.scroll_window.SetupScrolling(True, False)

        top_sizer = wx.BoxSizer(wx.VERTICAL)
        top_sizer.Add(self.scroll_window, 1, wx.EXPAND)
        self.SetSizer(top_sizer)
        self.Layout()

        self.Bind(wx.EVT_SIZE, self.on_size)

        normalization.parameter_change_listeners.append(self.update_parameters)

        # force an update of visible plots
        self.update_parameters()

    def update_plots(self):
        # XXX - should have some way to mark plots dirty
        for p in self.panels.values():
            p.draw()
            p.canvas.draw()
        self.Refresh()

    def update_parameters(self):
        self.panel_sizer.Show(self.panels['aligned platemaps'], self.normalization.align_when != ALIGN_NEVER)
        self.panel_sizer.Show(self.panels['merged before'], self.normalization.num_replicates > 1)
        self.panel_sizer.Show(self.panels['merged after'], self.normalization.num_replicates > 1)
        self.panel_sizer.Show(self.panels['transformed agreement'], self.normalization.num_replicates > 1)
        self.panel_sizer.Show(self.panels['cleaned agreement'], self.normalization.num_replicates > 1)
        self.on_size(None)

    def on_size(self, evt):
        self.Layout()  # allows the plots to draw when the window first appears
        height = self.scroll_window.ClientSize[1]
        width = int(height / np.sqrt(2))
        num_visible = len([p for p in self.panels.values() if p.IsShown()])
        self.subpanel.Size = (width * num_visible, height)
        self.scroll_window.VirtualSize = self.subpanel.Size

    def save_plots(self, filename):
        # XXX - progressbar
        pdfpages = wxplotpanel.start_pdf(filename)
        for p in self.panels.values():
            p.save_to_pdf(pdfpages)
        wxplotpanel.end_pdf(pdfpages)

class Frame(wx.Frame):
    def __init__(self, title, normalization):
        wx.Frame.__init__(self, None, title=title, size=(600, 600))
        self.normalization = normalization
        self.appname = title

        menuBar = wx.MenuBar()
        menu = wx.Menu()
        m_exit = menu.Append(wx.ID_EXIT, "E&xit\tAlt-X", "Close window and exit program.")
        self.Bind(wx.EVT_MENU, self.on_close, m_exit)
        menuBar.Append(menu, "&File")
        menu = wx.Menu()
        m_about = menu.Append(wx.ID_ABOUT, "&About", "Information about this program")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        menuBar.Append(menu, "&Help")
        self.SetMenuBar(menuBar)

        self.statusbar = self.CreateStatusBar()

        panel = wx.Panel(self)

        notebook = self.notebook = wx.Notebook(panel)

        self.plots = Plots(notebook, self.normalization)
        notebook.AddPage(PlateLayout(notebook, self.normalization), "File && Layout")
        notebook.AddPage(Controls(notebook, self.normalization), "Controls")
        notebook.AddPage(Feature(notebook, self.normalization), "Feature")
        notebook.AddPage(Parameters(notebook, self.normalization, self.plots.update_plots), "Parameters")
        notebook.AddPage(self.plots, "Plots")

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL | wx.EXPAND, 5)
        panel.SetSizer(sizer)

        self.Layout()
        self.Bind(wx.EVT_CLOSE, self.on_close)

    def on_close(self, event):
        dlg = wx.MessageDialog(self,
            "Do you really want to close this application?",
            "Confirm Exit", wx.OK | wx.CANCEL | wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK:
            self.Destroy()

    def on_about(self, event):
        dlg = AboutBox()
        dlg.ShowModal()
        dlg.Destroy()

    def update_title(self):
        title = self.appname
        if self.normalization.input_file != '':
            title += ' - %s' % (os.path.basename(self.normalization.input_file))
            if self.normalization.output_file != '':
                title += ' -> %s' % (os.path.basename(self.normalization.output_file))
        self.Title = title



def main():
    normalization = Normalization()
    app = wx.App(redirect=False)
    top = Frame(app_name, normalization)
    top.Centre()
    top.Show()
    if len(sys.argv) > 1:
        normalization.set_input_file(sys.argv[1])

    def show_errordialog(type, exc, tb):
        def doit():
            display_error_dialog(top, exc, None, tb=tb, continue_only=True,
                                 message="Exception in processing")
        # continue is really the only choice
        wx.CallAfter(doit)
    # replace default hook with error dialog
    orig_excepthook = sys.excepthook
    sys.excepthook = show_errordialog


    app.MainLoop()

if __name__ == "__main__":
    main()
