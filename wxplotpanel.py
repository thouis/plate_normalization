# from: http://www.scipy.org/Matplotlib_figure_in_a_wx_panel
import matplotlib
matplotlib.interactive(True)
matplotlib.use('WXAgg')
from mpl_toolkits.axes_grid1 import ImageGrid
import wx
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure, SubplotParams

class PlotPanel (wx.Panel):
    """The PlotPanel has a Figure and a Canvas"""
    def __init__(self, parent, color=(255, 255, 255), dpi=None, **kwargs):
        # initialize Panel
        if 'id' not in kwargs.keys():
            kwargs['id'] = wx.ID_ANY
        if 'style' not in kwargs.keys():
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__(self, parent, **kwargs)

        # subplotparams = SubplotParams(0.02, 0.02, 0.98, 0.98, 0.1, 0.1)
        # initialize matplotlib stuff
        self.figure = Figure((2.1, 2.97), dpi)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.SetColor(color)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND | wx.SHAPED, 0)
        self.SetSizer(sizer)

        self.Bind(wx.EVT_SIZE, self.set_size)
        self.draw()
        self.Refresh()

    def SetColor(self, rgbtuple=None):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE).Get()
        clr = [c/255. for c in rgbtuple]
        self.figure.set_facecolor(clr)
        self.figure.set_edgecolor(clr)
        self.canvas.SetBackgroundColour(wx.Colour(*rgbtuple))

    def set_size(self, evt=None):
        if self.ClientSize[0] > 0 and self.ClientSize[1] > 0:
            self.canvas.SetSize(self.ClientSize)
            self.canvas.draw()

    def draw(self):
        raise NoImplementedError # abstract, to be overridden by child classes

    def image_grid(self, num_rows, num_cols):
        return ImageGrid(self.figure, 111,
                         share_all = True,
                         nrows_ncols = (num_rows, num_cols),
                         cbar_size = "3%",
                         cbar_pad = 0.02,
                         cbar_mode = 'single')

    def get_norm(self, vmin, vmax):
        return matplotlib.colors.normalize(vmax=vmax, vmin=vmin)

    def save_to_pdf(self, pdfpages):
        old_fig = self.figure
        self.figure = Figure((8.5, 11), dpi=300)
        canvas = matplotlib.backends.backend_pdf.FigureCanvasPdf(self.figure)
        self.draw()
        pdfpages.savefig(self.figure)
        self.figure = old_fig

    def align_subplots(self):
        xmin = matplotlib.numpy.inf
        xmax = - matplotlib.numpy.inf
        for subplot in self.figure.get_axes():
            xmin = min(xmin, subplot.get_xlim()[0])
            xmax = max(xmax, subplot.get_xlim()[1])
        for subplot in self.figure.get_axes():
            subplot.set_xlim(xmin, xmax)



def start_pdf(filename):
    return PdfPages(filename)

def end_pdf(pdfpages):
    pdfpages.close()
