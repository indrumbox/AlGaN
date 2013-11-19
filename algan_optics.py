import wx
from wx.lib.stattext import GenStaticText 
import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
matplotlib.rcParams['legend.fontsize'] = 10
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import math
import AbsorptionCoefficient
import pylab

def calculatePolinom(a, b, hv):
    polinom_states = [1 for i in range(9)]
    for i in range(9):
        polinom_states[i] = polinom_states[i] * b[i] * ((hv - a[i]) ** 8)
        for j in range(9):
            if j != i:
                polinom_states[i] = polinom_states[i] / (a[j] - a[i])

    polinom = 0
    for i in range(9):
        polinom = polinom + polinom_states[i]
    return polinom

class AlGaN():
    def __init__(self, parent, d, x):
        self.d = d
        self.x = x
        self.setEg(x)

    def setEg(self, x):
        Eg_AlN = 6.13
        Eg_GaN = 3.42
        b = 1.08 # MBE material
        self.Eg = self.x * Eg_AlN + (1 - self.x) * Eg_GaN - b * self.x * (1 - self.x)
        

    def absorptionCoeff(self, Energy):
        alpha0 = 1e5 # cm-1*eV-1
        return alpha0 * math.sqrt(Energy - self.Eg)

    def refractiveIndex(self, wavelen, kind='o'):
        Eg_ = self.Eg
        B0k = [6.626, 7.042]
        B1k = [-0.934, -1.054]
        B2k = [0.0598, 0.0733]
        C0k = [396.8, 381.2]
        C1k = [-84.12, -76.68]
        C2k = [6.758, 6.068]
    
        if kind == 'o':
            B0 = B0k[0]
            B1 = B1k[0]
            B2 = B2k[0]
            C0 = C0k[0]
            C1 = C1k[0]
            C2 = C2k[0]
        else:
            B0 = B0k[1]
            B1 = B1k[1]
            B2 = B2k[1]
            C0 = C0k[1]
            C1 = C1k[1]
            C2 = C2k[1]
        
        A0 = B0 + B1 * Eg_ + B2 * Eg_ * Eg_
        lam0 = C0 + C1 * Eg_ + C2 * Eg_ * Eg_

        n_2 = 1.0 + A0 * (wavelen ** 2) / ((wavelen ** 2) - (lam0 ** 2))
        n = math.sqrt(abs(n_2))
        return n

    def transmissionBasis(self, wavelen, d, x):
        '''
        1.366309.pdf
        '''
        absorp = AbsorptionCoefficient.AbsorpKoeff(x, wavelen*1E-9)
        n_AlGaN = self.refractiveIndex(wavelen)
        n_air = 1
        n_Al2O3 = 1.77
        R1 = (n_air - n_AlGaN) / (n_air + n_AlGaN)
        R2 = (n_AlGaN - n_Al2O3) / (n_AlGaN + n_Al2O3)
        R3 = (n_Al2O3 - n_air) / (n_Al2O3 + n_air)
        R1 = abs(R1) ** 2
        R2 = abs(R2) ** 2
        R3 = abs(R3) ** 2
        Xi = R1 * R2 + R1 * R3 - 2 * R1 * R2 * R3
        A = - (1 - R1) * (1 - R2) * (1 + R3)

        T = A * math.exp(-absorp * d)
        zn = Xi * (math.exp(-absorp * d) ** 2) + R2 * R3 - 1
        T = T * 1.0 / zn
        return T


class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        global Sample
        super(MainWindow, self).__init__(parent, title=title, size=(750, 720))

        self.SetTitle('AlGaN optics')
        self.InitUI()

        Sample = AlGaN(None, 200, 0.45)

        self.editD.SetValue(str(Sample.d))
        self.editX.SetValue(str(Sample.x))
        self.drawGraphs()
        self.Centre()
        self.Show()

    def InitUI(self):
        panel = wx.Panel(self)
        sizer = wx.GridBagSizer(3, 3)
        
        # creating widgets on left panel
        dataPanel = wx.Panel(panel)
        dataSizer = wx.FlexGridSizer(3, 2, 2, 2)

        self.labelThickness = GenStaticText(dataPanel, label="Thickness of layer, nm:", style=wx.ALIGN_RIGHT)
        self.labelComposition = GenStaticText(dataPanel, label="Composition of layer, m.f.:", style=wx.ALIGN_RIGHT)
        self.editD = wx.TextCtrl(dataPanel, size=(90, -1), style=wx.TE_PROCESS_ENTER)
        self.editX = wx.TextCtrl(dataPanel, size=(90, -1), style=wx.TE_PROCESS_ENTER)
        
        self.buttonCalculate = wx.Button(dataPanel, label="Calculate", size=(90, 28))
        self.Bind(wx.EVT_BUTTON, self.buttonCalculatePress, self.buttonCalculate)

        dataSizer.Add(self.labelThickness, 1, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT)
        dataSizer.Add(self.editD, 1, wx.EXPAND | wx.ALL, 2)
        dataSizer.Add(self.labelComposition, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT)
        dataSizer.Add(self.editX, 1, wx.EXPAND | wx.ALL, 2)
        dataSizer.Add(wx.StaticText(dataPanel, label=''), wx.EXPAND) # space in sizer
        dataSizer.Add(self.buttonCalculate, 1, wx.EXPAND | wx.ALL, 3)

        dataPanel.SetSizer(dataSizer)

        
        sizer.Add(dataPanel, pos=(0, 0), flag=wx.EXPAND)
        # end working with left panel of controls

        graphsPanel = wx.Panel(panel)
        graphsPanel.SetBackgroundColour(wx.Colour(230,233,234))

        graphsSizer = wx.BoxSizer(wx.VERTICAL)

        self.dpi = 100
        
        self.figAbs = Figure((5.0, 2.0), dpi=self.dpi, tight_layout=True)
        self.canvas1 = FigureCanvas(graphsPanel, -1, self.figAbs)
        self.axesAbs = self.figAbs.add_subplot(111)
        graphsSizer.Add(self.canvas1, 1, wx.EXPAND | wx.ALL, 1)

        self.figRefr = Figure((5.0, 2.0), dpi=self.dpi, tight_layout=True)
        self.canvas2 = FigureCanvas(graphsPanel, -1, self.figRefr)
        self.axesRefr = self.figRefr.add_subplot(111)
        graphsSizer.Add(self.canvas2, 1, wx.EXPAND | wx.ALL, 1)

        self.figTrans = Figure((5.0, 2.0), dpi=self.dpi, tight_layout=True)
        self.canvas3 = FigureCanvas(graphsPanel, -1, self.figTrans)
        self.axesTrans = self.figTrans.add_subplot(111)
        graphsSizer.Add(self.canvas3, 1, wx.EXPAND | wx.ALL, 1)

        sizer.Add(graphsPanel, pos=(0, 1), span=(3, 2), flag=wx.EXPAND)
        sizer.AddGrowableCol(1)
        sizer.AddGrowableRow(1)

        # apply sizers        
        graphsPanel.SetSizerAndFit(graphsSizer)
        panel.SetSizerAndFit(sizer)

    def drawGraphs(self):
        # draws graphs in axesAbs, axesRefr, axesTrans
        global Sample
        Sample.d = float(self.editD.GetValue())*1e-9 # editD - in nanometers
        Sample.x = float(self.editX.GetValue())
        Sample.setEg(Sample.x)

        self.axesAbs.clear()        
        self.axesAbs.semilogy()  #logarithmic values of y axis
        self.axesAbs.set_title("Absorption Coefficient of AlGaN", fontsize=14)
        self.axesAbs.set_xlabel('wavelength, nm', fontsize=10)
        self.axesAbs.set_ylabel('absorption coefficient, sm-1', fontsize=10)
        self.axesAbs.grid(True)
        pylab.setp(self.axesAbs.get_xticklabels(), fontsize=10)
        pylab.setp(self.axesAbs.get_yticklabels(), fontsize=10)

        self.axesRefr.clear()        
        self.axesRefr.set_title("Refractive Index of AlGaN", fontsize=14)
        self.axesRefr.set_xlabel('wavelength, nm', fontsize=10)
        self.axesRefr.set_ylabel('refractive index', fontsize=10)
        self.axesRefr.grid(True)
        pylab.setp(self.axesRefr.get_xticklabels(), fontsize=10)
        pylab.setp(self.axesRefr.get_yticklabels(), fontsize=10)

        self.axesTrans.clear()
        self.axesTrans.set_title("Transmission of AlGaN", fontsize=14)
        self.axesTrans.set_xlabel('wavelength, nm', fontsize=10)
        self.axesTrans.set_ylabel('transmission, %', fontsize=10)
        self.axesTrans.grid(True)
        pylab.setp(self.axesTrans.get_xticklabels(), fontsize=10)
        pylab.setp(self.axesTrans.get_yticklabels(), fontsize=10)
        
        
        wavelengthSet = range(200, 801) # [200, ..., 800] - 601 position
        absorptionSet = []
        refractiveSet = []
        transmissiSet = []
        for waveLength in wavelengthSet:
            absorptionSet.append(AbsorptionCoefficient.AbsorpKoeff(Sample.x, waveLength*1e-9))
            refractiveSet.append(Sample.refractiveIndex(waveLength))
            transmissiSet.append(Sample.transmissionBasis(waveLength, Sample.d, Sample.x))
        self.axesAbs.plot(wavelengthSet, absorptionSet, label="x={0}".format(Sample.x), alpha=0.44, picker=5)
        self.axesRefr.plot(wavelengthSet, refractiveSet, label="x={0}".format(Sample.x), alpha=0.44, picker=5)
        self.axesTrans.plot(wavelengthSet, transmissiSet, label="x={0}".format(Sample.x), alpha=0.44, picker=5)

        #pylab.setp(self.axesAbs.legend(), size=8)
        self.axesAbs.legend()
        self.axesRefr.legend()
        self.axesTrans.legend()
        
        self.canvas1.draw()
        self.canvas2.draw()
        self.canvas3.draw()

    def buttonCalculatePress(self, event):
        self.drawGraphs()
        
        

        
        

# run the program
if __name__ == '__main__':
    app = wx.App()
    MainWindow(None, title='AlGaN optical properties')
    app.MainLoop()
