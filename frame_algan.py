import wx
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import math
import AbsorptionCoefficient

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
    
    def __init__(self, parent, x):
        self.x = x

        Eg_AlN = 6.13
        Eg_GaN = 3.42
        b = 1.08 # MBE material
        self.Eg = self.x * Eg_AlN + (1 - self.x) * Eg_GaN - b * self.x * (1 - self.x)

    def absorptionCoeff(self, Energy):
        alpha0 = 1e5 # cm-1*eV-1
        return alpha0 * math.sqrt(Energy - self.Eg)

    def refractiveIndex(self, wavelen, kind='o'):
        '''
        '''
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

    def transmission1(self, wavelen, d, x):
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

    def absorptionCoeff_my(self, wavelen):
        '''
        '''
        # wavelen [nm] -> hv [eV]
        hv = 1.2398619043e-6 / (wavelen * 1e-9)

        A1 = 0.465
        A2 = -9.189137089 * self.x + 7.669893919
        A3 = 0.2325
        C1 = -0.4570245289 * self.x + 1.151745336
        C2 = 27.26910615 * self.x - 23.03122697
        C3 = 4.141167261 - 0.0121604643 * self.x - 1.450673214 * (self.x ** 2) + 3.291528196 * (self.x ** 3) - 4.355129034 * (self.x ** 4)

        a = [3.683 + 0.02 * i - 4.676 * self.x for i in range(9)]
        b = [1 for i in range(9)]
        b[0] = A1 * a[0] + C1
        b[8] = A2 * a[8] + C2
        '''b[4] = 0.206 * b[0] + 0.794 * b[8]
        b[2] = 0.312 * b[0] + 0.688 * b[4]
        b[6] = 0.392 * b[4] + 0.608 * b[8]
        b[1] = 0.389 * b[0] + 0.611 * b[2]
        b[3] = 0.389 * b[2] + 0.611 * b[4]
        b[5] = 0.389 * b[4] + 0.611 * b[6]
        b[7] = 0.389 * b[6] + 0.611 * b[8]
        ''' #another var of coefficients b[]:
        b[4] = 0.8 * b[0] + 0.2 * b[8]
        b[2] = 0.7 * b[0] + 0.3 * b[4]
        b[6] = 0.6 * b[4] + 0.4 * b[8]
        b[1] = 0.6 * b[0] + 0.4 * b[2]
        b[3] = 0.6 * b[2] + 0.4 * b[4]
        b[5] = 0.6 * b[4] + 0.4 * b[6]
        b[7] = 0.6 * b[6] + 0.4 * b[8]
        

        L1 = calculatePolinom(a, b, hv)
        
        a = [4.067 + 0.01 * i - 4.870 * self.x for i in range(9)]
        b = [1 for i in range(9)]
        b[0] = A2 * a[0] + C2
        b[8] = A3 * a[8] + C3
        '''b[4] = 0.206 * b[0] + 0.794 * b[8]
        b[2] = 0.312 * b[0] + 0.688 * b[4]
        b[6] = 0.392 * b[4] + 0.608 * b[8]
        b[1] = 0.389 * b[0] + 0.611 * b[2]
        b[3] = 0.389 * b[2] + 0.611 * b[4]
        b[5] = 0.389 * b[4] + 0.611 * b[6]
        b[7] = 0.389 * b[6] + 0.611 * b[8]
        '''
        b[4] = 0.8 * b[0] + 0.2 * b[8]
        b[2] = 0.7 * b[0] + 0.3 * b[4]
        b[6] = 0.6 * b[4] + 0.4 * b[8]
        b[1] = 0.6 * b[0] + 0.4 * b[2]
        b[3] = 0.6 * b[2] + 0.4 * b[4]
        b[5] = 0.6 * b[4] + 0.4 * b[6]
        b[7] = 0.6 * b[6] + 0.4 * b[8]
        
        L2 = calculatePolinom(a, b, hv)

        linear = lambda hv_, A_, C_: hv_ * A_ + C_

        if hv <= 3.683 - 4.676 * self.x:
            return 10.0 ** linear(hv, A1, C1)
        elif 3.683 - 4.676 * self.x < hv <= 3.843 - 4.676 * self.x:
            return 10.0 ** 0
        elif 3.843 - 4.676 * self.x < hv <= 4.067 - 4.870 * self.x:
            return 10.0 ** linear(hv, A2, C2)
        elif 4.067 - 4.870 * self.x < hv <= 4.147 - 4.870 * self.x:
            return 10.0 ** 0
        else:
            if linear(hv, A3, C3) < 5.04:
                return 10.0 ** linear(hv, A3, C3)
            else:
                return 10.0 ** 5.04
            
            
class p1(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, size=(800, 400))

        # how to create canvas for plots... figsize=(inches, inches), dpi=pixels per inch (zoom)
        self.figure = matplotlib.figure.Figure(figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
        
        # here I changed visibility of my plots to look greater!
        self.figure.subplots_adjust(left=0.10, bottom=0.1, right=0.95, top=0.89, hspace=0.7, wspace=0.7)
        
        self.axes = self.figure.add_subplot(311)
        xval = []
        yval = []
        
        for ctr in range(200, 800):
            wavelen = ctr
            xval.append(wavelen)
            yval.append(AbsorptionCoefficient.AbsorpKoeff(sample.x, wavelen*1E-9))
            #print wavelen, ' ', AbsorptionCoefficient.AbsorpKoeff(sample.x, wavelen*1E-9)
            
        self.axes.semilogy()  #logarithmic values of y axis
        self.axes.plot(xval, yval)
        self.axes.set_title("Transmission, Absorption Coefficient and Refractive Index of AlGaN")
        self.axes.set_xlabel('wavelength, nm')
        self.axes.set_ylabel('absorption coefficient, sm-1')
        self.axes.grid(True)
        print self.figure.subplotpars

        self.axes2 = self.figure.add_subplot(312)
        xval = []
        yval = []
        
        for ctr in range(200, 800):
            wavelen = ctr
            xval.append(wavelen)
            yval.append(sample.refractiveIndex(wavelen))
            #print wavelen, ' ', sample.refractiveIndex(wavelen)

        #self.axes.semilogy()
        self.axes2.plot(xval, yval)
        self.axes2.set_xlabel('wavelength, nm')
        self.axes2.set_ylabel('refractive index')
        self.axes2.grid(True)


        self.axes3 = self.figure.add_subplot(313, axisbg='w', title='Transmission')
        xval = []
        yval = []
        
        for ctr in range(200, 800):
            wavelen = ctr
            xval.append(wavelen)
            yval.append(sample.transmission1(wavelen, 1*1e-6, 0.45))
            #print wavelen, ' ', sample.transmission1(wavelen, 1*1e-6, 0.45)

        #self.axes.semilogy()
        self.axes3.plot(xval, yval)
        self.axes3.set_xlabel('wavelength, nm')
        self.axes3.set_ylabel('transmission, %')
        self.axes3.grid(True)
        
        self.canvas = FigureCanvas(self, -1, self.figure)

class TestFrame(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title)

        self.Position = (200, 50)
        self.ClientSize = (900, 600)        # Set client area. Frame size is unimportant.

        frm_pnl = wx.Panel(self)
        frm_pnl.SetBackgroundColour((0, 179, 0))
        #self.sp = wx.SplitterWindow(frm_pnl)
        left_panel  = p1(frm_pnl)
        left_panel.SetBackgroundColour((179, 179, 179))
        right_panel = wx.Panel(frm_pnl)
        right_panel.SetBackgroundColour((255, 255, 255))
        
        #sp = wx.SplitterWindow(frm_pnl)
        #self.p1 = p1(self.sp)
        #self.p2 = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
        #self.sp.SplitVertically(left_panel, right_panel, 400)

        frmPnl_horzSizer = wx.BoxSizer(wx.HORIZONTAL)
        frmPnl_horzSizer.Add( left_panel,  proportion=1, flag=wx.EXPAND )
        frmPnl_horzSizer.Add( right_panel, proportion=2, flag=wx.EXPAND )

        #-----  Invoke the sizer via its container.
        
        frm_pnl.SetSizer( frmPnl_horzSizer )
        frm_pnl.Layout()



app = wx.App(redirect=False)
sample = AlGaN(None, 0.5)
frame = TestFrame(None, "AlGaN optical properties")
frame.Show()
app.MainLoop()
