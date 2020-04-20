import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
import numpy as np

# All the variable placeholders///////////////////////////////////////////////////////

NumberOfPoints = int(1e5)
spectrum = np.zeros(NumberOfPoints)
spectralMode = 'full'
concentration = 1 #molar
wlMin = 400 #nm
wlMax = 800
#wlrange = np.linspace(1,1000,100)

#////////////////////////////////////////////////////////////////////////////////////////
def wlRange(*args):
    try:
        global wlrange
        wlStart = float(wlMin.get())
        wlStop = float(wlMax.get())
        wlrange = np.linspace(wlStart,wlStop,NumberOfPoints)
    except ValueError:
        pass

def wlSingle(*args):
    try:
        value1 = float(wlMin.get())
        
    except ValueError:
        pass

    

def spectrumAbsorbance(*args):
    global spectrum 
    global photonsabsorbed
    global peakfreq
    global freqbw
    try:
        wavenumsource = 1/(wlrange/1e7)
        phosec = (wlrange*intensity)/(h*c) #photons per second per cm**2
        omega = wavenumsource
        gamma = HWHM
        Concentration = concentration
        concMolecules = Concentration*6.02e23
        Nmolecule = concMolecules/1000 #molecules per cm^3
        sigmap = CrossSection #cm^2
        peaknm = absorptionPeak #nm
        peak = peaknm*1e-9 #m
        peakwn = 1/(peaknm*1e-7)
        omegao = peakwn
        sigmaw = (sigmap*gamma**2)/((gamma**2)+((omegao-omega)**2))#normalized to peak area
        phosecs = phosec*np.exp((-(sigmaw)*path*Nmolecule)) #photons per second per cm^2 transmitted through 
        spectrum = np.log10(phosec/phosecs)
        plotSpectrum(spectrum)
    except ValueError:
        pass


def singleAbsorbance(*args):
    global count
    global photonsabsorbed
    global peakfreq
    global freqbw
    try:
        global value
        wavenumsource = 1/(wlrange/1e7)
        phosec = (wlrange*intensity)/(h*c) #photons per second per cm**2
        omega = wavenumsource
        gamma = HWHM
        Concentration = conc
        concMolecules = Concentration*6.02e23
        Nmolecule = concMolecules/1000 #molecules per cm^3
        sigmap = CrossSection #cm^2
        peaknm = wavelength #nm
        peak = peaknm*1e-9 #m
        peakwn = 1/(peaknm*1e-7)
        omegao = peakwn
        sigmaw = (sigmap*gamma**2)/((gamma**2)+((omegao-omega)**2))#normalized to peak area
        intensitys = intensity*np.exp((-(sigmaw)*path*Nmolecule))
        phosecs = phosec*np.exp((-(sigmaw)*path*N)) #photons per second per cm^2 transmitted through 
        value = np.log10(phosec/phosecs)
        
    except ValueError:
        pass


    
    
def giveResults(*args):
    if spectralMode is 'full':
        spectrumAbsorbance(*args)
    elif spectralMode is 'single':
        singleAbsorbance(*args)
    else:
        pass
    

option = tk.Tk()
option.wm_title("Spectrophotometry Parameters")
optionWindow = tk.ttk.Frame(option, padding="6 6 12 12")
optionWindow.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
optionWindow.columnconfigure(0, weight=10)
optionWindow.rowconfigure(0, weight=1)

#wavelength setting options//////////////////////////////////////////////////
wlMin = tk.StringVar()
wlMax = tk.StringVar()

wlMin_entry = tk.ttk.Entry(optionWindow, width=7, textvariable=wlMin)
wlMin_entry.grid(column=2, row=1, sticky=(tk.W, tk.E))

wlMax_entry = tk.ttk.Entry(optionWindow, width=7, textvariable=wlMax)
wlMax_entry.grid(column=2, row=2, sticky=(tk.W, tk.E))

ttk.Label(optionWindow, text="Minimum Wavelength").grid(column=1, row=1, sticky=tk.W)
ttk.Label(optionWindow, text="Maximum Wavelength").grid(column=1, row=2, sticky=tk.E)

ttk.Button(optionWindow, text="Set Range", command=wlRange).grid(column=3, row=1, sticky='w')

#///////////////////////////////////////////////////////////////////////////////


# Mode Selection////////////////////////////////////////////////
ttk.Radiobutton(optionWindow, text="Single Wavelength", variable = spectralMode, value = 'single').grid(column=5, row=5, sticky='w')
ttk.Radiobutton(optionWindow, text="Full Spectrum", variable = spectralMode, value = 'full').grid(column=5, row=6, sticky='w')
#///////////////////////////////////////////////////////////////////////

#User Parameter Input////////////////////////////////////////////////////////////
concentration_entry = tk.ttk.Entry(optionWindow, width=7, textvariable=concentration)
concentration_entry.grid(column=2, row=3, sticky=(tk.W, tk.E))

ttk.Label(optionWindow, text="concentration (M)").grid(column=1, row=3, sticky='w')

#////////////////////////////////////////////////////////////////





#Hidden Parameters///////////////////////////////////////////////////

PathLength = .1 #cm
LaserPower = 10e-3 #watts
BeamRadius= 12 #cm

path = PathLength
power = LaserPower
spotradius = BeamRadius
sampleradius= BeamRadius

spotarea = np.pi*(spotradius**2)
samplearea = np.pi*(sampleradius**2)
totalintensity = power/spotarea #watt/cm**2, assumes non-gaussian beam
intensity =  totalintensity*samplearea

HWHM = 1000 #wavenumbers
CrossSection = 1e-18 #cm^2
absorptionPeak = 500 #nanometers

h = 6.626e-34 #joule-seconds
c = 3e8 #m/s
e = 1.6e-19 #coulombs per electron


#/////////////////////////////////////////////////


# ANSWER OUTPUT////////////////////////////////////////////////////
#ttk.Label(optionWindow, textvariable=absorbance).grid(column=4, row=6, sticky=tk.W)

#/////////////////////////////////////////////////////////////////////


#Plotting/////////////////////////////////////////////////////////////////////////////
fig = Figure(figsize=(5, 4), dpi=100)
ax1 = fig.add_subplot(111)

def plotSpectrum(*args):    
    try:
        wavelength = wlrange
        ax1.clear()
        ax1.plot(wavelength, spectrum)
        ax1.autoscale(enable = True)
        
        root = tk.Toplevel()
        root.wm_title("Absorption Spectroscopy")
        
        canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        button = tk.Button(master=root, text="Quit", command=_quit)
        button.pack(side=tk.BOTTOM)
        
        canvas.mpl_connect("key_press_event", on_key_press)
    except ValueError:
        pass
    
"""def clearSpectrum(*args):    
    try:
        fig.clf()
        fig.add_subplot(111).plot()
    except ValueError:
        pass"""
    


ttk.Button(optionWindow, text="Get Results", command=giveResults).grid(column=3, row=6, sticky='w')
#ttk.Button(optionWindow, text="Clear Spectrum", command=clearSpectrum).grid(column=3, row=5, sticky='w')
        


def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, canvas, toolbar)





def _quit():
    option.quit()     # stops mainloop
    option.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

        

ttk.Button(optionWindow, text="QUIT", command=_quit).grid(column=6, row=1, sticky='w')


tk.mainloop()
# If you put root.destroy() here, it will cause an error if the window is
# closed with the window manager.