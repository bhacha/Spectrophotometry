import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
import numpy as np



option = tk.Tk()
option.wm_title("Spectrophotometry Parameters")
optionWindow = tk.ttk.Frame(option, padding="6 6 12 12")
optionWindow.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
optionWindow.columnconfigure(0, weight=10)
optionWindow.rowconfigure(0, weight=1)


# All the variable placeholders///////////////////////////////////////////////////////

NumberOfPoints = int(1e5)
spectrum = np.zeros(NumberOfPoints)
photonsabsorbed = spectrum
spectralMode = tk.StringVar()
spectralMode.set(0)
concentration = .001 #molar
concentrationString = tk.StringVar()
concentrationString.set(concentration)
absorbanceValue = tk.StringVar()

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
    global value1
    try:
        value1 = float(wlMin.get())
    except ValueError:
        pass

def Absorbance(wavelength,CrossSection,HWHM,conc):
    global spectrum 
    global photonsabsorbed
    global wlrange
    global phosec
    gamma = HWHM
    wavenumsource = 1/(wlrange/1e7)
    omega = wavenumsource
    Concentration = conc
    concMolecules = Concentration*6.02e23
    N = concMolecules/1000 #molecules per cm^3
    sigmap = CrossSection #cm^2
    peaknm = wavelength #nm
    peak = peaknm*1e-9#m
    peakwn = 1/(peaknm*1e-7)
    omegao = peakwn
    sigmaw = (sigmap*gamma**2)/((gamma**2)+((omegao-omega)**2))#normalized to peak area
    intensitys = intensity*np.exp((-(sigmaw)*path*N))
    phosecs = phosec*np.exp((-(sigmaw)*path*N)) #photons per second per cm^2 transmitted through 
    peakfreq = c/peak
    freqbw = peakfreq*((peak+(2e-9))-(peak-2e-9))/(peak)
    spectrum = spectrum+np.log10(phosec/phosecs)
    photonsabsorbed = photonsabsorbed + (phosec-phosecs)

def spectrumAbsorbance(*args):
    global spectrum 
    global photonsabsorbed
    global phosec
    global wlrange
    global concentration
    try:
        concentration = float(concentrationString.get())
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
#This line adds interfering species/////////////////////////////////////////////////////////////////
        #Absorbance(Peak Wavelength, Cross Section, HWHM in wavenumbers, Concentration)
        Absorbance(250,1e-20,20000,1)
#///////////////////////////////////////////////////////////////////////////////////////////////    
        plotSpectrum(spectrum)
    except ValueError:
        pass


def singleAbsorbance(*args):
    global concentration
    global count
    global photonsabsorbed
    global peakfreq
    global freqbw
    global CrossSection
    global absorbanceValue
    global value1
    wlSingle()
    try:
        concentration = float(concentrationString.get())
        wavenumsource = 1/(value1/1e7)
        phosec = (value1*intensity)/(h*c) #photons per second per cm**2
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
        absorbanceValue.set(np.log10(phosec/phosecs))
    except ValueError:
        pass


def enableEntry(*args):
    wlMax_entry.config(state = 'normal')
    wlMax_entry.update()

def disableEntry(*args):
    wlMax_entry.config(state = 'disabled')
    wlMax_entry.update()

def giveResults(*args):
    test = float(spectralMode.get())
    if test == 1:
        spectrumAbsorbance(*args)
    elif test == 0:
        singleAbsorbance(*args)
    else:
        pass
    


#wavelength setting options//////////////////////////////////////////////////
wlMin = tk.StringVar()
wlMax = tk.StringVar()

wlMin.set("300") #default minimum
wlMax.set("800") #default maximum

wlMin_entry = tk.ttk.Entry(optionWindow, width=7, textvariable=wlMin)
wlMin_entry.grid(column=2, row=1, sticky=(tk.W, tk.E))

wlMax_entry = tk.ttk.Entry(optionWindow, width=7, textvariable=wlMax)
wlMax_entry.grid(column=2, row=2, sticky=(tk.W, tk.E))

ttk.Label(optionWindow, text="Minimum Wavelength").grid(column=1, row=1, sticky=tk.W)
ttk.Label(optionWindow, text="Maximum Wavelength").grid(column=1, row=2, sticky=tk.E)

ttk.Button(optionWindow, text="Set Range", command=wlRange).grid(column=3, row=1, sticky='w')

#///////////////////////////////////////////////////////////////////////////////


# Mode Selection////////////////////////////////////////////////
ttk.Radiobutton(optionWindow, text="Single Wavelength", variable = spectralMode, value = 0, command = disableEntry).grid(column=5, row=5, sticky='w')
ttk.Radiobutton(optionWindow, text="Full Spectrum", variable = spectralMode, value = 1,command = enableEntry).grid(column=5, row=6, sticky='w')
#///////////////////////////////////////////////////////////////////////

#User Parameter Input////////////////////////////////////////////////////////////
concentration_entry = tk.ttk.Entry(optionWindow, width=7, textvariable=concentrationString)
concentration_entry.grid(column=2, row=3, sticky=(tk.W, tk.E))

ttk.Label(optionWindow, text="concentration (M)").grid(column=1, row=3, sticky='w')

#////////////////////////////////////////////////////////////////




#Hidden Parameters///////////////////////////////////////////////////

PathLength = 1 #cm
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

HWHM = 4000 #wavenumbers
CrossSection = 1e-19 #cm^2
absorptionPeak = 500 #nanometers

h = 6.626e-34 #joule-seconds
c = 3e8 #m/s


#/////////////////////////////////////////////////


#ANSWER OUTPUT////////////////////////////////////////////////////
ttk.Label(optionWindow, textvariable=absorbanceValue).grid(column=4, row=4, sticky=tk.W)



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