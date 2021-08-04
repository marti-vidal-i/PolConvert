#
# This file was generated using xslt from its XML file
#
# Copyright 2009, Associated Universities Inc., Washington DC
#
import sys
import os
from  casac import *
import string
from taskinit import casalog
from taskinit import xmlpath
#from taskmanager import tm
import task_polconvert
def polconvert(IDI='', OUTPUTIDI='', DiFXinput='', DiFXcalc='', doIF=[], linAntIdx=[1], Range=[], ALMAant='', spw=-1, calAPP='', calAPPTime=[0.,5.], APPrefant='', gains=[["NONE"]], interpolation=[], gainmode=[], XYavgTime=0.0, dterms=["NONE"], amp_norm=0.01, XYadd={}, XYdel={}, XYratio={}, usePcal=[], swapXY=[False], swapRL=True, feedRotation=[], correctParangle=False, IDI_conjugated=False, plotIF=-1, plotRange=[], plotAnt=-1, excludeAnts=[], excludeBaselines=[], doSolve=-1, solint=[1,1], doTest=True, npix=50, solveAmp=True, solveMethod='gradient', calstokes=[1.,0.,0.,0.], calfield=-1):

        """\n\nVersion 1.9.0 -- Converts VLBI visibilities polarization basis.

For more information about the internals of PolConvert, please read:

Marti-Vidal et al. 2016, Astronomy and Astrophysics, 587, 143 

PROCEDURE:

If a VLBI antenna used linear-feed receivers, PolConvert must first
estimate the X-Y cross gain for that antenna (phase, amplitude, 
and multi-band delay) before the final conversion. Use the plotting
option of PolConvert to plot a selected scan of a given baseline and 
that scan will be used to estimate the cross-polarization gains. 
PolConvert returns a list with the amplitude and phase cross-gains 
for all the antennas, as long as it is run in plotting mode. The user 
can then set these gain lists to XYadd and XYratio for a second run 
of PolConvert.

Given an FITS-IDI file (or set of SWIN DiFX files) with phased-ALMA 
observations, the user must first calibrate completely (i.e., in full 
polarization) the corresponding ALMA measurement set. It is important 
to create it using asis='CALAPPPHASE' in the importasdm task. 

If more than one asdm was created in the observations, the user must 
concatenate all the CALAPPPHASE tables of each asdm, for PolConvert to 
work properly. This can be done with the following commands (we assume 
here that MSLIST is a python list with the names of the measurement sets 
for each asdm):

for i,myms in enumerate(MSLIST):
  if i==0:
    os.system('cp -rf %s/ASDM_CALAPPPHASE ./CALAPPPHASE.tab'%myms)
  else:
    tb.open('%s/ASDM_CALAPPPHASE'%myms)
    tb.copyrows('./CALAPPPHASE.tab')
    tb.close()

These lines will create the table './CALAPPPHASE.tab', to be used by 
PolConvert (i.e., the table specified in the 'calAPP' keyword).

PolConvert can also do this for you, if you set calAPP = MSLIST. But
check carefully the information that it will print, regarding the 
time coverage of both the CALAPPPHASE table and the MSs.

Let us assume that the calibration tables are named 'gain.tb' and 
'bandpass.tb' (there can be many others, for delay, XY-phase, etc.) and 
the D-term table is called 'dterms.tb'. Then, with this assumption:

- If ALMA is the only station with linear receivers (let's say it is 
station number 1 in the FITS-IDI file), the call to PolConvert should be 
done with the following keyword values:

- linAntIdx = [1]

- Range = []    # i.e., all data will be converted

- ALMAvis = 'TheALMAvisibilities.ms' 
 
- spw = 0  # it may be a good idea to split first the science spw, 
           # before concatenating and calibrating the ms

- calAPP = './CALAPPPHASE.tab'

- calAPPTime = [0., 5.0]   # The CALAPP entries sometime start and end
                           # with time lags in between. The time tolerance
                           # of 5 seconds should avoid problems related
                           # to this.

- gains = [['gain.tb','bandpass.tb']]
- dterms = ['dterms.tb']

- doTest = False # to actually APPLY the changes, not only compute them!


###################
   
     SPECIAL CASE 1:

If there was a second antenna with linear-polarization receivers, 
it can also be converted, but has to be listed after ALMA. Let us assume
that this antenna has id=4 in the FITS-IDI file. Then:

- linAntIdx = [1,4]   # i.e., ALMA plus the other linear-pol. station.

- gains = [ ['gain.tb','bandpass.tb'] , ['NONE'] ]

     # Notice that 'NONE' can be used to tell PolConvert not to apply
     # any calibration to antenna #4 before conversion.

- dterms = [ 'dterms.dt' , 'NONE']


###################
   
     SPECIAL CASE 2:

If the user wants to check the conversion before applying it, PolConvert
can plot the fringes for a given IF, baseline and time range. Let us 
assume we want to plot the baseline to antenna 2 (i.e., baseline 1-2) in 
the time range 0-07:30:00 to 0-07:31:00 (AIPS format). Then:

- doTest = True  # i.e., do NOT write on the FITS-IDI file!

- plotIF = 1  # i.e., first IF in FITS-IDI file.

- plotRange = [0,7,30,0,0,7,31,0]

- Range = [0,7,30,0,0,7,31,0]  # i.e., to not waste resources computing
                               # things that we will not save nor plot.

- plotAnt = 2  


###################
   
     SPECIAL CASE 2:

For the two linear-pol antennas, use the pcal tones to correct the 
phase difference between X and Y. In addition to this, for the first 
antenna, the X/Y relative amplitude is estimated from the phasecal 
tones. For the second antenna, the X/Y amplitude spectrum is 
estimated from the autocorrelations, using a running median
filter of 51 channels (if the number of channels per IF is 510):


- usePcal = [True,True]
- XYratio = [-10., 0.0]


Notice that, if the GCPFF algorithm is used to solve for the X/Y gains,
these will be stored in the 'XYratio' and 'XYadd' keys of the returned
dictionary, whereas the a-priori gains computed from the pcals and the
autocorrelations will be stored as complex arrays in 'aPrioriXYGain'.



###################
   
     OTHER SPECIAL CASES (NOT FULLY TESTED):

1.- If two antennas have linear-pol receivers (i.e., ALMA plus another one)
and the second one was correlated with the pol. channels swapped, then:

- swapXY = [False, True]

If it was ALMA the antenna with swapped pol. channels, then:

- swapXY = [True, False]



2.- If the second antenna with linear-pol receivers had an offset of, say,
65 degrees between X and Y, this offset can be corrected before conversion:

- XYadd = [0.0, 65.]

  If there are 4 IFs and the X-Y phases for the second antenna differ among
  IFs, we can set them in the following way:

- XYadd = [[0.0, 0.0, 0.0, 0.0], [65., 30., 25., 10.]]


NOTICE THAT POLCONVERT CAN ESTIMATE THE XYADD AND XYRATIO, FROM 
THE SCAN THAT YOU ASK IT TO PLOT. IF IT IS A SCAN OF A STRONG CALIBRATOR,
POLCONVERT CAN ESTIMATE ALL THESE QUANTITIES FOR YOUR VLBI STATIONS WITH
LINEAR FEEDS.








#########################################
#########################################
##### A TYPICAL ALMA DATA REDUCTION #####

1.- Import all ASDMs into measurement sets. 
    BEWARE THAT YOU SET asis='CalAppPhase'

2.- Find out the spw that matches the VLBI frequency setup. Split it for
    each measurement set.

3.- Concatenate the splitted measurement sets. 

4.- Concatenate the 'ASDM_CALAPPPHASE' tables of the measurement sets
    into a new CALAPPPHASE TABLE (see help above).

5.- Calibrate the concatenated measurement set (full-polarization 
    calibration) using standard ALMA procedures.

6.- Execute polconvert, feeding it with the calibration tables and 
    the CALAPPPHASE table, in mode "doTest = True", and applying it
    only to a short (say, 1-2 minutes) scan of a strong source 
    (e.g., the phase calibrator). Select a given IF and antenna 
    to plot (it's better to select a short VLBI baseline to ALMA).

7.- The program will plot the fringes and print an estimate of the 
    extra X/Y phase that should be added to ALMA. This number should 
    be small, as long as the calibration is OK. If a large number 
    is found, and you see large cross-hand polarizations, add this 
    extra X-Y phase to polconvert, via the keyword "XYadd"

8.- Re-execute polconvert in test mode. Check whether the conversion 
    is satisfactory.

9.- Once satisfied with the conversion of the selected calibrator scan, 
    execute polconvert over the shole dataset with "doTest = False".

10.- It may be a good idea to do extra sanity checks, like the 
    possible dependence of XYadd with IF and/or its eventual time 
    evolution. All these effects should have been properly corrected 
    if the measurement set calibration was successful. Different 
    XYadd phases can be added to different IFs by converting each IF
    separately.


# END OF POLCONVERT DOCUMENTATION
#######################################



        """

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['IDI'] = IDI
        mytmp['OUTPUTIDI'] = OUTPUTIDI
        mytmp['DiFXinput'] = DiFXinput
        mytmp['DiFXcalc'] = DiFXcalc
        mytmp['doIF'] = doIF
        mytmp['linAntIdx'] = linAntIdx
        mytmp['Range'] = Range
        mytmp['ALMAant'] = ALMAant
        mytmp['spw'] = spw
        mytmp['calAPP'] = calAPP
        mytmp['calAPPTime'] = calAPPTime
        mytmp['APPrefant'] = APPrefant
        mytmp['gains'] = gains
        mytmp['interpolation'] = interpolation
        mytmp['gainmode'] = gainmode
        mytmp['XYavgTime'] = XYavgTime
        mytmp['dterms'] = dterms
        mytmp['amp_norm'] = amp_norm
        mytmp['XYadd'] = XYadd
        mytmp['XYdel'] = XYdel
        mytmp['XYratio'] = XYratio
        mytmp['usePcal'] = usePcal
        mytmp['swapXY'] = swapXY
        mytmp['swapRL'] = swapRL
        mytmp['feedRotation'] = feedRotation
        mytmp['correctParangle'] = correctParangle
        mytmp['IDI_conjugated'] = IDI_conjugated
        mytmp['plotIF'] = plotIF
        mytmp['plotRange'] = plotRange
        mytmp['plotAnt'] = plotAnt
        mytmp['excludeAnts'] = excludeAnts
        mytmp['excludeBaselines'] = excludeBaselines
        mytmp['doSolve'] = doSolve
        mytmp['solint'] = solint
        mytmp['doTest'] = doTest
        mytmp['npix'] = npix
        mytmp['solveAmp'] = solveAmp
        mytmp['solveMethod'] = solveMethod
        mytmp['calstokes'] = calstokes
        mytmp['calfield'] = calfield
        pathname="file:///home/marti/WORKAREA/GITHUB/PolConvert/"
        trec = casac.utils().torecord(pathname+'polconvert.xml')

        casalog.origin('polconvert')
        if trec.has_key('polconvert') and casac.utils().verify(mytmp, trec['polconvert']) :
            result = task_polconvert.polconvert(IDI, OUTPUTIDI, DiFXinput, DiFXcalc, doIF, linAntIdx, Range, ALMAant, spw, calAPP, calAPPTime, APPrefant, gains, interpolation, gainmode, XYavgTime, dterms, amp_norm, XYadd, XYdel, XYratio, usePcal, swapXY, swapRL, feedRotation, correctParangle, IDI_conjugated, plotIF, plotRange, plotAnt, excludeAnts, excludeBaselines, doSolve, solint, doTest, npix, solveAmp, solveMethod, calstokes, calfield)

        else :
          result = False
        return result
