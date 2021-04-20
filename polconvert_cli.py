#
# This file was generated using xslt from its XML file
#
# Copyright 2014, Associated Universities Inc., Washington DC
#
import sys
import os
import datetime
#from casac import *
import casac
import string
import time
import inspect
import numpy
from casa_stack_manip import stack_frame_find
from odict import odict
from types import *
from task_polconvert import polconvert
class polconvert_cli_:
    __name__ = "polconvert"
    rkey = None
    i_am_a_casapy_task = None
    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.

    def __init__(self) :
       self.__bases__ = (polconvert_cli_,)
       self.__doc__ = self.__call__.__doc__

       self.parameters={'IDI':None, 'OUTPUTIDI':None, 'DiFXinput':None, 'DiFXcalc':None, 'doIF':None, 'linAntIdx':None, 'Range':None, 'ALMAant':None, 'spw':None, 'calAPP':None, 'calAPPTime':None, 'APPrefant':None, 'gains':None, 'interpolation':None, 'gainmode':None, 'XYavgTime':None, 'dterms':None, 'amp_norm':None, 'XYadd':None, 'XYdel':None, 'XYratio':None, 'usePcal':None, 'swapXY':None, 'swapRL':None, 'feedRotation':None, 'correctParangle':None, 'IDI_conjugated':None, 'plotIF':None, 'plotRange':None, 'plotAnt':None, 'excludeAnts':None, 'excludeBaselines':None, 'doSolve':None, 'solint':None, 'doTest':None, 'npix':None, 'solveAmp':None, 'solveMethod':None, 'calstokes':None, 'calfield':None, }


    def result(self, key=None):
            #### and add any that have completed...
            return None


    def __call__(self, IDI=None, OUTPUTIDI=None, DiFXinput=None, DiFXcalc=None, doIF=None, linAntIdx=None, Range=None, ALMAant=None, spw=None, calAPP=None, calAPPTime=None, APPrefant=None, gains=None, interpolation=None, gainmode=None, XYavgTime=None, dterms=None, amp_norm=None, XYadd=None, XYdel=None, XYratio=None, usePcal=None, swapXY=None, swapRL=None, feedRotation=None, correctParangle=None, IDI_conjugated=None, plotIF=None, plotRange=None, plotAnt=None, excludeAnts=None, excludeBaselines=None, doSolve=None, solint=None, doTest=None, npix=None, solveAmp=None, solveMethod=None, calstokes=None, calfield=None, ):

        """\n\nVersion 1.8.2 -- Converts VLBI visibilities polarization basis.

        Detailed Description:
\n\nVersion 1.8.2 -- Converts VLBI visibilities from mixed-polarization (linear-circular) into circular basis. Works with single VLBI stations as well as with calibrated phased arrays (i.e., phased ALMA).\n\n
        Arguments :
                IDI: Input FITS-IDI file with VLBI visibilities. It can also be a direcotry containing SWIN files from DiFX.
                   Default Value: 

                OUTPUTIDI: Output FITS-IDI file (or SWIN directory). If equal to IDI, the file(s) will be overwritten
                   Default Value: 

                DiFXinput: If SWIN files are being converted, this must be the *.input file used by DiFX.
                   Default Value: 

                DiFXcalc: If SWIN files are being converted, this must be the *.calc file used by DiFX. This is optional, but if it is not provided, the cross-polarization gain estimates may be incorrect if doSolve>0.
                   Default Value: 

                doIF: List of IFs to convert. Default means all.
                   Default Value: []

                linAntIdx: List of indices of the linear-polarization antennas in the IDI file (lowest index starts with 1)
                   Default Value: [1]

                Range: Time range to convert (integer list; AIPS format). Default means all data
                   Default Value: []

                ALMAant: If ALMA has been used, this is the antenna table from the MS with the intra-ALMA visibilities.
                   Default Value: 

                spw: Spectral window in ALMAvis that contains the VLBI band. If negative, the program will derive it automatically.
                   Default Value: -1

                calAPP: If ALMA has been used, this is the combined ASDM_CALAPPPHASE table from the ASDM. The list of measurement sets can also be given (so the table is concatenated from all of them).
                   Default Value: 

                calAPPTime: Time shift and time tolerance (in sec) for the CALAPPPHASE table obtained from the ASDM.
                   Default Value: [0.,5.]

                APPrefant: If not empty, re-reference the TelCal phases, assuming that the X-Y phase-difference table provided in \'gains\' (see keyword below) uses APPrefant as the reference antenna. Notice that the name of the gain table with the X-Y phase differences has to contain the string \'.XY0\'.
                   Default Value: 

                gains: Gain tables to pre-calibrate the linear-pol VLBI stations (one list of gains per linear-pol station).
                   Default Value: [["NONE"]]

                interpolation:  Interpolation type to use (one per calibration table). Tells whether to apply linear or nearest interpolation. Default is to use linear for all tables.
                   Default Value: []

                gainmode:  Mode of gain calibration to impose (one per calibration table). Default is \'T\' for all tables, unless either the string \'XY0\', \'bandpass\' or \'Gxyamp\' appears in the table name. The gain types can be either \'G\' (i.e., split gains per polarization) or \'T\' (i.e., combine polarizations).
                   Default Value: []

                XYavgTime:  Re-compute the G-mode gains by adding a time smoothing of X-Y differences to the T-mode gains. Default is NOT to do this (i.e., use the G-mode gains as given). If positive, use a running time average with this size (in seconds).
                   Default Value: 0.0

                dterms: D-term tables to pre-calibrate the linear-pol VLBI stations (one table per linear-pol station).
                   Default Value: ["NONE"]

                amp_norm: If positive, normalize the amplitude correction to the X-Y average, and save the scaling factor (vs time) in an external (ASCII) file (ANTAB format, assuming a DPFU=amp_norm). If zero, or negative, apply the amplitude correction as is.
                   Default Value: 0.01

                XYadd: Add manually a phase between X and Y before conversion (in deg.). Either a list with one value per linear-pol station OR a list of lists (i.e., one value per IF for each antenna) OR a list of lists of lists (one value per channel, for each IF, for each linear-polarization antenna).
                   Default Value: {}

                XYdel: Add manually a multiband delay between X and Y before conversion (in deg./chanel). One value per linear-pol station.
                   Default Value: {}

                XYratio: Add manually an amplitude ratio between X and Y before conversion (R=X/Y). Follows the same format as XYadd. If a negative value is given for an antenna, the X/Y ratio will be estimated from its autocorrelations (the spectrum for antenna i will be computed using a running-median filter of width equal to -1/XYratio[i] of the IF bandwidth). If 0.0 is given for an antenna, the ratio will be estimated from the phasecal amplitudes (as long as usePcal is True).
                   Default Value: {}

                usePcal: List of booleans (one boolean per linear-polarization station). If True, use the X-Y difference of phasecal tones as an estimate of the X-Y cross-polarization phase. Default means to NOT use the phasecals.
                   Default Value: []

                swapXY: Swap X-Y before conversion. One boolean per linear-pol VLBI station.
                   Default Value: [False]

                swapRL: Swap R-L of the OTHER antenna(s) when plotting the fringes.
                   Default Value: True

                feedRotation: Rotation angle of the feed (one value per antenna, in degrees). Default means zero degrees (so that X is truly horizontal for the linear-pol. antennas). These angles are used in the gain-solver step.
                   Default Value: []

                correctParangle: If True, the correction for parallactic angle is applied to the converted antenna(s).
                   Default Value: False

                IDI_conjugated: Assume a swap in the baseline defintion (i.e., conjugation) of the FITS-IDI file. This has NO effect on SWIN files. 
                   Default Value: False

                plotIF: IF index(es) to plot. Default means to NOT plot. An empty list, [], means to plot ALL IFs being converted (but do not forget to set plotRange and plotAnt!).
                   Default Value: -1

                plotRange: Time range to plot (integer list; AIPS format). Default means to NOT plot
                   Default Value: []

                plotAnt: Index of the other antenna in the baseline to plot. Default means to NOT plot.
                   Default Value: -1

                excludeAnts: List of antennas (i.e., list of antenna codenames) to NOT use in the cross-polarization gain estimates.
                   Default Value: []

                excludeBaselines: List of baselines (i.e., a list of lists of two antenna codenames) to NOT use in the cross-polarization gain estimates.
                   Default Value: []

                doSolve: If negative, do not estimate the cross-polarization gains. If positive or zero, estimate the gains using a Global Cross-Pol Fringe Fitting (GCPFF). The gains are fitted with an error function (Chi Square) defined as:\n\n sum( doSolve*(RR/LL-1)^2 + (RL^2 + LR^2) ),\n\n so that doSolve=0 minimizes the cross-hand polarizations (so it assumes a small linear polarization of the source), whereas doSolve>>1 assumes a negligible Stokes V.
                   Default Value: -1

                solint: If solint[0] null or negative, solve the cross-polarization phase plus a multi-band delay (MBD). If not, solve in bandpass mode by averaging solint[0] channels per solution.\n Divide the solution time range (per scan) in solint[1] chunks (i.e., in solint[1] subscans). I.e., if solint[1]==1, the fringes are fully averaged in time for each scan (but corrected for the scan fringe rates) before the GPLFF condition is computed. solint[2] is the minimum time jump (in seconds) to split the data into different scans (default: 100 seconds).
                   Default Value: [1,1]

                doTest: If true, only compute (and eventually plot), the data, but leave OUTPUTIDI untouched.
                   Default Value: True

                npix: Number of pixels for the fringe plots (and fringe search).
                   Default Value: 50

                solveAmp: if the cross-polarization gains are being estimated, solve also for the X/Y amplitude ratios.
                   Default Value: True

                solveMethod: Method for the minimization of the Chi squared in the GCPFF. Can be \'gradient\', \'Levenberg-Marquardt\' or \'COBYLA\'.
                   Default Value: gradient

                calstokes: Stokes parameters, [I,Q,U,V] of the calibrator (of course, this is only used if doSolve is not negative). The total intensity is not needed in the calibration (i.e., calstokes[0] can be left to 1.0, so that the other parameters will correspond to fractional polarization). 
                   Default Value: [1.,0.,0.,0.]

                calfield: If not negative, field ID of the calibrator (useful if a time range covering several scans is being used in the GCPFF). If negative, use all data in the time range, regardless of the field ID.
                   Default Value: -1

        Returns: bool

        Example :


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
        if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
        #casac = self.__globals__['casac']
        casalog = self.__globals__['casalog']
        casa = self.__globals__['casa']
        #casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = 'polconvert'
        self.__globals__['taskname'] = 'polconvert'
        ###
        self.__globals__['update_params'](func=self.__globals__['taskname'],printtext=False,ipython_globals=self.__globals__)
        ###
        ###
        #Handle globals or user over-ride of arguments
        #
        if type(self.__call__.func_defaults) is NoneType:
            function_signature_defaults={}
        else:
            function_signature_defaults=dict(zip(self.__call__.func_code.co_varnames[1:],self.__call__.func_defaults))
        useLocalDefaults = False

        for item in function_signature_defaults.iteritems():
                key,val = item
                keyVal = eval(key)
                if (keyVal == None):
                        #user hasn't set it - use global/default
                        pass
                else:
                        #user has set it - use over-ride
                        if (key != 'self') :
                           useLocalDefaults = True

        myparams = {}
        if useLocalDefaults :
           for item in function_signature_defaults.iteritems():
               key,val = item
               keyVal = eval(key)
               exec('myparams[key] = keyVal')
               self.parameters[key] = keyVal
               if (keyVal == None):
                   exec('myparams[key] = '+ key + ' = self.itsdefault(key)')
                   keyVal = eval(key)
                   if(type(keyVal) == dict) :
                      if len(keyVal) > 0 :
                         exec('myparams[key] = ' + key + ' = keyVal[len(keyVal)-1][\'value\']')
                      else :
                         exec('myparams[key] = ' + key + ' = {}')

        else :
            print ''

            myparams['IDI'] = IDI = self.parameters['IDI']
            myparams['OUTPUTIDI'] = OUTPUTIDI = self.parameters['OUTPUTIDI']
            myparams['DiFXinput'] = DiFXinput = self.parameters['DiFXinput']
            myparams['DiFXcalc'] = DiFXcalc = self.parameters['DiFXcalc']
            myparams['doIF'] = doIF = self.parameters['doIF']
            myparams['linAntIdx'] = linAntIdx = self.parameters['linAntIdx']
            myparams['Range'] = Range = self.parameters['Range']
            myparams['ALMAant'] = ALMAant = self.parameters['ALMAant']
            myparams['spw'] = spw = self.parameters['spw']
            myparams['calAPP'] = calAPP = self.parameters['calAPP']
            myparams['calAPPTime'] = calAPPTime = self.parameters['calAPPTime']
            myparams['APPrefant'] = APPrefant = self.parameters['APPrefant']
            myparams['gains'] = gains = self.parameters['gains']
            myparams['interpolation'] = interpolation = self.parameters['interpolation']
            myparams['gainmode'] = gainmode = self.parameters['gainmode']
            myparams['XYavgTime'] = XYavgTime = self.parameters['XYavgTime']
            myparams['dterms'] = dterms = self.parameters['dterms']
            myparams['amp_norm'] = amp_norm = self.parameters['amp_norm']
            myparams['XYadd'] = XYadd = self.parameters['XYadd']
            myparams['XYdel'] = XYdel = self.parameters['XYdel']
            myparams['XYratio'] = XYratio = self.parameters['XYratio']
            myparams['usePcal'] = usePcal = self.parameters['usePcal']
            myparams['swapXY'] = swapXY = self.parameters['swapXY']
            myparams['swapRL'] = swapRL = self.parameters['swapRL']
            myparams['feedRotation'] = feedRotation = self.parameters['feedRotation']
            myparams['correctParangle'] = correctParangle = self.parameters['correctParangle']
            myparams['IDI_conjugated'] = IDI_conjugated = self.parameters['IDI_conjugated']
            myparams['plotIF'] = plotIF = self.parameters['plotIF']
            myparams['plotRange'] = plotRange = self.parameters['plotRange']
            myparams['plotAnt'] = plotAnt = self.parameters['plotAnt']
            myparams['excludeAnts'] = excludeAnts = self.parameters['excludeAnts']
            myparams['excludeBaselines'] = excludeBaselines = self.parameters['excludeBaselines']
            myparams['doSolve'] = doSolve = self.parameters['doSolve']
            myparams['solint'] = solint = self.parameters['solint']
            myparams['doTest'] = doTest = self.parameters['doTest']
            myparams['npix'] = npix = self.parameters['npix']
            myparams['solveAmp'] = solveAmp = self.parameters['solveAmp']
            myparams['solveMethod'] = solveMethod = self.parameters['solveMethod']
            myparams['calstokes'] = calstokes = self.parameters['calstokes']
            myparams['calfield'] = calfield = self.parameters['calfield']


        result = None

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
        pathname="file:///home/marti/WORKAREA/LAUNCHPAD/PolConvert/polconvertsd/"
        trec = casac.casac.utils().torecord(pathname+'polconvert.xml')

        casalog.origin('polconvert')
        try :
          #if not trec.has_key('polconvert') or not casac.casac.utils().verify(mytmp, trec['polconvert']) :
            #return False

          casac.casac.utils().verify(mytmp, trec['polconvert'], True)
          scriptstr=['']
          saveinputs = self.__globals__['saveinputs']

          # Save .last file for this task execution. MPI servers don't write it (CASR-329).
          from mpi4casa.MPIEnvironment import MPIEnvironment
          do_full_logging = MPIEnvironment.is_mpi_disabled_or_client()
          if type(self.__call__.func_defaults) is NoneType:
              saveinputs=''
          else:
              saveinputs('polconvert', 'polconvert.last', myparams, self.__globals__,scriptstr=scriptstr, do_save_inputs=do_full_logging)

          tname = 'polconvert'
          spaces = ' '*(18-len(tname))
          casalog.post('\n##########################################'+
                       '\n##### Begin Task: ' + tname + spaces + ' #####')
          # Don't do telemetry from MPI servers (CASR-329)
          if do_full_logging and casa['state']['telemetry-enabled']:
              #casalog.poststat('Begin Task: ' + tname)
              task_starttime = str(datetime.datetime.now())
          if type(self.__call__.func_defaults) is NoneType:
              casalog.post(scriptstr[0]+'\n', 'INFO')
          else:
              casalog.post(scriptstr[1][1:]+'\n', 'INFO')

          # Effective call to the task as defined in gcwrap/python/scripts/task_*
          result = polconvert(IDI, OUTPUTIDI, DiFXinput, DiFXcalc, doIF, linAntIdx, Range, ALMAant, spw, calAPP, calAPPTime, APPrefant, gains, interpolation, gainmode, XYavgTime, dterms, amp_norm, XYadd, XYdel, XYratio, usePcal, swapXY, swapRL, feedRotation, correctParangle, IDI_conjugated, plotIF, plotRange, plotAnt, excludeAnts, excludeBaselines, doSolve, solint, doTest, npix, solveAmp, solveMethod, calstokes, calfield)

          if do_full_logging and casa['state']['telemetry-enabled']:
              task_endtime = str(datetime.datetime.now())
              casalog.poststat( 'Task ' + tname + ' complete. Start time: ' + task_starttime + ' End time: ' + task_endtime )
          casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####'+
                       '\n##########################################')

        except Exception, instance:
          if(self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__['__rethrow_casa_exceptions']) :
             raise
          else :
             #print '**** Error **** ',instance
             tname = 'polconvert'
             casalog.post('An error occurred running task '+tname+'.', 'ERROR')
             pass
        casalog.origin('')

        return result
#
#
#
#    def paramgui(self, useGlobals=True, ipython_globals=None):
#        """
#        Opens a parameter GUI for this task.  If useGlobals is true, then any relevant global parameter settings are used.
#        """
#        import paramgui
#        if not hasattr(self, "__globals__") or self.__globals__ == None :
#           self.__globals__=stack_frame_find( )
#
#        if useGlobals:
#            if ipython_globals == None:
#                myf=self.__globals__
#            else:
#                myf=ipython_globals
#
#            paramgui.setGlobals(myf)
#        else:
#            paramgui.setGlobals({})
#
#        paramgui.runTask('polconvert', myf['_ip'])
#        paramgui.setGlobals({})
#
#
#
#
    def defaults(self, param=None, ipython_globals=None, paramvalue=None, subparam=None):
        if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
        if ipython_globals == None:
            myf=self.__globals__
        else:
            myf=ipython_globals

        a = odict()
        a['IDI']  = ''
        a['OUTPUTIDI']  = ''
        a['DiFXinput']  = ''
        a['DiFXcalc']  = ''
        a['doIF']  = []
        a['linAntIdx']  = [1]
        a['Range']  = []
        a['ALMAant']  = ''
        a['spw']  = -1
        a['calAPP']  = ''
        a['calAPPTime']  = [0.,5.]
        a['APPrefant']  = ''
        a['gains']  = [["NONE"]]
        a['interpolation']  = []
        a['gainmode']  = []
        a['XYavgTime']  = 0.0
        a['dterms']  = ["NONE"]
        a['amp_norm']  = 0.01
        a['XYadd']  = {}
        a['XYdel']  = {}
        a['XYratio']  = {}
        a['usePcal']  = []
        a['swapXY']  = [False]
        a['swapRL']  = True
        a['feedRotation']  = []
        a['correctParangle']  = False
        a['IDI_conjugated']  = False
        a['plotIF']  = -1
        a['plotRange']  = []
        a['plotAnt']  = -1
        a['excludeAnts']  = []
        a['excludeBaselines']  = []
        a['doSolve']  = -1
        a['solint']  = [1,1]
        a['doTest']  = True
        a['npix']  = 50
        a['solveAmp']  = True
        a['solveMethod']  = 'gradient'
        a['calstokes']  = [1.,0.,0.,0.]
        a['calfield']  = -1


### This function sets the default values but also will return the list of
### parameters or the default value of a given parameter
        if(param == None):
                myf['__set_default_parameters'](a)
        elif(param == 'paramkeys'):
                return a.keys()
        else:
            if(paramvalue==None and subparam==None):
               if(a.has_key(param)):
                  return a[param]
               else:
                  return self.itsdefault(param)
            else:
               retval=a[param]
               if(type(a[param])==dict):
                  for k in range(len(a[param])):
                     valornotval='value'
                     if(a[param][k].has_key('notvalue')):
                        valornotval='notvalue'
                     if((a[param][k][valornotval])==paramvalue):
                        retval=a[param][k].copy()
                        retval.pop(valornotval)
                        if(subparam != None):
                           if(retval.has_key(subparam)):
                              retval=retval[subparam]
                           else:
                              retval=self.itsdefault(subparam)
                     else:
                        retval=self.itsdefault(subparam)
               return retval


#
#
    def check_params(self, param=None, value=None, ipython_globals=None):
      if ipython_globals == None:
          myf=self.__globals__
      else:
          myf=ipython_globals
#      print 'param:', param, 'value:', value
      try :
         if str(type(value)) != "<type 'instance'>" :
            value0 = value
            value = myf['cu'].expandparam(param, value)
            matchtype = False
            if(type(value) == numpy.ndarray):
               if(type(value) == type(value0)):
                  myf[param] = value.tolist()
               else:
                  #print 'value:', value, 'value0:', value0
                  #print 'type(value):', type(value), 'type(value0):', type(value0)
                  myf[param] = value0
                  if type(value0) != list :
                     matchtype = True
            else :
               myf[param] = value
            value = myf['cu'].verifyparam({param:value})
            if matchtype:
               value = False
      except Exception, instance:
         #ignore the exception and just return it unchecked
         myf[param] = value
      return value
#
#
    def description(self, key='polconvert', subkey=None):
        desc={'polconvert': '\n\nVersion 1.8.2 -- Converts VLBI visibilities polarization basis.',
               'IDI': 'Input FITS-IDI file with VLBI visibilities. It can also be a direcotry containing SWIN files from DiFX.',
               'OUTPUTIDI': 'Output FITS-IDI file (or SWIN directory). If equal to IDI, the file(s) will be overwritten',
               'DiFXinput': 'If SWIN files are being converted, this must be the *.input file used by DiFX.',
               'DiFXcalc': 'If SWIN files are being converted, this must be the *.calc file used by DiFX. This is optional, but if it is not provided, the cross-polarization gain estimates may be incorrect if doSolve>0.',
               'doIF': 'List of IFs to convert. Default means all.',
               'linAntIdx': 'List of indices of the linear-polarization antennas in the IDI file (lowest index starts with 1)',
               'Range': 'Time range to convert (integer list; AIPS format). Default means all data',
               'ALMAant': 'If ALMA has been used, this is the antenna table from the MS with the intra-ALMA visibilities.',
               'spw': 'Spectral window in ALMAvis that contains the VLBI band. If negative, the program will derive it automatically.',
               'calAPP': 'If ALMA has been used, this is the combined ASDM_CALAPPPHASE table from the ASDM. The list of measurement sets can also be given (so the table is concatenated from all of them).',
               'calAPPTime': 'Time shift and time tolerance (in sec) for the CALAPPPHASE table obtained from the ASDM.',
               'APPrefant': 'If not empty, re-reference the TelCal phases, assuming that the X-Y phase-difference table provided in \'gains\' (see keyword below) uses APPrefant as the reference antenna. Notice that the name of the gain table with the X-Y phase differences has to contain the string \'.XY0\'.',
               'gains': 'Gain tables to pre-calibrate the linear-pol VLBI stations (one list of gains per linear-pol station).',
               'interpolation': ' Interpolation type to use (one per calibration table). Tells whether to apply linear or nearest interpolation. Default is to use linear for all tables.',
               'gainmode': ' Mode of gain calibration to impose (one per calibration table). Default is \'T\' for all tables, unless either the string \'XY0\', \'bandpass\' or \'Gxyamp\' appears in the table name. The gain types can be either \'G\' (i.e., split gains per polarization) or \'T\' (i.e., combine polarizations).',
               'XYavgTime': ' Re-compute the G-mode gains by adding a time smoothing of X-Y differences to the T-mode gains. Default is NOT to do this (i.e., use the G-mode gains as given). If positive, use a running time average with this size (in seconds).',
               'dterms': 'D-term tables to pre-calibrate the linear-pol VLBI stations (one table per linear-pol station).',
               'amp_norm': 'If positive, normalize the amplitude correction to the X-Y average, and save the scaling factor (vs time) in an external (ASCII) file (ANTAB format, assuming a DPFU=amp_norm). If zero, or negative, apply the amplitude correction as is.',
               'XYadd': 'Add manually a phase between X and Y before conversion (in deg.). Either a list with one value per linear-pol station OR a list of lists (i.e., one value per IF for each antenna) OR a list of lists of lists (one value per channel, for each IF, for each linear-polarization antenna).',
               'XYdel': 'Add manually a multiband delay between X and Y before conversion (in deg./chanel). One value per linear-pol station.',
               'XYratio': 'Add manually an amplitude ratio between X and Y before conversion (R=X/Y). Follows the same format as XYadd. If a negative value is given for an antenna, the X/Y ratio will be estimated from its autocorrelations (the spectrum for antenna i will be computed using a running-median filter of width equal to -1/XYratio[i] of the IF bandwidth). If 0.0 is given for an antenna, the ratio will be estimated from the phasecal amplitudes (as long as usePcal is True).',
               'usePcal': 'List of booleans (one boolean per linear-polarization station). If True, use the X-Y difference of phasecal tones as an estimate of the X-Y cross-polarization phase. Default means to NOT use the phasecals.',
               'swapXY': 'Swap X-Y before conversion. One boolean per linear-pol VLBI station.',
               'swapRL': 'Swap R-L of the OTHER antenna(s) when plotting the fringes.',
               'feedRotation': 'Rotation angle of the feed (one value per antenna, in degrees). Default means zero degrees (so that X is truly horizontal for the linear-pol. antennas). These angles are used in the gain-solver step.',
               'correctParangle': 'If True, the correction for parallactic angle is applied to the converted antenna(s).',
               'IDI_conjugated': 'Assume a swap in the baseline defintion (i.e., conjugation) of the FITS-IDI file. This has NO effect on SWIN files. ',
               'plotIF': 'IF index(es) to plot. Default means to NOT plot. An empty list, [], means to plot ALL IFs being converted (but do not forget to set plotRange and plotAnt!).',
               'plotRange': 'Time range to plot (integer list; AIPS format). Default means to NOT plot',
               'plotAnt': 'Index of the other antenna in the baseline to plot. Default means to NOT plot.',
               'excludeAnts': 'List of antennas (i.e., list of antenna codenames) to NOT use in the cross-polarization gain estimates.',
               'excludeBaselines': 'List of baselines (i.e., a list of lists of two antenna codenames) to NOT use in the cross-polarization gain estimates.',
               'doSolve': 'If negative, do not estimate the cross-polarization gains. If positive or zero, estimate the gains using a Global Cross-Pol Fringe Fitting (GCPFF). The gains are fitted with an error function (Chi Square) defined as:\n\n sum( doSolve*(RR/LL-1)^2 + (RL^2 + LR^2) ),\n\n so that doSolve=0 minimizes the cross-hand polarizations (so it assumes a small linear polarization of the source), whereas doSolve>>1 assumes a negligible Stokes V.',
               'solint': 'If solint[0] null or negative, solve the cross-polarization phase plus a multi-band delay (MBD). If not, solve in bandpass mode by averaging solint[0] channels per solution.\n Divide the solution time range (per scan) in solint[1] chunks (i.e., in solint[1] subscans). I.e., if solint[1]==1, the fringes are fully averaged in time for each scan (but corrected for the scan fringe rates) before the GPLFF condition is computed. solint[2] is the minimum time jump (in seconds) to split the data into different scans (default: 100 seconds).',
               'doTest': 'If true, only compute (and eventually plot), the data, but leave OUTPUTIDI untouched.',
               'npix': 'Number of pixels for the fringe plots (and fringe search).',
               'solveAmp': 'if the cross-polarization gains are being estimated, solve also for the X/Y amplitude ratios.',
               'solveMethod': 'Method for the minimization of the Chi squared in the GCPFF. Can be \'gradient\', \'Levenberg-Marquardt\' or \'COBYLA\'.',
               'calstokes': 'Stokes parameters, [I,Q,U,V] of the calibrator (of course, this is only used if doSolve is not negative). The total intensity is not needed in the calibration (i.e., calstokes[0] can be left to 1.0, so that the other parameters will correspond to fractional polarization). ',
               'calfield': 'If not negative, field ID of the calibrator (useful if a time range covering several scans is being used in the GCPFF). If negative, use all data in the time range, regardless of the field ID.',

              }

        if(desc.has_key(key)) :
           return desc[key]

    def itsdefault(self, paramname) :
        a = {}
        a['IDI']  = ''
        a['OUTPUTIDI']  = ''
        a['DiFXinput']  = ''
        a['DiFXcalc']  = ''
        a['doIF']  = []
        a['linAntIdx']  = [1]
        a['Range']  = []
        a['ALMAant']  = ''
        a['spw']  = -1
        a['calAPP']  = ''
        a['calAPPTime']  = [0.,5.]
        a['APPrefant']  = ''
        a['gains']  = [["NONE"]]
        a['interpolation']  = []
        a['gainmode']  = []
        a['XYavgTime']  = 0.0
        a['dterms']  = ["NONE"]
        a['amp_norm']  = 0.01
        a['XYadd']  = {}
        a['XYdel']  = {}
        a['XYratio']  = {}
        a['usePcal']  = []
        a['swapXY']  = [False]
        a['swapRL']  = True
        a['feedRotation']  = []
        a['correctParangle']  = False
        a['IDI_conjugated']  = False
        a['plotIF']  = -1
        a['plotRange']  = []
        a['plotAnt']  = -1
        a['excludeAnts']  = []
        a['excludeBaselines']  = []
        a['doSolve']  = -1
        a['solint']  = [1,1]
        a['doTest']  = True
        a['npix']  = 50
        a['solveAmp']  = True
        a['solveMethod']  = 'gradient'
        a['calstokes']  = [1.,0.,0.,0.]
        a['calfield']  = -1

        #a = sys._getframe(len(inspect.stack())-1).f_globals

        if a.has_key(paramname) :
              return a[paramname]
polconvert_cli = polconvert_cli_()
