#########################
### MASTER SCRIPT FOR EU-VGOS POLCONVERSION.
### VERSION 29/09/2021.
### I. MARTI-VIDAL & J. GONZALEZ.
#########################



### PROCEDURE (NOT AUTOMATED YET):

### 0.- COPY THIS SCRIPT IN A DEDICATED DIRECTORY (E.G., THE PARENT
###     DIRECTORY OF THE DiFX DATA DIRECTORY CORRESPONDING TO THE
###     SESSION TO BE POLCONVERTED).

### 1.- SET THE EXPERIMENT NAME (EXPNAME) AND THE NAME OF THE 
###     DIRECTORY WITH THE DIFX DATA (DIFX_DIR).
###
### 2.- RUN THIS SCRIPT BY SETTING "mysteps = [0]".
###     IT WILL GENERATE THE FILE "SOURCES_*.txt", WHERE YOU
###     CAN CHOOSE THE BEST CALIBRATOR SCANS. USE THIS 
###     INFORMATION TO SET THE "POLCAL_SCAN" LIST.
###     (BEWARE OF THE LEADING "0" IN THE SCAN NAMES).
###
### 3.- RUN THIS SCRIPT BY SETTING "mysteps = [1]" AND "DOIF".
###     YOU MAY ALSO NEED TO SET CHAVG, SOLVE_AMP, etc.
###
### 4.- IF THE ASSESSMENT PLOTS ARE OK, SET "mysteps = [2]" 
###     AND RUN THIS SCRIPT. IT WILL POLCONVERT THE WHOLE 
###     SESSION.




######################################
#####################################################

thesteps = {0: 'Source scanner',
            1: 'Estimate cross-polarization gains',
            2: 'PolConvert the whole experiment'}

PYTHON_CALL = 'python3 %s.py'

#####################################################
######################################

# List of steps to be performed:
mysteps = [2]





######################################


# Number of processes allowed to run simultaneously:
NCPU = 11


# Params about the dataset to process (USED BY ALL STEPS):
EXPNAME = 'ev0287'

# Directory with the ORIGINAL data:
DIFX_DIR = 'DiFX_ORIG'

# Destionation directory of the POLCONVERTED data (difx + metadata):
PCONV_DIR = 'DiFX'


# LIST OF IFS (STARTING FROM 1 !!!):
DOIF = list(range(1,33))

# IF SOME AUTOCORRS ARE STORED IN DIFFERENT IFs, SET THIS TO THE OFFSET IF NUMBER:
IF_OFFSET = 32  # (CASE OF YJ, WHERE I --->  I + 32

# FREQUENCY AVERAGING FOR CROSS-POL GAINS:
CHAVG = 8

# WIDTH OF MEDIAN WINDOW FOR AUTOCORRELATIONS 
# (USED TO REMOVE PHASECALS):
PCAL_MED_WINDOW = 4

# PRE-AVERAGING TIME OF THE DATA (SECONDS):
INTTIME = 20.

# DO WE ALSO SOLVE FOR AMPLITUDE GAINS ???
SOLVE_AMP = True

# DO WE APPLY THE AMPLITUDES??
APPLY_AMP = True

# WHICH ALGORITHM DO WE USE??

#'Nelder-Mead', #'COBYLA', #"Levenberg-Marquardt"

#SOLVER = 'BFGS'                 ## NO SO GOOD
#SOLVER = 'Levenberg-Marquardt'  ## BROKEN
#SOLVER = 'gradient'             ## BROKEN
SOLVER = 'COBYLA'               ## GOOD
#SOLVER = 'Newton-CG'            ## BAD
#SOLVER = 'SLSQP'                ## BAD

## AS NEW ANTENNAS ARE ADDED, SET HERE WHETHER THEIR PCALS
## ARE TO BE USED:
USE_PCAL = {'OE':True, 'OW':True, 'YJ':True, 'IS':True,'WS':True}


## LIST OF FREQUENCY RANGES TO ZERO FOR THE PCALS (in MHz):
ZERO_PCALS = {'YJ':[[1000.,4000.]]}




## REFANT (FOR ASSESSMENT PLOTS):
REFANT = 'YJ'

## LIST OF ANTENNAS AND BASELINES TO EXCLUDE FROM THE
## GAIN SOLUTIONS. IN PRINCIPLE, ALL THE INTRA-SITE BASELINES
## MAY NOT BE USABLE:
EXCLUDE_ANTENNA = []
EXCLUDE_BASELINE = [['OE','OW']]

# SCANS TO USE FOR THE POL. CALIBRATION. IDEALLY, SHOULD BE SCANS 
# COVERING A RANGE OF PARALLACTIC ANGLES, WITH GOOD SNRs 
# ***AND*** WITH THE SAME ANTENNA CONFIGURATION!!!!
POLCAL_SCAN = ['030','031','079'] # Little difference if we add more scans, e.g.:
#POLCAL_SCAN = ['030','031','079','112','129','146']

#POLCAL_SCAN = ['267','266','268', '237', '166', '080']

## SUFFIX TO ADD TO THE ORIGINAL COPY OF THE DATA:
SUFFIX = ''




######################################
#####################################################


#################
# SCRIPT STARTS #
#################


import os
import pickle as pk
import numpy as np
import glob
import matplotlib.pyplot as pl
import multiprocessing




# All auxiliary scripts will start with these lines:
Start = 'import pickle as pk\n'
Start += 'from EUVGOS_PY3 import SOURCE_SCANNER as Sscan\n'
Start += 'from EUVGOS_PY3 import POL_CALIBRATE_MP as PolCal\n'
Start += 'from EUVGOS_PY3 import POLCONVERTER as PConv\n'
Start += 'from PolConvert import polconvert_standalone as PC\n'

# Loads the keywords for the execution of the corresponding task:
Start += 'IFF = open(\'keywords_%s.dat\',\'rb\'); kww = pk.load(IFF); IFF.close()\n'

# Place to save all intermediate CASA log files (can get VERY large!!)
if not os.path.exists('LOGS'):
  os.system('mkdir LOGS')  
currlogs = [] #glob.glob('*.log')








# STEP 0: SCAN INFO (SOURCES AND SNR).
if 0 in mysteps:

  if len(list(filter(lambda x: 'SOURCE_SCANNER' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      


  SCRIPT_NAME = 'STEP0'

  keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':DIFX_DIR}
  keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys, protocol=0); keys.close()

  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print(Start%SCRIPT_NAME,file=OFF)
  print('Sscan.SOURCE_SCANNER(**kww)',file=OFF)
  OFF.close()

  os.system(PYTHON_CALL%SCRIPT_NAME) 

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)      


  if os.path.exists('SOURCE_SCANNER.FAILED'):
     raise Exception('STEP 0 FAILED!') 











# STEP 1: FIND CROSS-POL GAINS FROM A CALIBRATOR SCAN.
if 1 in mysteps:

  os.system('rm -rf POLCAL_OUTPUT_SCAN-*_IF-*.dat')
  os.system('rm -rf FRINGE.PEAKS FRINGE.PLOTS POLCONVERT.FRINGE')

  CALDIR = os.path.join(PCONV_DIR,'%s_PC_CALIB'%EXPNAME)
  if not os.path.exists(PCONV_DIR):
    os.system('mkdir %s'%PCONV_DIR)
  os.system('rm -rf %s'%CALDIR)
  os.system('mkdir %s'%CALDIR)

  for SI in POLCAL_SCAN:
    os.system('cp -r %s %s/.'%('%s_%s.difx'%(os.path.join(DIFX_DIR,EXPNAME),SI),  CALDIR))  
    os.system('cp -r %s %s/.'%('%s_%s.calc'%(os.path.join(DIFX_DIR,EXPNAME),SI),  CALDIR))  
    os.system('cp -r %s %s/.'%('%s_%s.input'%(os.path.join(DIFX_DIR,EXPNAME),SI), CALDIR))
    
    pcals = glob.glob(os.path.join(CALDIR,'%s_%s.difx/PCAL*'%(EXPNAME,SI )))
  #  print('pcals: ',pcals)
    for pcal in pcals:
      for CURRIF in DOIF:
        os.system('cp %s %s_IF%i'%(pcal,pcal,CURRIF))

  NCPU = int(NCPU)
  if NCPU < 1:
    NCPU = multiprocessing.cpu_count() - 1


  if len(list(filter(lambda x: 'POL_CALIBRATE' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      


  SCRIPT_NAMES = []
  for CURRIF in DOIF:

    SCRIPT_NAME = 'STEP1_%i'%CURRIF
    keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':CALDIR, 'DOSCAN':POLCAL_SCAN, 'CHANSOL': CHAVG, 
            'USE_PCAL':USE_PCAL, 'INTTIME':INTTIME, 'EXCLUDE_BASELINES':EXCLUDE_BASELINE, 
            'DOAMP':SOLVE_AMP,'DOIF':[CURRIF],'PLOTANT':REFANT, 'APPLY_AMP':APPLY_AMP, 
            'SOLVER':SOLVER, 'APPLY_POLCAL':False, 'PCAL_SUFFIX':'_IF%i'%CURRIF, 'IF_OFFSET':IF_OFFSET}

    keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys,protocol=0); keys.close()

    OFF = open('%s.py'%SCRIPT_NAME,'w')
    print(Start%SCRIPT_NAME,file=OFF)
    print('PolCal.POL_CALIBRATE(**kww)',file=OFF)
    OFF.close()
    SCRIPT_NAMES.append(SCRIPT_NAME)

  def DO_PARALLEL(filename):
    print('GOING TO RUN %s'%filename)
    os.system(PYTHON_CALL%filename) 
 

  if NCPU>1:
    pool = multiprocessing.Pool(processes=NCPU)
    pool.map(DO_PARALLEL,SCRIPT_NAMES)
    pool.close()
    pool.join()
  else:
    for filename in SCRIPT_NAMES:
      DO_PARALLEL(filename)


 # input("HOLD")

  CALGAINS = {'XYratio':{}, 'XYadd':{},'PARAMETERS':{},'Frequency':{},'XYratioOriginal':{}}

  for CURRIF in DOIF:
    IFF=open('POLCAL_OUTPUT_SCAN-%s_IF-%i.dat'%(POLCAL_SCAN[0],CURRIF),'rb')
    IFGAIN = pk.load(IFF)    
    IFF.close()
    for key0 in CALGAINS.keys():
      for key1 in IFGAIN[key0].keys():
        if key1 not in CALGAINS[key0].keys():
          CALGAINS[key0][key1] = {}
        if type(IFGAIN[key0][key1]) is dict:
          for key2 in IFGAIN[key0][key1].keys():
            CALGAINS[key0][key1][key2] = IFGAIN[key0][key1][key2]
        else:
            CALGAINS[key0][key1] = IFGAIN[key0][key1]

  OFF = open('POLCAL_GAINS_%s.dat'%(EXPNAME),'wb')
  pk.dump(CALGAINS,OFF,protocol=0)
  OFF.close()


  os.system('rm -rf FRINGE.PLOTS/ALL_IFs*.png')
  os.system('rm -rf FRINGE.PLOTS/RL_LR*.png')

  os.system('rm -rf FRINGE.PEAKS')
  os.system('rm -rf FRINGE.PLOTS_NOCALIB')
  os.system('mv FRINGE.PLOTS FRINGE.PLOTS_NOCALIB')




## SECOND ROUND!!








#  input('HOLD')
  
  SCRIPT_NAME = 'STEP1B'
  XYG = 'POLCAL_GAINS_%s.dat'%(EXPNAME)
  keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':CALDIR, 'XYGAINS':XYG, 'SUFFIX': SUFFIX, 'IF_OFFSET':IF_OFFSET, 
          'USE_PCAL':USE_PCAL, 'DOPLOT':True, 'SCAN_LIST':[POLCAL_SCAN[0]], 'REFANT':REFANT}
  keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys); keys.close()
  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print(Start%SCRIPT_NAME,file=OFF)
  print('PConv.POLCONVERTER(**kww)',file=OFF)
  OFF.close()
  os.system(PYTHON_CALL%SCRIPT_NAME) 









  fig = pl.figure()
  MaxG = 0.0; NuMin = 1.0e20 ; NuMax = 0.0
  color = ['r','g','b','k','m','y','c']
  symbol = ['o','+','x']

  sub1 = fig.add_subplot(211)
  sub2 = fig.add_subplot(212,sharex=sub1)
 
  for CURRIF in DOIF:
    NUS = CALGAINS['Frequency'][CURRIF]/1.e9
    NuMin = np.min([NuMin,np.min(NUS)])
    NuMax = np.max([NuMax,np.max(NUS)])
    for antii,anti in enumerate(sorted(CALGAINS['XYadd'].keys())):
      toplotPh = np.array(CALGAINS['XYadd'][anti][CURRIF])
      toplotAp = np.array(CALGAINS['XYratioOriginal'][anti][CURRIF])
      MaxG = np.max([MaxG,np.max(toplotAp)])
      if CURRIF==DOIF[0]:
        sub1.plot(np.array(NUS),toplotPh,
                symbol[int(((antii)//len(color))%len(symbol))]+color[int((antii)%len(color))],
                label='ANT. '+str(anti))
        sub2.plot(np.array(NUS),toplotAp,
            symbol[int(((antii)//len(color))%len(symbol))]+color[int((antii)%len(color))],
            label='ANT. '+str(anti))
      else:
        sub1.plot(np.array(NUS),toplotPh,
                symbol[int(((antii)//len(color))%len(symbol))]+color[int((antii)%len(color))])
        sub2.plot(np.array(NUS),toplotAp,
            symbol[int(((antii)//len(color))%len(symbol))]+color[int((antii)%len(color))])


  sub1.set_ylim((-180.,180.))
#  sub2.set_ylim((0.,1.1*MaxG))
  pl.setp(sub1.get_xticklabels(),'visible',False)
  Dnu = NuMax-NuMin
  sub1.set_xlim((NuMin - Dnu*0.1, NuMax + Dnu*0.45))
  sub2.set_xlim((NuMin - Dnu*0.1, NuMax + Dnu*0.45))
  sub2.set_ylim((0.,2.5))


  sub1.legend(numpoints=1)
  sub1.set_ylabel('Cross-Phase (deg.)')
  sub2.set_ylabel('Cross-Amp (Norm.)')
  sub2.set_xlabel('Frequency (GHz)')
 
  fig.suptitle('XPOL GAINS %s'%EXPNAME)
  pl.savefig('Cross-Gains_%s_JOINED.png'%(EXPNAME))


  os.system('rm -rf FRINGE.PLOTS/ALL_IFs*.png')
  os.system('rm -rf FRINGE.PLOTS/RL_LR*.png')







  ANTS = [ai for ai in sorted(CALGAINS['XYadd'].keys()) if ai!=REFANT]
  NROW = 3
  Ncol = len(ANTS)//NROW


  fig = pl.figure(figsize=(8,3*Ncol))
  fig.subplots_adjust(wspace=0.01,hspace=0.01,right=0.98,top=0.90)
  fig.suptitle('FRINGE PEAKS. SCAN %s'%POLCAL_SCAN[0])

  for i in range(len(ANTS)):
    coli = (i+1)%NROW
    sub = fig.add_subplot(Ncol,NROW,i+1)
    ALLIF = glob.glob('FRINGE.PEAKS/FRINGE.PEAKS_IF*_%s*.dat'%ANTS[i])
    toplot = []
    MAX = 0.0
    for iffile in ALLIF:
      IFF = open(iffile)
      lines = IFF.readlines()
      IFF.close()
      REFANT = lines[0].split()[-1]
      currIF = float(lines[1].split()[2].replace('#','').replace('.',''))
      PEAK = float(lines[1].split()[-2])
      RR = float(lines[2].split()[1])*PEAK
      LL = float(lines[3].split()[1])*PEAK
      RL = float(lines[4].split()[1])*PEAK
      LR = float(lines[5].split()[1])*PEAK
      MAX = np.max([MAX,RR,LL,RL,LR])
      toplot.append([currIF,RR,LL,RL,LR])
    toplot = np.array(toplot)
    MAXIF = np.max(toplot[:,0])
    sub.plot(toplot[:,0],toplot[:,1]/MAX,'or')
    sub.plot(toplot[:,0],toplot[:,2]/MAX,'sb')
    sub.plot(toplot[:,0],toplot[:,3]/MAX,'xg')
    sub.plot(toplot[:,0],toplot[:,4]/MAX,'+k')
    sub.set_ylim((0.,1.2))
    sub.set_xlim((0.,MAXIF*1.2))
    pl.setp(sub.get_yticklabels(),'visible',False)
    if i < Ncol*(NROW-1):
      pl.setp(sub.get_xticklabels(),'visible',False)
    else:
      sub.set_xlabel('IF Number')

    if i==0:
      sub.plot([],[],'or',label='RR')
      sub.plot([],[],'sb',label='LL')
      sub.plot([],[],'xg',label='RL')
      sub.plot([],[],'+k',label='LR')
      pl.legend(numpoints=1,loc=1)

    pl.text(MAXIF*0.4,1.1,'%s-%s'%(ANTS[i],REFANT))
    
  pl.savefig('FRINGE.PLOTS/ALL_IF_PEAKS_SCAN_%s.png'%POLCAL_SCAN[0])




  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)    



  os.system('rm -rf keywords_STEP1_*.dat')
  os.system('rm -rf POLCAL_OUTPUT_SCAN*.dat')
  os.system('rm -rf STEP1_*.py')
  os.system('rm -rf STEP1B*.py')
  os.system('rm -rf PolConvert.XYGains_IF*.dat')
  os.system('rm -rf Cross-Gains_*_CALIB_IF*.png')
  


  if os.path.exists('POL_CALIBRATE.FAILED'):
     raise Exception('STEP 1 FAILED!') 








# STEP 2: POL-CONVERT THE WHOLE EXPERIMENT.
if 2 in mysteps:

  if len(list(filter(lambda x: 'POLCONVERTER' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      


  if not os.path.exists(PCONV_DIR):
    os.system('mkdir %s'%PCONV_DIR)


  NCPU = int(NCPU)
  if NCPU < 1:
    NCPU = multiprocessing.cpu_count() - 1
    

  SCRIPT_NAME = 'STEP2'
  XYG = 'POLCAL_GAINS_%s.dat'%(EXPNAME)


  if not APPLY_AMP:
    IFF = open(XYG,'rb')
    TEMP = pk.load(IFF)
    IFF.close()
    for anti in TEMP['XYratio'].keys():
      for ki in TEMP['XYratio'][anti].keys():
        TEMP['XYratio'][anti][ki][:] = 1.0
    OFF = open('POLCAL_GAINS_NOAMP_%s.dat'%(EXPNAME),'wb')
    pk.dump(TEMP,OFF)
    OFF.close()
    XYG = 'POLCAL_GAINS_NOAMP_%s.dat'%(EXPNAME)


  IFF = open('SOURCES_%s.txt'%EXPNAME)
 
  lines = IFF.readlines()
  SCANS = []
  REFANTS = []
 
  IFF.close()
 
  for li,line in enumerate(lines):
    if line.startswith(EXPNAME):
      TEMP = line.split()[0][:-1]
      SCANS.append(TEMP.split('_')[1])
      foundRef = False
      i = li+1
      while not foundRef:
        if lines[i].startswith(EXPNAME):
          REFANTS.append(-1)
          foundRef = True
        elif '+' in lines[i].split()[-1]:
          foundRef = True
          REFANTS.append(int(lines[i].split()[1][:-1]))
        else:
          i += 1
  
  for sci in range(len(SCANS)):
    if REFANTS[sci]<0:
      print('WARNING! SCAN %s DOES NOT HAVE ANY VALID ANTENNA!'%SCANS[sci])


  SCRIPT_NAMES = []


  for SCAN in SCANS:

    SCRIPT_NAME = 'STEP2_%s'%SCAN

   # os.system('cp -r %s %s'%(os.path.join(DIFX_DIR,'%s_%s*'%(EXPNAME,SCAN)), PCONV_DIR))

    keyw = {'EXPNAME':EXPNAME, 'ORIG_DIR':DIFX_DIR, 'DIFX_DIR':PCONV_DIR, 'XYGAINS':XYG, 'SUFFIX': SUFFIX, 'USE_PCAL':USE_PCAL,'SCAN_LIST':[SCAN],'ZERO_PCALS':ZERO_PCALS, 'IF_OFFSET':int(IF_OFFSET), 'AC_WINDOW':int(PCAL_MED_WINDOW)}
    keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys); keys.close()

    OFF = open('%s.py'%SCRIPT_NAME,'w')
    print(Start%SCRIPT_NAME,file=OFF)
    print('PConv.POLCONVERTER(**kww)',file=OFF)
    OFF.close()
    SCRIPT_NAMES.append(SCRIPT_NAME)


  def DO_PARALLEL(filename):

    print('GOING TO RUN %s'%filename)
    os.system(PYTHON_CALL%filename) 

  if NCPU>1:
    pool = multiprocessing.Pool(processes=NCPU)
    pool.map(DO_PARALLEL,SCRIPT_NAMES)
    pool.close()
    pool.join()
  else:
    for filename in SCRIPT_NAMES:
      DO_PARALLEL(filename)

  for filename in SCRIPT_NAMES:
    os.system('rm -rf %s.py'%filename)


  os.system('rm -rf keywords_STEP2_*.dat')

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)      

  if os.path.exists('POLCONVERTER.FAILED'):
     raise Exception('STEP 2 FAILED!') 





