#########################
### MASTER SCRIPT FOR EU-VGOS POLCONVERSION.
### VERSION 23/07/2021.
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
mysteps = [0]





######################################

# Params about the dataset to process (USED BY ALL STEPS):
EXPNAME = 'ev0009'
DIFX_DIR = 'DiFX'

# LIST OF IFS (STARTING FROM 1 !!!):
DOIF = list(range(1,33))

# FREQUENCY AVERAGING FOR CROSS-POL GAINS:
CHAVG = 32

# DO WE ALSO SOLVE FOR AMPLITUDE GAINS ???
SOLVE_AMP = False


## AS NEW ANTENNAS ARE ADDED, SET HERE WHETHER THEIR PCALS
## ARE TO BE USED.
USE_PCAL = {'OE':True, 'OW':True, 'YJ':True, 'IS':True}
## NOTE: IF THERE ARE PROBLEMS WITH YJ (MULTI-FILE STILL
## UNDER ASSESSMENT), JUST SET ITS PCAL TO "False".


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
POLCAL_SCAN = ['38','17','24','31','45']

## SUFFIX TO ADD TO THE POLCONVERTED DIFX FILES:
SUFFIX = '_PC_v0'




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



# All auxiliary scripts will start with these lines:
Start = 'import pickle as pk\n'
Start += 'from EUVGOS_PY3 import SOURCE_SCANNER as Sscan\n'
Start += 'from EUVGOS_PY3 import POL_CALIBRATE as PolCal\n'
Start += 'from EUVGOS_PY3 import POLCONVERTER as PConv\n'
Start += 'from PolConvert import polconvert_standalone as PC\n'

# Loads the keywords for the execution of the corresponding task:
Start += 'IFF = open(\'keywords_%s.dat\',\'rb\'); kww = pk.load(IFF); IFF.close()\n'

# Place to save all intermediate CASA log files (can get VERY large!!)
if not os.path.exists('LOGS'):
  os.system('mkdir LOGS')  
currlogs = glob.glob('*.log')








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

  if len(list(filter(lambda x: 'POL_CALIBRATE' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      




  SCRIPT_NAME = 'STEP1'
  keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':DIFX_DIR, 'DOSCAN':POLCAL_SCAN, 'CHANSOL': CHAVG, 'USE_PCAL':USE_PCAL,'EXCLUDE_BASELINES':EXCLUDE_BASELINE,'DOAMP':SOLVE_AMP,'DOIF':DOIF,'PLOTANT':REFANT}

  keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys,protocol=0); keys.close()

  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print(Start%SCRIPT_NAME,file=OFF)
  print('PolCal.POL_CALIBRATE(**kww)',file=OFF)
  OFF.close()

  os.system(PYTHON_CALL%SCRIPT_NAME) 

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)      

  if os.path.exists('POL_CALIBRATE.FAILED'):
     raise Exception('STEP 1 FAILED!') 




# STEP 2: POL-CONVERT THE WHOLE EXPERIMENT.
if 2 in mysteps:

  if len(list(filter(lambda x: 'POLCONVERTER' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

    

  SCRIPT_NAME = 'STEP2'
  if type(POLCAL_SCAN) is list:
    SC0 = POLCAL_SCAN[0]
  else:
    SC0 = POLCAL_SCAN

  XYG = 'POLCAL_OUTPUT_SCAN-%s.dat'%SC0
  keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':DIFX_DIR, 'XYGAINS':XYG, 'SUFFIX': SUFFIX, 'USE_PCAL':USE_PCAL}
  keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys); keys.close()

  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print(Start%SCRIPT_NAME,file=OFF)
  print('PConv.POLCONVERTER(**kww)',file=OFF)
  OFF.close()

  os.system(PYTHON_CALL%SCRIPT_NAME) 

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)      

  if os.path.exists('POLCONVERTER.FAILED'):
     raise Exception('STEP 2 FAILED!') 





