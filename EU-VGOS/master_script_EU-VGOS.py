## EU-VGOS POLCONVERT MASTER SCRIPT.
#
# Copyright (c) Ivan Marti-Vidal 2018-2021 
#               Yebes Observatory (IGN, Spain)
#               Universitat de Valencia (Spain)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# a. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# b. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
# c. Neither the name of the author nor the names of contributors may 
#    be used to endorse or promote products derived from this software 
#    without specific prior written permission.
#
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#


###############################################
#### ANALYSIS PROCEDURE (TO BE IMPROVED/AUTOMATED):

## 0. Run the "source scanner" (i.e., mysteps=[0]) over the whole experiment.
#     Check the "SOURCES_*.txt" file, to select a good calibrator 
#     (and a set of calibrator scans where all antennas are participating).
#
## 1. Configure this script, i.e. set the variables "POLCAL_SCAN", "CHAVG",...
#     and "CALSOU" (according to the results from step 0). Then, run 
#     PolConvert in "solving mode" (i.e., mysteps=[1]). 
#     Check the quality of the solution (i.e., Fringe plots after conversion 
#     (in the *CHECK directories) and the "Cross-Gains*.png" figure.
#
## 2. Configure the script (variable "SUFFIX") and convert the whole
#     experiment (i.e., mysteps=[2]).
#
## 3. Correct the DiFX metadata (run the script "prepolconvert.py",
#     located in the "PP" directory). This step can also be done before
#     step 0, if desired.
###############################################




import os
import pickle as pk
import numpy as np
import glob
import matplotlib.pyplot as pl




######################################
#####################################################

thesteps = {0: 'Source scanner',
            1: 'Estimate cross-polarization gains',
            2: 'PolConvert the whole experiment'}

CASA_CALL = 'casa --nologger -c %s.py'

#####################################################
######################################

# List of steps to be performed:
mysteps = [0]





######################################
# Setup parameters (USED BY ALL STEPS):

# Location of the "EUVGOS_CASA" library:
EUVGOS_SCRIPT_DIR = '/home/marti/WORKAREA/POLCONVERT/EU-VGOS/'

# Experiment code (as it appears in the DiFX-related filenames):
EXPNAME = 'ev9203'

# Name of the data directory:
DIFX_DIR = 'DiFX'

## NOTE: THE DATA (i.e., DiFX directories and all the metadata, 
## like, e.g., the "*.input" files) MUST BE LOCATED IN "EXPNAME/DIFX_DIR"



# Whether to use (or not) the pulse-cal phase tones at each antenna:
USE_PCAL = [True, True, True, True]

# List of antenna (names) to be excluded from the solution estimate:
EXCLUDE_ANTENNA = []

# List of baselines to be excluded from the solution estimate 
# (normally, the intra-sites):
EXCLUDE_BASELINE = [ ['OE','OW'] ]


#########################################
# Params about the pol. calibrator (USED BY STEP 1):
# ALL scans here must have ALL the solvable antennas observing!!
# Otherwise, PolConvert will FAIL!
# Give either one number or a list of numbers (ideally, of the same source).
POLCAL_SCAN = ['074','113','148'] ## List of scans to use by the solver.
CHAVG = 128  # Channels to average (for each IF) in the solution.
CALSOU = 'OJ287'  # Calibrator source (i.e., the source observed in the "POLCAL_SCAN" scans).


#########################################
# Params about polconversion + export to CASA (USED BY STEP 2):
SUFFIX = '_PC_v0'  # Suffix to be added to all the polconverted DiFX files.



######################################
#####################################################


#################
# SCRIPT STARTS #
#################



# All auxiliary scripts will start with these lines:
Start = 'import sys\n'
Start += 'import pickle as pk\n'
Start += 'sys.path.append(\'%s\')\n'%EUVGOS_SCRIPT_DIR
Start += 'from EUVGOS_CASA import SOURCE_SCANNER as Sscan\n'
Start += 'from EUVGOS_CASA import POL_CALIBRATE as PolCal\n'
Start += 'from EUVGOS_CASA import POLCONVERTER as PConv\n'

# Loads the keywords for the execution of the corresponding task:
Start += 'IFF = open(\'keywords_%s.dat\'); kww = pk.load(IFF); IFF.close()\n'


# Place to save all intermediate CASA log files (can get VERY large!!)
if not os.path.exists('LOGS'):
  os.system('mkdir LOGS')  
currlogs = glob.glob('*.log')








# STEP 0: SCAN INFO (SOURCES AND SNR).
if 0 in mysteps:

  if len(filter(lambda x: 'SOURCE_SCANNER' not in x, glob.glob('*.FAILED')))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      


  SCRIPT_NAME = 'STEP0'
  keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':DIFX_DIR}
  keys = open('keywords_%s.dat'%SCRIPT_NAME,'w'); pk.dump(keyw, keys); keys.close()

  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print >> OFF, Start%SCRIPT_NAME
  print >> OFF, 'Sscan.SOURCE_SCANNER(**kww)'
  OFF.close()

  os.system(CASA_CALL%SCRIPT_NAME) 

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)      


  if os.path.exists('SOURCE_SCANNER.FAILED'):
     raise Exception('STEP 0 FAILED!') 



# STEP 1: FIND CROSS-POL GAINS FROM A CALIBRATOR SCAN.
if 1 in mysteps:

  if len(filter(lambda x: 'POL_CALIBRATE' not in x, glob.glob('*.FAILED')))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      




  SCRIPT_NAME = 'STEP1'
  keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':DIFX_DIR, 'DOSCAN':POLCAL_SCAN, 'CHANSOL': CHAVG, 'USE_PCAL':USE_PCAL,'EXCLUDE_BASELINES':EXCLUDE_BASELINE}
  keys = open('keywords_%s.dat'%SCRIPT_NAME,'w'); pk.dump(keyw, keys); keys.close()

  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print >> OFF, Start%SCRIPT_NAME
  print >> OFF, 'PolCal.POL_CALIBRATE(**kww)'
  OFF.close()

  os.system(CASA_CALL%SCRIPT_NAME) 

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)      

  if os.path.exists('POL_CALIBRATE.FAILED'):
     raise Exception('STEP 1 FAILED!') 




# STEP 2: POL-CONVERT THE WHOLE EXPERIMENT.
if 2 in mysteps:

  if len(filter(lambda x: 'POLCONVERTER' not in x, glob.glob('*.FAILED')))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

    

  SCRIPT_NAME = 'STEP2'
  if type(POLCAL_SCAN) is list:
    SC0 = POLCAL_SCAN[0]
  else:
    SC0 = POLCAL_SCAN

  XYG = 'POLCAL_OUTPUT_SCAN-%s.dat'%SC0
  keyw = {'EXPNAME':EXPNAME, 'DIFX_DIR':DIFX_DIR, 'XYGAINS':XYG, 'SUFFIX': SUFFIX, 'USE_PCAL':USE_PCAL}
  keys = open('keywords_%s.dat'%SCRIPT_NAME,'w'); pk.dump(keyw, keys); keys.close()

  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print >> OFF, Start%SCRIPT_NAME
  print >> OFF, 'PConv.POLCONVERTER(**kww)'
  OFF.close()

  os.system(CASA_CALL%SCRIPT_NAME) 

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)      

  if os.path.exists('POLCONVERTER.FAILED'):
     raise Exception('STEP 2 FAILED!') 





