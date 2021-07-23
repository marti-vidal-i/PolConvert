from __future__ import print_function
import numpy as np
import pylab as pl
import struct as stk
import glob
#from PolConvert import polconvert_standalone as PC
#from PolConvert import _XPCal as XP
import pickle as pk
import os, sys



if __name__ == '__main__':
    
 EXPNAME = 'vgt274'
 XYGAINS = 'PC_OUT_WPCAL_52.dat'
 SUFFIX = '_PC_v1'
 DIFX_DIR = 'DiFX'



#################################
# COMMENT OUT THIS LINE WHEN DEBUGGING WITH execfile(...)
def POLCONVERTER(EXPNAME, XYGAINS, DIFX_DIR, SUFFIX = '_PC', USE_PCAL=True):
 """ Converts all scans in a SWIN directory, using a cross-polarization
 gain file computed by PolConvert from a calibrator scan. It saves the new
 SWIN files in the same directory, adding SUFFIX to the *.difx subdirecotries."""
#################################

# try:
 if True:

######################
# Get scan info:
  IFF = open('SOURCES_%s.txt'%EXPNAME)
 
  lines = IFF.readlines()
  SCANS = []
  REFANTS = []
 
  IFF.close()
 
  OFF = open('%s.POLCONVERT.msg'%EXPNAME,'w')
 
  for li,line in enumerate(lines):
   if line.startswith(EXPNAME):
     SCANS.append(line.split()[0][:-1])
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
######################


# Get number of antennas and IFs:
  IFF = open(XYGAINS,'rb')
  XYG = pk.load(IFF)
  NAMS = list(XYG['XYadd'].keys())
  NANT = len(NAMS)
  At = list(XYG['XYadd'].keys())[0]
  NIF = len(XYG['XYadd'][At])
  IFF.close()

  doIF = range(NIF)

  if type(USE_PCAL) is bool:
    temp = bool(USE_PCAL)
    USE_PCAL = {}
    for ANT in NAMS:
      USE_PCAL[ANT] = bool(temp)
  

# Convert if there is valid data:
  for sci, scn in enumerate(SCANS):

   if REFANTS[sci]<0:
    print('WARNING! SCAN %s DOES NOT HAVE ANY VALID ANTENNA!'%scn,file=OFF)
  
   else:
    print('POLCONVERTING SCAN %s'%scn,file=OFF)

    DIFX = './%s/%s.difx'%(DIFX_DIR, scn)
    OUTPUT = './%s/%s.difx%s'%(DIFX_DIR, scn,SUFFIX)

    print('\nDOING SCAN: %s\n'%scn)

    INP = './%s/%s.input'%(DIFX_DIR, scn)
    CAL = './%s/%s.calc'%(DIFX_DIR, scn)
    OUT = './%s/PC_OUTPUT_%s.dat'%(DIFX_DIR, scn)

    command =  'import pickle as pk\n'
    command += 'import numpy as np\n\n'
 #   command += 'import sys\n'
 #   command += 'sys.path.append(\'/home/marti/WORKAREA/GITHUB\')\n'
    command += 'from PolConvert import polconvert_standalone as PC\n'
    command += 'IFF = open(\'%s\',\'rb\') ; XYG = pk.load(IFF); IFF.close()\n\n'%XYGAINS
#    command += 'TotalPhase = []; TotalAmp = []\n'
#    command += 'for i in range(%i):\n'%NANT
#    command += '  TotalPhase.append([]); TotalAmp.append([])\n'
#    command += '  for j in range(%i):\n'%NIF
#    command += '    TotalPhase[i].append(np.array(XYG[\'XYadd\'][i+1][j]))\n'
#    command += '    TotalAmp[i].append(np.abs(XYG[\'XYratio\'][i+1][j]))\n\n\n'

    USE_PCAL_STR = [] #"{" #+",".join(map(str,USE_PCAL))+"]" 
    for ANT in USE_PCAL.keys():
       USE_PCAL_STR.append('\'%s\':%s'%(ANT,USE_PCAL[ANT]))
    
    command +="MY_PCONV = PC.polconvert(IDI=\'%s\',\n"%DIFX
    command +="  OUTPUTIDI = \'%s\',\n"%OUTPUT 
    command +="  DiFXinput = \'%s\',\n"%INP
    command +="  DiFXcalc = \'%s\',\n"%CAL
    command +="  doIF = list(range(1,%i)),solveMethod = \'COBYLA\', solveAmp = False,\n"%(NIF+1)
    command +="  linAntIdx = [%s],swapXY=[%s],XYadd=XYG[\'XYadd\'],\n"%(','.join(['\'%s\''%NAM for NAM in NAMS]), ','.join(['False' for i in range(NANT)]))
    command +="  plotIF = [], Range=[],XYratio=XYG[\'XYratio\'], usePcal = {%s},\n"%(','.join(USE_PCAL_STR))
    command +="  correctParangle=True,doSolve=-1,doTest=False)\n"
  #  command +="  dterms = [%s],amp_norm=1.0,plotAnt=%i)\n"%(','.join(['\'NONE\'' for i in range(NANT)]), REFANTS[sci]+1)

    command +="OFF=open(\'%s\',\'wb\')\n"%OUT       #; DONE = True"
    command +="pk.dump(MY_PCONV,OFF,protocol=0) ; OFF.close()\n"


    comfile = open('./%s/AUXSCRIPT_%s.py'%(DIFX_DIR,scn),'w')
    print(command,file=comfile)
    comfile.close()
    
    os.system('python3 ./%s/AUXSCRIPT_%s.py'%(DIFX_DIR,scn))
    
   # IFF = open(OUT)
   # HOLA = pk.load(IFF)
   # IFF.close()

  OFF.close()



  if os.path.exists('POLCONVERTER.FAILED'):   
    os.system('rm -rf POLCONVERTER.FAILED')   

# except:

#  e = sys.exc_info()[0]  
#  OFF = open('POLCONVERTER.FAILED','w')
#  print(e,file=OFF)
#  OFF.close()



