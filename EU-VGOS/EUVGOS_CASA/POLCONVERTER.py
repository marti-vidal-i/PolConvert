## EU-VGOS POLCONVERT HELPER SCRIPT.
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



import numpy as np
import pylab as pl
import struct as stk
import glob
import _XPCal as XP
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

 try:

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
  IFF = open(XYGAINS)
  XYG = pk.load(IFF)
  NANT = len(XYG['XYadd'].keys())
  A1 = XYG['XYadd'].keys()[0]
  NIF = len(XYG['XYadd'][A1])
  IFF.close()

  doIF = range(NIF)

  if type(USE_PCAL) is bool:
      temp = bool(USE_PCAL)
      USE_PCAL = [temp for i in range(NANT)]
  

# Convert if there is valid data:
  for sci, scn in enumerate(SCANS):

   if REFANTS[sci]<0:
    print >> OFF, 'WARNING! SCAN %s DOES NOT HAVE ANY VALID ANTENNA!'%scn
  
   else:
    print >> OFF, 'POLCONVERTING SCAN %s'%scn

    DIFX = './%s/%s.difx'%(DIFX_DIR, scn)
    OUTPUT = './%s/%s.difx%s'%(DIFX_DIR, scn,SUFFIX)

    print '\nDOING SCAN: %s\n'%scn

    INP = './%s/%s.input'%(DIFX_DIR, scn)
    CAL = './%s/%s.calc'%(DIFX_DIR, scn)
    OUT = './%s/PC_OUTPUT_%s.dat'%(DIFX_DIR, scn)

    command =  'import pickle as pk\n'
    command += 'import numpy as np\n\n'
    command += 'IFF = open(\'%s\') ; XYG = pk.load(IFF); IFF.close()\n\n'%XYGAINS
#    command += 'TotalPhase = []; TotalAmp = []\n'
#    command += 'for i in range(%i):\n'%NANT
#    command += '  TotalPhase.append([]); TotalAmp.append([])\n'
#    command += '  for j in range(%i):\n'%NIF
#    command += '    TotalPhase[i].append(np.array(XYG[\'XYadd\'][i+1][j]))\n'
#    command += '    TotalAmp[i].append(np.abs(XYG[\'XYratio\'][i+1][j]))\n\n\n'

    USE_PCAL_STR = "["+",".join(map(str,USE_PCAL))+"]" 
    command +="MY_PCONV = polconvert(IDI=\'%s\',\n"%DIFX
    command +="  OUTPUTIDI = \'%s\',\n"%OUTPUT 
    command +="  DiFXinput = \'%s\',\n"%INP
    command +="  DiFXcalc = \'%s\',\n"%CAL
    command +="  doIF = range(1,%i),solveMethod = \'COBYLA\', solveAmp = False,\n"%(NIF+1)
    command +="  linAntIdx = [%s],swapXY=[%s],XYadd=XYG[\'XYadd\'],\n"%(','.join(map(str,range(1,NANT+1))), ','.join(['False' for i in range(NANT)]),  )
    command +="  plotIF = [], Range=[],XYratio=XYG[\'XYratio\'], usePcal = %s,\n"%(USE_PCAL_STR)
    command +="  ALMAant = \'\',spw=-1,calAPP=\'\',XYdel={},correctParangle=True,\n"
    command +="  gains = [%s], solint = [1,1],\n"%(','.join(['[\'NONE\']' for i in range(NANT)]))
    command +="  gainmode = [%s],XYavgTime=0.0,doSolve=-1,doTest=False,\n"%tuple(','.join(['[]' for i in range(NANT)]))
    command +="  dterms = [%s],amp_norm=1.0,plotAnt=%i)\n"%(','.join(['\'NONE\'' for i in range(NANT)]), REFANTS[sci]+1)

    command +="OFF=open(\'%s\',\'w\')\n"%OUT       #; DONE = True"
    command +="pk.dump(MY_PCONV,OFF) ; OFF.close()\n"
 
    comfile = open('./%s/AUXSCRIPT_%s.py'%(DIFX_DIR,scn),'w')
    print >> comfile, command
    comfile.close()
    
    os.system('casa --nologger -c ./%s/AUXSCRIPT_%s.py'%(DIFX_DIR,scn))
    
   # IFF = open(OUT)
   # HOLA = pk.load(IFF)
   # IFF.close()

  OFF.close()



  if os.path.exists('POLCONVERTER.FAILED'):   
    os.system('rm -rf POLCONVERTER.FAILED')   

 except:

  e = sys.exc_info()[0]  
  OFF = open('POLCONVERTER.FAILED','w')
  print >> OFF, e
  OFF.close()



