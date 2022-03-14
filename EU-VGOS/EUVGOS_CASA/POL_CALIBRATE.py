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

import pickle as pk
import os, sys
from polconvert_cli import polconvert_cli as polconvert



if __name__=='__main__':
 EXPNAME = 'vgt274'
 DIFX_DIR = 'DiFX'
 DOSCAN = 52
 NIF = 32
 CHANSOL = 32

 DO_RAW = False
 DO_WPCAL = True
 DO_FINAL = False



########################
# MODEL OF X-Y DIFFERENCE WITH A MB-DELAY:
 DO_WITH_DELAYS = False
 Delays = [[-1.095e-10, 9.40041e9], [-5.9e-11, -0.20843e9], [5.229e-10, 5.69512e9]]
 Phases= [0.,0.,0.]
########################




#################################
# COMMENT OUT THIS LINE WHEN DEBUGGING WITH execfile(...)
def POL_CALIBRATE(EXPNAME='',DIFX_DIR = '', DOSCAN=-1,CHANSOL=32,USE_PCAL=True,EXCLUDE_BASELINES=[],DOIF=[],DOAMP=True):
 """ Estimates cross-polarization gains using a scan in a SWIN directory.
  The channel resolution is set to CHANSOL. Saves the gains in a dictionary
  that can ge used by PolConvert."""
#################################

# try:
 if True:

  EXP = EXPNAME
  if DOSCAN < 0:
    raise Exception("POL_CALIBRATE ERROR! SCAN(S) NOT SPECIFIED!\n")

  DIFX = DIFX_DIR

  if not (type(DOSCAN) is list):
    DOSCAN = [DOSCAN]


# Figure out number of antennas, IFs and channels per IF:
  #temp = open('%s_%03i.input'%(os.path.join(DIFX,EXP),DOSCAN[0]))  
  temp = open('%s_%s.input'%(os.path.join(DIFX,EXP),DOSCAN[0]))
  lines = temp.readlines()
  for li in lines:
    if li.startswith('TELESCOPE ENTRIES:'):
      Nants = int(li.split()[-1])

    if li.startswith('FREQ ENTRIES:'):
      NIF = int(li.split()[-1])
    if li.startswith('NUM CHANNELS 0:'):
      Nchan = int(li.split()[-1])

  temp.close()
  del lines

  if type(USE_PCAL) is bool:
      temp = bool(USE_PCAL)
      USE_PCAL = [temp for i in range(Nants)]

  if len(DOIF)==0:
    DOIF = range(1,NIF+1)


############################################################
# X-Y cross-pol gain estimate (with phasecal correction):

  os.system('rm -rf %s'%os.path.join(DIFX,'POLCONVERT_CALIB_SCANS'))
  os.system('mkdir %s'%os.path.join(DIFX,'POLCONVERT_CALIB_SCANS'))
  for SI in DOSCAN:
    os.system('cp -r %s %s/.'%('%s_%s.difx'%(os.path.join(DIFX,EXP),SI), os.path.join(DIFX,'POLCONVERT_CALIB_SCANS')))  

  
  WITH_PCAL = polconvert(IDI='%s/POLCONVERT_CALIB_SCANS'%DIFX,
    OUTPUTIDI = '%s_POL_CALIBRATE_BLIND'%(os.path.join(DIFX,EXP)),
    DiFXinput = '%s_%s.input'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    DiFXcalc = '%s_%s.calc'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    doIF = DOIF, plotIF = DOIF,
    plotRange=[0,0,0,0,2,0,0,0], plotAnt=1,excludeBaselines=EXCLUDE_BASELINES,
    linAntIdx = range(1,Nants+1),swapXY=[False for i in range(Nants)], 
    usePcal = USE_PCAL, XYadd={}, XYratio={}, XYdel={},
# ALMA-specific:
    Range = [], ALMAant = '',spw=-1,calAPP='',
    gains = [['NONE'] for i in range(Nants)], 
    dterms = ['NONE' for i in range(Nants)],
    gainmode = [[] for i in range(Nants)],XYavgTime=0.0,amp_norm=1.0,
# Gain-solver configuration:    
    solveAmp = DOAMP, solveMethod = 'COBYLA', #'Nelder-Mead', #'COBYLA', #"Levenberg-Marquardt",
    doSolve=0.0,doTest=False, solint = [CHANSOL, 1])

 # raw_input('HOLD')

  FPK = 'FRINGE.PEAKS'
  FPL = 'FRINGE.PLOTS'
  PCF = 'POLCONVERT.FRINGE'
  Plot = 'Cross-Gains.png'

  os.system('rm -rf %s_POLCAL_%s %s_POLCAL_%s %s_POLCAL_%s'%(FPK, DOSCAN[0], FPL, DOSCAN[0], PCF, DOSCAN[0]))
  os.system('mv %s %s_POLCAL_%s'%(FPK, FPK, DOSCAN[0]))
  os.system('mv %s %s_POLCAL_%s'%(FPL, FPL, DOSCAN[0]))
  os.system('mv %s %s_POLCAL_%s'%(PCF, PCF, DOSCAN[0]))

  os.system('mv %s %s_POLCAL_%s.png'%(Plot, Plot, DOSCAN[0]))

  OFF=open('POLCAL_OUTPUT_SCAN-%s.dat'%DOSCAN[0],'w')
  pk.dump(WITH_PCAL,OFF) ; OFF.close()


  #IFF = open('POLCAL_OUTPUT_SCAN-%03i.dat'%DOSCAN)
  #WITH_PCAL = pk.load(IFF)
  #IFF.close()
  
  FINAL = polconvert(IDI='%s/POLCONVERT_CALIB_SCANS'%DIFX,
    OUTPUTIDI = '%s_POL_CALIBRATE_RESULTS'%(os.path.join(DIFX,EXP)),          
    DiFXinput = '%s_%s.input'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    DiFXcalc = '%s_%s.calc'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    doIF = DOIF, plotIF = DOIF,
    plotRange=[0,0,0,0,2,0,0,0], plotAnt=1,excludeBaselines=EXCLUDE_BASELINES,
    linAntIdx = range(1,Nants+1),swapXY=[False for i in range(Nants)], 
    usePcal = USE_PCAL, 
    XYadd=WITH_PCAL['XYadd'], XYratio=WITH_PCAL['XYratio'], XYdel={},
# ALMA-specific:
    Range = [], ALMAant = '',spw=-1,calAPP='',
    gains = [['NONE'] for i in range(Nants)], 
    dterms = ['NONE' for i in range(Nants)],
    gainmode = [[] for i in range(Nants)],XYavgTime=0.0,amp_norm=1.0,
# Gain-solver configuration:    
    doSolve=-1,doTest=False)

  os.system('rm -rf %s_CHECK_%s %s_CHECK_%s %s_CHECK_%s'%(FPK, DOSCAN[0], FPL, DOSCAN[0], PCF, DOSCAN[0]))
  os.system('mv %s %s_CHECK_%s'%(FPK, FPK, DOSCAN[0]))
  os.system('mv %s %s_CHECK_%s'%(FPL, FPL, DOSCAN[0]))
  os.system('mv %s %s_CHECK_%s'%(PCF, PCF, DOSCAN[0]))


  if os.path.exists('POL_CALIBRATE.FAILED'):   
    os.system('rm -rf POL_CALIBRATE.FAILED')   

# except:
 else:

  e = sys.exc_info()[0]  
  OFF = open('POL_CALIBRATE.FAILED','w')
  print >> OFF, e
  OFF.close()


