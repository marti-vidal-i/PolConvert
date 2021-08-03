# Copyright (c) Ivan Marti-Vidal 2012-2021 
#               EU ALMA Regional Center. Nordic node.
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

from __future__ import absolute_import
from __future__ import print_function
__version__ = "1.9.0  "  # 7 characters
date = 'Mar 31 2021'     


################
# Import all necessary modules. 

########
# Execute twice, to avoid the silly (and harmless) 
# error regarding the different API versions of 
# numpy between the local system and CASA:
import os, sys
try: 
# mypath = os.path.dirname(os.path.realpath('__file__'))
# sys.path.append(mypath)
# import _PolConvert as PC
# goodclib = True
# print('\nC++ shared library loaded successfully (first try)\n')
#except:
# goodclib=False
# print('\n There has been an error related to the numpy') 
# print(' API version used by CASA. This is related to PolConvert')
# print(' (which uses the API version of your system) and should') 
# print(' be *harmless*.\n')

#mypath = os.path.dirname(os.path.realpath('__file__'))
#print mypath

#if not goodclib:
#  try: 
#   mypath = os.path.dirname(os.path.realpath('__file__'))
#   sys.path.append(mypath)
#   import _PolConvert as PC
#   goodclib = True
#   print('\nC++ shared library loaded successfully (2nd try)\n')
#  except:
#   goodclib=False
#   print('\nNo, the shared library did not load successfully\n')
#############

#if not goodclib:
#  try:
#   import _PolConvert as PC
#   goodclib = True
#   print('\nC++ shared library loaded successfully (3rd try)\n')
#  except:
#   goodclib=False
#   print('\nNo, the shared library still did not load successfully\n')
#############

### ... and we could keep on trying! (with an infinite loop! XD ).



import os,sys,shutil,re
import gc
import time
import struct as stk
import scipy.optimize as spopt
import numpy as np
import pylab as pl
import datetime as dt
import sys
import pickle as pk
import polconvert_standalone as PCONV

if sys.version_info.major < 3:
  try:
    from taskinit import *
    ms = gentools(['ms'])[0]
    tb = gentools(['tb'])[0]
  except Exception as ex:
    print('unable to load casa tools in python2\n\n')
    raise ex
else:
  try:
    # from taskutil import *
    from casatools import ms as ms_casa
    from casatools import table as tb_casa
    ms = ms_casa()
    tb = tb_casa()
  except Exception as ex:
    print('unable to load casa tools in python3\n\n')
    raise ex





def polconvert(IDI, OUTPUTIDI, DiFXinput, DiFXcalc, doIF, linAntIdx, Range, ALMAant, spw, calAPP, calAPPTime, APPrefant, gains, interpolation, gainmode, XYavgTime, dterms, amp_norm, XYadd, XYdel, XYratio, usePcal, swapXY, swapRL, feedRotation, correctParangle, IDI_conjugated, plotIF, plotRange, plotAnt,excludeAnts,excludeBaselines,doSolve,solint,doTest,npix,solveAmp,solveMethod, calstokes, calfield):


############################################

# this turns into the verbosity argument of _PolConvert.so
  print('Entered task_polconvert::polconvert()')
  DEBUG = False

  if 'POLCONVERTDEBUG' in os.environ:
     if os.environ['POLCONVERTDEBUG'] == 'True': DEBUG = True
     else:                                       DEBUG = False
  print('DEBUG setting is ' + str(DEBUG))
  print('__name__ is ' + __name__)


# Auxiliary function: print error and raise exception:
  def printError(msg):
    print(msg,'\n') 
    lfile = open("PolConvert.log","a")
    print('\n'+msg+'\n', file=lfile)
    lfile.close()
    sys.stdout.flush()
    raise Exception(msg)




# Auxiliary function: print message (terminal and log):
  def printMsg(msg, doterm=True, dolog=True):
    if doterm:
      print(msg)
    if dolog:
      lfile = open("PolConvert.log","a")
      print(msg,file=lfile)
      lfile.close()





# Auxiliary function: Geometric Median of a complex number:
  def geoMedian(Window, method='phasor'):

    WinData = np.array(Window)

 # Simplest approach (Amp & Phase separate). Assume NO WRAPS in phase:
    pAvg = np.median(np.abs(WinData))*np.exp(1.j*np.median(np.angle(WinData)))
    if method=='phasor': 
      return pAvg


 # A bit more complicated (point of minimum distance).
    elif method=='Torricelli':

      def L1Dist(p):
        return np.sum(np.abs(WinData-p[0]-1.j*p[1]))

      Torr = spopt.minimize(L1Dist, [pAvg.real, pAvg.imag], method='COBYLA')
      return Torr.x[0] + 1.j*Torr.x[1]








# Auxiliary function: Smooth the X-Y difference of G-mode gains
# using a running average:
  def XYsmooth(GAINTABLE, DTIME, SPW, IANT=-1):

    os.system('rm -rf %s.XYsmooth.PolConvert'%GAINTABLE)
    try:
      os.system('cp -r %s %s.XYsmooth.PolConvert'%(GAINTABLE,GAINTABLE))
      tb.open('%s.XYsmooth.PolConvert/ANTENNA'%GAINTABLE)
      GallAnts = tb.getcol('NAME')
      tb.close()

      if IANT>=0:
        allAnts = [IANT]
      else:
        allAnts = list(range(len(GallAnts)))

      tb.open('%s.XYsmooth.PolConvert'%GAINTABLE,nomodify=False)
      Gtimes = tb.getcol('TIME')
      Gants = tb.getcol('ANTENNA1')
      Gspws = tb.getcol('SPECTRAL_WINDOW_ID')
      Gflg = np.logical_not(tb.getcol('FLAG'))
      if np.shape(Gflg)[0] > 1:
        isT = 1
      else: 
        isT = 0
  #    print isT, np.shape(Gflg)
      Ggood = np.logical_and(Gflg[0,0,:],Gflg[isT,0,:])
      Mask = np.logical_and(Gspws == SPW, Ggood)
      Ggains = tb.getcol('CPARAM')

    except:
      printError('ERROR: Bad gain table %s!'%GAINTABLE)

   # Get the X-Y cross-phase gain:
    GDiff = Ggains[0,0,:]/Ggains[isT,0,:]  

   # Get the polarization-independent gain:
    TMode = Ggains[0,0,:]*Ggains[isT,0,:]



  # Smooth them:
    for iant in allAnts:
      Mask2 = np.where(np.logical_and(Mask,Gants==iant))[0]
      sys.stdout.write('\rSmoothing X-Y difference for antenna %s (%i of %i)   '%(GallAnts[iant],iant+1,len(GallAnts)))
      sys.stdout.flush()

      for tii,ti in enumerate(Mask2):
        Window = []
        tij = tii
        while tij>=0 and np.abs(Gtimes[ti]-Gtimes[Mask2[tij]])<DTIME/2.:
          Window.append(GDiff[Mask2[tij]])
          tij -= 1
        tij = tii+1
        while tij<len(Mask2) and np.abs(Gtimes[ti]-Gtimes[Mask2[tij]])<DTIME/2.:
          Window.append(GDiff[Mask2[tij]])
          tij += 1


  # Median (normalized):
        AvgDiff = geoMedian(Window)
        AvgDiff /= np.abs(AvgDiff)


        Ggains[0,0,ti] = np.sqrt(TMode[ti]*AvgDiff)
        Ggains[isT,0,ti] = np.sqrt(TMode[ti]/AvgDiff)


  # Update the table:
    tb.putcol('CPARAM',Ggains) 
    tb.close()         
    print('\nDONE!\n\n') 











# Auxiliary function: Re-reference XY-phases to 
# another ALMA refant, using the CalAPPPhase table:
  def ReReference(CALAPPPHASE,XY0,SPW,REFANT):

    printMsg("\n\n  GOING TO RE-REFERENCE X-Y PHASES TO %s.\n\n"%REFANT)

    DTMax = 180. # Minimum time gap (sec) to assume that TelCal has been reset.

# Figure out the baseband:
    tb.open('%s/SPECTRAL_WINDOW'%XY0)
    try:
      spname = tb.getcol('NAME')[SPW]
      BB =  [ii for ii in spname.split('#') if "BB_" in ii][0]
    except:
      printError("\n ERROR: BAD NAME FOR SPW %i IN TABLE %s\n"%(SPW,XY0))
    tb.close()


# Read TelCal's phases:
    tb.open(CALAPPPHASE)
    IMAX = np.argmax(tb.getcol('numPhasedAntennas'))
    ANT = list(tb.getcell('phasedAntennas',IMAX))

    try:
      REFIDX = ANT.index(REFANT)
    except:
      printError("\n\n ERROR: ANTENNA %s IS NOT IN CALAPP-PHASE TABLE!"%REFANT)

    ti = tb.getcol('startValidTime')
    bb = tb.getcol('basebandName')

    Tsort = np.argsort(ti)
    Torder = np.copy(ti[Tsort])
    UT = 24.*(Torder/86400.-int(Torder[0]/86400.))

# Arrange phases of the REFANT:
    Gains = []
    for i in range(len(ti)):
      aux = list(tb.getcell('phasedAntennas', rownr = i))
      try:
        REFI = aux.index(REFANT)
        aux2 = tb.getcell('phaseValues', rownr = i)
        NIF = tb.getcell('numChannels', rownr = i)
        Gains.append([aux2[NIF*REFI:NIF*(REFI+1)],aux2[NIF*(REFI+len(aux)):NIF*(REFI+len(aux)+1)]])
      except:
        printMsg("WARNING: ANTENNA %s NOT IN LIST OF PHASE ANTENNAS AT TIME %.1f.\n      THE RESULTING X-Y PHASES MAY BE *WRONG* AT THIS TIME!"%(REFANT,ti))
        Gains.append([np.zeros(NIF),np.zeros(NIF)])

    GainsA = np.array(Gains)


# Filter phases for the SPW:
    b1 = bb[Tsort] == BB

    Nchan = np.shape(GainsA)[-1]
    
    IX = np.copy(GainsA[Tsort,:][b1,0,:])
    IY = np.copy(GainsA[Tsort,:][b1,1,:])

    tb.close()
    Mask = np.where(np.logical_and(Torder[b1][1:]-Torder[b1][:-1]>DTMax,IX[:-1,0]!=0.0))

# Plot TelCal phases for REFANT:
    plfig = pl.figure()
    plsub = plfig.add_subplot(111)
    symb = ['or','og','ob','ok','om','oy','^r','^g','^b','^k','^m','^y']
    XYDiff = []
    AvXY = []
    for ni in range(Nchan):
      IntX = IX[1:,ni][Mask]
      IntY = IY[1:,ni][Mask]
      XYDiff.append(np.mod(IntX-IntY,2.*np.pi))

# We take the time MEDIAN as the best XY-phase estimate for the REFANT:
      AvXY.append(np.median(XYDiff[-1]))

      plsub.plot(UT[b1][1:][Mask],XYDiff[-1]*180./np.pi,symb[ni],label='CH %i'%(ni+1))
      printMsg('For antenna %s, Chan %i, found TelCal median X-Y phase of %.1f deg.'%(REFANT,ni+1,180./np.pi*AvXY[-1]))

    pl.legend(numpoints=1)
    plsub.set_xlabel('UT (h)')
    plsub.set_ylabel('X-Y phase (deg.)')
    plsub.set_xlim((0,24))
    pl.title('TelCal X-Y phases for new REFANT: %s'%REFANT)
    pl.savefig('%s.RE-REFERENCING.png'%CALAPPPHASE)


# Correct XY=phase table:
    printMsg("\n\n ADDING PHASES TO NEW XY0 TABLE\n\n")
    os.system("rm -rf %s.REFANT_%s"%(XY0,REFANT))
    os.system("cp -r %s %s.REFANT_%s"%(XY0,XY0,REFANT))

    tb.open("%s.REFANT_%s"%(XY0,REFANT),nomodify=False)
    spwi = tb.getcol('SPECTRAL_WINDOW_ID')
    XYData = []
    Mask = spwi==SPW

    for si in np.where(Mask)[0]:
      XYData.append(tb.getcell("CPARAM",si))
    XYData = np.array(XYData)

#    XYData = tb.getcol("CPARAM")
#    Mask = spwi==SPW

    NNu = np.shape(XYData)[1]
    NuPerChan = NNu/Nchan
    for ni in range(Nchan):
      XYData[0,ni*NuPerChan:(ni+1)*NuPerChan,Mask] *= np.exp(-1.j*(AvXY[ni]))

    for sii,si in enumerate(np.where(Mask)[0]):
      tb.putcell("CPARAM",si,XYData[sii])
   
#    tb.putcol("CPARAM",XYData)
    tb.close()

# Plot new vs. old XY-phases:
    Idx = np.where(Mask)[0][0]
    tb.open(XY0)
    XOrig = tb.getcell("CPARAM",Idx)[0,:]
    YOrig = tb.getcell("CPARAM",Idx)[1,:]
    tb.close()
    plfig.clf()
    plsub = plfig.add_subplot(111)
    plsub.plot(np.angle(XOrig/YOrig)*180./np.pi,'or',label='Original')

    tb.open("%s.REFANT_%s"%(XY0,REFANT))
    XOrig = tb.getcell("CPARAM",Idx)[0,:]
    YOrig = tb.getcell("CPARAM",Idx)[1,:]
    tb.close()
    plsub.plot(np.angle(XOrig/YOrig)*180./np.pi,'ob',label='Re-Referenced')

    pl.legend(loc=0,numpoints=1)
    plsub.set_ylim((-250,250))
    pl.savefig('%s_RE-REF_XYPHASES.png'%XY0)













# Auxiliary function: derive job label from DiFXinput
  def jobLabel(inputname):
    label = inputname
    try:
      label = re.sub('.input','', os.path.basename(label))
    except:
      pass
    return label





# Auxiliary function: unwrap phases for time interpolation
  def unwrap(phases, check=False):

    dims = np.shape(phases)
 #   print dims
    if dims[1]==0:  # Bandpass type
     for i in range(len(phases)-1):
      if phases[i+1]-phases[i] > np.pi:
        phases[i+1,:] -= 2.*np.pi
      elif phases[i+1]-phases[i] < -np.pi:
        phases[i+1,:] += 2.*np.pi
     if check:
       pl.figure()
       pl.plot(180./np.pi*phases)
    #   pl.show()
    #   raw_input('CHECK')


    elif dims[0]>1:  # Bandpass-gain type
 #    pl.plot(phases[:,0])
     for j in range(dims[1]):
      for i in range(dims[0]-1):
       if phases[i+1,j]-phases[i,j] > np.pi:
        phases[i+1:,j] -= 2.*np.pi
        printMsg('Adding phase wrap to gain at channel %i'%i)
       elif phases[i+1,j]-phases[i,j] < -np.pi:
        phases[i+1:,j] += 2.*np.pi
        printMsg('Removing phase wrap to gain at channel %i'%i)
     if check:
       pl.figure()
       pl.plot(180./np.pi*phases[:,0])
    #   pl.show()
    #   raw_input('CHECK')




# Auxiliary function: prepare the CALAPPPHASE, from a set of MSs:
  def makeCalAPP(mslist):
    printMsg('Will create CALAPPPHASE table from measurement sets.')
    os.system('rm -rf ./CALAPPPHASE.tab')
    for i,asd in enumerate(mslist):
      printMsg('Working out MS #%i - %s'%(i+1,asd))
      if i==0:
       os.system('cp -rf %s/ASDM_CALAPPPHASE ./CALAPPPHASE.tab'%asd)
      else:
       tb.open(os.path.join(asd,'ASDM_CALAPPPHASE'))
       tb.copyrows('./CALAPPPHASE.tab')
       tb.close()

    tb.open(os.path.join(asd,'ASDM_CALAPPPHASE'))
    time0 = tb.getcol('startValidTime')
    time1 = tb.getcol('endValidTime')
    tb.close()
    tb.open(asd)
    mstime = tb.getcol('TIME')
    tb.close()
    if len(time0)==0:
      printMsg('WARNING! NO APPCAL DATA FOR %s\n'%asd)
    else:
      print('\n\nMEAS SET %s:'%asd)
      tmin = np.min(time0)/86400.
      tmax = np.max(time1)/86400.
    mstmin = np.min(mstime)/86400.
    mstmax = np.max(mstime)/86400.
    printMsg('MS TIME RUNS FROM %8.5f TO %8.5f DAYS'%(mstmin,mstmax))
    if len(time0)>0:
     printMsg('FOR APP TABLE, TIME RUNS FROM %8.5f TO %8.5f DAYS'%(tmin,tmax))


    return timeranges


  tic = time.time()



#  greet  = ''' 
  ##########################################################################'
  # POLCONVERT --  version.                                                #'
  #       Please, add the POLCONVERT reference to your publications:       #'
  #                                                                        #'
  #          Marti-Vidal, Roy, Conway & Zensus 2016, A&A, 587, 143         #'
  #                                                                        #'
  ##########################################################################'
#  '''   

#  greetings = re.sub('version', __version__, greet)
#  printMsg(greetings,dolog=False)
#  printMsg('\n\nPOLCONVERT - VERSION %s'%__version__, doterm=False)



#########################################
# DO SOME SANITY CHECKS OF PARAMETERS

  try:
    doSolve = float(doSolve)
  except:
    printError("ERROR! doSolve should be a float!")
  

#  allMethods = ['gradient','Levenberg-Marquardt','COBYLA','Nelder-Mead']
#  scipyMethods = ['COBYLA','Nelder-Mead']
#  if solveMethod not in allMethods:
#    printError("ERROR! \'solveMethod\' must be any of: %s"%(', '.join(allMethods)))


#  if type(calstokes) is not list:
#    printError("ERROR! Wrong calstokes!")
#  elif len(calstokes)!=4 and doSolve>0.0:
#    printError("ERROR! calstokes should have 4 elements, not %d (%s)!" % (
#        len(calstokes), str(calstokes)))
#  for item in calstokes:
#    if type(item) is not float:
#      printError("ERROR! calstokes should only have float elements; got %s!" %
#         str(type(item)))

#  Stokes = list(calstokes)

#  if doSolve >= 0.0 and (Stokes[0]<=0. or 
#      Stokes[0]<np.sqrt(Stokes[1]**2.+Stokes[2]**2.+Stokes[3]**2.)):
#      printError("ERROR! Inconsistent Stokes parameters!")


# Will implement solveQU soon!
#  solveQU = False
#  calfield = -1

#  if type(solveQU) is not bool:
#    printError("ERROR! Wrong solveQU!")
#  if type(calfield) is not int:
#    printError("ERROR! Wrong calfield!")

#  if calfield>=0:
#    printMsg("Will use field #%i\n"%calfield)


#  doConj = True
#  if type(IDI_conjugated) is not bool:
#    printError("ERROR! IDI_cojugated should be a boolean!")
#  else:
#    doConj = IDI_conjugated

#  if type(swapRL) is not bool:
#    printError("ERROR! swapRL should be a boolean!")


#  if type(doIF) is not list:
#    printError("ERROR! doIF should be a list of integers!")
#  else:
#    for elem in doIF:
#      if type(elem) is not int:
#        printError("ERROR! doIF should be a list of integers!")
#      if elem==0:
#        printError("ERROR! IF numbers are given in AIPS (and FITS-IDI) convention\n(i.e., starting from 1; not from 0).\n")  

#  if type(doTest) is not bool:
#    printError("ERROR! doTest should be a boolean!")

#  if doTest:
#    printMsg("Will only compute, but not update the output file(s)",doterm=False)

#  if type(linAntIdx) is not list:
#    printError("ERROR! linAntIdx should be a list of antenna indices or names!")
#  for elem in linAntIdx:
#    if type(elem) not in [int, str]:
#      printError("ERROR! linAntIdx should be a list of antenna indices or names!")


#  if type(solint) is not list or (len(solint) not in [2,3]):
#    printError("ERROR! solint (%s) must be a list of two/three numbers!"%str(solint))
#  else:
#    try:
#      solint[0] = int(solint[0])
#      solint[1] = int(solint[1])
#    except:
#      printError("ERROR! solint (%s) must be a list of two/three numbers!"%str(solint))
#
#  if len(solint)==3:
#    solint[2] = float(solint[2])
#  else: # Default dt
#    solint.append(100.)
#
#  nALMA = len(linAntIdx)


#  if not os.path.exists(IDI):
#    printError("ERROR! IDI file (or SWIN folder) does not exist!")


  if type(calAPP) is list:
    try:
      makeCalAPP(calAPP)
      calAPP = './CALAPPPHASE.tab'
    except:
      printError('ERROR: Could not create nor interprete the CALAPPPHASE table!')

  c0 = len(calAPP) == 0 ; c1 = len(ALMAant) == 0
  if c0 != c1: 
    printError("ERROR! either both or none of calAPP and ALMAant need to be set!")

  if c0:
    printMsg("WARNING: calAPP and ALMAant are not set.")
    printMsg("ALMA CAL TABLES DO NOT SEEM TO BE USED, ACCORDING TO USER.")
    printMsg("This is not a problem in some cases...(isPhased now False).")
    isPhased=False
  else:
    if not os.path.exists(calAPP) or not os.path.exists(ALMAant):
      printError("ERROR! calAPP and/or ALMAant WERE NOT FOUND!")
    else:
      isPhased = True


#  if type(OUTPUTIDI) is not str:
#    printError("ERROR! OUTPUTIDI should be a string!")

#  if type(plotIF) is int:
#    if plotIF >0:
#      plotIF = [plotIF]
#    else:
#      plotIF = []
#  for pli in plotIF:
#    if type(pli) is not int:
#      printError("ERROR! plotIF should be an integer or a list of integers!")
#
#  for pli in plotIF:
#    if pli not in doIF:
#      printError("ERROR! Only converted IFs can be plotted!")



  try:
    spw = int(spw)
  except:
    printError("ERROR! spw should be an integer!")



  if len(gains)!= nALMA or len(dterms)!= nALMA:
    printError("Invalid format for gains and/or dterms!\n Should be lists as large as the number of linear-polarization VLBI stations!")
  

# Sanity check for interpolation:
  if type(interpolation) is not list:
    printError("Interpolation must be a list (or a list of lists)!")

  if len(interpolation)==0:
    interpolation = [[] for i in range(nALMA)]

  if len(interpolation)!= nALMA:
    printError("Wrong length for interpolation!")

  for i,intype in enumerate(interpolation):
    if type(intype) is not list:
      printError("Interpolation must be a list (or a list of of lists)!")
    if len(intype)==0:
      interpolation[i] = ['linear' for g in gains[i]]
      intype = interpolation[i]
    if len(intype) != len(gains[i]):
      printError("Interpolation must have the same dimensions as gains!")
    for ints in intype:
      if ints not in ['linear','nearest']:
        printMsg("integration type " + ints + " requested.")
        printError("Only \'linear\' and \'nearest\' interpolations are supported!")


  try:
    XYavgTime = float(XYavgTime)

  except:
    printError("XYavgTime must be a positive (or zero) double!")

  if XYavgTime < 0.0:
    printError("XYavgTime must be positive or zero!")



  if type(gainmode) is not list:
    printError("gainmode must be a list (or a list of lists)!")

  if len(gainmode)==0:
    gainmode = [[] for i in range(nALMA)]

  if len(gainmode)!= nALMA:
    printError("Wrong length for gainmode!")

  for i,gtype in enumerate(gainmode):
    if type(gtype) is not list:
      printError("gainmode must be a list (or a list of of lists)!")
    if len(gtype)==0:
      gainmode[i] = [{True:'G',False:'T'}['XY0' in g or 'bandpass' in g or 'Gxyamp' in g] for g in gains[i]]
      gtype = gainmode[i]
    if len(gtype) != len(gains[i]):
      printError("gainmode must have the same dimensions as gains!")
    for ints in gtype:
      if ints not in ['G','T','S']:
        printMsg("Gain type " + ints + " requested.")
        printError("Only \'G\', \'S\' and \'T\' interpolations are supported!")

  for gi in range(len(gains)):
    for gii in range(len(gains[gi])):
      printMsg('Will calibrate with table %s in %s mode.'%(os.path.basename(gains[gi][gii]),gainmode[gi][gii]))

# Sanity check for gains and dterms:
  for g in gains:
    if type(g) is not list:
      printError("Invalid format for gains!\n Each element should be a list with the names of the gain tables for the ith linear-pol. VLBI station") 
    else:
      for elem in g:
        if type(elem) is not str:
          printError("Invalid format for gains!\n Each element (of each element) should be a string (the name of a calibration table)") 

  for elem in dterms:
    if type(elem) is not str:
      printError("Invalid format for dterms!\n Each element should be a string (the name of a calibration table)") 



#  if type(usePcal) is not list:
#    printError("Invalid format for usePcal! Should be a list of booleans!\n")
#  elif len(usePcal)==0:
#    usePcal = [False for i in range(nALMA)]
#  else:
#    for pi in usePcal:
#      if type(pi) is not bool:
#        printError("Invalid format for usePcal! " + 
#            "It should be a list of booleans!\n")

#  if len(np.where(usePcal)[0]) > 0:
#    isPcalUsed = True
#    printMsg("Info: Pcal used in %s" % str(np.where(usePcal)[0])) 
#  else:
#    isPcalUsed = False
#    printMsg("Info: Pcal is not in use")


#  if len(swapXY) != nALMA:
#    printError("Invalid format for swapXY!\n Should be a list of booleans, as large as the number of linear-polarization VLBI stations!")
#  for sxy in swapXY:
#      if type(sxy) is not bool:
#       printError("Invalid format for swapXY!\n Should be a list of booleans, as large as the number of linear-polarization VLBI stations!")


#  if len(Range) not in [0,8]:
#    printError("Invalid format for Range! Should be either an empty list or a list of 8 integers!")
#  for rr in Range:
#    if type(rr) is not int:
#      printError("Invalid format for Range! Should be either an empty list or a list of 8 integers!")

#  if len(plotRange) not in [0,8]:
#    printError("Invalid format for Range! Should be either an empty list or a list of 8 integers!")
#  for rr in plotRange:
#    if type(rr) is not int:
#      printError("Invalid format for Range! Should be either an empty list or a list of 8 integers!")


  if len(calAPPTime) != 2:
    printError("Bad format for calAPPTime. Should be a list of 2 floats!")

  try:
    CALAPPTSHIFT, CALAPPDT = [float(cc) for cc in calAPPTime]
  except:
    printError("Bad format for calAPPTime. Should be a list of 2 floats!")


#########################################




######################
# READ IDs OF ALMA ANTENNAS USED IN THE PHASING:

  if isPhased:

#   print calAPP
#   os.system('ls %s'%calAPP)

   success = tb.open(calAPP)

   if not success:
    printError('ERROR: INVALID calAPP TABLE!')
    
   try:

    time0 = tb.getcol('startValidTime')-CALAPPDT+CALAPPTSHIFT
    time1 = tb.getcol('endValidTime')+CALAPPDT+CALAPPTSHIFT



#########
# Figure out times with unphased data:
    ADJ = tb.getcol('adjustToken')
    BBN = tb.getcol('basebandName')
    BBC = ['BB_%i'%i for i in range(1,5)]
    sc2flag = [[] for i in range(4)] # Currently, we assume 1 spw per BBC in the APP mode.
    scgood = [[] for i in range(4)] # Currently, we assume 1 spw per BBC in the APP mode.
    timeranges = [[] for i in range(4)] # Currently, we assume 1 spw per BBC in the APP mode.

    for j in range(4):
      sc2flag[j] = [i for i in range(len(ADJ)) if BBN[i]==BBC[j] and ('PHASE_UPDATED' not in ADJ[i]) and ('PHASE_NOCHANGE' not in ADJ[i])]
      scgood[j] = [i for i in range(len(ADJ)) if BBN[i]==BBC[j] and ('PHASE_UPDATED' in ADJ[i] or 'PHASE_NOCHANGE' in ADJ[i])]

    SUBSCANDUR = time1[sc2flag[0][0]] - time0[sc2flag[0][0]]
    sec = 1.

    for j in range(4):
     if len(sc2flag[j])>0:
      for si in sc2flag[j]:
        timerange = [time1[si]-SUBSCANDUR-sec, time1[si]+sec]
        if timerange not in timeranges[j]:
          timeranges[j].append(timerange)

     for si in scgood[j]:
      timerange = [time1[si]-SUBSCANDUR-sec,time0[si]+sec]
      if si-1 in sc2flag[j]:  
        if timerange not in timeranges[j]:
          timeranges[j].append(timerange)

# For testing:
#    timeranges.append([4998635280., 4998635350.])
    timerangesArr = [np.array(ti) for ti in timeranges]
#########


    nphant = tb.getcol('numPhasedAntennas')
    refs = tb.getcol('refAntennaIndex')
    asdmtimes = [time0,time1]
    phants = []
    for i in range(len(nphant)):
      phants.append(list(tb.getcell('phasedAntennas', rownr = i)))

    tb.close()

   except:
#   else:
    printError('ERROR: INVALID calAPP TABLE CONTENT!')

   success = tb.open(ALMAant)
   if not success:
    printError('ERROR: NO VALID ANTENNA TABLE FOUND!')

   allants = list(tb.getcol('NAME'))
   allantidx = list(range(len(allants)))
   tb.close()

   refants = np.array([allants.index(phants[t][refs[t]]) for t in range(len(phants))])

   antimes = [[[0.,0.]] for ian in allants]
   for ian,antenna in enumerate(allants):
    for j,phan in enumerate(phants):
      if antenna in phan:
        antimes[ian].append([time0[j]-CALAPPDT+CALAPPTSHIFT,time1[j]+CALAPPDT+CALAPPTSHIFT])

    auxarr = np.zeros((len(antimes[ian]),2))
    auxarr[:] = antimes[ian]
    antimes[ian] = auxarr
    
   nphtimes = [len(anti) for anti in antimes]


  else:

# Not a phased array. Dummy values for phasing info:
   allants = [1]
   allantidx = [0]
   antimes = [np.array([0.,1.e20])]
   nphtimes = [1]
   refants = np.array([0],dtype=np.int)
   asdmtimes = [np.array([0.]),np.array([1.e20])]
   timerangesArr = [np.array([0.,0.]) for i in range(4)]

######################



#######
# CHECK IF THIS IS A FITS-IDI FILE OR A SWIN DIR.:
  if os.path.isdir(IDI):
    isSWIN = True
#    printMsg('\n\nYou have asked to convert a set of SWIN files.')
#    if len(DiFXinput)==0 or not os.path.exists(DiFXinput) or not os.path.isfile(DiFXinput):
#      printError("Invalid DiFX input file! %s"%DiFXinput)
#    printMsg('Opening calc file... %s' % DiFXcalc)
#    try:
#      printMsg('Opening "%s"' % (DiFXcalc))
#      antcoords = []  ; soucoords = [[],[]]; antmounts = []; antcodes = [];
#      calc = open(DiFXcalc)
#      lines = calc.readlines()
#      calc.close()
#      printMsg('Read %d lines from %s' % (len(lines), DiFXcalc))
#      for ii, line in enumerate(lines):
#        if 'TELESCOPE' in line and 'NAME' in line:
#          antcodes.append(line.split()[-1])
#        if 'TELESCOPE' in line and 'X (m):' in line:
#    #      printMsg(line.rstrip())
#          antcoords.append(list(map(float,[ll.split()[-1] for ll in lines[ii:ii+3]])))
#          printMsg('TELESCOPE %s AT X: %.3f ; Y: %.3f ; Z: %.3f'%tuple([antcodes[-1]] + antcoords[-1]))
## CURRENTLY, ONLY ALTAZ MOUNTS SUPPORTED FOR SWIN FILES:
#          antmounts.append(0)
#
#        if 'SOURCE' in line and ' NAME: ' in line:
#          SNAM = line.split()[-1]
#          soucoords[0].append(float(lines[ii+1].split()[-1]))
#          soucoords[1].append(float(lines[ii+2].split()[-1]))
#          printMsg('SOURCE %s AT RA: %.8f rad, DEC: %.8f rad'%(SNAM,soucoords[0][-1],soucoords[1][-1]))
#      antcoords = np.array(antcoords,dtype=np.float)
#      antmounts = np.array(antmounts)
#      soucoords[0] = np.array(soucoords[0],dtype=np.float)
#      soucoords[1] = np.array(soucoords[1],dtype=np.float)
#      printMsg('done parsing calc')
#    except Exception as ex:
#      printMsg(str(ex))
#      printMsg(("WARNING! Invalid DiFX calc file '%s'!\n" + 
#        "PolConvert may not calibrate properly.") % DiFXcalc)
#    if len(antmounts)==0:
#      printError('ERROR! NO ANTENNAS FOUND IN CALC FILE!')
#    else:
#      printMsg('There are %i antennas.'%len(antmounts))
  elif os.path.isfile(IDI):
    isSWIN = False
#    printMsg('\n\nYou have asked to convert a FITS-IDI file.')
#    printMsg('Reading array geometry...')
#    try:
#      import pyfits as pf
#      ffile = pf.open(IDI)
#      for ii,group in enumerate(ffile):
#        if group.name == 'ARRAY_GEOMETRY':
#          grarr = ii
#        elif group.name == 'SOURCE':
#          grsou = ii
#
#
#      raappUnit = ffile[grsou].data.columns[ [coln.name for coln in ffile[grsou].data.columns].index('RAAPP') ].unit
#      decappUnit = ffile[grsou].data.columns[ [coln.name for coln in ffile[grsou].data.columns].index('DECAPP') ].unit
#
#      soucoords = [np.array(ffile[grsou].data['RAAPP'],dtype=np.float),np.array(ffile[grsou].data['DECAPP'],dtype=np.float)]
#
#      if raappUnit=='DEGREES':
#        soucoords[0] *= np.pi/180.
#
#      if decappUnit=='DEGREES':
#        soucoords[1] *= np.pi/180.
#
#      antcodes = [ff[:2] for ff in ffile['ANTENNA'].data['ANNAME']]
#      ffile.close()
#
## THESE LINES FAIL IF ORBPARM IS PRESENT IN ARRAY GEOMETRY!
##      antcoords = np.array(ffile[grarr].data['STABXYZ'],dtype=np.float)
##      antmounts = np.array(ffile[grarr].data['MNTSTA'],dtype=np.float)
#      import _getAntInfo as gA
#      
#      success = gA.getAntInfo(IDI)
#      if success != 0:
#        printError("ERROR GETTING FITS-IDI METADATA! ERR: %i"%success)
#      else:
#        antcoords = gA.getCoords()
#        antmounts = gA.getMounts()
#    # FIXME: need a way to find antenna codes if above fails:   
#    #        antcodes = ['%02i'%i for i in range(1,len(antmounts)+1)]
#    except:
#      printMsg('WARNING! This FITS-IDI file has missing information!\nPolConvert may not calibrate properly.')
#  else:
#    printError("Invalid input data!") 



######




######
# IF THIS IS A SWIN DATASET, READ THE INPUT FILE INTO 
# A METADATA LIST:

  if isSWIN:
    printMsg('Reading the DiFX input file\n')
    ifile = open(DiFXinput)
    inputlines = ifile.readlines()
    ifile.close()
    FreqL = [inputlines.index(l) for l in inputlines if 'FREQ TABLE' in l]


#    if len(antcoords) == 0:
#      Nteles = [inputlines.index(l) for l in inputlines if 'TELESCOPE ENTRIES' in l][0]
#      antcoords = np.ones((Nteles,3),dtype=np.float)
#      antmounts = np.zeros(Nteles,dtype=np.int)
#    if len(soucoords[0])==0:
#      soucoords = [np.zeros(1,dtype=np.float),np.zeros(1,dtype=np.float)]

# ONLY ONE FREQ TABLE IS ALLOWED:
    try:
      fr = FreqL[0]
      Nfreq = int(list(filter(
        lambda x: 'FREQ ENTRIES' in x, inputlines[fr+1:]))[0].split()[-1])
      Nr = list(range(Nfreq))
    except Exception as ex:
      printMsg(str(ex))
      printError("BAD input file!")

    FrInfo = {'FREQ (MHZ)':[0. for i in Nr], 'BW (MHZ)':[0. for i in Nr], 
              'SIDEBAND':['U' for i in Nr], 'NUM CHANNELS':[1 for i in Nr],
              'CHANS TO AVG':[1 for i in Nr], 'OVERSAMPLE FAC.':[1 for i in Nr], 
              'DECIMATION FAC.':[1 for i in Nr], 'SIGN' :[1. for i in Nr]}

# READ METADATA FOR ALL FREQUENCIES:
    for entry in FrInfo.keys():
     for line in inputlines[fr+1:]:
       if entry in line:
         index = int((line.split(':')[0]).split()[-1])
         FrInfo[entry][index] = type(FrInfo[entry][0])(line.split()[-1])

# SORT OUT THE CHANNEL FREQUENCIES:

    if len(doIF)==0:
      doIF = list(range(1,len(Nr)+1))

    metadata = []
    IFchan = 0
    for nu in Nr:
      nu0 = FrInfo['FREQ (MHZ)'][nu] 
      bw = FrInfo['BW (MHZ)'][nu]
      nchan = FrInfo['NUM CHANNELS'][nu]
# MAX. NUMBER OF CHANNELS:
      chav = FrInfo['CHANS TO AVG'][nu]
      if nu in doIF:
        IFchan = max([IFchan,int(nchan/chav)])
      sb = {True: 1.0 , False: -1.0}[FrInfo['SIDEBAND'][nu] == 'U']
      FrInfo['SIGN'][nu] = float(sb)
      freqs = (nu0 + np.linspace((sb-1.)/2.,(sb+1.)/2.,nchan/chav,endpoint=False)*bw)*1.e6
      if float(nchan//chav) != nchan/chav:
        printMsg("linspace check chan: %d / %d = %f" %
            (nchan, chav, nchan/chav))
      freqs = (nu0 + np.linspace((sb-1.)/2.,(sb+1.)/2.,
        nchan//chav,    # should be exactly divisible
        endpoint=False)*bw)*1.e6
      metadata.append(freqs)



#####
 

  else:

# READ FREQUENCY INFO TO HELP SELECTING THE SPW AUTOMATICALLY:
    import pyfits as pf
    fitsf = pf.open(IDI)
    nu0 = fitsf['FREQUENCY'].header['REF_FREQ']
    bw = fitsf['FREQUENCY'].header['CHAN_BW']
    nch = fitsf['FREQUENCY'].header['NO_CHAN'] #*fitsf['FREQUENCY'].header['NO_BAND']
    IFchan = nch
    Nr = fitsf['FREQUENCY'].header['NO_BAND']
    sgn = {True:1.0,False:-1.0}[bw>0.0]
    FrInfo = {'FREQ (MHZ)':[], 'BW (MHZ)':[], 'SIGN':[], 'NUM CHANNELS':[]}
    if sgn:
      FrInfo['SIDEBAND'] = ['U' for i in range(Nr)]
    else:
      FrInfo['SIDEBAND'] = ['L' for i in range(Nr)]

    metadata = []
    for i in range(Nr):
      FrInfo['FREQ (MHZ)'] += [(nu0 + i*bw*nch)/1.e6]
      FrInfo['BW (MHZ)'] += [bw*nch/1.e6]
      FrInfo['SIGN'] += [sgn]
      FrInfo['NUM CHANNELS'] += [int(nch)]
      freqs = nu0 + np.linspace((sgn-1.)/2.,(sgn+1.)/2.,nch,endpoint=False)*bw
      metadata.append(freqs)

    FrInfo['CHANS TO AVG'] = [1 for i in range(Nr)]
    FrInfo['OVERSAMPLE FAC.'] = [1 for i in range(Nr)] 
    FrInfo['DECIMATION FAC.']=[1 for i in range(Nr)]


    if len(doIF)==0:
     doIF = list(range(1,1+fitsf['FREQUENCY'].header['NO_BAND']))

    fitsf.close()




# ANTENNAS TO PARTICIPATE IN THE GAIN ESTIMATES:
#  nTotAnt = len(antcoords)
#
#  calAnts = []
#  for exA in antcodes:
#    if exA not in excludeAnts:
#      calAnts.append(antcodes.index(exA)+1)
#    else:
#      printMsg("Excluding antenna %s from the solution."%str(exA))

##z following section
#  try:
#    plotAnt = int(plotAnt)
#  except:
#    if plotAnt not in antcodes:
#      printError("Reference antenna %s is not found in metadata!"%str(plotAnt))
#    else:
#      plotAnt = antcodes.index(plotAnt)+1 
#
#  for i in range(len(linAntIdx)):
#    try:  
#      linAntIdx[i] = int(linAntIdx[i])
#    except:
#      if linAntIdx[i] not in antcodes:
#        linAntIdx[i] = antcodes.index(linAntIdx[i])+1
#        
#
#  if plotAnt in linAntIdx:
#    printMsg(
#      "WARNING: Plotting will involve autocorrelations. \nThis has not been fully tested!") 
#

## In some cases, these arrays are REALLY long (i.e., one value per channel, per IF).
  #printMsg("XYadd is %s"%str(XYadd))
  #printMsg("XYdel is %s"%str(XYdel))
  #printMsg("XYratio is %s"%str(XYratio))





#  FlagBas1 = []
#  FlagBas2 = []
#  for fbi in excludeBaselines:
#    printMsg("Excluding baseline %s from solution."%str(fbi))
#    if fbi[0] in antcodes and fbi[1] in antcodes:
#      FlagBas1.append(antcodes.index(fbi[0])+1) ### = np.array([int(i[0]+1) for i in excludeBaselines])
#      FlagBas2.append(antcodes.index(fbi[1])+1) ### = np.array([int(i[1]+1) for i in excludeBaselines])
#    else:
#      printError('Unknown antenna(s) %s and/or %s in excludeBaselines!\n'%(fbi[0],fbi[1]))
#
#  FlagBas1 = np.array(FlagBas1); FlagBas2 = np.array(FlagBas2)
# 
#
#  if plotAnt not in calAnts:
#    if (doSolve>=0):
#      printError("ERROR! plotAnt/Reference antenna is NOT in list of calibratable antennas!")
#    else:
#      printMsg("plotAnt (%d) is not in antenna list, so plots will be missing" % plotAnt)
#
#
#  if type(feedRotation) is not list:
#    printError("feedRotation must be a list of numbers")
#  elif len(feedRotation)==0:
#    feedRot = np.zeros(nTotAnt) 
#  elif len(feedRotation)!= nTotAnt:
#    printError("feedRotation must have %i entries!"%nTotAnt)
#  else:
#    feedRot = np.pi/180.*np.array(feedRotation, dtype=np.float)  







#######################
##### GET SPECTRAL WINDOW AUTOMATICALLY:
  if isPhased and spw < 0:
    printMsg("Deducing SPW...\n")
 
    from collections import Counter

    BPidx = -1
    for i in range(len(gains[0])):
      tb.open(os.path.join(gains[0][i],'SPECTRAL_WINDOW'))
      if tb.getcol('NUM_CHAN')[0]>1:
        BPidx = i
        tb.close()
        break
      tb.close()

    if BPidx >= 0:
      printMsg("Using BP table %d\n" % BPidx)

    if BPidx == -1:
      printError("I cannot guess the right spw if there are no BP-like tables!") 
      
    tb.open(os.path.join(gains[0][BPidx],'SPECTRAL_WINDOW'))

    calfreqs = []
    calfreqs2 = []
    for nis in range(len(tb.getcol('NAME'))):
      calfreqs.append(tb.getcell('CHAN_FREQ',nis)[0]/1.e6)
      calfreqs2.append(tb.getcell('CHAN_FREQ',nis)[-1]/1.e6)

  #  calfreqs = tb.getcol('CHAN_FREQ')[0,:]/1.e6
  #  nchansp = np.shape(tb.getcol('CHAN_FREQ'))[0]
  #  calfreqs2 = calfreqs + tb.getcol('CHAN_WIDTH')[0,:]*nchansp/1.e6

    tb.close()
    nurange = [[np.min([calfreqs[i],calfreqs2[i]]),np.max([calfreqs[i],calfreqs2[i]])] for i in range(len(calfreqs))]
    spwsel = -np.ones(len(doIF),dtype=np.int)   #[-1 for nu in doIF]
    slop = 5.0 # MHz
    for nui,nu in enumerate(doIF):
      for spwi in range(len(calfreqs)):
       try:
        nu0 = FrInfo['FREQ (MHZ)'][nu-1]
        nu1 = FrInfo['FREQ (MHZ)'][nu-1] + FrInfo['BW (MHZ)'][nu-1]*FrInfo['SIGN'][nu-1]
        nus = [np.min([nu0,nu1]),np.max([nu0,nu1])]
        print(nu, ':', nurange[spwi][0], '<', nus[0], nus[1], '<', nurange[spwi][1], end=' ')
        if (nurange[spwi][0] - slop) < nus[0] and nus[1] < (nurange[spwi][1] + slop):
          spwsel[nui] = spwi
          print(' pass')
        else:
          print(' fail')
       except:
        printMsg("WARNING! spw %i is NOT in SWIN file! Will skip it!"%nu)
        spwsel[nui] = -2
    errmsg = []
    isErr = False
    for i,spws in enumerate(spwsel):
       if spws == -1:
         isErr = True
         errmsg += [str(doIF[i])]

    if isErr:
         printMsg("WARNING! There is no spw that covers all the IF frequencies!\n" +
            "Problematic IFs are:  %s"%(','.join(errmsg)))

         doIF = [doIF[i] for i in range(len(doIF)) if i in list(np.where(spwsel>=0)[0])]
         printMsg('\n\n  NEW LIST OF IFs: '+','.join(list(map(str,doIF))))

    spwsel = list(set(spwsel[spwsel>=0]))
    if len(spwsel)>1:
       printError("There is more than one possible spw for some IFs!")

    spw = spwsel[0]
    printMsg('Selected spw: %d\n' % spw)
########################




# Get the REAL number (and names) of linear-pol antennas in this dataset:
#  nALMATrue = 0
#  linAntIdxTrue = []
#  linAntNamTrue = []
#  OrigLinIdx = []
#  for i,idd in enumerate(linAntIdx):
#    if type(idd) is int and idd<=len(antcodes):
#      nALMATrue += 1
#      linAntIdxTrue.append(idd)
#      linAntNamTrue.append(antcodes[idd-1])
#      OrigLinIdx.append(i)
#    elif idd in antcodes:
#      nALMATrue += 1
#      linAntNamTrue.append(idd)
#      linAntIdxTrue.append(antcodes.index(idd)+1)
#      OrigLinIdx.append(i)
#
#      
#  printMsg("There are %i linear-polarization antennas in THIS dataset"%nALMATrue)

  



# COMPUTE XY delays:


#  XYdelF = [[0.0,0.0] for i in range(nALMATrue)]
#  if type(XYdel) is not dict: #or len(XYdel) != nALMA:
#    printError("Invalid format for XYdel!\n") # Should be a LIST of numbers, as large as the number of linear-polarization VLBI stations!")
#
#
#  for i,doant in enumerate(linAntNamTrue):
#    if doant in XYdel.keys():
#      if type(XYdel[doant]) is list:
#        try:
#          XYdelF[i] = list(map(float,XYdel[doant])) #float(XYdel[i]*np.pi/180.)
#        except:
#          printError("Invalid format for XYdel!\n Should be a dictionary with LISTS of numbers!")
#      else:
#        try:
#          XYdelF[i] = [float(XYdel[doant]),0.0] #float(XYdel[i]*np.pi/180.)
#        except:
#          printError("Invalid format for XYdel!\n Should be a dictionary with LISTS of numbers!")
#




#  XYaddF = [[] for i in range(nALMATrue)]
#
#  for i in range(nALMATrue):
#    for j in doIF: # range(len(FrInfo['SIGN'])):
#      sgn = FrInfo['SIGN'][j-1]
#
#      if (float(FrInfo['NUM CHANNELS'][j-1]//FrInfo['CHANS TO AVG'][j-1]) !=
#          FrInfo['NUM CHANNELS'][j-1]/FrInfo['CHANS TO AVG'][j-1]):
#            printMsg("linspace check freq: %d / %d = %f" % (
#              FrInfo['NUM CHANNELS'][j-1],FrInfo['CHANS TO AVG'][j-1],
#              FrInfo['NUM CHANNELS'][j-1]/FrInfo['CHANS TO AVG'][j-1]))
#      if isSWIN:
#        NuChan = np.linspace((sgn-1.)/2.,(sgn+1.)/2.,
#            FrInfo['NUM CHANNELS'][j-1]//FrInfo['CHANS TO AVG'][j-1],
#            endpoint=False)
#      else:
#        NuChan = np.linspace(0.,sgn,
#            FrInfo['NUM CHANNELS'][j-1]//FrInfo['CHANS TO AVG'][j-1],
#            endpoint=False)
#      Nus = 1.e6*np.array(
#        FrInfo['FREQ (MHZ)'][j-1] + FrInfo['BW (MHZ)'][j-1]*NuChan,
#            dtype=np.float)
#      XYaddF[i].append(2.*np.pi*(Nus-XYdelF[i][1])*XYdelF[i][0])






# Prepare memory of XY amplitude ratios:

#  if type(XYratio) is not dict: # or len(XYratio) != nALMA:
#    printError("Invalid format for XYratio!") #\n Should be a list as large as the number of linear-polarization VLBI stations!\nThe elements of that list shall be either numbers or lists as large as the number of IFs")

#  XYratioF = [[] for i in range(nALMATrue)]
#  for i in range(nALMATrue):
#    for j in doIF:
#   #   XYratioF[i].append(np.ones(FrInfo['NUM CHANNELS'][j-1]/FrInfo['CHANS TO AVG'][j-1],dtype=np.float))
#      if (float(FrInfo['NUM CHANNELS'][j-1]//FrInfo['CHANS TO AVG'][j-1]) !=
#          FrInfo['NUM CHANNELS'][j-1]/FrInfo['CHANS TO AVG'][j-1]):
#            printMsg("linspace check freq: %d / %d = %f" % (
#              FrInfo['NUM CHANNELS'][j-1],FrInfo['CHANS TO AVG'][j-1],
#              FrInfo['NUM CHANNELS'][j-1]/FrInfo['CHANS TO AVG'][j-1]))
#      XYratioF[i].append(np.ones(
#        FrInfo['NUM CHANNELS'][j-1]//FrInfo['CHANS TO AVG'][j-1],
#        dtype=np.float))
#




# COMPUTE TIME RANGES:


#  if len(plotRange)==0:
#    plRan = np.array([0.,0.])
#    plotFringe = False
#  else:
#   try:
#     plRan = np.array([plotRange[0]+plotRange[1]/24.+plotRange[2]/1440.+plotRange[3]/86400.,plotRange[4]+plotRange[5]/24.+plotRange[6]/1440.+plotRange[7]/86400.])
#     plotFringe = True
#     if len(plotIF)==0: plotIF = list(doIF)
#   except:
#     printError("Bad time range format for plotRange!")
#
#
#  if len(Range) == 0:
#    Ran = np.array([0.,1.e20])
#  else:
#   try:
#     Ran = np.array([Range[0]+Range[1]/24.+Range[2]/1440.+Range[3]/86400.,Range[4]+Range[5]/24.+Range[6]/1440.+Range[7]/86400.])
#   except:
#     printError("Bad time range format for Range!")




#######
# WARNING! UNCOMMENT THIS IF NOT DEBUGGING!
#  if os.path.exists(OUTPUTIDI) and IDI != OUTPUTIDI:
#    printMsg('Will REMOVE the existing OUTPUT file (or directory)!\n')
#    printMsg('Copying IDI to OUTPUTIDI!\n')
#    os.system('rm -rf %s'%OUTPUTIDI)
#    os.system('cp -rf %s %s'%(IDI,OUTPUTIDI))
#  elif not os.path.exists(OUTPUTIDI):
#    printMsg('Copying IDI to OUTPUTIDI!\n')
#    os.system('cp -rf %s %s'%(IDI,OUTPUTIDI))
#     
#######


#  if isSWIN:
#    OUTPUT = []
#    PHASECALS = []
#    walk = [f for f in os.walk(OUTPUTIDI)]
#    for subd in walk:
#      OUTPUT += [os.path.join(subd[0],fi) for fi in filter(lambda x: x.startswith('DIFX_'),subd[2])]
#      PHASECALS += [os.path.join(subd[0],fi) for fi in filter(lambda x: x.startswith('PCAL_'),subd[2])]
#
#    if len(OUTPUT) == 0:
#      printError("No *.difx files found in directory!")
#
## Derive the days of observation of each difx file:
#    mjd = np.zeros(len(OUTPUT))
#    mjs = np.zeros(len(OUTPUT))
#    for i,fi in enumerate(OUTPUT):
#      mjd[i],mjs[i] = list(map(float,((os.path.basename(fi)).split('.')[0]).split('_')[1:3]))
#    mjd0 = np.min(mjd)
#    mjp = Ran + mjd0
#    metadata.append(mjd0)
#
## Filter out files outside the computing time window:
#    t0 = mjd + mjs/86400.
#    i0 = np.logical_and(t0<=mjp[1],t0>=mjp[0])
#    OUTPUT = [OUTPUT[i] for i in range(len(i0)) if i0[i]]
#
#  else:
#    metadata = []
#    OUTPUT = OUTPUTIDI




# Set XYadd and XYratio:
#
#  if type(XYadd) is not dict: # or len(XYadd.keys) != nALMA:
#    printError("Invalid format for XYadd!\n") # Should be a list as large as the number of 
#                                              #linear-polarization VLBI stations!
#                                              #The elements of that list shall be either 
#                                              #numbers or lists as large as the number of IFs
#
#
#
#  if isPcalUsed:
#    import _XPCal as XP
#    import scipy.interpolate as spint
#    printMsg('keys of XYadd:' + str(XYadd.keys()))
#    printMsg('keys of XYratio:' + str(XYratio.keys()))
#
#
#  for i,doant in enumerate(linAntNamTrue):
#      
#        
##########################
##### CORRECTIONS BASED ON PHASECAL TONES:
#
#    if usePcal[i]:
#        printMsg("Using Pcal for %d" % i)
#        PCFile = filter(lambda x: x.endswith(doant),PHASECALS)
#        if len(PCFile)==0:
#          printError("\n\n SANITY-TEST FAILURE! NO PHASECAL FILE FOR %i\n"%doant)
#        tempArr = XP.XPCal(PCFile[0],0,0.0,len(Nr))
#        print('\nDONE READING PCAL %s!\n'%doant)
#        # Update pcal files (if not doing a test):
#        if not doTest:
#          ErrCode = XP.XPConvert(PCFile[0])  
#          if ErrCode != 0:
#            printError("\n\n ERROR Converting phasecal file %s\n"%os.path.basename(PCFile[0]))
#
#
#
#        if len(tempArr[0])==0:
#          printError("\n\n ERROR! No phasecal information for antenna %i\n Will NOT convert!\n"%i)
#        else:
#
#          CPhase = spint.interp1d(tempArr[0],-tempArr[1],bounds_error=False,fill_value = 'extrapolate')
#          CAmpl = spint.interp1d(tempArr[0],tempArr[4],bounds_error=False,fill_value = 1.0)
#
#          for ji,j in enumerate(doIF):
#            sgn = FrInfo['SIGN'][j-1]  
#            if isSWIN:
#              NuChan = np.linspace((sgn-1.)/2.,(sgn+1)/2.,FrInfo['NUM CHANNELS'][j-1]/FrInfo['CHANS TO AVG'][j-1],endpoint=False)
#            else:
#              NuChan = np.linspace(0.,sgn,FrInfo['NUM CHANNELS'][j-1]/FrInfo['CHANS TO AVG'][j-1],endpoint=False)
#            Nus = np.array(FrInfo['FREQ (MHZ)'][j-1] + FrInfo['BW (MHZ)'][j-1]*NuChan,dtype=np.float)
#            XYaddF[i][ji] += (CPhase(Nus))*np.pi/180.
#
#            if doant in XYratio.keys() and XYratio[doant] == 0.0:
#              XYratioF[i][ji] *= CAmpl(Nus)          
#          del CPhase
#          Nelem = len(tempArr)
#          for j in range(Nelem-1,-1,-1):
#              del tempArr[j]
#          del tempArr
###################################
#    if doant not in XYadd.keys():
#       printMsg('WARNING IN APPLYING XYadd! ANTENNA %s IS NOT FOUND IN DICT!'%doant)
#
#    else:
#        
#      if type(XYadd[doant]) is list:
#        for j in range(len(XYadd[doant])):
#         if type(XYadd[doant][j]) in [list,np.ndarray]:
#           arrSize = np.shape(np.array(XYadd[doant][j]))
#           if len(arrSize)!=1 or arrSize[0] != IFchan:
#             printError("Shape of XYadd array(s) does not coincide with number of IF channels\n")
#           else:
#             XYaddF[i][j] += np.array(XYadd[doant][j])*np.pi/180.
#         else:
#           try:
#             XYaddF[i][j] += np.ones(IFchan)*float(XYadd[doant][j])*np.pi/180.
#           except Exception as ex:
#             printMsg(str(ex))
#             printError("Invalid format for XYadd!\nShould be a LIST of numbers (or a list of lists of numbers),\nor a list of lists of arrays")
#      else:
#        try:
#          for j in range(len(doIF)):           
#            XYaddF[i][j] += np.ones(len(XYaddF[i][j]))*float(XYadd[doant])*np.pi/180.
#        except Exception as ex:
#          printMsg(str(ex))
#          printError("Invalid format for XYadd!\n Should be a LIST of numbers (or a list of lists),\n as large as the number of linear-polarization VLBI stations!")
#
#
#
##  raw_input('HOLD')
#
#
#  UseAutoCorrs = np.zeros(nALMATrue,dtype=np.int32)
#  
#  
#  
#  for i,doant in enumerate(linAntNamTrue):
#
#    try:  
#  
#      if doant not in XYratio.keys():
#        printMsg('WARNING IN APPLYING XYratio! ANTENNA %s IS NOT FOUND IN DICT!'%doant)
#
#      else:
#          
#        if type(XYratio[doant]) not in [list,np.ndarray]:
#
#          if float(XYratio[doant]) < 0.0:
#            NchanAutos = max([1,int(float(IFchan)/float(-XYratio[doant]))])
#            UseAutoCorrs[i] = NchanAutos
#            printMsg("Will correct Antenna %i with auto-correlations, applying a median filter of %i channels.\n"%(linAntIdxTrue[i],NchanAutos))
#          elif float(XYratio[doant]) > 0.0:
#            for j in range(len(doIF)): 
#              XYratioF[i][j] *= float(XYratio[doant])
#
#        else:
#          
#          for j in range(len(XYratio[doant])):
#            if type(XYratio[doant][j]) in [list,np.ndarray]:
#              tempArr = np.array(XYratio[doant][j])
#              arrSize = np.shape(tempArr)
#              if len(arrSize)!=1 or arrSize[0] != IFchan:
#                printError("Shape of XYratio array(s) does not coincide with number of IF channels\n")
#              else:
#                XYratioF[i][j] *= tempArr
#            else:
#              try:
#                XYratioF[i][j] *= float(XYratio[doant][j])
#              except Exception as ex:
#                printMsg(str(ex))
#                printError("Invalid format for XYratio!\nShould be a list (or a list of lists,\nor a list of lists of arrays))")
#
#
#    except Exception as ex:
#      printMsg(str(ex))
#      printError("Invalid format for XYratio!\n Should be a LIST of numbers (or a list of lists),\n as large as the number of linear-polarization VLBI stations!")




# A-PRIORI GAINS:
#  PrioriGains = []
#  for i in range(len(linAntIdxTrue)):
#    PrioriGains.append([])
#    for j in range(len(XYaddF[i])):
#       PrioriGains[i].append(np.array(XYratioF[i][j]*np.exp(1.j*XYaddF[i][j]),dtype=np.complex64))






# TEMPORARY FILE TO STORE FRINGE FOR PLOTTING:
#  if os.path.exists('POLCONVERT.FRINGE'):
#    os.system('rm -rf POLCONVERT.FRINGE')



  ngain = [len(g) for g in gains]
  kind = []


# ONLY FIRST ANTENNA IN LIN-POL LIST IS ALLOWED TO BE A PHASED ARRAY (ALMA):
  NSUM = [1 for i in range(nALMA)]
  NSUM[0] = len(allants)



########################################
# Re-reference phases, if needed:
  if len(APPrefant)>0:
    if APPrefant not in allants:
      printError("\n\nERROR: Antenna %s NOT FOUND in the table!"%APPrefant)

# Get the first gain table in list that has the ".XY0" string:
    XY0 = [gi for gi in range(ngain[0]) if '.XY0' in gains[0][gi]][0]

# Re-reference and re-assign gain table:
    ReReference(calAPP, gains[0][XY0], int(spw), APPrefant)
    gains[0][XY0] = "%s.REFANT_%s"%(gains[0][XY0],APPrefant) 

########################################









####################################
# READ CALIBRATION TABLES:
  gaindata = []
  dtdata = []
  isLinear = []
  for i in OrigLinIdx:
   isLinear.append(np.zeros(len(gains[i]),dtype=np.bool))
   gaindata.append([])
   kind.append([])
   dtdata.append([])
   if dterms[i]=="NONE":
     nchan = 1
     ntime = 1
     dtdata[-1].append(np.ones(nchan).astype(np.float64))
     for ant in range(NSUM[i]):
       dtdata[-1].append([])
       dtdata[-1][-1].append(np.zeros((nchan,ntime)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((nchan,ntime)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((nchan,ntime)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((nchan,ntime)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((nchan,ntime)).astype(np.bool))
   else:
    success = tb.open(os.path.join(dterms[i],'SPECTRAL_WINDOW'))
    if not success:
      printError("ERROR READING TABLE %s"%dterms[i])

    
    dtfreqs = tb.getcell('CHAN_FREQ',int(spw))   #[:,int(spw)]
    nchan = len(dtfreqs)
    tb.close()
    success = tb.open(os.path.join(dterms[i],'ANTENNA'))
    tabants = tb.getcol('NAME')
    tb.close()
    print('Reading ',dterms[i])
    tb.open(dterms[i])
    spmask = tb.getcol('SPECTRAL_WINDOW_ID')==int(spw)
    data = []
    flagrow = []
    for di in np.where(spmask)[0]: 
      data.append(tb.getcell('CPARAM',di))  #[:,:,spmask]
      flagrow.append(tb.getcell('FLAG',di))  #[:,:,spmask]

    data = np.array(data).transpose(1,2,0)
    flagrow = np.array(flagrow).transpose(1,2,0)
 #   print np.shape(data)
 #   raw_input('HOLD')
    antrow = []
    for ai in tb.getcol('ANTENNA1')[spmask]:
      if tabants[ai] in allants: 
        antrow.append(allants.index(tabants[ai]))
      else:
        antrow.append(-1)
    antrow = np.array(antrow)


    trow = tb.getcol('TIME')[spmask]
#    flagrow = tb.getcol('FLAG')[:,:,spmask]
    flagsf = np.logical_or(flagrow[0,:,:],flagrow[1,:,:])
    flagsd = np.logical_or(np.abs(data[0,:,:])==0.0,np.abs(data[1,:,:])==0.0)
    flags = np.logical_or(flagsf,flagsd)
    tb.close()
    dtdata[-1].append(np.zeros(nchan).astype(np.float64))
    dtdata[-1][0][:] = dtfreqs
   
    for ant in range(NSUM[i]):
      dtdata[-1].append([])
      dd0 = data[0,:,:]
      dd1 = data[1,:,:]
      dims = np.shape(dd0[:,antrow==ant])
      if dims[1]>0:
       dtdata[-1][-1].append(np.zeros(dims).astype(np.float64))
       dtdata[-1][-1].append(np.zeros(dims).astype(np.float64))
       dtdata[-1][-1].append(np.zeros(dims).astype(np.float64))
       dtdata[-1][-1].append(np.zeros(dims).astype(np.float64))
       dtdata[-1][-1].append(np.zeros(dims).astype(np.bool))
       dtdata[-1][-1][0][:] = (dd0[:,antrow==ant]).real
       dtdata[-1][-1][1][:] = (dd0[:,antrow==ant]).imag
 #     unwrap(dtdata[i][-1][1][:])
       dtdata[-1][-1][2][:] = (dd1[:,antrow==ant]).real
       dtdata[-1][-1][3][:] = (dd1[:,antrow==ant]).imag
 #     unwrap(dtdata[i][-1][3][:])
       dtdata[-1][-1][4][:] = flags[:,antrow==ant]
      else:
       dtdata[-1][-1].append(np.zeros((dims[0],1)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((dims[0],1)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((dims[0],1)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((dims[0],1)).astype(np.float64))
       dtdata[-1][-1].append(np.zeros((dims[0],1)).astype(np.bool))
 
 #  print np.where(np.isnan((dd0[:,antrow==ant]).real)) 
 #  print np.where(np.isnan((dd1[:,antrow==ant]).real)) 
  
   for j,gain in enumerate(gains[i]):
     gaindata[-1].append([])
     print('Reading ',gain)
     isLinear[-1][j] = interpolation[i][j]=='linear'
     if gain=="NONE":
      nchan = 1
      ntime = 1
      kind[-1].append(0)
      gaindata[-1][j].append(np.ones(nchan).astype(np.float64))
      for ant in range(NSUM[i]):
       gaindata[-1][j].append([])
       gaindata[-1][j][-1].append(np.ones(ntime).astype(np.float64))
       gaindata[-1][j][-1].append(np.ones((nchan,ntime)).astype(np.float64))
       gaindata[-1][j][-1].append(np.zeros((nchan,ntime)).astype(np.float64))
       gaindata[-1][j][-1].append(np.ones((nchan,ntime)).astype(np.float64))
       gaindata[-1][j][-1].append(np.zeros((nchan,ntime)).astype(np.float64))
       gaindata[-1][j][-1].append(np.zeros((nchan,ntime)).astype(np.bool))
     else:

# Smooth X-Y differences:
      if gainmode[i][j]=='S' and XYavgTime>0.0:
        printMsg("Will average X-Y phase differences over %.1f seconds"%XYavgTime)
        XYsmooth(gain, XYavgTime, int(spw)) 
        gain = gain+'.XYsmooth.PolConvert'

# Read data and metadata:
      sucess = tb.open(os.path.join(gain,'SPECTRAL_WINDOW'))

      if not success:
        printError("ERROR READING TABLE %s"%gain)

      gfreqs = tb.getcell('CHAN_FREQ',int(spw))  #[:,int(spw)]
      nchan = len(gfreqs)
      tb.close()
      tb.open(os.path.join(gain,'ANTENNA'))
      tabants = tb.getcol('NAME')
      tb.close()
      tb.open(gain)
      spmask = tb.getcol('SPECTRAL_WINDOW_ID')==int(spw)
      trowns = tb.getcol('TIME')[spmask]
      tsort = np.argsort(trowns)
      trow = trowns[tsort]
      data = []
      flagrow = []

      if 'CPARAM' in tb.colnames():
        kind[-1].append(0)
        for di in np.where(spmask)[0]:
          data.append(tb.getcell('CPARAM',di))
          flagrow.append(tb.getcell('FLAG',di))
      else:
        if tb.info()['subType'] == 'B TSYS':
          kind[-1].append(2)
        else:
          kind[-1].append(1)

        for di in np.where(spmask)[0]:
          data.append(tb.getcell('FPARAM',di)) #[:,:,spmask])
          flagrow.append(tb.getcell('FLAG',di))

      data = (np.array(data)).transpose(1,2,0)[:,:,tsort]
      flagrow = (np.array(flagrow)).transpose(1,2,0)[:,:,tsort]


   #   antrow = np.array([allants.index(tabants[ai]) for ai in tb.getcol('ANTENNA1')[spmask]])[tsort]

      antrow = []
      for ai in tb.getcol('ANTENNA1')[spmask]:
        if tabants[ai] in allants: 
          antrow.append(allants.index(tabants[ai]))
        else:
          antrow.append(-1)
      antrow = np.array(antrow)[tsort]


 #     flagrow = (tb.getcol('FLAG')[:,:,spmask])[:,:,tsort]
      if np.shape(data)[0] == 2:  # A DUAL-POL GAIN (i.e., mode 'G')
        flagsf = np.logical_or(flagrow[0,:],flagrow[1,:])
        flagsd = np.logical_or(np.abs(data[0,:])==0.0,np.abs(data[1,:])==0.0)
      else: # A GAIN IN MODE 'T'
        flagsf = np.copy(flagrow[0,:])
        flagsd = np.abs(data[0,:])==0.0
      flags = np.logical_or(flagsf,flagsd)
      tb.close()
      gaindata[-1][j].append(np.zeros(nchan).astype(np.float64))
      gaindata[-1][j][0][:] = gfreqs
 #     print np.where(np.isnan(gfreqs))
      for ant in range(NSUM[i]):
        gaindata[-1][j].append([])
        if np.shape(data)[0] == 2:  # A DUAL-POL GAIN (i.e., mode 'G')
          if gainmode[i][j] in ['G','S']:
            dd0 = data[0,:,:]
            dd1 = data[1,:,:]
          else:  # DUAL-POL GAIN FORCED TO 'T' MODE:
            Aux = np.sqrt(data[0,:,:]*data[1,:,:])
            dd0 = Aux
            dd1 = Aux
        else:  # A GAIN ALREADY IN MODE 'T'
          dd0 = data[0,:,:]
          dd1 = data[0,:,:]
        antrowant = antrow==ant
        dims = np.shape(dd0[:,antrowant])
        isFlagged=False
   # All antennas MUST have the re-ref XY0 phase, even if not used 
   # in the pol. calibration!
        if ".XY0" in gain:
          if dims[1]==0:
            antrowant = antrow==refants[0]
            dims = np.shape(dd0[:,antrowant])
        else:
          if dims[1]==0:
            dims = (dims[0],1)
            isFlagged=True
            antrowant = antrow==refants[0]
       #     ant=refants[0]
        gaindata[-1][j][-1].append(np.zeros(np.shape(trow[antrowant])).astype(np.float64))
        gaindata[-1][j][-1].append(np.ones(dims).astype(np.float64))
        gaindata[-1][j][-1].append(np.zeros(dims).astype(np.float64))
        gaindata[-1][j][-1].append(np.ones(dims).astype(np.float64))
        gaindata[-1][j][-1].append(np.zeros(dims).astype(np.float64))
        gaindata[-1][j][-1].append(np.zeros(dims).astype(np.bool))
        if not isFlagged:
         gaindata[-1][j][-1][0][:] = trow[antrowant]
     #    print np.where(np.isnan(trow[antrowant]))

         if j==0:
          gaindata[-1][j][-1][1][:] = np.abs(dd0[:,antrowant])
         else: # CHANGE TO = 1.0 FOR TESTING:
          gaindata[-1][j][-1][1][:] = np.abs(dd0[:,antrowant])
         gaindata[-1][j][-1][2][:] = np.angle(dd0[:,antrowant])
         unwrap(gaindata[-1][j][-1][2]) #, check=ant<3)
         gaindata[-1][j][-1][3][:] = np.abs(dd1[:,antrowant])
         gaindata[-1][j][-1][4][:] = np.angle(dd1[:,antrowant])
         unwrap(gaindata[-1][j][-1][4]) #, check=ant<3)
         gaindata[-1][j][-1][5][:] = flags[:,antrowant]

    #     print np.where(np.isnan(np.abs(dd0[:,antrowant])))
    #     print np.where(np.isnan(np.abs(dd1[:,antrowant])))
    #     print np.where(np.isnan(flags[:,antrowant]))

#  if plotFringe and len(plotIF)==0:
#    plotIF = list(map(int,doIF))

# CALL POLCONVERT. THE LENGTH OF THE "metadata" LIST WILL TELL
# POLCONVERT WHETHER THIS IS A FITS-IDI OR A SWIN DATASET:



#  if amp_norm>0.0:
#    os.system('rm -rf POLCONVERT.GAINS')
#
#  if len(plotIF)>0:
#    os.system('rm -rf CONVERSION.MATRIX; mkdir CONVERSION.MATRIX')
#    os.system('rm -rf FRINGE.PEAKS; mkdir FRINGE.PEAKS')
#    os.system('rm -rf FRINGE.PLOTS; mkdir FRINGE.PLOTS')
#  os.system('rm -rf POLCONVERT.FRINGE; mkdir POLCONVERT.FRINGE')
#
#
#  printMsg("\n###\n### Going to PolConvert\n###")

 # print PrioriGains
#  raw_input('HOLD!')

  doAmpNorm = amp_norm>0.0


#  ALMAstuff = {'allants':allants, 'ngain':ngain, 'NSUM':NSUM, 'kind':kind,
#               'gaindata':gaindata, 'dtdata':dtdata, 'allantidx':allantidx,
#               'nphtimes':nphtimes, 'antimes':antimes, 'refants':refants,
#               'asdmtimes':asdmtimes, 'timeranges':timerangesArr[int(spw)]}

  ALMAstuff = [ngain, NSUM, kind,
               gaindata, dtdata, allantidx,
               nphtimes, antimes, refants,
               asdmtimes, timerangesArr[int(spw)]]

#  if DEBUG:
#    PC_Params = [nALMATrue, plotIF, plotAnt, len(allants), doIF, 
#        swapXY, ngain, NSUM, kind, len(gaindata), len(dtdata), OUTPUT, 
#        linAntIdxTrue, plRan, Ran, allantidx, len(nphtimes), len(antimes), 
#        len(refants), len(asdmtimes),  doTest, doSolve, doConj, doAmpNorm, 
#        np.shape(PrioriGains), len(metadata), soucoords, antcoords, antmounts, 
#        isLinear,calfield,timerangesArr[int(spw)],UseAutoCorrs,correctParangle,DEBUG]
#    printMsg("POLCONVERT CALLED WITH: %s"%str(PC_Params))

  # plotAnt is no longer used by PC.PolConvert(), but is required by doSolve
  # the second argument is "PC:PolConvert::plIF" and controls whether the huge binary fringe files are written.
#  try:
  PCONV.PolConvert(IDI=IDI, OUTPUTIDI=OUTPUTIDI, DiFXinput = DiFXinput, DiFXcalc, doIF = doIF,
                   linAntIdx=linAntIdx, Range=Range, XYadd=XYadd, XYdel=XYdel, XYratio=XYratio,
                   usePcal=usePcal, swapXY=swpaXY, swapRL=swapRL, feedRotation=feedRotation,
                   correctParangle=correctParangle, IDI_conjugated=IDI_conjugated, plotIF=plotIF,
                   plotRange=plotRange, plotAnt=plotAnt, excludeAnts=excludeAnts, 
                   excludeBaselines=excludeBaselines, doSolve=doSolve, solint=solint, 
                   doTest=doTest, npix=npix, solveAmp=solveAmp, solveMethod=solveMethod,
                   calstokes=calstokes, calfield=calfield, ALMAstuff=ALMAstuff)




#        nALMATrue, plotIF, plotAnt, len(allants), doIF, 
#        swapXY, ngain, NSUM, kind, gaindata, dtdata, OUTPUT, linAntIdxTrue, 
#        plRan, Ran, allantidx, nphtimes, antimes, refants, asdmtimes, 
#        doTest, doSolve, doConj, doAmpNorm, PrioriGains, metadata, 
#        soucoords, antcoords, antmounts, isLinear,calfield,
#        timerangesArr[int(spw)],UseAutoCorrs,bool(correctParangle),DEBUG)


#  except Exception as ex:
#    printMsg(str(ex))
#    printMsg("Continuing despite the exception, just for the fun of it")
   # didit = 0

#  printMsg("\n###\n### Done with PolConvert (status %d).\n###" % (didit))






