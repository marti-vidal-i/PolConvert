# PY_PHASES: PYTHON-BASED GLOBAL FRINGE FITTER OF SWIN FILES FOR (EU-)VGOS.
#            IT WRITES FOURFIT CONFIG FILES WITH THE PHASE SOLUTIONS.
#
#             Copyright (C) 2021  Ivan Marti-Vidal
#             University of Valencia (Spain)  
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#  
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#  
#You should have received a copy of the GNU General Public License   
#along with this program.  If not, see <http://www.gnu.org/licenses/>
#  
#



from __future__ import print_function
import pylab as pl
#import matplotlib as pl
from PolConvert import _XPCalMF as XP
import numpy as np
import struct as stk
import os, sys, glob
from multiprocessing import Pool
import gc
import scipy.interpolate as spint
import scipy.optimize as spopt
import datetime as dt
from ftplib import FTP_TLS
import scipy.interpolate as spint
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


__version__ = '0.9b (18 Nov 2021)'


def getTEC(SCAN, fig, SET='jplg', LFACT=1.0):

  K = list(filter(lambda x: "EU-VGOS" in x, sys.path))
  for ki in K:
    fname = os.path.join(ki,'EUVGOS_PY3','World_Map.png')
    if os.path.exists(fname):
       break

  if os.path.exists(fname):
    World = mpimg.imread(fname)
  else:
    raise Exception("EU-VGOS library path not set correctly")

#  SET = 'jplg'
  REARTH = 6371.e3
  SPEED_OF_LIGHT = 2.99458792e8
#  LFACT = 1.0

## Parse CALC file:
  TELS = {}
  CALC = SCAN[:-4]+'calc'
  IFF = open(CALC)
  for line in IFF.readlines():
    if line.startswith('START MJD'):
      MJD = float(line.split(':')[1])
    if line.startswith('START YEAR'):
      YY = int(line.split(':')[1])
    if line.startswith('START MONTH'):
      MM = int(line.split(':')[1])
    if line.startswith('START DAY'):
      DD = int(line.split(':')[1])
    if line.startswith('START HOUR'):
      hh = int(line.split(':')[1])
    if line.startswith('START MINUTE'):
      mm = int(line.split(':')[1])
    if line.startswith('START SECOND'):
      ss = int(line.split(':')[1])
    if line.startswith('TELESCOPE'):
      temp = line.split()
      if temp[2]=='NAME:':
        TELS[int(temp[1])] = [temp[-1].replace(' ',''),0.,0.,0.]
      if temp[2]=='X':
        TELS[int(temp[1])][1] = float(temp[-1])
      if temp[2]=='Y':
        TELS[int(temp[1])][2] = float(temp[-1])
      if temp[2]=='Z':
        TELS[int(temp[1])][3] = float(temp[-1])
    if line.startswith('SOURCE'):
      temp = line.split()
      if temp[2]=='RA:':
        RA = float(temp[-1])
      if temp[2]=='DEC:':
        DEC = float(temp[-1])

  IFF.close()


## Get Observing date and download IONEX maps:
  d0 = dt.date(YY,1,1)
  d1 = dt.date(YY,MM,DD)
  DOY = int((d1-d0).days+1)

  if not os.path.exists('TEC_ARCHIVE_%i_%i'%(YY,DOY)):
    os.system('rm TEC_ARCHIVE_%i_%i.gz'%(YY,DOY))
    destFile = '%s%3i0.%si.Z'%(SET,DOY,str(YY)[-2:])

    ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
    ftps.login(user='anonymous', passwd='imarvi2@uv.es')
    ftps.prot_p()
    directory='gps/products/ionex/%4i/%3i'%(YY,DOY)
    ftps.cwd(directory)
    ftps.retrbinary("RETR " + destFile, open('TEC_ARCHIVE_%i_%i.gz'%(YY,DOY), 'wb').write)

    os.system('gunzip TEC_ARCHIVE_%i_%i.gz'%(YY,DOY))


# Parse IONEX file to get the maps:
  EXPO = -1.0
  SCALING = 1.0
  IFF = open('TEC_ARCHIVE_%i_%i'%(YY,DOY),'r')
  lines = IFF.readlines()
  TIMES = []
  for li,line in enumerate(lines):
    if line[60:79]==  '# OF MAPS IN FILE  ':
      NMAP = int(line.split()[0])
    elif line[60:79]=='MAP DIMENSION      ':
      NDIM = int(line.split()[0])
    elif line[60:79]=='HGT1 / HGT2 / DHGT ':
      HGT1,HGT2,DHGT = map(float,line.split()[:3])
      HGT1 *= 1.e3 ; HGT2 *= 1.e3 ; DHGT *= 1.e3 ;
    elif line[60:79]=='LAT1 / LAT2 / DLAT ':
      LAT1,LAT2,DLAT = map(float,line.split()[:3])
    elif line[60:79]=='LON1 / LON2 / DLON ':
      LON1,LON2,DLON = map(float,line.split()[:3])
    elif line[60:79]=='EXPONENT           ':
      EXPO = float(line.split()[0])
    elif line[60:79]=='START OF TEC MAP   ':
      hour = list(map(int,lines[li+1].split()[:6]))
      d2 = dt.date(hour[0],hour[1],hour[2])
      dday = (d2-d1).days
      mhour = hour[3]+hour[4]/60.+hour[5]/3600.+24.*dday
      TIMES.append([int(line.split()[0]),mhour,li+2])
    elif 'COMMENT' in line[60:79] and 'TECU;' in line.split():
      temp = line.split()
      SCALING = float(temp[temp.index('TECU;')-1])
  IFF.close()

  TS = np.argsort([ti[1] for ti in TIMES])

# Get the map times bracketing the scan:
  found = False
  for ti in range(len(TIMES)-1):
    if TIMES[TS[ti]][1]<hh+mm/60. and TIMES[TS[ti+1]][1]>hh+mm/60.:
      found = True
      break

  if not found: 
    raise Exception("ERROR! IONEX MAP DOES NOT CONTAIN OBSERVING TIME!")

# Interpolation times:
  #DT = (TIMES[TS[ti+1]][1]-(hh+mm/60.))/(TIMES[TS[ti+1]][1]-TIMES[TS[ti]][1])
  DT1 = (hh+mm/60.) - TIMES[TS[ti]][1]
  DT2 = TIMES[TS[ti+1]][1]-(hh+mm/60.)

#  DT2 = 0.
#  DT1 = 1.0

## Prepare memory for maps:
  NLAT = int((LAT2-LAT1)/DLAT) ; NLON = int((LON2-LON1)/DLON)
  MAP1 = np.zeros((NLAT,NLON+1),dtype=np.float32)
  MAP2 = np.zeros((NLAT,NLON+1),dtype=np.float32)

  LATGRID = np.linspace(LAT2,LAT1,NLAT)
  LONGRID = np.linspace(LON1,LON2,NLON+1)

# Read maps:
  rLat = 0
  lread = TIMES[TS[ti]][2]
  for i in range(NLAT):
    nlonRead = 0
    lread += 1
    while nlonRead<NLON:
      line = list(map(float,lines[lread][:-1].split()))
      nlon = len(line)
      MAP1[i,nlonRead:nlonRead+nlon]=line
      nlonRead += nlon
      lread += 1

  rLat = 0
  lread = TIMES[TS[ti+1]][2]
  for i in range(NLAT):
    nlonRead = 0
    lread += 1
    while nlonRead<NLON:
      line = list(map(float,lines[lread][:-1].split()))
      nlon = len(line)
      MAP2[i,nlonRead:nlonRead+nlon]=line
      nlonRead += nlon
      lread += 1

# Build map Interpolations:
  MAP1 *= SCALING
  MAP2 *= SCALING
  MapInterp1 = spint.RectBivariateSpline(LATGRID,LONGRID,MAP1[::-1,:],kx=1,ky=1)
  MapInterp2 = spint.RectBivariateSpline(LATGRID,LONGRID,MAP2[::-1,:],kx=1,ky=1)


  
## Get GMST:
  t = (MJD-51544.0)/36525.
  Hh = (MJD - np.floor(MJD))
  GMsec = 24110.54841 + 8640184.812866*t + 0.093104*t*t - 0.0000062*t*t*t
  GMST = (GMsec/86400. + Hh)*2.*np.pi #- 0.047*np.pi/12.

  CosDec = np.cos(DEC)
  SinDec = np.sin(DEC)

  TECORR = {}

  TELCOORDS = {}

  for ant in TELS.keys():

## Get Antenna pointing direction and intersection with Ionosphere:

    TNAM = TELS[ant][0]
    LAT = np.arctan2(TELS[ant][3],np.sqrt(TELS[ant][2]**2.+TELS[ant][1]**2.))
    LON = np.arctan2(TELS[ant][2],TELS[ant][1])

    TELCOORDS[TNAM] = [LAT*180./np.pi,LON*180./np.pi]

    HANG = (GMST - RA)%(2.*np.pi) + LON

    ELEV = np.arcsin(SinDec*np.sin(LAT)+np.cos(LAT)*CosDec*np.cos(HANG))
    ZANG = np.pi/2. - ELEV

    if np.cos(ELEV)!=0.0:
      AZIM = np.arctan2(-CosDec*np.sin(HANG),np.cos(LAT)*SinDec - np.sin(LAT)*CosDec*np.cos(HANG))
    else:
      AZIM = 0.0

    if AZIM<0.0:
      AZIM += 2.*np.pi



    ZAION = np.arcsin(REARTH*np.sin(ZANG)/(REARTH+HGT1))
    THETA = ZANG - ZAION
    LATION = np.arcsin( np.sin(LAT)*np.cos(THETA)+np.cos(LAT)*np.sin(THETA)*np.cos(AZIM))
    DLATI = LATION - LAT

    SAZION = np.sin(AZIM)*np.cos(LAT)/np.cos(LATION)
    if SAZION >= 1.0: AZION = np.pi/2.
    elif SAZION <= -1.0: AZION = -np.pi/2.
    else: AZION = np.arcsin(SAZION)

    DLONG = np.arcsin(np.sin(AZIM)*np.sin(THETA)/np.cos(LATION))

    if np.abs(AZIM) > np.pi/2.:
      if AZION > 0.0:
         AZION = np.pi - AZION
      else:
         AZION = -np.pi - AZION

    IONLON = LON + DLONG
    IONLAT = LAT + DLATI


## Apply Ionosphere rotation:
    TLONG1 = IONLON*180./np.pi + 360.0/24.*DT1*LFACT
    TLONG2 = IONLON*180./np.pi - 360.0/24.*DT2*LFACT

    TLAT = IONLAT*180./np.pi

    if TLONG1 < -180.:  TLONG1 += 360.
    elif TLONG1 > 180.: TLONG1 -= 360.

    if TLONG2 < -180.:  TLONG2 += 360.
    elif TLONG2 > 180.: TLONG2 -= 360.


## Estimate TEC:
    TEC1 = MapInterp1(TLAT,TLONG1)[0][0]
    TEC2 = MapInterp2(TLAT,TLONG2)[0][0]

    TEC = (DT2*TEC1 + DT1*TEC2)/(DT1+DT2)

    TEPATH = TEC/np.cos(ZAION)

    TECORR[TNAM] = [TEC,-40.28*TEPATH*1.e16/2.99458792e8,  DLONG*180./np.pi,  DLATI*180./np.pi]

    print('%s: LAT: %.3f (%.3f)  LON: %.3f (%.3f) | EL: %.3f  AZ: %.3f | TEC: %.3f '%(TNAM,LAT*180./np.pi,DLATI*180./np.pi,LON*180./np.pi,DLONG*180./np.pi,ELEV*180./np.pi, AZIM*180/np.pi,TEPATH))


  PLOTLON1 = LONGRID + 360.0/24.*DT1*LFACT
  PLOTLON2 = LONGRID - 360.0/24.*DT2*LFACT
  PLOTLON1[PLOTLON1>180.] -= 360.
  PLOTLON1[PLOTLON1<-180.] += 360.
  PLOTLON2[PLOTLON2>180.] -= 360.
  PLOTLON2[PLOTLON2<-180.] += 360.

  MAP2PLOT = np.zeros(np.shape(MAP1))
  for i,li in enumerate(PLOTLON1):
    MAP2PLOT[:,i] = ((MapInterp1(LATGRID,li)*DT2 + MapInterp2(LATGRID,PLOTLON2[i])*DT1)/(DT1+DT2))[:,0]


 # fig = pl.figure(figsize=(20,5))
  sub= fig.add_subplot(221)  
  sub.imshow(World,extent=[-180,180,-90,90])
  cbp = sub.imshow(MAP2PLOT[:,:],origin='lower',extent=[-180,180,LAT2,LAT1],alpha=0.5)
  sub.set_title('IONEX MAP (INTERPOLATED)')
  cb = pl.colorbar(cbp)
  cb.set_label('TECU')

  sub2 = fig.add_subplot(223)
  NuFreq = np.linspace(2.0,12.0,512)
  telnam = sorted(TECORR.keys())

  fig.subplots_adjust(left=0.04,right=0.99,wspace=0.126)

  for t1 in range(len(telnam)):
    for t2 in range(len(telnam)):
      if t1>t2:
        DTEC = TECORR[telnam[t1]][1]-TECORR[telnam[t2]][1]
        Phases = DTEC/(NuFreq*1.e9)
     #   Phases = np.mod(Phases,360.)
        sub2.plot(NuFreq,Phases,'.',label='%s-%s'%(telnam[t1],telnam[t2]))
        DELAY = (Phases[1]-Phases[0])/((NuFreq[1]-NuFreq[0])*1.e9)/360.
  sub2.legend(numpoints=1)
  sub2.set_xlabel('Frequency (GHz)')
  sub2.set_ylabel('TEC Phase (cycles)')
  sub2.set_title('IONEX PREDICTION')
  #print(np.shape(MAP1))
 # print(np.shape(MAP2PLOT))




  for tel in TELCOORDS.keys():
    sub.plot(TELCOORDS[tel][1],TELCOORDS[tel][0],'.w')
    sub.text(TELCOORDS[tel][1]+2.,TELCOORDS[tel][0]+2.,tel,color='w')
    sub.plot(np.array([TELCOORDS[tel][1],TELCOORDS[tel][1]+TECORR[tel][2]]),np.array([TELCOORDS[tel][0],TELCOORDS[tel][0]+TECORR[tel][3]]),'-w')

  fig.suptitle(os.path.basename(SCAN),fontsize=20)

  return TECORR
















TWOPI = 2.*np.pi

def QuinnTau(FRN):
     return 0.25*np.log1p(3.*FRN*FRN+6.*FRN)-np.sqrt(6.)/24.*np.log1p(-2.*np.sqrt(2./3.)/(FRN+1.+np.sqrt(2./3.)))

def Quinn(FFT):
    Denom = FFT[1].real*FFT[1].real+FFT[1].imag*FFT[1].imag
    AP = (FFT[2].real*FFT[1].real + FFT[2].imag*FFT[1].imag)/Denom
    AM = (FFT[0].real*FFT[1].real + FFT[0].imag*FFT[1].imag)/Denom
    DP = -AP/(1.-AP)
    DM = AM/(1.-AM)
    return (DP+DM)/2. + QuinnTau(DP*DP) - QuinnTau(DM*DM)




#filename = 'DiFX/ev0287_030.difx'
#FLAGBAS = [['OE','OW']]

#HOPSNAMES = {'OE':'S','OW':'T','WS':'V','YJ':'Y'}
#IFNAMES = 'abcdefghijklmnopqrstuvwxyzABCDEF'
#PCALDELAYS = {'YJ':[1.e9,4.e9,'abcdefgh']}

#REFANT = 'WS'

#PCALPLOT = [['OW','YJ']]

def GET_FOURFIT_PHASES(SCAN='',HOPSNAMES={'OE':'S','OW':'T','WS':'V','YJ':'Y'}, 
                       IFNAMES = 'abcdefghijklmnopqrstuvwxyzABCDEF', FLAGBAS = [['OE','OW']],
                       PCALDELAYS = {'YJ':[1000.,4000.,'abcdefgh']}, REFANT='WS'):



  if os.path.isdir(SCAN):
    print('Path is a directory. Will look for SWIN files.')
    calcFile = '.'.join(SCAN.split('.')[:-1])+'.calc'
    inputFile = '.'.join(SCAN.split('.')[:-1])+'.input'
  else:
    raise Exception("Argument must be a difx directory")



## READ BANDWIDTH AND (CENTER) FREQUENCY OF EACH IF:
  INPF = open(inputFile,'r')
  BWs = {} ; NUs = {}; SBs = {}; CHANFREQ = {}; NCHAN = {}

  for line in INPF.readlines():
    if line.startswith('BW (MHZ)'):
      BWs[int(line.split(':')[0].split()[-1])] = float(line.split(':')[-1][:-1])*1.e6      
    if line.startswith('FREQ (MHZ)'):
      isIn = False
      for key in NUs.keys():
        if NUs[key] == float(line.split(':')[-1][:-1])*1.e6:
          isIn = True
          break
      if not isIn:
        NUs[int(line.split(':')[0].split()[-1])] = float(line.split(':')[-1][:-1])*1.e6      
    if line.startswith('SIDEBAND'):
      side = line.split()[-1]
      SBs[int(line.split(':')[0].split()[-1])] = {True: 1.0 , False: -1.0}['U' in side] 
    if line.startswith('NUM CHANNELS'):
      NCHAN[int(line.split(':')[0].split()[-1])] = int(line.split(':')[-1][:-1])    

  INPF.close()

  for IF in NUs.keys():
    CHANFREQ[IF] = NUs[IF] + BWs[IF]*np.linspace((SBs[IF]-1.)/2.,(SBs[IF]+1.)/2.,NCHAN[IF])
    NUs[IF] += BWs[IF]/2.*SBs[IF] 


  FREQORDER = []
  for nui in sorted(NUs.keys()):
    FREQORDER.append(NUs[nui])
  FREQORDER = np.array(FREQORDER)
  SORTEDIF = np.argsort(FREQORDER)





## Get Ionospheric Electron Content:
  fig = pl.figure(figsize=(15,10))
  TECOR = getTEC(SCAN,fig)





## READ ANTENNA NAMES AND PHASECALS:    
  REFID = -1
  ANAMES = {}
  PCALS = {}
  PCALADDITIVE = {}
  if os.path.isdir(SCAN):
    if os.path.exists(calcFile):
      IFF = open(calcFile)
      for line in IFF.readlines():
        if line.startswith('TELESCOPE ') and ' NAME: ' in line:
          ID = int(line.split()[1])+1
          ANAMES[ID] = line.split()[-1]
          if ANAMES[ID]==REFANT:
            REFID = int(ID)
          pcalFile = glob.glob(os.path.join(SCAN,'PCAL_*_%s'%ANAMES[ID]))
          PCALS[ID] = [[0., 0., 0.] for fri in NUs.keys()]
          if len(pcalFile)>0:
            fName = XP.XPCalMF(pcalFile[0],[],0,1)
            IFFCP = open(fName,'r')
            tempArr = []
            for line in IFFCP.readlines():
              temp = line.split()
              tempArr.append([float(temp[0])*1.e6, float(temp[1])*np.pi/180.,float(temp[2])])
            IFFCP.close()

            tempArr = np.array(tempArr)

            if ANAMES[ID] in PCALDELAYS:
              PCALADDITIVE[ANAMES[ID]] = np.copy(tempArr)

# Fit tone delay for each IF:
            for spi in CHANFREQ.keys():
              Nu0 = np.min(CHANFREQ[spi]) ; Nu1 = np.max(CHANFREQ[spi])
              NuAv = (Nu1+Nu0)/2.
              mask = np.logical_and(tempArr[:,0]>=Nu0,tempArr[:,0]<=Nu1)
              Ntone = int(np.sum(mask))

              Phases = tempArr[mask,1] 
          #    for tii in range(Ntone-1):
          #      if Phases[tii+1]-Phases[tii]>np.pi:
          #        Phases[tii+1:] -= 2.*np.pi
          #      if Phases[tii+1]-Phases[tii]<-np.pi:
          #        Phases[tii+1:] += 2.*np.pi

              Freqs = tempArr[mask,0]
              Weights = tempArr[mask,2]
              ToneDel = np.sum((Freqs-np.average(Freqs))*(Phases-np.average(Phases)))/np.sum(np.power(Freqs-np.average(Freqs),2.))/(2.*np.pi)
              Phasors = np.exp(1.j*tempArr[mask,1])

          #    TECDel = TECOR[ANAMES[ID]][1]/NuAv**2.
          #    TECPhas = -TWOPI*TECOR[ANAMES[ID]][1]/NuAv

              ResPhase = np.angle(np.sum(Phasors*np.exp(-1.j*TWOPI*(ToneDel)*(Freqs-NuAv))))

              PCALS[ID][spi] = [ToneDel,ResPhase,NuAv] #,TECDel,TECPhas]

          PCALS[ID] = np.array(PCALS[ID])


#            Dnu = np.abs(tempArr[1][0]-tempArr[0][0])
#
#            AddPcals = []
#            AddNus = []
#       ## Fill phasecals at edge frequencies:
#            for spi in CHANFREQ.keys():
#              tini = -1; tend = -1
#              Nu0 = np.min(CHANFREQ[spi]) ; Nu1 = np.max(CHANFREQ[spi])
#              if ID==4: print('FREQS for IF %i: %.1f - %.1f'%(spi,Nu0,Nu1))
#
#              if Nu1 > tempArr[-1][0]: 
#                 if Nu1 not in AddNus:
#                   AddPcals.append([Nu1,tempArr[-1][1]+(tempArr[-1][1]-tempArr[-2][1])/Dnu*(Nu1-tempArr[-1][0])])
#                   AddNus.append(Nu1)
#                   if ID==4: print('END: ', AddPcals[-1], 'TOT')
#              if Nu0 < tempArr[0][0]: 
#                 if Nu0 not in AddNus:
#                   AddPcals.append([Nu0,tempArr[0][1]-(tempArr[1][1]-tempArr[0][1])/Dnu*(tempArr[0][0]-Nu0)])
#                   AddNus.append(Nu0)
#                   if ID==4: print('INI: ',AddPcals[-1], 'TOT')
#              
#              for ti in range(1,len(tempArr)):
#                if tempArr[ti-1][0]<Nu0 and tempArr[ti][0]>Nu0:
#                  tini = ti
#                  break
#              for ti in range(1,len(tempArr)):
#                if tempArr[ti-1][0]<Nu1 and tempArr[ti][0]>Nu1:
#                  tend = ti
#                  break
#
#              if tini>=0 and Nu0 not in AddNus:
#                AddPcals.append([Nu0,tempArr[tini][1]-(tempArr[tini+1][1]-tempArr[tini][1])/Dnu*(tempArr[tini][0]-Nu0)])
#                AddNus.append(Nu0)
#                if ID==4: print('INI: ',AddPcals[-1],tini)
#              if tend>=0 and Nu1 not in AddNus:
#                AddPcals.append([Nu1,tempArr[tend-1][1]+(tempArr[tend-1][1]-tempArr[tend-2][1])/Dnu*(Nu1 - tempArr[tend-1][0])])
#                AddNus.append(Nu1)
#                if ID==4: print('END: ',AddPcals[-1],tend)
#
#
#            tempArr2 = np.array(tempArr+AddPcals)
#            SORT = np.argsort(tempArr2[:,0])
#            tempArr3 = np.copy(tempArr2[SORT,:])
#
#          #  for ti in range(len(tempArr3)-1):
#          #    if tempArr3[ti+1,1]-tempArr3[ti,1]>np.pi:
#          #      tempArr3[ti+1:,1] -= 2.*np.pi
#          #    if tempArr3[ti+1,1]-tempArr3[ti,1]<-np.pi:
#          #      tempArr3[ti+1:,1] += 2.*np.pi
#
#
#            PCALS[ID] = spint.interp1d(tempArr3[:,0],tempArr3[:,1],kind='linear',bounds_error=False,fill_value = 'extrapolate')
#          else:
#            PCALS[ID] = lambda x: 0.0




    filename = glob.glob(os.path.join(SCAN,'DIFX_*'))
    if len(filename)==0:
      raise Exception('File (or dir) not found.\n')
    else:
      filename = filename[0]



  if REFID < 0:
    raise Exception("Reference antenna %s not found!\n"%REFANT)


  pcalNu = np.zeros(2*len(NUs.keys())) #len(NUs.keys()))
  pcalPlot1 = np.zeros(2*len(NUs.keys())) #len(NUs.keys()))
  pcalPlot2 = np.zeros(2*len(NUs.keys())) #len(NUs.keys()))

  for i in range(len(NUs.keys())):   
    pcalNu[2*i] = SORTEDIF[i]
    pcalNu[2*i+1] = SORTEDIF[i]+1


  if False:
   for plpair in PCALPLOT:
    pc1 = -1; pc2 = -1
    for aid in ANAMES.keys():
      if ANAMES[aid]==plpair[0]:
        pc1 = int(aid)
      if ANAMES[aid]==plpair[1]:
        pc2 = int(aid)

    fig2 = pl.figure(figsize=(15,3))
    sub = fig2.add_subplot(111)
    fig2.subplots_adjust(left=0.05,right=0.983)
    for i in range(len(NUs.keys())):   
      pcalPlot1[2*i:2*i+2] = PCALS[pc1][i][1]*180./np.pi
      pcalPlot2[2*i:2*i+2] = PCALS[pc2][i][1]*180./np.pi


    for i in range(len(NUs.keys())):
      sub.plot(pcalNu[2*i:2*i+2],-pcalPlot1[2*i:2*i+2],'-g')
      sub.plot(pcalNu[2*i:2*i+2],-pcalPlot2[2*i:2*i+2],'-m')
      sub.plot(np.array([pcalNu[2*i],pcalNu[2*i]]),np.array([-180.,180.]),':k')
      sub.plot(np.array([pcalNu[2*i+1],pcalNu[2*i+1]]),np.array([-180.,180.]),':k')
      if pc1 > 0:
        sub.text(pcalNu[2*i],-240.,'%.1f'%(PCALS[pc1][int(i)][0]*1.e9),color='g')
      if pc2 > 0:
        sub.text(pcalNu[2*i],-220.,'%.1f'%(PCALS[pc2][int(i)][0]*1.e9),color='m')
    sub.set_ylim((-250,185))
    pl.text(10,200,plpair[0],color='g')
    pl.text(12,200,plpair[1],color='m')

    pl.savefig('PCALS_%s-%s.png'%(plpair[0],plpair[1]))
 # pl.show()

 





## READ VISIBILTIES:

  frfile = open(filename,"rb")

  WORD = b'\x00\xff\x00\xff\x01\x00\x00\x00'


## Figure out number of channels:
  temp = frfile.read(8+4+4+8+4+4+4+2+4+8+8*3)
  for i in range(4096):
    test = frfile.read(8)
    if test==WORD:
      break

  NCHAN = i
  print("There seem to be %i channels.\n"%i)
  frfile.close()

## Read data:
  frfile = open(filename,"rb")
  alldats = frfile.read(8)
  i=0
  DATA = {'ANTS':[],'JDT':[],'IF':[],'POL':[],'VIS':[],'UVW':[]}
  ALLANTS = []
  while True:
    if i%1024==0:
      sys.stdout.write('\r Reading VIS %i'%i)
      sys.stdout.flush()
    alldats = frfile.read(4+4+8+4+4+4)
    if not alldats: break
    BASEL,MJD,SEC,CFI,SI,SPI = stk.unpack("iidiii",alldats)
    A1 = BASEL//256
    A2 = BASEL%256

    if i==0:
      MJD0 = MJD

    if A1 not in ALLANTS:
      ALLANTS.append(A1)
    if A2 not in ALLANTS:
      ALLANTS.append(A2)

    alldats = frfile.read(2)
    P1,P2 = stk.unpack("cc",alldats)
    alldats = frfile.read(4)
    PB = stk.unpack("i",alldats)
    alldats = frfile.read(8)
    PW = stk.unpack("d",alldats)
    alldats = frfile.read(8*3)
    U,V,W = stk.unpack("ddd",alldats)
    visib = np.zeros(NCHAN,dtype=np.complex64)
    for chi in range(NCHAN):
      alldats = frfile.read(8)
      Re,Im = stk.unpack("ff",alldats)
      visib[chi] = Re + 1.j*Im
    hola = frfile.read(8)

## Only store cross-correlations:
    if A1!=A2:
      i += 1
      DATA['ANTS'].append([A1,A2])
      DATA['JDT'].append(86400.*(MJD-MJD0)+SEC)
      DATA['IF'].append(SPI)
      PC1 = str(P1.decode('utf-8')); PC2 = str(P2.decode('utf-8'))
      if PC1 == 'X':
         PC1 = 'R'
      if PC1 == 'Y':
         PC1 = 'L'
      if PC2 == 'X':
         PC2 = 'R'
      if PC2 == 'Y':
         PC2 = 'L'
      DATA['POL'].append(PC1+PC2)
      DATA['UVW'].append([U,V,W])
      DATA['VIS'].append(visib)



  print('\n Done read')
  NVIS = i
  frfile.close()
# Convert to arrays:
  for key in DATA.keys():
    DATA[key] = np.array(DATA[key])




# DETERMINE SET OF ANTENNAS, BASELINES AND IFs:

  ALLBASNOFLAG = np.unique(DATA['ANTS'],axis=0)
  ALLBAS = []

  for bif in ALLBASNOFLAG:
    good = True
    for bi in FLAGBAS:
      if ANAMES[int(bif[0])] in bi and ANAMES[int(bif[1])] in bi:
        good = False ; break
    if good:
      ALLBAS.append(bif)




  ALLANTS = list(np.unique(ALLBAS))
  ALLSPW = np.unique(DATA['IF'])


## Apply PCALS and Ionosphere corrections!

#  for bi in ALLBAS:
#    mask = np.logical_and(DATA['ANTS'][:,0]==bi[0],DATA['ANTS'][:,1]==bi[1])
#    for spi in ALLSPW:
#      mask2 = np.logical_and(mask,DATA['IF']==int(spi))
#      DATA['VIS'][mask2,:] *= np.exp(1.j*(PCALS[bi[0]](CHANFREQ[spi])-PCALS[bi[1]](CHANFREQ[spi])))
#      del mask2
#      gc.collect()
#    del mask
#    gc.collect()

  for bi in ALLBAS:
    mask = np.logical_and(DATA['ANTS'][:,0]==bi[0],DATA['ANTS'][:,1]==bi[1])
    for spi in ALLSPW:
      mask2 = np.logical_and(mask,DATA['IF']==int(spi))
      DATA['VIS'][mask2,:] *= np.exp(1.j*( TWOPI*(PCALS[bi[0]][spi][0]-PCALS[bi[1]][spi][0])*(CHANFREQ[spi]-np.average(CHANFREQ[spi])) + PCALS[bi[0]][spi][1] - PCALS[bi[1]][spi][1] + TWOPI*(TECOR[ANAMES[bi[0]]][1]-TECOR[ANAMES[bi[1]]][1])/CHANFREQ[spi]  ))
      del mask2
      gc.collect()
    del mask
    gc.collect()




## ESTIMATE TEC OVER EACH STATION:
#  TECS = getTEC(filename)




## REMOVE DIFFERENTIAL TEC FROM ALL BASELINES:






### PERFORM A GLOBAL FRINGE FITTING (PER IF):


  maskRR = DATA['POL']=='RR'
  maskLL = DATA['POL']=='LL'



  RESIDUALS = {}

  OBSTIMES = {}

  for bi in ALLBAS:
    mask = np.logical_and(DATA['ANTS'][:,0]==bi[0],DATA['ANTS'][:,1]==bi[1])
    bistr = '%i-%i'%(bi[0],bi[1])
    RESIDUALS[bistr] = {}
    print('FFT visibilities for baseline %s-%s\n'%(ANAMES[bi[0]],ANAMES[bi[1]]))
    SCTIMES = DATA['JDT'][mask]
    SCDUR = np.max(SCTIMES) - np.min(SCTIMES)
    for si in ALLSPW:
   #   sys.stdout.write('\r IF %02i'%si); sys.stdout.flush()
      mask2 = np.logical_and(mask,DATA['IF']==int(si))
      VIS = np.copy(DATA['VIS'][mask2*maskRR,:]+DATA['VIS'][mask2*maskLL,:])
      FR = np.fft.fftshift(np.fft.fft2(VIS))
      FRA = np.abs(FR)
      Nt,Nch = np.shape(VIS)
    #  print(Nt,Nch)
      if si==ALLSPW[0]:
        OBSTIMES[bistr] = np.linspace(-SCDUR/2.,SCDUR/2.,Nt)
        
      Ntot = Nt*Nch
      PEAK = np.unravel_index(np.argmax(FRA),np.shape(FR))
      RMS = np.sqrt(((np.std(FRA)**2. + np.mean(FRA)**2.) - np.sum(np.power(FRA[PEAK[0]-1:PEAK[0]+2,PEAK[1]-1:PEAK[1]+2],2.))/Ntot))
      SNR= FRA[PEAK[0],PEAK[1]]/RMS

      Taround = [PEAK[0]-1,PEAK[0],PEAK[0]+1]
      Faround = [PEAK[1]-1,PEAK[1],PEAK[1]+1]
      if Taround[1]==0: Taround[0]=Nt-1
      if Taround[1]==Nt-1: Taround[2]=0
      if Faround[1]==0: Faround[0]=Nch-1
      if Faround[1]==Nch-1: Faround[2]=0

      BLDelay = (PEAK[1] - Nch/2. + Quinn(FR[PEAK[0],Faround]))
      if bool(Nt%2):
        BLRate = (PEAK[0] - (Nt-1)/2. + Quinn(FR[Taround,PEAK[1]]))
      else:
        BLRate = (PEAK[0] - Nt/2. + Quinn(FR[Taround,PEAK[1]]))

      if BLDelay > Nch/2.0: BLDelay -= Nch
      if BLRate > Nt/2.0: BLRate -= Nt

      BLDelay /= BWs[int(si)]*Nch/(Nch-1.)
      BLRate /= SCDUR

      RESIDUALS[bistr][si] = [BLDelay,BLRate,1., SNR]

      del VIS,FR,FRA,Taround,Faround,PEAK,mask2
      gc.collect()
    del SCTIMES,mask
    gc.collect()

# Code for globalization:
  HESSIAN = np.zeros((len(ALLANTS)-1,len(ALLANTS)-1))
  DELRES = np.zeros(len(ALLANTS)-1)
  RATRES = np.copy(DELRES)
  COVMAT = np.copy(HESSIAN)
  DELAYS = np.copy(DELRES)
  RATES = np.copy(RATRES)
  GAINS = {}

  for i in ALLANTS:
    GAINS[i] = [np.zeros(len(ALLSPW)),np.zeros(len(ALLSPW))]

  for si in ALLSPW:
    HESSIAN[:] = 0.0 ; DELRES[:] = 0.0 ; RATRES[:] = 0.0
    for bi in ALLBAS:
      bistr = '%i-%i'%(bi[0],bi[1])
      a1o = int(bi[0]) ## ALLANTS.index(bi[0])+1
      a2o = int(bi[1]) ## ALLANTS.index(bi[1])+1

      if a1o<REFID:
        a1 = a1o-1
      elif a1o>REFID:
        a1 = a1o-2

      if a2o<REFID:
        a2 = a2o-1
      elif a2o>REFID:
        a2 = a2o-2

      if a1o!=REFID:
        RATRES[a1] += RESIDUALS[bistr][si][2]*RESIDUALS[bistr][si][1]
        DELRES[a1] += RESIDUALS[bistr][si][2]*RESIDUALS[bistr][si][0]
        HESSIAN[a1,a1] += RESIDUALS[bistr][si][2]
      if a2o!= REFID:
        RATRES[a2] -= RESIDUALS[bistr][si][2]*RESIDUALS[bistr][si][1]
        DELRES[a2] -= RESIDUALS[bistr][si][2]*RESIDUALS[bistr][si][0]
        HESSIAN[a2,a2] += RESIDUALS[bistr][si][2]
      if a1o!=REFID and a2o!=REFID:
        HESSIAN[a1,a2] -= RESIDUALS[bistr][si][2]
        HESSIAN[a2,a1] -= RESIDUALS[bistr][si][2]

    COVMAT[:] = np.linalg.inv(HESSIAN)
    DELAYS[:] = COVMAT.dot(DELRES)
    RATES[:] = COVMAT.dot(RATRES)
    for ai in ALLANTS:
      if ai<REFID:
        GAINS[ai][0][si] = -DELAYS[ai-1]
        GAINS[ai][1][si] = -RATES[ai-1]
      elif ai>REFID:
        GAINS[ai][0][si] = -DELAYS[ai-2]
        GAINS[ai][1][si] = -RATES[ai-2]



  NuPlot = np.zeros(len(ALLSPW))
  for i in range(len(ALLSPW)):
    NuPlot[i] = NUs[i]



## DETERMINE PHASES OF EACH IF AND BASELINE (AFTER DELAY-RATE CORRECTIONS):
### El delay the antenna 1 se suma para corregir la franja
    




  AUXVIS = {}
  AUXVIS2 = {}
  AUXVIS3 = {}

  PHASES = {}
  WEIGHTS = {}
  colors = {'1-3':'r','1-4':'b','2-3':'g','3-4':'k','2-4':'c'}

  if False:
    fig3 = pl.figure(figsize=(15,5))
    sub = fig3.add_subplot(111)

  for bi in ALLBAS:
    bistr = '%i-%i'%(bi[0],bi[1])
    mask = np.logical_and(DATA['ANTS'][:,0]==bi[0],DATA['ANTS'][:,1]==bi[1])
    PHASES[bistr] = []
    for si in ALLSPW:
      mask2 = np.logical_and(mask,DATA['IF']==int(si))
      VIS = np.copy(DATA['VIS'][mask2*maskRR,:]+DATA['VIS'][mask2*maskLL,:])
    #  if si==20:
    #    AUXVIS[bistr] = np.copy(VIS)
      VIS *= np.exp(1.j*TWOPI*(GAINS[bi[0]][0][si]-GAINS[bi[1]][0][si])*(CHANFREQ[si]-NUs[si]))[np.newaxis,:]
    #  if si==20:
    #    AUXVIS3[bistr] = np.copy(VIS)
      VIS *= np.exp(1.j*TWOPI*(GAINS[bi[0]][1][si]-GAINS[bi[1]][1][si])*OBSTIMES[bistr])[:,np.newaxis]
      PHASES[bistr].append(np.angle(np.average(VIS)))
      CHPLOT = si*NCHAN + NCHAN*np.linspace(0,1,NCHAN)
   #   sub.plot(CHPLOT,180./np.pi*np.angle(np.average(VIS,axis=0)),'.%s'%colors[bistr])

      if si==20:
        AUXVIS2[bistr] = np.copy(VIS)

    PHASES[bistr] = np.array(PHASES[bistr])

  #pl.show()



  def globPhases(p,IF):
    ChiSq = 0.0
    resid = []
    for bistr in PHASES.keys():
      bi0,bi1=map(int,bistr.split('-'))

      if bi0<REFID:
        p1 = bi0-1
      elif bi0>REFID:
        p1 = bi0-2

      if bi1<REFID:
        p2 = bi1-1
      elif bi1>REFID:
        p2 = bi1-2


      if bi0 == REFID:
        GainPh = -p[p2]
      elif bi1 == REFID:
        GainPh = p[p1]
      else:
        GainPh = p[p1] - p[p2]
      ObsPh = PHASES[bistr][IF]
      resid += [np.cos(ObsPh)-np.cos(GainPh), np.sin(ObsPh)-np.sin(GainPh)]
      ChiSq += resid[-2]**2. + resid[-1]**2. #(np.cos(ObsPh)-np.cos(GainPh))**2. + (np.sin(ObsPh)-np.sin(GainPh))**2.
    return resid

  GlobalPhases = {}
  for ai in ALLANTS:
      GlobalPhases[ai] = np.zeros(len(ALLSPW)) #[0. for spi in range(len(ALLSPW))] 


  for spi in ALLSPW:
#    myfit = spopt.fmin(globPhases,[0. for i in range(len(ALLANTS)-1)],args=(spi,))
    pini = [0. for i in range(len(ALLANTS)-1)]
    for i in ALLANTS:
      bistr = '%i-%i'%(REFID,i)
      if bistr in PHASES.keys():

        if i<REFID:
          pi = i-1
        elif i>REFID:
          pi = i-2
        pini[pi] = PHASES[bistr][spi]

    myfit = spopt.leastsq(globPhases,pini,args=(spi,))[0]
    for ai in ALLANTS:
      if ai<REFID:
        aig = ai-1
      elif ai>REFID:
        aig = ai-2
      if ai!=REFID:
        GlobalPhases[ai][spi] = myfit[aig]*180./np.pi
    #  GlobalPhases[ai+2][spi] += 360.*np.median(GAINS[ai+2][0])*(NUs[spi]-6.e9)



  if False:
    fig4 = pl.figure(figsize=(15,5))
    sub = fig4.add_subplot(111)

  for bi in ALLBAS:
    bistr = '%i-%i'%(bi[0],bi[1])
    mask = np.logical_and(DATA['ANTS'][:,0]==bi[0],DATA['ANTS'][:,1]==bi[1])
    WEIGHTS[bistr] = []
    for si in ALLSPW:
      mask2 = np.logical_and(mask,DATA['IF']==int(si))
      VIS = np.copy(DATA['VIS'][mask2*maskRR,:]+DATA['VIS'][mask2*maskLL,:])
      VIS *= np.exp(1.j*TWOPI*(GAINS[bi[0]][0][si]-GAINS[bi[1]][0][si])*(CHANFREQ[si]-NUs[si])-1.j*(GlobalPhases[bi[0]][si]-GlobalPhases[bi[1]][si])*np.pi/180.)[np.newaxis,:]
      VIS *= np.exp(1.j*TWOPI*(GAINS[bi[0]][1][si]-GAINS[bi[1]][1][si])*OBSTIMES[bistr])[:,np.newaxis]
      WEIGHTS[bistr].append(1./np.std(np.angle(np.average(VIS,axis=0)))**2.)

    WEIGHTS[bistr] = np.array(WEIGHTS[bistr])

    #  CHPLOT = si*NCHAN + NCHAN*np.linspace(0,1,NCHAN)
    #  if si==ALLSPW[0]:
    #    sub.plot(CHPLOT,180./np.pi*np.angle(np.average(VIS,axis=0)),'.%s'%colors[bistr],label='%s-%s'%(ANAMES[bi[0]],ANAMES[bi[1]]))
    #  else:
    #    sub.plot(CHPLOT,180./np.pi*np.angle(np.average(VIS,axis=0)),'.%s'%colors[bistr])
    

  if False:
    sub.set_xlabel('Channel #')
    sub.set_ylabel('Residual Phase (deg)')
    sub.set_ylim((-181.,181.))
    pl.legend(numpoints=1)
 # pl.show()





  K = list(filter(lambda x: "EU-VGOS" in x, sys.path))

  for ki in K:
    fname = os.path.join(ki,'EUVGOS_PY3','cf_TEMPLATE')
    if os.path.exists(fname):
      break

  if os.path.exists(fname):
    os.system('cp %s cf_PyPhases'%(fname))
  else:
    raise Exception("EU-VGOS library path not properly set")

  IFF = open('cf_PyPhases','a')

  print('\n\n  *** ADDITIVE PHASES ESTIMATED WITH PyPhases. VERSION %s ***\n\n'%__version__, file=IFF)

  for ant in PCALDELAYS.keys():
    Iant = -1
    for ID in ANAMES.keys():
      if ANAMES[ID]==ant:
        Iant = int(ID)
        break
    if Iant <0:
      raise Exception("Unknown antenna %s"%ant)

    msg = '\nif station %s\n  delay_offs %s'%(HOPSNAMES[ant],PCALDELAYS[ant][2])
    mask = np.logical_and(FREQORDER>=PCALDELAYS[ant][0]*1.e6,FREQORDER<=PCALDELAYS[ant][1]*1.e6)
    nomask = np.logical_not(mask)
    AVDEL1 = np.median(GAINS[Iant][0][mask])*1.e9
    AVDEL2 = np.median(GAINS[Iant][0][nomask])*1.e9
    print('\n* Delays estimated from pcals/nopcal fits (in ns): \n*   %s delays: %.3e (ALL)  %.3e (NO PCAL)  %.3e (PCAL)\n'%(ant,np.median(GAINS[Iant][0])*1.e9,AVDEL1,AVDEL2),file=IFF)

    INST_DELAY = AVDEL1-AVDEL2
    GAINS[Iant][0][mask] -= INST_DELAY*1.e-9

    msg += (' %.2f '%(INST_DELAY))*np.sum(mask)
    print(msg,file=IFF);
    print('\n',file=IFF)
    for spi in np.where(mask)[0]:
      Freqs = PCALADDITIVE[ant][np.logical_and(PCALADDITIVE[ant]>=np.min(CHANFREQ[spi]),PCALADDITIVE[ant]<=np.max(CHANFREQ[spi]))]
      AvFreq = (np.min(CHANFREQ[spi])+np.max(CHANFREQ[spi]))/2.

      AddPhase =  np.angle(np.sum(np.exp(1.j*TWOPI*INST_DELAY*(1.e-9)*(Freqs-AvFreq))))*180/np.pi
      GlobalPhases[Iant][spi] -= AddPhase
#      print('Pcal phase for %s %s: %.2f deg.'%(ant,IFNAMES[SORTEDIF[spi]],AddPhase))


  for spi in ALLSPW:
    for ai in range(len(ALLANTS)):
      GlobalPhases[ai+1][spi] = np.mod(GlobalPhases[ai+1][spi],360.)
      if GlobalPhases[ai+1][spi] > 180:
        GlobalPhases[ai+1][spi] -= 360.
      if GlobalPhases[ai+1][spi] < -180:
        GlobalPhases[ai+1][spi] += 360.




  for ai in GlobalPhases.keys():
    msg = '\nif station %s\n   pc_phases %s '%(HOPSNAMES[ANAMES[ai]],IFNAMES)
    SortedPhases = GlobalPhases[ai][SORTEDIF]
    msg += ('%.1f '*len(ALLSPW))%tuple(SortedPhases)
    print(msg,file=IFF);
    print('\n',file=IFF)






## Plot SBD vs. frequency (i.e., try to get signal from residual TEC):
  sub = fig.add_subplot(122)
  TFAC = 2.99792458e8/(40.26e16)
  for i in ALLANTS:

    ONEOVERNUSQ = 1./(NuPlot**2.)
    AVNUSQ = np.average(ONEOVERNUSQ)
    AVGAIN = np.average(GAINS[i][0])

    SXX = np.sum(ONEOVERNUSQ**2.)
    XXN = np.sum(ONEOVERNUSQ)**2./len(NuPlot)
    SXY = np.sum(ONEOVERNUSQ*GAINS[i][0])
    XYN = np.sum(ONEOVERNUSQ)*np.sum(GAINS[i][0])/len(NuPlot)
    
    SLOPE = (SXY - XYN)/(SXX - XXN)
    SLOPE_ERR = np.sqrt(np.sum((GAINS[i][0]-AVGAIN)**2.)/(len(NuPlot)-2.))/np.sqrt(np.sum((ONEOVERNUSQ-AVNUSQ)**2.))
    sub.plot(NuPlot/1.e9,(GAINS[i][0] - np.median(GAINS[i][0]))*1.e9,'o',label='%s: %.2f +/- %.2f TECU'%(ANAMES[i],SLOPE*TFAC,SLOPE_ERR*TFAC))  
    sub.set_title('RESIDUAL SBDs vs. FREQUENCY')

  sub.set_xlabel('Frequency (GHz)')
  sub.set_ylabel('Residual Antenna SBD (ns)')
  pl.sca(sub)
  sub.set_ylim((-10,10))
  pl.legend(numpoints=1)

#  NUMI2 = 1./(NuPlot**2.)
#  for bi in ALLBAS:
#    bistr = '%i-%i'%(bi[0],bi[1])
#    TOFIT = np.zeros(len(ALLSPW))
#    for si in ALLSPW:
#      TOFIT[si] = RESIDUALS[bistr][si][0]
#    SUMWGT = np.sum(WEIGHTS[bistr])
#    SUMXY = np.sum(WEIGHTS[bistr]*TOFIT*NUMI2)
#    SUMX = np.sum(WEIGHTS[bistr]*NUMI2)
#    SUMY = np.sum(WEIGHTS[bistr]*TOFIT)
#    SUMXX = np.sum(WEIGHTS[bistr]*NUMI2*NUMI2)
#
#    DET = SUMWGT*SUMXX - SUMX**2.
#
#    SLOPE = (SUMWGT*SUMXY - SUMX*SUMY)/DET
#    SBD0 =  (SUMXX*SUMY-SUMX*SUMXY)/DET
# 
#    CHI2 = np.sum(WEIGHTS[bistr]*np.power(TOFIT - (SLOPE*NUMI2+SBD0),2.))
#    ERROR = np.sqrt(SUMWGT*CHI2/DET)
#
#    sub.errorbar(NuPlot/1.e9,(TOFIT-SBD0)*1.e9,1./np.sqrt(WEIGHTS[bistr]),fmt='k',linestyle='none')  
#
#    sub.plot(NuPlot/1.e9,(TOFIT-SBD0)*1.e9,'o',label='%s: %.2f +/- %.2f TECU'%(ANAMES[i],SLOPE*TFAC,ERROR*TFAC))  
#    sub.set_title('RESIDUAL SBDs vs. FREQUENCY')
#
#  sub.set_xlabel('Frequency (GHz)')
#  sub.set_ylabel('Residual Antenna SBD (ns)')
#  pl.sca(sub)
#  sub.set_ylim((-10,10))
#  pl.legend(numpoints=1)
#




## DETERMINE ADDITIVE PHASES FOR IF ALIGNMENT:
#TOPLOT = '2-4'

#pl.plot(180./np.pi*np.angle(np.average(AUXVIS[TOPLOT],axis=0)),'ok'); pl.plot(180./np.pi*np.angle(np.average(AUXVIS2[TOPLOT],axis=0)),'or'); pl.show()

#pl.plot(180./np.pi*np.angle(np.average(AUXVIS3[TOPLOT],axis=1)),'ok'); pl.plot(180./np.pi*np.angle(np.average(AUXVIS2[TOPLOT],axis=1)),'or'); pl.show()

#pl.plot(NuPlot,PHASES[TOPLOT],'ok'); 
  pl.savefig(os.path.basename(SCAN)[:-5]+'_TEC.png')


  IFF.close()

