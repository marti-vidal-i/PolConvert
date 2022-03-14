# SOURCE_SCANNER: SIMPLE CHECK OF SNR FROM A SET OF SWIN FILES.
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



import numpy as np
import pylab as pl
import os, sys, glob

if __name__=='__main__':
 DIRE = 'DiFX'
 EXP  = 'vgt274'
 Nchan = 128; NIF = 32

 SNRBas = [[0,1],[0,2],[1,2]]
 SNRCut = 10.0


#######################################
# COMMENT THIS LINE OUT WHEN DEBUGGING AS execfile(...)
def SOURCE_SCANNER(EXPNAME='', DIFX_DIR='',SNRCut = 10.0):
 """ Reads all the scans in a SWIN directory and creates an ASCII file
  with information about the sources, participating antennas, and SNRs 
  of the correlation products. If the minimum SNR is equal or higher than 
  SNRCut, the scan is marked as good.""" 
#######################################



 try:

  EXP = EXPNAME; DIRE = DIFX_DIR

  print(EXP,DIFX_DIR)

# Figure out number of antennas, IFs and channels per IF:
  inps = [f[:-1] for f in os.popen('find ./%s -name \"*.input\" -print'%(DIRE))]
  Nants = 0; NIF = 0; Nchan=0
  for i,inp in enumerate(inps):
    temp = open(inp)
    lines = temp.readlines()
    for li in lines:
        if li.startswith('TELESCOPE ENTRIES:'):
          Nants = max([Nants,int(li.split()[-1])])

        if i==0 and li.startswith('FREQ ENTRIES:'):
            NIF = int(li.split()[-1])
        if i==0 and li.startswith('NUM CHANNELS 0:'):
            Nchan = int(li.split()[-1])
    temp.close()
    del lines


  print('There are %i IFs with %i channels each.'%(NIF,Nchan))



# ID list of baselines:
  SNRBas = []
  for i in range(Nants-1):
    for j in range(i+1,Nants):
      SNRBas.append([i,j])



# Read visibilities for SNR estimate:
  OUTPUT = open('SOURCES_%s.txt'%EXP,'w')
  calcs = [f[:-1] for f in os.popen('find ./%s -name \"*.calc\" -print'%(DIRE))]



  fmt = "  ANT %i: X = [%.2f  %.2f] | Y = [ %.2f  %.2f]  %s  \n" * Nants
  fmt += "  SNR PASS: %s\n\n"

  for dd in sorted(calcs):

# Read visibilities:
    DIFXFile = glob.glob('%s.difx/DIFX*'%dd[:-5])[0]
    frfile = open(DIFXFile,"rb")
    alldats = frfile.read(8)

    dtype2 = np.dtype([("BAS",np.int32),("MJD",np.int32),("SEC",np.float64),
                ("AUX1",np.int32),("AUX2",np.int32),("IF",np.int32),
                ("POL1",np.dtype('b')),("POL2",np.dtype('b')),("AUX3",np.int32),
                ("WEIGHT",np.float64), ("U",np.float64),("V",np.float64),
                ("W",np.float64),("VISIB",np.complex64,Nchan+1)])

    # Stupid numpy error:
    try:
      fringe2 = np.fromfile(frfile,dtype=dtype2)
    except:
      fringe2 = np.fromfile(frfile,dtype=dtype2)

    frfile.close()

    SNR_X = [[1.e18,0.] for i in range(Nants)]  
    SNR_Y = [[1.e18,0.] for i in range(Nants)]  

    for bb in SNRBas:
      BSel = (bb[0]+1)*256 + (bb[1]+1)
      X1m = 1.e18 ; X2m = 1.e18
      Y1m = 1.e18 ; Y2m = 1.e18
    
      X1M = 0. ; X2M = 0.
      Y1M = 0. ; Y2M = 0.

      for IFp in range(NIF): #doIF:

        MASK = np.logical_and(fringe2["BAS"]==BSel,fringe2["IF"]==IFp)

        if np.sum(MASK)>0:
          P1 = np.array([chr(x) for x in fringe2["POL1"][MASK]])
          P2 = np.array([chr(x) for x in fringe2["POL2"][MASK]])
          XXp = np.logical_and(np.logical_or(P1=='R',P1=='X'),np.logical_or(P2=='R',P2=='X'))
          XYp = np.logical_and(np.logical_or(P1=='R',P1=='X'),np.logical_or(P2=='L',P2=='Y'))
          YXp = np.logical_and(np.logical_or(P1=='L',P1=='Y'),np.logical_or(P2=='R',P2=='X'))
          YYp = np.logical_and(np.logical_or(P1=='L',P1=='Y'),np.logical_or(P2=='L',P2=='Y'))

          XXmatrix = np.fft.fftshift(np.abs(np.fft.fft2(fringe2["VISIB"][MASK,:][XXp,:-1])))
          XYmatrix = np.fft.fftshift(np.abs(np.fft.fft2(fringe2["VISIB"][MASK,:][XYp,:-1])))
          YXmatrix = np.fft.fftshift(np.abs(np.fft.fft2(fringe2["VISIB"][MASK,:][YXp,:-1])))
          YYmatrix = np.fft.fftshift(np.abs(np.fft.fft2(fringe2["VISIB"][MASK,:][YYp,:-1])))

          Peak = np.unravel_index(np.argmax(XXmatrix),np.shape(XXmatrix))
          SNR_XX = XXmatrix[Peak[0],Peak[1]]
          XXmatrix[Peak[0]-1:Peak[0]+1,Peak[1]-1:Peak[1]+1] = 0.0; 

          Peak = np.unravel_index(np.argmax(XXmatrix),np.shape(XYmatrix))
          SNR_XY = XYmatrix[Peak[0],Peak[1]]
          XYmatrix[Peak[0]-1:Peak[0]+1,Peak[1]-1:Peak[1]+1] = 0.0; 

          Peak = np.unravel_index(np.argmax(XXmatrix),np.shape(YXmatrix))
          SNR_YX = YXmatrix[Peak[0],Peak[1]]
          YXmatrix[Peak[0]-1:Peak[0]+1,Peak[1]-1:Peak[1]+1] = 0.0; 

          Peak = np.unravel_index(np.argmax(XXmatrix),np.shape(YYmatrix))
          SNR_YY = YYmatrix[Peak[0],Peak[1]]
          YYmatrix[Peak[0]-1:Peak[0]+1,Peak[1]-1:Peak[1]+1] = 0.0; 

# Compute the std:
          SNR_XX /= np.sqrt(np.var(XXmatrix) + np.average(XXmatrix)**2.)
          SNR_XY /= np.sqrt(np.var(XYmatrix) + np.average(XYmatrix)**2.)
          SNR_YX /= np.sqrt(np.var(YXmatrix) + np.average(YXmatrix)**2.)
          SNR_YY /= np.sqrt(np.var(YYmatrix) + np.average(YYmatrix)**2.)

          X1m = np.min([SNR_XX, SNR_XY, X1m]) 
          X2m = np.min([SNR_XX, SNR_YX, X2m]) 
          Y1m = np.min([SNR_YY, SNR_YX, Y1m]) 
          Y2m = np.min([SNR_YY, SNR_XY, Y2m]) 

          X1M = np.max([SNR_XX, SNR_XY, X1M]) 
          X2M = np.max([SNR_XX, SNR_YX, X2M]) 
          Y1M = np.max([SNR_YY, SNR_YX, Y1M]) 
          Y2M = np.max([SNR_YY, SNR_XY, Y2M]) 

          SNR_X[bb[0]][0] = np.min([SNR_X[bb[0]][0],X1m])
          SNR_X[bb[0]][1] = np.max([SNR_X[bb[1]][1],X1M])
          SNR_X[bb[1]][0] = np.min([SNR_X[bb[1]][0],X2m])
          SNR_X[bb[1]][1] = np.max([SNR_X[bb[1]][1],X2M])

          SNR_Y[bb[0]][0] = np.min([SNR_Y[bb[0]][0],Y1m])
          SNR_Y[bb[0]][1] = np.max([SNR_Y[bb[1]][1],Y1M])
          SNR_Y[bb[1]][0] = np.min([SNR_Y[bb[1]][0],Y2m])
          SNR_Y[bb[1]][1] = np.max([SNR_Y[bb[1]][1],Y2M])

          del P1, P2, XXp, YYp, XYp, YXp, XXmatrix, YYmatrix, XYmatrix, YXmatrix

        #else:
        #  print("Problem with baseline %i-%i, IF %i"%(bb[0],bb[1],IFp))
          
        del MASK

    del fringe2 
  
    IFF = open(dd)
    lines = IFF.readlines()
    IFF.close()
    for li, line in enumerate(lines):
      if line.startswith('NUM SOURCES:'):
        lsou = li
        Nsou = int(line.split()[-1])
        break
  
    for j in range(Nsou):
      Snam = lines[lsou+1+j*5].split()[-1]
      print('%s:   %s'%(os.path.basename(dd).split('.')[0], Snam), file=OUTPUT)
    
    IsGood = True    
    SNROut = []  
    for j in range(Nants):
      Observes = '+'
      if SNR_X[j][1]<SNRCut or SNR_Y[j][1]<SNRCut:
          IsGood = False
          if (SNR_X[j][1]==0.00 and SNR_Y[j][1]==0.0): Observes = '-'
      SNROut += [j, SNR_X[j][0], SNR_X[j][1], SNR_Y[j][0], SNR_Y[j][1], Observes]
    SNROut += [{True:'Y',False:'N'}[IsGood]]
  
    print(fmt%tuple(SNROut),file=OUTPUT)
    
  OUTPUT.close()  




  if os.path.exists('SOURCE_SCANNER.FAILED'):   
    os.system('rm -rf SOURCE_SCANNER.FAILED')   

 except:

  e = sys.exc_info()[0]  
  OFF = open('SOURCE_SCANNER.FAILED','w')
  print(e,file=OFF)
  OFF.close()

  
