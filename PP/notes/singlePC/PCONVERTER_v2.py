from PolConvert import polconvert_standalone as PCONV
import os
import sys

#FILE = '../b1bsnip/ta037b-2-b1_3060'
FILE = 'ta037b-2-b1_3060'

doSolve = False
doApply = False
if len(sys.argv) == 1:
    print('one argument: solve or apply')
    sys.exit(0)
elif sys.argv[1] == 'solve':
    doSolve = True
elif sys.argv[1] == 'apply':
    doApply = True
else:
    print('one argument: solve or apply')
    sys.exit(0)

if doSolve:
  os.system('cp -r %s.calc kakalamaka.calc'%FILE)
  RESULT = PCONV.polconvert(IDI = FILE+'.difx',
               OUTPUTIDI = 'kakalamaka.difx',
               DiFXinput = FILE + '.input',
               DiFXcalc = FILE + '.calc',
               doIF = list(range(37,43)),
               linAntIdx = [1],
               plotIF = list(range(37,43)),
               plotRange = [0,0,0,0,2,0,0,0],
               plotAnt = 8,
## This number works well. Higher values make the AA gains noisier
## While lower values mya introduce noisier gains for the other 
## antennas and spurious 180 deg jumps for AA:
               doSolve = 0.1,
## First number is channel averaging, second number is pre-averaging
## time (in seconds). IF playing with continuum data, the first number
## should be MUCH smaller:
               solint = [64,20],
               useDelays=True,
               doTest = True)



if doApply:
  import pickle as pk
  IFF = open("PolConvert.XYGains.dat","rb")
  RESULT = pk.load(IFF)
  IFF.close()
  RESULT = PCONV.polconvert(IDI = FILE+'.difx',
               OUTPUTIDI = 'kakalamaka.difx',
               DiFXinput = FILE + '.input',
               DiFXcalc = FILE + '.calc',
               doIF = list(range(37,43)),
               linAntIdx = [1],
               XYadd=RESULT["XYadd"],
               XYratio=RESULT["XYratio"],
               plotIF = list(range(37,43)),
               plotRange = [0,0,0,0,2,0,0,0],
               plotAnt = 8,
               doSolve = -1,
               doTest = False)

