import pickle as pk
import os, sys
from PolConvert import polconvert_standalone as PC



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
def POL_CALIBRATE(EXPNAME='',DIFX_DIR = '', DOSCAN=-1,CHANSOL=32,USE_PCAL=True,EXCLUDE_BASELINES=[],DOIF=[],DOAMP=True,PLOTANT=1,APPLY_AMP=True):
 """ Estimates cross-polarization gains using a scan in a SWIN directory.
  The channel resolution is set to CHANSOL. Saves the gains in a dictionary
  that can ge used by PolConvert."""
#################################

# try:
 if True:

  EXP = EXPNAME

  if not (type(DOSCAN) is list):
    DOSCAN = [DOSCAN]

  if len(DOSCAN)==0 or int(DOSCAN[0]) < 0:
    raise Exception("POL_CALIBRATE ERROR! SCAN(S) NOT SPECIFIED!\n")

  DIFX = DIFX_DIR



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
    DOIF = list(range(1,NIF+1))


############################################################
# X-Y cross-pol gain estimate (with phasecal correction):

  os.system('rm -rf %s'%os.path.join(DIFX,'%s_PC_CALIB'%EXP))
  os.system('mkdir %s'%os.path.join(DIFX,'%s_PC_CALIB'%EXP))
  for SI in DOSCAN:
    os.system('cp -r %s %s/.'%('%s_%s.difx'%(os.path.join(DIFX,EXP),SI), os.path.join(DIFX,'%s_PC_CALIB'%EXP)))  
    os.system('cp -r %s %s/.'%('%s_%s.calc'%(os.path.join(DIFX,EXP),SI), os.path.join(DIFX,'%s_PC_CALIB'%EXP)))  

  
  WITH_PCAL = PC.polconvert(IDI='%s/%s_PC_CALIB'%(DIFX,EXP),
    OUTPUTIDI = '%s_POL_CALIBRATE_BLIND'%(os.path.join(DIFX,EXP)),
    DiFXinput = '%s_%s.input'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    DiFXcalc = '%s_%s.calc'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    doIF = DOIF, plotIF = DOIF,
    plotRange=[0,0,0,0,2,0,0,0], plotAnt=PLOTANT,excludeBaselines=EXCLUDE_BASELINES,
    linAntIdx = list(range(1,Nants+1)),swapXY=[False for i in range(Nants)], 
    usePcal = USE_PCAL, XYadd={}, XYratio={}, XYdel={},
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

  if not APPLY_AMP:
    for anti in WITH_PCAL['XYratio'].keys():
      NIF = len(WITH_PCAL['XYratio'][anti])
      for ki in range(NIF):
        WITH_PCAL['XYratio'][anti][ki][:] = 1.0

  OFF=open('POLCAL_OUTPUT_SCAN-%s.dat'%DOSCAN[0],'wb')
  pk.dump(WITH_PCAL,OFF,protocol=0) ; OFF.close()


  #IFF = open('POLCAL_OUTPUT_SCAN-%03i.dat'%DOSCAN)
  #WITH_PCAL = pk.load(IFF)
  #IFF.close()
  
  FINAL = PC.polconvert(IDI='%s/%s_PC_CALIB'%(DIFX,EXP),
    OUTPUTIDI = '%s_POL_CALIBRATE_RESULTS'%(os.path.join(DIFX,EXP)),          
    DiFXinput = '%s_%s.input'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    DiFXcalc = '%s_%s.calc'%(os.path.join(DIFX,EXP),DOSCAN[0]),
    doIF = DOIF, plotIF = DOIF,
    plotRange=[0,0,0,0,2,0,0,0], plotAnt=PLOTANT,excludeBaselines=EXCLUDE_BASELINES,
    linAntIdx = list(range(1,Nants+1)),swapXY=[False for i in range(Nants)], 
    usePcal = USE_PCAL, 
    XYadd=WITH_PCAL['XYadd'], XYratio=WITH_PCAL['XYratio'], XYdel={},
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


