
######################################
# SCRIPT TO POLCONVERT IN NON-ALMA (A.K.A. EVN) MODE.
# I. Marti-Vidal, June 2019
######################################


# NAME OF DIFX FILE (CAL SCAN):
CALIDI = 'e18c21-0-b2_1022.difx'

# NAMES OF ALL DIFX TO CONVERT:
IDIs = ['e18c21-0-b2_1022.difx']


# RANGE OF IFs TO CONVERT:
IF = range(36,68)


# METHOD TO MINIMIZE THE ERROR FUNCTION (EXPERIMENTAL BEING ADDED):
METHOD = "COBYLA"



# INDEX (INDICES) OF THE LINEAR-POLARIZATION ANTENNA(S).
# THE FIRST ANTENNA HAS INDEX 1:
LINANT = [1]


# ANTENNA WHOSE BASELINES ARE TO BE PLOTTED FOR ASSESSMENT:
PLOTANT = 2


# GCPFF LAGRANGE MULTIPLIER (ZERO WILL NOT USE RR/LL, WHICH
# MAY INTRODUCE PI AMBIGUITIES IN THE XY PHASE):
SOLVEMOD = 0.2

# NUMBER OF CHANNELS TO AVERAGE WITHIN EACH IF FOR THE SOLUTION:
NCHAV = 8

# SEARCH WINDOW OF THE FRINGES IN DELAY-RATE MATRIC (PIXEL UNITS):
NPIX = 100

#################
# SCRIPT STARTS #
#################


# Remove de "difx" extension:
BASENAM = CALIDI[:-5]




# STEP 1. ESTIMATE THE CROSS-POL GAINS.
# THIS WILL RETURN A DICTIONARY (GAIN ENTRIES PER ANTENNA):
# RESULTS OF THE CALIBRATION IS PLACED IN RESULT WHICH
# IS A PYTHON DICTIONARY WHICH CONTAINS THE NUMPY ARRAYS
# THIS MUST BE PICKLED TO BREAK INTO TWO STEPS
# THE POLCONVERT OPTIONS MAY NEED TWEAKING BETWEEN
# ALMA AND EVN OR VGOS OR ... CASES
RESULT = polconvert(IDI  =  CALIDI,
  OUTPUTIDI          =  "TEST_EVN",
  DiFXinput          =  "%s.input"%BASENAM,
  DiFXcalc           =  "%s.calc"%BASENAM,
  doIF               =  list(IF),
  linAntIdx          =  LINANT,
  Range              =  [],
  ALMAant            =  "",
  spw                =  -1,
  calAPP             =  "",
  calAPPTime         =  [0.0, 8.0],
  APPrefant          =  "",
  gains              =  [['NONE']],
  interpolation      =  [[]],
  gainmode           =  [[]],
  XYavgTime          =  0.0,
  dterms             =  ['NONE'],
  amp_norm           =  1.0,
  XYadd              =  {},
  XYdel              =  {},
  XYratio            =  {},
  swapXY             =  [False],
  swapRL             =  False,
  IDI_conjugated     =  True,
  plotIF             =  list(IF),
  plotRange          =  [0, 0, 0, 0, 14, 0, 0, 0],
  plotAnt            =  PLOTANT,
  excludeAnts        =  [],
  doSolve            =  SOLVEMOD,
  solint             =  [NCHAV, 1],
  doTest             =  True,
  npix               =  NPIX,
  solveAmp           =  True,
  solveMethod        =  METHOD,
  calstokes          =  [1.0, 0.0, 0.0, 0.0],
  calfield           =  -1)


# STEP 2: CALIBRATE THE WHOLE EXPERIMENT:
# NOTE THAT AT THE MOMENT ONLY RESULT['XYadd'] IS USED
# WE NEED TO EXPERIMENT MORE HERE...
for DIFX in IDIs:
  BASENAM = DIFX[:-5]

  polconvert(IDI  =  DIFX,
    OUTPUTIDI          = '%s.difx_polconverted'%BASENAM ,
    DiFXinput          =  "%s.input"%BASENAM,
    DiFXcalc           =  "%s.calc"%BASENAM,
    doIF               =  list(IF),
  linAntIdx          = LINANT,
  Range              =  [],
  ALMAant            =  "",
  spw                =  -1,
  calAPP             =  "",
  calAPPTime         =  [0.0, 8.0],
  APPrefant          =  "",
  gains              =  [['NONE']],
  interpolation      =  [[]],
  gainmode           =  [[]],
  XYavgTime          =  0.0,
  dterms             =  ['NONE'],
  amp_norm           =  1.0,
  XYadd              =  RESULT['XYadd'],
  XYdel              =  {},
  XYratio            =  RESULT['XYratio'],
  swapXY             =  [False],
  swapRL             =  False,
  IDI_conjugated     =  True,
  plotIF             =  list(IF),
  plotRange          =  [0, 0, 0, 0, 14, 0, 0, 0],
  plotAnt            =  PLOTANT,
  excludeAnts        =  [],
  doSolve            =  -1.0,
  solint             =  [1, 1],
  doTest             =  False,
  npix               =  100,
  solveAmp           =  False,
  calfield           =  -1)





