
######################################
# SCRIPT TO POLCONVERT IN NON-ALMA (A.K.A. EVN) MODE.
# I. Marti-Vidal, June 2019
######################################






from PolConvert import polconvert_CASA as PCONV
import pickle as pk


REFANT = 'O8' # Antenna to which refer the conversion gain solution (O8)
LINANT = 'EF' # Antenna with linear feed (EF)

REF_IDI = 'eo014_1_1.IDI6' # IDI with calibrator scan
CALRANGE = [0,23,28,0,0,23,39,45] # time range of calibrator scan (J0927+3902)
NCHAN = 32 # Number of channels per IF (to compute the multi-band delay)
NIF = 8  # Number of IFs.

NITER = 1  # Number of X-Y phase-estimate iterations (just 1 should suffice)

# List with the names of all FITS-IDIs to convert:
ALL_IDIs = ['eo014_1_1.IDI6']

# Add this suffix to the converted IDIs (empty string -> overwrite!)
SUFFIX = '.POLCONVERT'

import os
import numpy as np


# Initial gain estimates (dummy gains):
#EndGainsAmp = [1.0 for i in range(NIF)]
#EndGainsPhase = [0.0 for i in range(NIF)]
#TotGains = []


# Estimate cross-gains with PolConvert:
for i in range(NITER):


##################################
# Convert XYadd and XYratio to lists, in order
# to avoid corruption of the *.last file
  if i==0:
    XYadd = {}
    XYratio = {}
  else:
    IFF = open('PolConvert.XYGains.dat','rb')
    XYG = pk.load(IFF)
    XYadd = XYG['XYadd']
    XYratio = XYG['XYratio']
##################################


## WORKS: doSolve = 0.01  solint = [1,30] (OLD??)

  GainsIt = PCONV.polconvert(IDI=REF_IDI,
             OUTPUTIDI=REF_IDI,
             linAntIdx=[LINANT],
             plotIF = [],
             correctParangle = False,
             doIF = list(range(1,NIF+1)),
             XYadd = XYadd,
             XYratio = XYratio,
             Range = CALRANGE,
             plotRange = CALRANGE,
             IDI_conjugated = False,
             doSolve = 0.1,   #4,
             solint = [1,30],  #BP MODE
             plotAnt = REFANT,
             amp_norm = 0.0,
             excludeAnts = ["HH","T6"] , #8,9],
             doTest=True)
  
  os.system('rm -rf FRINGE.PLOTS.ITER%i'%i)
  os.system('mv FRINGE.PLOTS FRINGE.PLOTS.ITER%i'%i)
  os.system('mv Cross-Gains.png Cross-Gains.ITER%i.png'%i)

# HERE WE CAN CONVERT ALL IDIs:
if False:
 for IDI in ALL_IDIs:
  polconvert(IDI=IDI,
             OUTPUTIDI=IDI+SUFFIX,
             linAntIdx=[LINANT],
             correctParangle = False,
             doIF = list(range(1,NIF+1)),
             doSolve = -1,
             plotIF = [],
             XYadd = GainsIt['XYadd'],
             XYratio = GainsIt['XYratio'],
             IDI_conjugated = False,
             amp_norm = 0.0,
             doTest=False)











