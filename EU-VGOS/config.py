#########################
### CONFIG FILE FOR EU-VGOS POLCONVERSION PIPELINE (MASTER SCRIPT)
### VERSION 24/04/2026.
### E. ALBENTOSA-RUIZ
#########################
###
### PROCEDURE (PARTIALLY AUTOMATED):
### 0.- COPY MASTER SCRIPT (AND CONFIG FILE) INTO A DEDICATED DIRECTORY
###     (e.g., PARENT DIRECTORY OF DiFX DATA DIRECTORY FOR THE SESSION).
### 1.- SET THE EXPERIMENT PARAMETERS:
###     - EXPNAME, EXP_CODE, EXP_YEAR
###     - DIFX_DIR (ORIGINAL SWIN DATA), PCONV_DIR (POLCONVERT OUTPUT)
###     - DOIF, IF_OFFSET
###     - POLCAL_SCAN, REFANT
###     - ANTENNA / PCAL / FLAGGING CONTROLS
###     - (OPTIONAL) ANT_COORDS (IF NOT PROVIDED, THEY CAN BE AUTO-FILLED
###                              FROM VEX, IF AVAILABLE ON THE CDDIS GNSS PRODUCTS)
###     - (OPTIONAL) XPOL_DELAYS (ONLY IF AUTO_XPOL_DELAY = False, OR FOR BACKUP A PRIORI)
### RUN MASTER SCRIPT (e.g. python master_script_python3_MP.py config_myexp.py)
###
### Recommended run sequence:
### 2.- RUN MASTER SCRIPT WITH mysteps = [0].
###     RUNS SOURCE_SCANNER AND WRITES:
###       - SOURCES_<EXPNAME>.txt
###       - <EXPNAME>_ScanStats.dat/.txt
###       - <EXPNAME>_PipelineAuto.dat
###     THEN:
###       - IF POLCAL_SCAN IS EMPTY, IT AUTO-FILLS POLCAL_SCAN FROM <EXPNAME>_ScanStats.dat
###         AND ALSO SETS ADDITIVE_PHASE_SCANS = POLCAL_SCAN.
###       - AUTO-CLEANS ANTENNA / BASELINE / PCAL SETTINGS.
###
### 3.- RUN MASTER SCRIPT WITH mysteps = [1].
###     ESTIMATES CROSS-POL GAINS FROM POLCAL_SCAN:
###       - CONCATENATES POLCAL_SCAN SWIN DATA INTO: <PCONV_DIR>/<EXPNAME>_PC_CALIB/
###       - RUNS POL_CALIBRATE PER IF
###       - WRITES POLCAL_GAINS_<EXPNAME>.dat
###       - IF TWO_ITER_STEP1=True:
###             FIRST RUN DIAGNOSES POORLY CONSTRAINED ANTENNAS FROM FRINGE PEAKS,
###             THEN SOLVES WELL-BEHAVED ANTENNAS AND RECOVERS PROBLEMATIC ANTENNAS.
###       - IF AUTO_XPOL_DELAY=True:
###             ESTIMATES THE BEST A PRIORI XPOL_DELAYS FROM THE DATA
###             AND WRITES XPOL_DELAYS_<EXPNAME>.txt/.dat
###       - UPDATES <EXPNAME>_PipelineAuto.dat.
###
### IF THE ASSESSMENT PLOTS ARE OK:
### 4.- RUN MASTER SCRIPT WITH mysteps = [2].
###     POLCONVERTS THE WHOLE EXPERIMENT, USING:
###       - POLCAL_GAINS_<EXPNAME>.dat
###       - XPOL_DELAYS FROM <EXPNAME>_PipelineAuto.dat OR USER CONFIG
###     OUTPUT DATA ARE WRITTEN TO:
###       - PCONV_DIR/
###
### 5.- RUN mysteps = [3] TO PREPARE A CF FILE FOR ALIGNMENT.
### 6.- RUN mysteps = [4] TO REMOVE IONEX-BASED TEC / OPTIONAL BANDPASS CALIBRATION.
### 7.- RUN mysteps = [5] TO PERFORM BROADBAND GLOBAL FRINGE FITTING.
### 8.- RUN mysteps = [6] TO WRITE FULLY CALIBRATED SWIN FILES.
###
#########################
### WHAT EACH STEP READS/WRITES QUICK REFERENCE
### Step 0 SOURCE_SCANNER:
###   reads:  DIFX_DIR/
###   writes: SOURCES_<EXPNAME>.txt, <EXPNAME>_ScanStats.dat/.txt,
###           <EXPNAME>_PipelineAuto.dat
###
### Step 1 POL_CALIBRATE + AUTO_XPOL_DELAY:
###   reads:  DIFX_DIR/, selected POLCAL_SCAN, <EXPNAME>_PipelineAuto.dat
###   writes: POLCAL_GAINS_<EXPNAME>.dat,
###           XPOL_DELAYS_<EXPNAME>.txt/.dat if AUTO_XPOL_DELAY,
###           updated <EXPNAME>_PipelineAuto.dat,
###           diagnostic plots/logs
###
### Step 2 POLCONVERTER:
###   reads:  DIFX_DIR/, POLCAL_GAINS_<EXPNAME>.dat, XPOL_DELAYS
###   writes: PCONV_DIR/
###
### Step 3 GET_FOURFIT_PHASES -> CF file:
###   reads:  PCONV_DIR/, ADDITIVE_PHASE_SCANS
###   writes: CF_FILENAME, PYPHASES.PLOTS/
###
### Step 4 removeTEC / BP+TEC calibration:
###   reads:  PCONV_DIR/
###   writes: BPCAL_DIR/
###
### Step 5 DO_GFF:
###   reads:  BPCAL_DIR/, CF_FILENAME
###   writes: GFF products/logs
###
### Step 6 WRITE_CALIBRATED export:
###   reads:  PCONV_DIR/, CF_FILENAME and optionally GFF products
###   writes: FINAL_DIR/
#########################

# =============================================================
# 1. STEPS / RESOURCES
# =============================================================

# List of steps to be performed, e.g. [0], [1], [2], [0,1,2], ...
# NOTE:
#   - mysteps=[1] expects <EXPNAME>_PipelineAuto.dat from Step 0,
#     unless Step 0 is included in the same run.
#   - Steps 2/4/5/6 expect SOURCES_<EXPNAME>.txt to exist from Step 0.
mysteps = [0,1]

# Number of processes allowed to run simultaneously:
NCPU = 15


# =============================================================
# 2. EXPERIMENT IDENTITY AND PATHS
# =============================================================

# Experiment identifiers:
# - EXPNAME: DiFX prefix for scan folders/files
# - EXP_CODE/EXP_YEAR: used for remote metadata helpers
EXPNAME = 'vo5120'
EXP_CODE = 'vo5120'
EXP_YEAR = 2025

# Directory containing ORIGINAL DiFX SWIN data:
DIFX_DIR = 'swin'

# Output directory for POLCONVERTED data:
PCONV_DIR = 'DiFX_PCONV'

# Destination directory of BP+IONEX corrected data (step 4):
BPCAL_DIR = 'DiFX_TECOR'

# Output directory for fully calibrated data products (step 6):
FINAL_DIR = "DiFX_FINAL"


# =============================================================
# 3. FREQUENCY / IF SETUP
# =============================================================

# List of IFs to process, starting from 1.
# VGOS typical: range(1,33). Adjust for selected IFs / setup.
DOIF = range(1,33) # [5, 13, 21, 28] # range(1,33)

# If autocorrelations are stored in different IFs, set IF offset:
# some setups use IF -> IF + IF_OFFSET mapping for certain products.
# Typical VGOS dual-pol setups: 32
IF_OFFSET = 32

# List of IFs to explore in SOURCE_SCANNER.
# Include at least one per band.
SCAN_IF = [5, 13, 21, 28]

# Number of IFs per band. If scanning only selected IFs, this is used
# to propagate detected problems to the full band.
NIF_PER_BAND = 8

# IF labels used by fourfit, string in frequency order.
IFNAMES = "abcdefghijklmnopqrstuvwxyzABCDEF"


# =============================================================
# 4. STEP 0 SCAN SELECTION / POLCAL_SCAN
# =============================================================

# POLCAL scans used for polarization calibration.
# Leave empty to auto-fill from <EXPNAME>_ScanStats.dat produced by Step 0.
# Scan names are strings with leading zeros, e.g. "0002".
POLCAL_SCAN = []
# Example:
# POLCAL_SCAN = ["1045", "1079", "1121", "1129", "1509", "1992", "1995"]

# Optional suffix used by POLCONVERTER for outputs / file naming.
SUFFIX = ''

# Main scan-selection parameter: minimum per-antenna SNR required.
MIN_SCAN_SNR = 15.

# Minimum elevation allowed, degrees:
MIN_ELEV_DEG = 10.

# Minimum uv-span for a scan, meters:
MIN_UV_SPAN = 1e6

# Minimum number of antennas for a scan:
MIN_ANTS_PRESENT = 5

# Minimum number of unique baselines required for a scan:
MIN_UNIQUE_BASELINES = MIN_ANTS_PRESENT * (MIN_ANTS_PRESENT - 1) // 2

# Target parallactic-angle coverage over selected calibration scans, degrees:
MIN_TOTAL_PA_COVERAGE_DEG = 30.

# Target number of PA-usable calibration scans per antenna:
PA_SCANS_REQUIRED = 2


# =============================================================
# 5. ANTENNA METADATA
# =============================================================

# Antenna coordinates in meters, ECEF: antenna name -> (X, Y, Z)
# NOTE:
#  - If empty {}, master script will try to download/read the VEX file.
#  - Any antenna missing from ANT_COORDS may be added from the VEX file.
ANT_COORDS = {
    "GS": (1130729.877, -4831245.972, 3994228.300),   # GGAO12M
    "HV": (5085502.673, 2668368.171, -2768493.227),   # HARTVGS
    "HB": (-3949991.094, 2522421.259, -4311707.721),  # HOBART12
    "IS": (-3959636.203, 3296825.448, 3747042.571),   # ISHIOKA
    "KE": (-4147354.925, 4581542.266, -1573302.777),  # KATH12M
    "K2": (-5543831.745, -2054585.590, 2387828.974),  # KOKEE12M
    "MG": (-1330788.462, -5328106.593, 3236427.492),  # MACGO12M
    "NN": (1200992.670, 252098.630, 6238038.646),     # NYALE13N
    "NS": (1200992.670, 252098.630, 6238038.646),     # NYALE13S
    "MB": (4641823.841, 1393011.179, 4133445.299),    # MATVGOS
    "SA": (4618524.302, -2166020.720, 3816270.345),   # RAEGSMAR
    "OE": (3370889.298, 711571.199, 5349692.048),     # ONSA13NE
    "OW": (3370946.779, 711534.507, 5349660.925),     # ONSA13SW
    "YJ": (4848831.021, -261629.388, 4122976.576),    # RAEGYEB
    "S6": (-2831646.947, 4675729.541, 3275365.102),   # SESHAN13
    "T1": (-2826801.860, 4679253.933, 3274516.038),   # TIANMA13
    "UM": (228671.597, 4631855.932, 4367130.506),     # URUMQI13
    "WF": (1492206.223, -4458130.552, 4296015.629),   # WESTFORD
    "WN": (4075627.541, 931774.394, 4801552.451),     # WETTZ13N
    "WS": (4075658.807, 931824.883, 4801516.273),     # WETTZ13S
    "YG": (-2388896.500, 5043350.051, -3078590.462),  # YARRA12M
}

# Antenna name translation SWIN -> fourfit/HOPS naming.
# If .corr can be downloaded, master script may fill/update this.
HOPSNAMES = {}

# Sampler delays in X, ordered by band from low to high.
# If .corr can be downloaded, master script may fill/update this.
SAMP_DELAYS = {}


# =============================================================
# 6. ANTENNA / BASELINE / PHASECAL FLAGGING CONTROLS
# =============================================================

# Reference antenna.
# Used for diagnostic plots of POL_CALIBRATE/POLCONVERTER.
# If REFANT is set to a station name, it is always protected and forced:
#   Step 0 will use it as the only REFANT candidate,
#     (after ensuring it participates in the experiment
#      and is not excluded / BAD_IF-flagged)
# If REFANT = None or REFANT = "": 
#   Step 0 selects REFANT automatically according to REFANT_SELECTION_MODE
REFANT = ''
#
# REFANT auto-selection mode:
#   "max_connectivity": choose REFANT that connects most stations;
#                       tie-break by tone/SNR score.
#   "highest_snr": choose clean station with highest SNR,
#                  then add scans needed to connect to other stations.
REFANT_SELECTION_MODE = "highest_snr"

# List of possible reference antennas for GFF, in order of preference.
# If empty, Step 0 auto-fills using best SNR / valid antennas.
GFF_REFANTS = []

# Antennas globally excluded from solving.
# Step 0 will automatically add stations with no usable scan / bad diagnostics.
EXCLUDE_ANTENNA = []

# Baselines globally excluded from solving.
# Step 0 will clean invalid baselines and can auto-add very short intra-site baselines.
EXCLUDE_BASELINE = [['WN', 'WS'], ['OE', 'OW']]

# Individual bad IFs per antenna, dropped from Step 1 solving.
# Format: {"ANT": [if1, if2, ...]} or {"ANT": range(...)}
# Step 0 may automatically fill/update this.
BAD_IF = {}

# Minimum number of pcal tones per IF required to consider an antenna tone-rich.
TONE_POOR_PCAL_PER_IF = 5

# Make optional pcal tone inspection plots in Step 1:
DO_PCAL_TONE_PLOTS = True

# Per-antenna whether phasecals are used in solving/application logic.
# Step 0 cleans this to participating antennas and defaults missing antennas to True.
USE_PCAL = {
    'OE': True, 'OW': True, 'YJ': True, 'IS': True, 'WS': True,
    'WF': True, 'GS': True, 'K2': True, 'MG': True, 'NN': True,
    'HB': True, 'HV': True, 'KE': True, 'WN': True, 'S6': True,
    'SA': True, 'UM': True, 'T1': True, 'YG': True,
}

# TODO: may be possible to fill automatically both ZERO_PCALS and FLAG_PCALS
# Frequency ranges MHz where pcals are zeroed, e.g. RFI / known bad ranges.
# Format: {'ANT': [[fmin_MHz, fmax_MHz], ...], ...}
ZERO_PCALS = {}
# Flag specific problematic pcal tones per antenna.
# Content is tone MHz frequencies; project-dependent.
FLAG_PCALS = {
    # 'GS': [3090, 3070, 3075],
    # 'OE': [6425, 6430, 5445],
    # 'MG': [5715, 5720, 5615, 5620],
    # 'OW': [10265, 10270, 3215, 10255, 10625],
    # 'WF': [3215, 3310, 3375, 3355],
    # 'YJ': [6650],
}
# TODO: FLAG_PCALS is supposedly set before step 3, but only after step 3 we can plot residual phasecal

# Per-antenna GFF weights. Lower = downweight bad stations. Default is 1.0
# Master script will set weight to 0 if station is in EXCLUDE_ANTENNA.
ANT_WEIGHTS = {}
# TODO: COULD BE SET AUTOMATICALLY IN SOME WAY? We can identify problematic stations based on final fringe peaks, or there may be an even better way to evaluate. Only needed for WBGFF, not PolConvert 

# =============================================================
# 7. STEP 1 POLCAL / XPOL / PCAL SETTINGS
# =============================================================

# Run XPOL_DELAY auto-estimation in Step 1?
AUTO_XPOL_DELAY = False # True

# XPOL_DELAYS: per antenna [tau_s, f_zero_hz]
#   tau_s: cross-pol delay in seconds
#   f_zero_hz: frequency Hz where fitted cross-pol phase is 0
#
# Usage rules:
# (1) To use trusted XPOL delays:
#     AUTO_XPOL_DELAY = False
#     XPOL_DELAYS = {'ANT': [tau_s, f_zero_hz], ...}
# (2) To auto-estimate XPOL delays:
#     AUTO_XPOL_DELAY = True
#     XPOL_DELAYS will be ignored/reset by the master script for Step 1.
# (3) Step 1 writes:
#     XPOL_DELAYS_<EXPNAME>.txt/.dat and updates <EXPNAME>_PipelineAuto.dat.
XPOL_DELAYS = {
    'GS': [1.554685872057e-09, 7.080205000744e+09],
    'HV': [6.511773286101e-11, -1.591050369428e+08],
    'IS': [3.255886643051e-11, 1.503498469265e+10],
    'K2': [1.465148989373e-10, 6.870641300314e+09],
    'KE': [-8.139716607627e-12, 6.566930935278e+10],
    'NN': [-9.767659929152e-11, 5.265099514968e+09],
    'SA': [-8.139716607627e-12, 5.639814523177e+10],
    'T1': [-3.744269639508e-10, 5.740519122339e+09],
    'WF': [1.457009272765e-09, 6.845755053340e+09],
    'WN': [-8.139716607627e-12, 6.245897449949e+10],
    'WS': [5.697801625339e-11, 2.281217393587e+09],
}
#    'GS': [1.554685872057e-09, 7.076680440611e+09],
#    'IS': [3.255886643051e-11, 1.461937733008e+10],
#    'K2': [1.465148989373e-10, 6.982025447199e+09],
#    'OE': [-1.139560325068e-10, 9.087069826677e+09],
#    'OW': [2.197723484059e-10, 6.735863392582e+09],
#    'S6': [-4.346608668473e-09, 6.899552301369e+09],
#    'T1': [-4.069858303813e-10, 6.047530363822e+09],
#    'UM': [-3.255886643051e-11, 8.516617006594e+09],
#    'WF': [1.457009272765e-09, 6.850618288019e+09],
    #'WN': [-8.139716608e-12, 6.495866856e+10]
#    'WS': [4.883829964576e-11, -1.733838703760e+08],
#    'YG': [3.988461137737e-10, 5.588542748356e+09],
#    'YJ': [3.744269639508e-10, 6.277601567673e+09],

# Run Step 1 in two stages:
#   RUN1: diagnose problematic stations from fringe peaks.
#   SET1: solve well-behaved stations while excluding problematic ones.
#   SET2: recover problematic stations with well-behaved stations fixed.
TWO_ITER_STEP1 = True

# If True, allow poorly constrained/problematic antennas to use PCAL again in SET2.
# If False, keep PCAL disabled for antennas classified as poorly constrained.
REENABLE_BAD_PCAL_USE = True

# Conditions to detect problematic stations from RUN1 fringe peaks.
# Maximum cross-hand fraction relative to strongest parallel-hand term:
CROSS_HAND_MAX_FRAC = 0.475
# Maximum fraction of failed IFs tolerated before diagnosing station as problematic:
MAX_BAD_IF_FRAC = 0.475

# Frequency averaging in channels used for solving cross-pol gains.
# Larger => smoother/more robust, but less spectral detail.
CHAVG = 12

# Pre-averaging time in seconds.
# Larger => higher SNR, but risks decorrelation if phase changes quickly.
INTTIME = 20.0

# Apply amplitude gains in POLCONVERTER?
# If False, POLCAL gains are forced to unity amplitude for application.
APPLY_AMP = True

# Solve for amplitude gains in POL_CALIBRATE?
SOLVE_AMP = True

# If TWO_ITER_STEP1=True, solve amplitude on SET2?
SET2DOAMP = True

# How to interpolate X-Y pcal differences:
# "multitone" typical VGOS, or "bandpass".
XYPCALMODE = "multitone"

# Solver backend used inside POL_CALIBRATE.
# Others exist but may be unstable/bad: BFGS, Newton-CG, SLSQP.
SOLVER = "COBYLA"

# Use delay and/or rate to help gain estimates?
USE_DELAY = False # True # False
USE_RATE  = False

# UV taper in meters for X-pol gain estimator.
# Positive value lowers weight of long baselines.
# Set as negative value to flag long baselines.
UVTAPER = 1e9

# Minimum gap Hz between VGOS bands for band-aware logic:
BAND_GAP_SPLIT_HZ = 0.5e9

# Width of median window for autocorrelations, used to remove phasecals:
PCAL_MED_WINDOW = 4


# =============================================================
# 8. STEP 3 / STEP 4 / STEP 5 / STEP 6 CALIBRATION SETTINGS
# =============================================================

# Scans used to align phases among IFs.
# Usually POLCAL_SCAN; Step 0 auto-sets this when POLCAL_SCAN is auto-filled.
ADDITIVE_PHASE_SCANS = POLCAL_SCAN

# Output CF filename:
CF_FILENAME = "cf_PyPhases"

# Calibrate complex bandpass during processing?
CALIB_BPASS = False

# Max residual phasecal RMS in degrees used for automatic IF flagging:
MAX_PCAL_RMS = 10.0

# TODO: COULD BE SET UP FROM PCAL INVESTIGATION RESULTS
# Frequency range MHz + IF labels where instrumental delay is extrapolated.
# Used when pcals are unusable, e.g. RFI wipes a band.
# Format: {"ANT": [fmin_MHz, fmax_MHz, "iflabels"]}
PCAL_DELAYS = {}

# Ionosphere sun follow-up factor.
# Should be between 0 and 1; 1.0 is the sensible value.
LDFAC = 1.0

# IONEX model tag.
# The master script will try to download the IONEX map starting from this model;
# if not available, it may try others.
IONEX = "codg"

# Apply phasecals to new data in TEC/bandpass steps?
APPLY_PHASECAL = False
