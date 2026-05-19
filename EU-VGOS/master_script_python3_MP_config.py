#########################
### MASTER SCRIPT FOR EU-VGOS POLCONVERSION.
### VERSION 15/05/2026.
### E. ALBENTOSA-RUIZ & I. MARTI-VIDAL & J. GONZALEZ
#########################

######################################
#####################################################
thesteps = {0: 'Source scanner',
            1: 'Estimate cross-polarization gains',
            2: 'PolConvert the whole experiment',
            3: 'Estimate additive phases & write CF file',
            4: 'Calibrate bandpass and remove IONEX TEC',
            5: 'Perform broadband Global Fringe Fitting',
            6: "Write fully calibrated SWIN files (for FITS/MS export only!)"}

PYTHON_CALL = 'python3 %s.py'
#####################################################
######################################

#################
# SCRIPT STARTS #
#################

import os
import pickle as pk
import numpy as np
import glob
import matplotlib.pyplot as pl
import multiprocessing
import sys
import pathlib
import importlib.util


# Allow: python master_script.py config_myexp.py
# Default config if none is provided:
DEFAULT_CONFIG = "config.py"
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CONFIG
cfg_path = pathlib.Path(CONFIG_FILE)
if not cfg_path.exists():
    raise FileNotFoundError(f"Config file not found: {CONFIG_FILE}")
spec = importlib.util.spec_from_file_location("cfg", str(cfg_path))
cfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cfg)
# Export cfg vars into globals
globals().update({k: getattr(cfg, k) for k in dir(cfg) if not k.startswith("_")})


########### +++++++++++++++++++++++++++++++++++++++++++++  ###########
import logging

MASTER_LOG = "%s_master_script_steps_%s.log_masterscript_run" % ( EXPNAME, "-".join(str(s) for s in mysteps) )

logger = logging.getLogger("MASTER")
logger.setLevel(logging.INFO)
logger.handlers[:] = []

_fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

_fh = logging.FileHandler(MASTER_LOG, mode="a")
_fh.setFormatter(_fmt)
logger.addHandler(_fh)

_ch = logging.StreamHandler(sys.stdout)
_ch.setFormatter(logging.Formatter("[MASTER] %(message)s"))
logger.addHandler(_ch)

def log_header(msg):
    logger.info("=" * 72)
    logger.info(msg)
    logger.info("=" * 72)

def log_kv(title, dct):
    logger.info(title)
    for k in sorted(dct):
        logger.info("  %s = %s", k, dct[k])

def log_and_print(msg, level="info"):
    print(msg)
    if level == "warning":
        logger.warning(msg)
    elif level == "error":
        logger.error(msg)
    else:
        logger.info(msg)

def unique_log_copy_name(base_name):
    if not os.path.exists(base_name):
        return base_name
    root, ext = os.path.splitext(base_name)
    i = 1
    while True:
        candidate = "%s_try%i%s" % (root, i, ext)
        if not os.path.exists(candidate):
            return candidate
        i += 1

########### +++++++++++++++++++++++++++++++++++++++++++++  ###########

PIPELINE_AUTO_FILE = "%s_PipelineAuto.dat" % EXPNAME

from EUVGOS_PY3.DOWNLOAD_EXP_INFO import (
    ensure_ant_coords_from_vex,
    ensure_hopsnames_from_corr,
    ensure_samp_delays_from_corr,
)
ANT_COORDS = ensure_ant_coords_from_vex(
    exp_code=EXP_CODE,
    exp_year=EXP_YEAR,
    user_ant_coords=ANT_COORDS,
)
HOPSNAMES  = ensure_hopsnames_from_corr(
    exp_code=EXP_CODE,
    exp_year=EXP_YEAR,
    user_hopsnames=HOPSNAMES,
)
SAMP_DELAYS = ensure_samp_delays_from_corr(
    exp_code=EXP_CODE,
    exp_year=EXP_YEAR,
    user_samp_delays=SAMP_DELAYS,
)

if 'XPOL_DELAYS' not in globals() or XPOL_DELAYS is None:
    XPOL_DELAYS = {}
# XPOL workflow rule
if AUTO_XPOL_DELAY and (1 in mysteps):
    XPOL_DELAYS = {}


########### +++++++++++++++++++++++++++++++++++++++++++++  ###########
log_header("EU-VGOS PolConversion pipeline started")
log_kv("Initial experiment configuration", {
    "EXPNAME": EXPNAME,
    "EXP_CODE": EXP_CODE,
    "EXP_YEAR": EXP_YEAR,
    "DIFX_DIR": DIFX_DIR,
    "PCONV_DIR": PCONV_DIR,
    "mysteps": mysteps,
    "NCPU": NCPU,
    "DOIF": list(DOIF),
    "SCAN_IF": SCAN_IF,
    "IF_OFFSET": IF_OFFSET,
    "AUTO_XPOL_DELAY": AUTO_XPOL_DELAY,
    "TWO_ITER_STEP1": TWO_ITER_STEP1,
    "REFANT": REFANT,
    "REFANT_SELECTION_MODE": globals().get("REFANT_SELECTION_MODE", "highest_snr"),
})
logger.info("Loaded metadata: ANT_COORDS=%d, HOPSNAMES=%d, SAMP_DELAYS=%d", len(ANT_COORDS), len(HOPSNAMES), len(SAMP_DELAYS))
########### +++++++++++++++++++++++++++++++++++++++++++++  ###########

# Keywords used in all steps will be saved here:
if not os.path.exists('STEP_KEYWORDS'):
  os.system('mkdir STEP_KEYWORDS')


# All auxiliary scripts will start with these lines:
Start = 'import pickle as pk\n'
Start += 'from EUVGOS_PY3 import SOURCE_SCANNER_MP as Sscan\n'
Start += 'from EUVGOS_PY3 import SWIN_CONCAT as SwConcat\n'
Start += 'from EUVGOS_PY3 import POL_CALIBRATE_MP as PolCal\n'
Start += 'from EUVGOS_PY3 import POLCONVERTER as PConv\n'
Start += 'from EUVGOS_PY3 import PY_PHASES as PYF\n'
Start += 'from PolConvert import polconvert_standalone as PC\n'

# Loads the keywords for the execution of the corresponding task:
Start += 'IFF = open(\'keywords_%s.dat\',\'rb\'); kww = pk.load(IFF); IFF.close()\n'

# Place to save all intermediate CASA log files (can get VERY large!!)
if not os.path.exists('LOGS'):
  os.system('mkdir LOGS')  
currlogs = [] #glob.glob('*.log')




# STEP 0: SCAN INFO (SOURCES AND SNR).
if 0 in mysteps:

  if len(list(filter(lambda x: 'SOURCE_SCANNER' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

  SCRIPT_NAME = 'STEP0'

  log_header("STEP 0: Source scanner")
  logger.info("Running SOURCE_SCANNER with SNRCut=%s, MIN_ELEV_DEG=%s, MIN_UV_SPAN=%s", MIN_SCAN_SNR, MIN_ELEV_DEG, MIN_UV_SPAN)

  keyw = {
          'EXPNAME':EXPNAME, 'DIFX_DIR':DIFX_DIR, 'SNRCut': MIN_SCAN_SNR,
          'SCAN_IF':SCAN_IF, 'NIF':len(DOIF), 'IF_OFFSET':IF_OFFSET,
          'MIN_ELEV_DEG':MIN_ELEV_DEG,
          'REQUIRE_CROSSHANDS': True,
          'MIN_UNIQUE_BASELINES': MIN_UNIQUE_BASELINES,
          'MIN_UV_SPAN': MIN_UV_SPAN,
          'ANT_COORDS': ANT_COORDS,
          'MIN_TOTAL_PA_COVERAGE_DEG':MIN_TOTAL_PA_COVERAGE_DEG,
          'PA_SCANS_REQUIRED':PA_SCANS_REQUIRED,
          'NIF_PER_BAND': NIF_PER_BAND,
          'NCPU':NCPU}
  keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys, protocol=0); keys.close()

  OFF = open('%s.py'%SCRIPT_NAME,'w')
  print(Start%SCRIPT_NAME,file=OFF)
  print('Sscan.SOURCE_SCANNER(**kww)',file=OFF)
  OFF.close()

  os.system(PYTHON_CALL%SCRIPT_NAME) 

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)

  os.system('mv keywords_STEP0*.dat STEP0*.py STEP_KEYWORDS/.')

  if os.path.exists('SOURCE_SCANNER.FAILED'):
     raise Exception('STEP 0 FAILED!') 

  # Auto-set POLCAL_SCAN for later steps if user did not specify it.
  # Uses CALIB_SCANS from <EXPNAME>_ScanStats.dat written by STEP 0.
  scanstat_file = "%s_ScanStats.dat" % EXPNAME
  polcal_unset = ('POLCAL_SCAN' not in globals()) or (POLCAL_SCAN is None) or (len(POLCAL_SCAN) == 0)

  if os.path.exists(scanstat_file):
    with open(scanstat_file, "rb") as f:
      scanstats = pk.load(f)

    if len(scanstats) == 5:
      BEST_SCAN, CALIB_SCANS, FALLBACK_ANTS, BAD_ANTENNAS, BAD_IF_STEP0 = scanstats
      ALL_ANTENNAS_STEP0 = sorted(BEST_SCAN.keys())
      SCAN_INFO = []
      print("WARNING: old ScanStats format detected; SCAN_INFO unavailable. REFANT-baseline connectivity cannot be checked.")
      logger.warning("Old ScanStats format detected; SCAN_INFO unavailable. REFANT-baseline connectivity cannot be checked.")
    else:
      BEST_SCAN, CALIB_SCANS, FALLBACK_ANTS, BAD_ANTENNAS, BAD_IF_STEP0, ALL_ANTENNAS_STEP0, SCAN_INFO = scanstats

    ALL_ANTENNAS_STEP0 = sorted(set(ALL_ANTENNAS_STEP0))

    if polcal_unset:
      POLCAL_SCAN = list(CALIB_SCANS)
      if len(POLCAL_SCAN) == 0:
        raise RuntimeError("POLCAL_SCAN is empty after STEP0. \n  Select manually using SOURCES_{EXPERIMENT}.txt file or soften SNR requirement.")
      ADDITIVE_PHASE_SCANS = POLCAL_SCAN
      print("POLCAL_SCAN auto-filled from STEP0: %s" % str(POLCAL_SCAN))
      logger.info("POLCAL_SCAN auto-filled from STEP0: %s", POLCAL_SCAN)
    else:
      print("Using user-provided POLCAL_SCAN before connectivity enrichment: %s" % str(POLCAL_SCAN))
      logger.info("Using user-provided POLCAL_SCAN before connectivity enrichment: %s", POLCAL_SCAN)

    # Clean EXCLUDE_ANTENNA: remove stations not in experiment
    old_excl = list(EXCLUDE_ANTENNA)
    EXCLUDE_ANTENNA = sorted(set(EXCLUDE_ANTENNA) & set(ALL_ANTENNAS_STEP0))
    removed_excl = sorted(set(old_excl) - set(EXCLUDE_ANTENNA))
    if removed_excl:
      print("WARNING: removed non-participating antennas from EXCLUDE_ANTENNA: %s" %
            ",".join(removed_excl))
      logger.warning("Removed non-participating antennas from EXCLUDE_ANTENNA: %s", ",".join(removed_excl))

    # Add antennas that STEP0 found globally bad
    EXCLUDE_ANTENNA = sorted(set(EXCLUDE_ANTENNA) | set(BAD_ANTENNAS))

    # Merge STEP0 BAD_IF into existing BAD_IF
    if 'BAD_IF' not in globals() or BAD_IF is None:
      BAD_IF = {}
    for ant, ifs in BAD_IF_STEP0.items():
      if ant not in ALL_ANTENNAS_STEP0:
        continue
      old_ifs = list(BAD_IF.get(ant, []))
      BAD_IF[ant] = sorted(set(old_ifs) | set(ifs))

    # Clean BAD_IF for non-participating antennas
    BAD_IF = {ant: sorted(set(ifs)) for ant, ifs in BAD_IF.items()
              if ant in ALL_ANTENNAS_STEP0}

    # Fill / clean ANT_WEIGHTS
    if 'ANT_WEIGHTS' not in globals() or ANT_WEIGHTS is None:
      ANT_WEIGHTS = {}
    ANT_WEIGHTS = {ant: w for ant, w in ANT_WEIGHTS.items()
                   if ant in ALL_ANTENNAS_STEP0}
    for ant in EXCLUDE_ANTENNA:
      ANT_WEIGHTS[ant] = 0.0

    # Fill USE_PCAL
    if 'USE_PCAL' not in globals() or USE_PCAL is None:
      USE_PCAL = {}
    USE_PCAL_ORIG = dict(USE_PCAL)

    # Remove non-participating stations from USE_PCAL
    for ant in list(USE_PCAL.keys()):
      if ant not in ALL_ANTENNAS_STEP0:
        USE_PCAL.pop(ant, None)

    # Default all stations to True unless user already set False
    for ant in ALL_ANTENNAS_STEP0:
      if ant in EXCLUDE_ANTENNA:
        USE_PCAL[ant] = False
      elif ant not in USE_PCAL:
        USE_PCAL[ant] = True
      elif USE_PCAL[ant] is False:
        print("WARNING: USE_PCAL[%s] was set False by user; keeping False." % ant)
        logger.warning("USE_PCAL[%s] was set False by user; keeping False.", ant)

  else:
    raise Exception("STEP 0 did not create %s. Cannot build pipeline auto settings." % scanstat_file)

  # Auto-fill / clean EXCLUDE_BASELINE from close antenna pairs
  INTRASITE_BASELINE_MAX_M = 500.0

  # Antennas participating in automatic intra-site / very short baselines.
  # These will be avoided as REFANT candidates.
  SHORT_BASELINE_ANTS = set()

  if 'EXCLUDE_BASELINE' not in globals() or EXCLUDE_BASELINE is None:
    EXCLUDE_BASELINE = []

  # keep only valid user baselines
  clean_baselines = set()
  for bl in EXCLUDE_BASELINE:
    if len(bl) != 2:
      continue
    a1, a2 = bl[0], bl[1]
    if a1 == a2:
      continue
    if a1 not in ALL_ANTENNAS_STEP0 or a2 not in ALL_ANTENNAS_STEP0:
      print("WARNING: removing EXCLUDE_BASELINE %s-%s: not in experiment." % (a1, a2))
      logger.warning("Removing EXCLUDE_BASELINE %s-%s: not in experiment.", a1, a2)
      continue
    if a1 in EXCLUDE_ANTENNA or a2 in EXCLUDE_ANTENNA:
      continue
    clean_baselines.add(tuple(sorted([a1, a2])))

  # add automatic intra-site / very short baselines
  ants_for_baselines = [
    a for a in ALL_ANTENNAS_STEP0
    if a not in EXCLUDE_ANTENNA and a in ANT_COORDS
  ]
  for i, a1 in enumerate(ants_for_baselines):
    x1, y1, z1 = ANT_COORDS[a1]
    for a2 in ants_for_baselines[i+1:]:
      x2, y2, z2 = ANT_COORDS[a2]
      dist = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
      if dist <= INTRASITE_BASELINE_MAX_M:
        clean_baselines.add(tuple(sorted([a1, a2])))
        SHORT_BASELINE_ANTS.add(a1)
        SHORT_BASELINE_ANTS.add(a2)
        print("WARNING: auto-excluding short baseline %s-%s: %.1f m" %
              (a1, a2, dist))
        logger.warning("Auto-excluding short baseline %s-%s: %.1f m", a1, a2, dist)

  EXCLUDE_BASELINE = [list(bl) for bl in sorted(clean_baselines)]

  logger.info("EXCLUDE_BASELINE after automatic short-baseline filtering: %s", EXCLUDE_BASELINE)
  logger.info("SHORT_BASELINE_ANTS avoided for REFANT: %s", sorted(SHORT_BASELINE_ANTS))

  def _ant_has_bad_if(ant):
    return ant in BAD_IF and len(BAD_IF.get(ant, [])) > 0

  def _baseline_key(a1, a2):
    return tuple(sorted([a1, a2]))

  def _baseline_is_excluded(a1, a2):
    pair = _baseline_key(a1, a2)
    excluded = set(
        tuple(sorted(bl)) for bl in EXCLUDE_BASELINE
        if len(bl) == 2
    )
    return pair in excluded

  def _ant_snr_score(ant):
    try:
      return float(BEST_SCAN.get(ant, {}).get("ant_snr", -1.0))
    except Exception:
      return -1.0

  def _scan_has_ant_pair(s, a1, a2):
    """
    True if scan has the direct a1-a2 baseline.

    Uses SOURCE_SCANNER's metrics['baselines']['pairs_seen'] when available.
    Falls back to antenna co-presence only for old ScanStats.dat files.
    """
    pair_txt = "%s-%s" % tuple(sorted([a1, a2]))

    try:
      pairs_seen = s.get("metrics", {}).get("baselines", {}).get("pairs_seen", [])
      if pairs_seen:
        return pair_txt in pairs_seen
    except Exception:
      pass

    ants = s.get("ants", set())
    return (a1 in ants) and (a2 in ants)

  def _inspect_pcal_tones_for_scans(polcal_scans, make_plots=True):
    from EUVGOS_PY3.INSPECT_PCAL_TONES import INSPECT_PCAL_TONES

    tone_stats_global_local = {}

    for SI in polcal_scans:
      INPUT_PATH = os.path.join(DIFX_DIR, "%s_%s.input" % (EXPNAME, SI))
      if not os.path.exists(INPUT_PATH):
        print("WARNING: cannot find DiFX .input file for PCAL tone count: %s" % INPUT_PATH)
        logger.warning("Cannot find DiFX .input file for PCAL tone count: %s", INPUT_PATH)
        continue

      scan_dir = os.path.join(DIFX_DIR, "%s_%s.difx" % (EXPNAME, SI))
      if not os.path.isdir(scan_dir):
        print("WARNING: scan dir not found for POLCAL_SCAN %s: %s" % (SI, scan_dir))
        logger.warning("Scan dir not found for POLCAL_SCAN %s: %s", SI, scan_dir)
        continue

      plot_outdir = os.path.join("PCAL_TONE_PLOTS", "SCAN_%s" % SI)

      scan_res = INSPECT_PCAL_TONES(
          input_path=INPUT_PATH,
          scan_difx_dir=scan_dir,
          MAX_NTONES_PER_IF=TONE_POOR_PCAL_PER_IF,
          make_plots=make_plots,
          plot_outdir=plot_outdir,
          band_gap_hz=BAND_GAP_SPLIT_HZ,
      )

      for ant, info in scan_res.items():
        if ant not in ALL_ANTENNAS_STEP0:
          continue
        g = tone_stats_global_local.setdefault(
          ant, {"max_ntones": 0, "per_if_counts": {}, "scans": set()}
        )
        g["max_ntones"] = max(g["max_ntones"], int(info.get("max_ntones", 0)))
        for ifnum, ntones in info.get("per_if_counts", {}).items():
          g["per_if_counts"][ifnum] = max(g["per_if_counts"].get(ifnum, 0), int(ntones))
        g["scans"].add(SI)

    tone_rich_ants_local = set()
    tone_poor_ants_local = set()

    for ant in ALL_ANTENNAS_STEP0:
      info = tone_stats_global_local.get(ant, {})
      max_nt = int(info.get("max_ntones", 0))
      if max_nt >= TONE_POOR_PCAL_PER_IF:
        tone_rich_ants_local.add(ant)
      else:
        tone_poor_ants_local.add(ant)

    tone_stats_dump_local = {}
    for ant, info in tone_stats_global_local.items():
      tone_stats_dump_local[ant] = {
        "max_ntones": int(info.get("max_ntones", 0)),
        "per_if_counts": dict(info.get("per_if_counts", {})),
        "scans": sorted(info.get("scans", set())),
      }

    return tone_stats_global_local, tone_stats_dump_local, tone_rich_ants_local, tone_poor_ants_local

  # First PCAL-tone inspection on the initial POLCAL_SCAN.
  # This is used only to help rank REFANT candidates.
  logger.info("Inspecting PCAL tones on initial POLCAL_SCAN for REFANT scoring: %s", POLCAL_SCAN)

  tone_stats_global, tone_stats_dump, tone_rich_ants, tone_poor_ants = _inspect_pcal_tones_for_scans(
      POLCAL_SCAN,
      make_plots=False,
  )

  def _ref_score(ant):
    """
    REFANT intrinsic score.
    Used only after maximizing the number of connected stations.
    """
    snr = _ant_snr_score(ant)
    tone_bonus = 1e6 if ant in tone_rich_ants else 0.0
    return tone_bonus + snr

  # ------------------------------------------------------------------
  # REFANT + POLCAL_SCAN selection for inspection
  #
  # Philosophy:
  #   - Keep the original POLCAL_SCAN from SOURCE_SCANNER.
  #     These are the best individual scans per station.
  #   - Select REFANT among the cleanest possible candidates:
  #       * not excluded
  #       * no BAD_IF
  #       * not involved in automatic intra-site / short baselines
  #       * preferably tone-rich
  #       * high SNR
  #   - For each REFANT candidate, look for REFANT-station baselines
  #     using only scans that at least passed the SNR test.
  #   - Choose the REFANT that maximizes the number of connected stations.
  #   - Add only the best REFANT-connection scan per station.
  # ------------------------------------------------------------------

  def _snr_pass_scan(s):
    """
    True if scan passed the original SNR test.

    In SCAN_INFO, scan_list is the original SNROut list from getFringeSNR.
    The last element of that list is the original SNR PASS flag: "Y" or "N".
    """
    try:
      return s.get("scan_list", [])[-1] == "Y"
    except Exception:
      return False

  def _inspection_scan_ok(s):
    """
    For REFANT inspection connectivity we require at least SNR PASS.
    We do not require FINAL_PASS here, because FINAL_PASS may include
    extra cuts like elevation, uv span, etc.
    """
    return _snr_pass_scan(s)

  def _connection_scan_score(s, refant, ant):
    """
    Score the scan used to connect ant to REFANT.
    This keeps the scan that best fits the individual station.
    """
    score = 0.0
    # SNR of target antenna matters most for inspecting that station.
    try:
      score += 1000.0 * float(s.get("ant_snr", {}).get(ant, -1.0))
    except Exception:
      pass
    # REFANT SNR is useful, but secondary.
    try:
      score += 700.0 * float(s.get("ant_snr", {}).get(refant, -1.0))
    except Exception:
      pass
    # Prefer globally good scans if available, but do not require them.
    if s.get("final_pass") == "Y":
      score += 1e5
    return score

  def _evaluate_refant_candidate(refant, initial_polcal_scans):
    """
    For one REFANT candidate:
      - keep the initial POLCAL_SCAN untouched
      - find the best SNR-passing direct-baseline scan to each usable station
      - return coverage and scans to add
    """
    usable_ants = [
        ant for ant in ALL_ANTENNAS_STEP0
        if ant not in EXCLUDE_ANTENNA
    ]
    connected = {}
    missing = []
    added = []
    calib_set = set(initial_polcal_scans)

    for ant in usable_ants:

      if ant == refant:
        continue
      if _baseline_is_excluded(refant, ant):
        missing.append({
            "ant": ant,
            "reason": "REFANT baseline is excluded",
        })
        continue

      # First check if the original POLCAL_SCAN already contains
      # a usable direct REFANT-ant baseline.
      already_candidates = [
          s for s in SCAN_INFO
          if s.get("code") in calib_set
          and _inspection_scan_ok(s)
          and _scan_has_ant_pair(s, refant, ant)
      ]
      if len(already_candidates) > 0:
        best_s = max(
            already_candidates,
            key=lambda s: _connection_scan_score(s, refant, ant)
        )
        connected[ant] = {
            "scan": best_s.get("code"),
            "already_selected": True,
            "final_pass": best_s.get("final_pass", "N"),
            "snr_pass": _snr_pass_scan(best_s),
            "ant_snr": best_s.get("ant_snr", {}).get(ant, None),
        }
        continue

      # Otherwise search all scans for the best SNR-passing direct baseline.
      candidates = [
          s for s in SCAN_INFO
          if _inspection_scan_ok(s)
          and _scan_has_ant_pair(s, refant, ant)
      ]
      if len(candidates) == 0:
        missing.append({
            "ant": ant,
            "reason": "no SNR-passing scan with direct REFANT baseline",
        })
        continue

      best_s = max(
          candidates,
          key=lambda s: _connection_scan_score(s, refant, ant)
      )
      calib_set.add(best_s.get("code"))
      connected[ant] = {
          "scan": best_s.get("code"),
          "already_selected": False,
          "final_pass": best_s.get("final_pass", "N"),
          "snr_pass": _snr_pass_scan(best_s),
          "ant_snr": best_s.get("ant_snr", {}).get(ant, None),
      }
      added.append({
          "ant": ant,
          "scan": best_s.get("code"),
          "final_pass": best_s.get("final_pass", "N"),
          "ant_snr": best_s.get("ant_snr", {}).get(ant, None),
      })

    return {
        "refant": refant,
        "polcal_scan": sorted(calib_set),
        "connected": connected,
        "missing": missing,
        "added": added,
        "n_connected": len(connected),
        "n_missing": len(missing),
        "score": _ref_score(refant),
    }

  # Build strict REFANT candidates.
  usable_refants_all = [
      ant for ant in ALL_ANTENNAS_STEP0
      if ant not in EXCLUDE_ANTENNA
      and not _ant_has_bad_if(ant)
      and ant not in SHORT_BASELINE_ANTS
  ]
  # Prefer tone-rich REFANTs. If no tone-rich clean candidate exists,
  # fall back to all clean candidates.
  usable_refants_tone_rich = [
      ant for ant in usable_refants_all
      if ant in tone_rich_ants
  ]
  if len(usable_refants_tone_rich) > 0:
    refant_candidates = sorted(
        usable_refants_tone_rich,
        key=_ref_score,
        reverse=True,
    )
    print("REFANT candidates restricted to clean tone-rich antennas: %s" % str(refant_candidates))
    logger.info("REFANT candidates restricted to clean tone-rich antennas: %s", refant_candidates)
  else:
    refant_candidates = sorted(
        usable_refants_all,
        key=_ref_score,
        reverse=True,
    )
    print("WARNING: no clean tone-rich REFANT candidate found; using clean high-SNR candidates.")
    print("REFANT candidates: %s" % str(refant_candidates))
    logger.warning("No clean tone-rich REFANT candidate found; using clean high-SNR candidates.")
    logger.info("REFANT candidates: %s", refant_candidates)

  # If user supplied REFANT and it is clean, use ONLY that antenna.
  # If it is invalid, discard it and continue with automatic selection.
  USER_REFANT = REFANT
  if USER_REFANT is not None and USER_REFANT != "":
    if (USER_REFANT in ALL_ANTENNAS_STEP0 and USER_REFANT not in EXCLUDE_ANTENNA and not _ant_has_bad_if(USER_REFANT) ):
      refant_candidates = [USER_REFANT]
      print("User REFANT=%s accepted as the only REFANT candidate." % USER_REFANT)
      logger.info("User REFANT=%s accepted as the only REFANT candidate.", USER_REFANT)
    else:
      print("WARNING: user REFANT=%s is invalid, excluded, non-participating, or has BAD_IF. Discarding it and using automatic REFANT selection." % USER_REFANT)
      logger.warning("User REFANT=%s is invalid, excluded, non-participating, or has BAD_IF. Discarding it and using automatic REFANT selection.", USER_REFANT)
      USER_REFANT = None
      REFANT = None

  if len(refant_candidates) == 0:
    raise RuntimeError(
        "No valid REFANT candidates available. "
        "No user REFANT was accepted, and all automatic candidates are excluded, have BAD_IF, or are involved in short/intra-site baselines."
    )

  initial_polcal_scans = list(POLCAL_SCAN)
  REFANT_SELECTION_MODE = globals().get("REFANT_SELECTION_MODE", "highest_snr")

  if len(SCAN_INFO) == 0:
    print("WARNING: SCAN_INFO unavailable. REFANT baseline coverage cannot be optimized; selecting best clean REFANT by score.")
    logger.warning("SCAN_INFO unavailable. REFANT baseline coverage cannot be optimized; selecting best clean REFANT by score.")
    REFANT = refant_candidates[0]
    final_solution = {
      "refant": REFANT,
      "polcal_scan": sorted(set(initial_polcal_scans)),
      "connected": {},
      "missing": [],
      "added": [],
      "n_connected": 0,
      "n_missing": 0,
      "score": _ref_score(REFANT),
    }

  else:
    solutions = []

    if REFANT_SELECTION_MODE == "highest_snr" and (USER_REFANT is None or USER_REFANT == ""):
      best_snr_refant = max(refant_candidates, key=_ant_snr_score)
      refant_candidates = [best_snr_refant]
      print("REFANT selection mode highest_snr: using highest-SNR clean antenna %s as only REFANT candidate." % best_snr_refant)
      logger.info("REFANT selection mode highest_snr: using highest-SNR clean antenna %s as only REFANT candidate.", best_snr_refant)

    for cand_refant in refant_candidates:
      sol = _evaluate_refant_candidate(cand_refant, initial_polcal_scans)
      solutions.append(sol)
      print("REFANT candidate %s: connects %d stations; missing %d; score %.3f" % (cand_refant, sol["n_connected"], sol["n_missing"], sol["score"]))
      logger.info("REFANT candidate %s: connects %d stations; missing %d; score %.3f", cand_refant, sol["n_connected"], sol["n_missing"], sol["score"])

    if REFANT_SELECTION_MODE == "highest_snr":
      final_solution = solutions[0]
    else:
      # Select the REFANT that maximizes direct-baseline coverage.
      # Tie-breaker: intrinsic REFANT score, mostly tone-rich + SNR.
      final_solution = max(solutions, key=lambda sol: (sol["n_connected"], sol["score"]))
    REFANT = final_solution["refant"]

  POLCAL_SCAN = final_solution["polcal_scan"]
  ADDITIVE_PHASE_SCANS = list(POLCAL_SCAN)

  if USER_REFANT is not None and USER_REFANT != "":
    print("REFANT selected from user configuration: %s" % REFANT)
    logger.info("REFANT selected from user configuration: %s", REFANT)
  elif REFANT_SELECTION_MODE == "highest_snr":
    print("REFANT selected by highest clean station SNR: %s" % REFANT)
    logger.info("REFANT selected by highest clean station SNR: %s", REFANT)
  else:
    print("REFANT selected by maximum SNR-passing direct-baseline coverage: %s" % REFANT)
    logger.info("REFANT selected by maximum SNR-passing direct-baseline coverage: %s", REFANT)

  print("REFANT=%s connects %d stations; missing %d." % (REFANT, final_solution["n_connected"], final_solution["n_missing"]) )
  print("POLCAL_SCAN after preserving original scans and adding REFANT-connection scans: %s" % str(POLCAL_SCAN))
  logger.info("REFANT=%s connects %d stations; missing %d.", REFANT, final_solution["n_connected"], final_solution["n_missing"] )
  logger.info("POLCAL_SCAN after preserving original scans and adding REFANT-connection scans: %s", POLCAL_SCAN)

  if final_solution["connected"]:
    print("REFANT direct-baseline coverage:")
    logger.info("REFANT direct-baseline coverage:")

    for ant in sorted(final_solution["connected"]):
      item = final_solution["connected"][ant]
      print("  %s -- %s: scan %s  ant_snr=%s  FINAL_PASS=%s  already_selected=%s" % (
              REFANT,
              ant,
              item["scan"],
              str(item["ant_snr"]),
              item["final_pass"],
              str(item["already_selected"]),
          )
      )
      logger.info("  %s -- %s: scan %s  ant_snr=%s  FINAL_PASS=%s  already_selected=%s",
          REFANT,
          ant,
          item["scan"],
          str(item["ant_snr"]),
          item["final_pass"],
          str(item["already_selected"]),
      )

  if final_solution["added"]:
    print("Added scans for REFANT-station inspection:")
    logger.info("Added scans for REFANT-station inspection:")

    for item in final_solution["added"]:
      print(
          "  %s -- %s in scan %s  ant_snr=%s  FINAL_PASS=%s"
          % (
              REFANT,
              item["ant"],
              item["scan"],
              str(item["ant_snr"]),
              item["final_pass"],
          )
      )
      logger.info(
          "  %s -- %s in scan %s  ant_snr=%s  FINAL_PASS=%s",
          REFANT,
          item["ant"],
          item["scan"],
          str(item["ant_snr"]),
          item["final_pass"],
      )

  if final_solution["missing"]:
    print("WARNING: these antennas still lack an SNR-passing direct baseline to REFANT=%s:" % REFANT)
    logger.warning("These antennas still lack an SNR-passing direct baseline to REFANT=%s:", REFANT)

    for item in final_solution["missing"]:
      print("  %s: %s" % (item["ant"], item["reason"]))
      logger.warning("  %s: %s", item["ant"], item["reason"])

  # Re-inspect PCAL tones on the final POLCAL_SCAN.
  final_polcal_changed = sorted(set(POLCAL_SCAN)) != sorted(set(initial_polcal_scans))
  if final_polcal_changed or DO_PCAL_TONE_PLOTS:
    if final_polcal_changed:
      print("Re-inspecting PCAL tones on final POLCAL_SCAN after ensuring REFANT-connectivity with the full VGOS.")
      logger.info("Re-inspecting PCAL tones on final POLCAL_SCAN after ensuring REFANT-connectivity with the full VGOS array.")
    else:
      print("Inspecting PCAL tones on final POLCAL_SCAN for final statistics/plots.")
      logger.info("Inspecting PCAL tones on final POLCAL_SCAN for final statistics/plots.")

    tone_stats_global, tone_stats_dump, tone_rich_ants, tone_poor_ants = _inspect_pcal_tones_for_scans(
        POLCAL_SCAN,
        make_plots=DO_PCAL_TONE_PLOTS,
    )

  # GFF_REFANTS after connectivity-aware REFANT selection
  def _gff_ref_score(ant):
    snr = _ant_snr_score(ant)
    tone_bonus = 1e6 if ant in tone_rich_ants else 0.0
    bad_penalty = -1e9 if _ant_has_bad_if(ant) else 0.0
    excl_penalty = -1e12 if ant in EXCLUDE_ANTENNA else 0.0
    short_penalty = -1e8 if ant in SHORT_BASELINE_ANTS else 0.0
    return tone_bonus + snr + bad_penalty + excl_penalty + short_penalty

  valid_refants = [
      ant for ant in ALL_ANTENNAS_STEP0
      if ant not in EXCLUDE_ANTENNA
      and not _ant_has_bad_if(ant)
  ]

  valid_refants = sorted(valid_refants, key=_gff_ref_score, reverse=True)

  # Ensure selected REFANT is first.
  if REFANT in valid_refants:
    valid_refants.remove(REFANT)
  valid_refants = [REFANT] + valid_refants

  if 'GFF_REFANTS' not in globals() or GFF_REFANTS is None:
    GFF_REFANTS = []

  kept_gff = []

  for ant in GFF_REFANTS:
    if ant not in ALL_ANTENNAS_STEP0:
      print("WARNING: removing GFF_REFANT %s: not in experiment." % ant)
      logger.warning("Removing GFF_REFANT %s: not in experiment.", ant)
      continue
    if ant in EXCLUDE_ANTENNA:
      print("WARNING: removing GFF_REFANT %s: excluded." % ant)
      logger.warning("Removing GFF_REFANT %s: excluded.", ant)
      continue
    if _ant_has_bad_if(ant):
      print("WARNING: removing GFF_REFANT %s: has BAD_IF entries." % ant)
      logger.warning("Removing GFF_REFANT %s: has BAD_IF entries.", ant)
      continue
    kept_gff.append(ant)

  for ant in valid_refants:
    if ant not in kept_gff:
      kept_gff.append(ant)

  GFF_REFANTS = kept_gff
  print("GFF_REFANTS final list: %s" % str(GFF_REFANTS))

  logger.info("GFF_REFANTS final list: %s", GFF_REFANTS)

  with open(PIPELINE_AUTO_FILE, "wb") as f:
    pk.dump(
        {
          "ALL_ANTENNAS": ALL_ANTENNAS_STEP0,
          "USE_PCAL": USE_PCAL,
          "REFANT": REFANT,
          "GFF_REFANTS": GFF_REFANTS,
          "ANT_WEIGHTS": ANT_WEIGHTS,
          "BAD_IF": BAD_IF,
          "EXCLUDE_ANTENNA": EXCLUDE_ANTENNA,
          "EXCLUDE_BASELINE": EXCLUDE_BASELINE,
          "tone_stats_global": tone_stats_dump,
          "tone_rich_ants": sorted(tone_rich_ants),
          "tone_poor_ants": sorted(tone_poor_ants),
          "POLCAL_SCAN": POLCAL_SCAN,
          "ADDITIVE_PHASE_SCANS": ADDITIVE_PHASE_SCANS,
          "XPOL_DELAYS":XPOL_DELAYS,
        },
        f,
    )

  logger.info("POLCAL_SCAN selected: %s", POLCAL_SCAN)
  logger.info("ALL_ANTENNAS: %s", ALL_ANTENNAS_STEP0)
  logger.info("BAD_ANTENNAS from Step 0: %s", BAD_ANTENNAS)
  logger.info("BAD_IF after Step 0 merge: %s", BAD_IF)
  logger.info("EXCLUDE_ANTENNA final after Step 0: %s", EXCLUDE_ANTENNA)
  logger.info("EXCLUDE_BASELINE final after Step 0: %s", EXCLUDE_BASELINE)
  logger.info("SHORT_BASELINE_ANTS avoided for REFANT: %s", sorted(SHORT_BASELINE_ANTS))
  logger.info("REFANT selected: %s", REFANT)
  logger.info("GFF_REFANTS selected: %s", GFF_REFANTS)
  logger.info("Tone-rich antennas: %s", sorted(tone_rich_ants))
  logger.info("Tone-poor antennas: %s", sorted(tone_poor_ants))
  logger.info("Wrote pipeline auto file: %s", PIPELINE_AUTO_FILE)

  # Free STEP0-only bulky objects before later pipeline steps.
  # SCAN_INFO is only needed for REFANT/POLCAL_SCAN selection above.
  try:
    del SCAN_INFO
  except Exception:
    pass
  try:
    del scanstats
  except Exception:
    pass
  try:
    del tone_stats_global
  except Exception:
    pass

need_pipeline_auto = (len(mysteps) > 0) and (max(mysteps) > 0) ### NOTE: SET TO False IF TESTING MANUAL CONFIGURATION SKIPPING STEP 0!!!
if os.path.exists(PIPELINE_AUTO_FILE):
  with open(PIPELINE_AUTO_FILE, "rb") as f:
    auto = pk.load(f)

  BAD_IF = auto.get("BAD_IF", BAD_IF)
  EXCLUDE_ANTENNA = auto.get("EXCLUDE_ANTENNA", EXCLUDE_ANTENNA)
  EXCLUDE_BASELINE = auto.get("EXCLUDE_BASELINE", EXCLUDE_BASELINE)
  ANT_WEIGHTS = auto.get("ANT_WEIGHTS", ANT_WEIGHTS)
  USE_PCAL = auto.get("USE_PCAL", USE_PCAL)
  REFANT = auto.get("REFANT", REFANT)
  GFF_REFANTS = auto.get("GFF_REFANTS", GFF_REFANTS)
  POLCAL_SCAN = auto.get("POLCAL_SCAN", POLCAL_SCAN)
  ADDITIVE_PHASE_SCANS = auto.get("ADDITIVE_PHASE_SCANS", ADDITIVE_PHASE_SCANS)
  XPOL_DELAYS = auto.get("XPOL_DELAYS", XPOL_DELAYS)

  logger.info("Loaded pipeline auto file: %s", PIPELINE_AUTO_FILE)
  logger.info("Auto settings loaded: POLCAL_SCAN=%s REFANT=%s EXCLUDE_ANTENNA=%s BAD_IF=%s", POLCAL_SCAN, REFANT, EXCLUDE_ANTENNA, BAD_IF)

elif need_pipeline_auto:
  raise Exception(
    "%s not found. Run STEP 0 first so the pipeline can set POLCAL_SCAN, "
    "BAD_IF, USE_PCAL, REFANT, GFF_REFANTS, ANT_WEIGHTS, EXCLUDE_ANTENNA and EXCLUDE_BASELINE."
    % PIPELINE_AUTO_FILE
  )





# STEP 1: FIND CROSS-POL GAINS FROM A CALIBRATOR SCAN.
if 1 in mysteps:
  os.system('rm -rf POLCAL_OUTPUT_SCAN-*_IF-*.dat')
  os.system('rm -rf FRINGE.PEAKS FRINGE.PLOTS POLCONVERT.FRINGE')

  poorly_constrained_ants = []
  final_use_pcal = dict(USE_PCAL)

  if len(POLCAL_SCAN) == 0:
    raise RuntimeError("POLCAL_SCAN is empty. Run STEP 0 or set POLCAL_SCAN manually.")

  # ---------------------------------------
  # 1) SwConcat + replicate PCAL files per IF
  CALDIR = os.path.join(PCONV_DIR, '%s_PC_CALIB' % EXPNAME)
  if not os.path.exists(PCONV_DIR):
    os.system('mkdir %s' % PCONV_DIR)
  os.system('rm -rf %s' % CALDIR)
  os.system('mkdir %s' % CALDIR)

  log_header("STEP 1: SWIN_CONCAT and POL_CALIBRATE")
  logger.info("POLCAL_SCAN: %s", POLCAL_SCAN)
  logger.info("CALDIR: %s", CALDIR)
  logger.info("USE_PCAL initial: %s", USE_PCAL)
  logger.info("BAD_IF used in Step 1: %s", BAD_IF)
  logger.info("EXCLUDE_ANTENNA used in Step 1: %s", EXCLUDE_ANTENNA)
  logger.info("EXCLUDE_BASELINE used in Step 1: %s", EXCLUDE_BASELINE)
  logger.info("Substep 1A: concatenating POLCAL scans")

  swinToCal = ['%s_%s.difx' % (os.path.join(DIFX_DIR, EXPNAME), SI) for SI in POLCAL_SCAN]

  SCRIPT_NAME = 'STEP1_CONCAT'
  keyw = {'SWINs': swinToCal, 'concatName': CALDIR}
  with open('keywords_%s.dat' % SCRIPT_NAME, 'wb') as keys:
    pk.dump(keyw, keys, protocol=0)

  with open('%s.py' % SCRIPT_NAME, 'w') as OFF:
    print(Start % SCRIPT_NAME, file=OFF)
    print('SwConcat.swinConcat(**kww)', file=OFF)

  os.system(PYTHON_CALL % SCRIPT_NAME)

  logger.info("Substep 1A completed: concatenated %d POLCAL scan(s) into %s", len(POLCAL_SCAN), CALDIR)

  for SI in POLCAL_SCAN:
    pcals = glob.glob(os.path.join(CALDIR, '%s_%s.difx/PCAL*' % (EXPNAME, SI)))
    for pcal in pcals:
      for CURRIF in DOIF:
        os.system('cp -p "%s" "%s_IF%i"' % (pcal, pcal, CURRIF))

  logger.info("Replicated PCAL files with per-IF suffixes for %d IFs", len(list(DOIF)))

  NCPU = int(NCPU)
  if NCPU < 1:
    NCPU = max(1, multiprocessing.cpu_count() - 1)

  if len(list(filter(lambda x: 'POL_CALIBRATE' not in x, glob.glob('*.FAILED')))) > 0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')

  print("Finished swin concat.\nWill start the polarization calibration\n")

  # ---------------------------------------
  def RUN_POLCAL_AND_WRITE_GAINS(curr_xpol_delays, tag='', fit_antenna=None, exclude_antenna=None, use_pcal=None, solve_amp=True, use_rate = False, use_delay = False):
    """
    Runs POL_CALIBRATE over DOIF using curr_xpol_delays, then writes
    POLCAL_GAINS_<EXPNAME><tag>.dat and returns that filename.
    """
    if fit_antenna is None:
      fit_antenna = []
    if exclude_antenna is None:
      exclude_antenna = []
    if use_pcal is None:
      use_pcal = {}

    logger.info("Running POL_CALIBRATE tag=%s fit_antenna=%s exclude_antenna=%s solve_amp=%s", tag, fit_antenna, exclude_antenna, solve_amp)
    logger.info("XPOL delays supplied for %d antennas", len(curr_xpol_delays))

    os.system('rm -rf POLCAL_OUTPUT_SCAN-*_IF-*.dat')

    SCRIPT_NAMES = []
    for CURRIF in DOIF:
      BADANTS = [str(ai) for ai in exclude_antenna]
      for ant in use_pcal.keys():
        if (ant not in BADANTS) and (ant in BAD_IF.keys()) and CURRIF in BAD_IF[ant]:
          BADANTS.append(str(ant))

      SCRIPT_NAME = 'STEP1_%i' % CURRIF
      keyw = {
        'EXPNAME': EXPNAME,
        'DIFX_DIR': CALDIR,
        'DOSCAN': POLCAL_SCAN,
        'CHANSOL': CHAVG,
        'USE_PCAL': use_pcal,
        'INTTIME': INTTIME,
        'EXCLUDE_BASELINES': EXCLUDE_BASELINE,
        'FIT_ANTENNA': fit_antenna,
        'EXCLUDE_ANTENNA': BADANTS,
        'DOAMP': solve_amp,
        'DOIF': [CURRIF],
        'PLOTANT': REFANT,
        'APPLY_AMP': APPLY_AMP,
        'SOLVER': SOLVER,
        'APPLY_POLCAL': False,
        'XPOL_DELAYS': curr_xpol_delays,
        'USE_DELAY': use_delay,
        'USE_RATE': use_rate,
        'PCAL_SUFFIX': '_IF%i' % CURRIF,
        'IF_OFFSET': IF_OFFSET,
        'XYPCALMODE': XYPCALMODE,
        'UVTAPER': UVTAPER,
      }

      with open('keywords_%s.dat' % SCRIPT_NAME, 'wb') as keys:
        pk.dump(keyw, keys, protocol=0)

      with open('%s.py' % SCRIPT_NAME, 'w') as OFF:
        print(Start % SCRIPT_NAME, file=OFF)
        print('PolCal.POL_CALIBRATE(**kww)', file=OFF)

      SCRIPT_NAMES.append(SCRIPT_NAME)

    cmds = [PYTHON_CALL % name for name in SCRIPT_NAMES]

    if NCPU > 1:
      pool = multiprocessing.Pool(processes=NCPU)
      rcodes = pool.map(os.system, cmds)
      pool.close()
      pool.join()
    else:
      rcodes = [os.system(cmd) for cmd in cmds]

    bad = [(SCRIPT_NAMES[i], rcodes[i]) for i in range(len(rcodes)) if rcodes[i] != 0]
    if bad:
      for name, rc in bad:
        #print("  - %s (rcode=%s)" % (name, rc))
        log_and_print("  - %s (rcode=%s)" % (name, rc), level="warning")
      #print("WARNING: POL_CALIBRATE job(s) returned non-zero")
      log_and_print("WARNING: POL_CALIBRATE job(s) returned non-zero", level="warning")

    missing = []
    for CURRIF in DOIF:
      fn = 'POLCAL_OUTPUT_SCAN-%s_IF-%i.dat' % (POLCAL_SCAN[0], CURRIF)
      if not os.path.exists(fn):
        missing.append(fn)

    if missing:
      raise RuntimeError("Missing POLCAL output files: %s" % missing)

    CALGAINS = {
      'XYratio': {},
      'XYadd': {},
      'PARAMETERS': {},
      'Frequency': {},
      'XYratioOriginal': {},
    }

    for CURRIF in DOIF:
      with open('POLCAL_OUTPUT_SCAN-%s_IF-%i.dat' % (POLCAL_SCAN[0], CURRIF), 'rb') as IFF:
        IFGAIN = pk.load(IFF)

      for key0 in CALGAINS.keys():
        for key1 in IFGAIN[key0].keys():
          if key1 not in CALGAINS[key0].keys():
            CALGAINS[key0][key1] = {}
          if type(IFGAIN[key0][key1]) is dict:
            for key2 in IFGAIN[key0][key1].keys():
              CALGAINS[key0][key1][key2] = IFGAIN[key0][key1][key2]
          else:
            CALGAINS[key0][key1] = IFGAIN[key0][key1]

    outname = 'POLCAL_GAINS_%s%s.dat' % (EXPNAME, tag)
    with open(outname, 'wb') as OFF:
      pk.dump(CALGAINS, OFF, protocol=0)

    logger.info("Wrote POLCAL gains: %s", outname)

    return outname

  # ---------------------------------------
  def RUN_XPOL_SCANNER(gains_file, explore_candidates, out_subdir):
    from EUVGOS_PY3 import CROSS_POL_DELAY_SCANNER as XYDel

    logger.info("Running XPOL delay scanner on %s; explore_candidates=%s; out=%s", gains_file, explore_candidates, out_subdir)

    outdir = os.path.join(os.getcwd(), out_subdir)
    return XYDel.CROSS_POL_DELAY_SCANNER(
      gains_file,
      out_dir=outdir,
      EXPLORE_CANDIDATES=explore_candidates,
      BAND_GAP_HZ=BAND_GAP_SPLIT_HZ,
    )

  # ---------------------------------------
  def RUN_STEP1B_POLCONVERT(gains_file, script_name, suffix='', use_pcal_for_run=None):
    if use_pcal_for_run is None:
      use_pcal_for_run = USE_PCAL

    keyw = {
      'EXPNAME': EXPNAME,
      'DIFX_DIR': CALDIR,
      'XYGAINS': gains_file,
      'SUFFIX': suffix,
      'IF_OFFSET': IF_OFFSET,
      'XPOL_DELAYS': XPOL_DELAYS,
      'USE_PCAL': use_pcal_for_run,
      'DOPLOT': True,
      'SCAN_LIST': POLCAL_SCAN,
      'REFANT': REFANT,
      'XYPCALMODE': XYPCALMODE,
    }

    with open('keywords_%s.dat' % script_name, 'wb') as keys:
      pk.dump(keyw, keys)

    with open('%s.py' % script_name, 'w') as OFF:
      print(Start % script_name, file=OFF)
      print('PConv.POLCONVERTER(**kww)', file=OFF)

    os.system(PYTHON_CALL % script_name)

  # ---------------------------------------
  def PLOT_CROSS_GAINS(CALGAINS, tag=''):
    fig = pl.figure()
    MaxG = 0.0; NuMin = 1.0e20; NuMax = 0.0

    color = ['r', 'g', 'b', 'k', 'm', 'y', 'c']
    symbol = ['o', '+', 'x']

    sub1 = fig.add_subplot(211)
    sub2 = fig.add_subplot(212, sharex=sub1)

    for CURRIF in DOIF:
      NUS = CALGAINS['Frequency'][CURRIF] / 1.e9
      NuMin = np.min([NuMin, np.min(NUS)])
      NuMax = np.max([NuMax, np.max(NUS)])

      for antii, anti in enumerate(sorted(CALGAINS['XYadd'].keys())):
        toplotPh = np.array(CALGAINS['XYadd'][anti][CURRIF])
        toplotAp = np.array(CALGAINS['XYratioOriginal'][anti][CURRIF])
        MaxG = np.max([MaxG, np.max(toplotAp)])

        style = symbol[int(((antii) // len(color)) % len(symbol))] + color[int((antii) % len(color))]

        if CURRIF == DOIF[0]:
          sub1.plot(np.array(NUS), toplotPh, style, label='ANT. ' + str(anti), ms=2)
          sub2.plot(np.array(NUS), toplotAp, style, label='ANT. ' + str(anti), ms=2)
        else:
          sub1.plot(np.array(NUS), toplotPh, style, ms=2)
          sub2.plot(np.array(NUS), toplotAp, style, ms=2)

    sub1.set_ylim((-180., 180.))
    pl.setp(sub1.get_xticklabels(), 'visible', False)

    Dnu = NuMax - NuMin
    sub1.set_xlim((NuMin - Dnu * 0.1, NuMax + Dnu * 0.45))
    sub2.set_xlim((NuMin - Dnu * 0.1, NuMax + Dnu * 0.45))
    sub2.set_ylim((0., 2.5))

    sub1.legend(numpoints=1)
    sub1.set_ylabel('Cross-Phase (deg.)')
    sub2.set_ylabel('Cross-Amp (Norm.)')
    sub2.set_xlabel('Frequency (GHz)')

    fig.suptitle('XPOL GAINS %s' % EXPNAME)
    pl.savefig('Cross-Gains_%s_JOINED%s.png' % (EXPNAME, tag))
    pl.close(fig)

  # ---------------------------------------
  def PLOT_FRINGE_PEAKS(CALGAINS, tag='', scaled=False, save_to_fringe_dir=False):
    ANTS = [ai for ai in sorted(CALGAINS['XYadd'].keys()) if ai != REFANT]
    NROW = 3
    Ncol = max(1, int(np.ceil(float(len(ANTS)) / NROW)))

    title = 'FRINGE PEAKS'
    suffix = ''
    if scaled:
      title = 'FRINGE PEAKS (SCALED)'
      suffix = '_scaled'

    fig = pl.figure(figsize=(8, 3 * Ncol))
    fig.subplots_adjust(wspace=0.01, hspace=0.01, right=0.98, top=0.90)
    fig.suptitle('%s. SCAN %s' % (title, POLCAL_SCAN[0]))

    for i in range(len(ANTS)):
      plot_refant = REFANT
      sub = fig.add_subplot(Ncol, NROW, i + 1)

      ALLIF = sorted(glob.glob('FRINGE.PEAKS/FRINGE.PEAKS_IF*_%s*.dat' % ANTS[i]))
      toplot = []
      MAX = 0.0
      MAXIF = 1.0
      ifRead = []
      for iffile in ALLIF:
        try:
          with open(iffile) as IFF:
            lines = IFF.readlines()
          plot_refant = lines[0].split()[-1]
          currIF = float(lines[1].split()[2].replace('#', '').replace('.', ''))
          if currIF in ifRead:
            continue
          ifRead.append(currIF)

          PEAK = float(lines[1].split()[-2])
          RR = float(lines[2].split()[1]) * PEAK
          LL = float(lines[3].split()[1]) * PEAK
          RL = float(lines[4].split()[1]) * PEAK
          LR = float(lines[5].split()[1]) * PEAK

          if scaled:
            SCALE = np.max([RR, LL])
            if SCALE <= 0.0:
              SCALE = 1.0
            toplot.append([currIF, RR / SCALE, LL / SCALE, RL / SCALE, LR / SCALE])
          else:
            MAX = np.max([MAX, RR, LL, RL, LR])
            toplot.append([currIF, RR, LL, RL, LR])

        except Exception:
          if not scaled:
            toplot.append([0., 0., 0., 0., 0.])

      if len(toplot) > 0:
        toplot = np.array(sorted(toplot, key=lambda x: x[0]))
        MAXIF = np.max(toplot[:, 0])

        if scaled:
          sub.plot(toplot[:, 0], toplot[:, 1], 'or')
          sub.plot(toplot[:, 0], toplot[:, 2], 'sb')
          sub.plot(toplot[:, 0], toplot[:, 3], 'xg')
          sub.plot(toplot[:, 0], toplot[:, 4], '+k')
        else:
          if MAX == 0.0:
            MAX = 1.0
          sub.plot(toplot[:, 0], toplot[:, 1] / MAX, 'or')
          sub.plot(toplot[:, 0], toplot[:, 2] / MAX, 'sb')
          sub.plot(toplot[:, 0], toplot[:, 3] / MAX, 'xg')
          sub.plot(toplot[:, 0], toplot[:, 4] / MAX, '+k')

      sub.set_ylim((0., 1.2))
      sub.set_xlim((0., MAXIF * 1.2))
      pl.setp(sub.get_yticklabels(), visible=False)

      if i < Ncol * (NROW - 1):
        pl.setp(sub.get_xticklabels(), visible=False)
      else:
        sub.set_xlabel('IF Number')

      if i == 0:
        sub.plot([], [], 'or', label='RR')
        sub.plot([], [], 'sb', label='LL')
        sub.plot([], [], 'xg', label='RL')
        sub.plot([], [], '+k', label='LR')
        pl.legend(numpoints=1, loc=1)

      pl.text(MAXIF * 0.4, 1.1, '%s-%s' % (ANTS[i], plot_refant))

    outname = 'ALL_IF_PEAKS_SCAN_%s%s%s.png' % (POLCAL_SCAN[0], suffix, tag)
    pl.savefig(outname)

    if save_to_fringe_dir and os.path.isdir('FRINGE.PLOTS'):
      pl.savefig(os.path.join('FRINGE.PLOTS', outname))

    pl.close(fig)

  # ---------------------------------------
  def CLEAN_STEP1_OUTPUTS():
    newlogs = glob.glob('*.log')
    for log in newlogs:
      if log not in currlogs:
        os.system('mv %s LOGS/.' % log)

    os.system('mv keywords_STEP1_*.dat keywords_STEP1B*.dat STEP1B*.py STEP1_*.py STEP_KEYWORDS/.')

    if not os.path.exists('CROSS-POL_GAINS.PLOTS'):
      os.system('mkdir CROSS-POL_GAINS.PLOTS')
    os.system('rm -f CROSS-POL_GAINS.PLOTS/*.png')
    os.system('mv Cross-Gains_*_CALIB_IF*.png CROSS-POL_GAINS.PLOTS/.')

    if not os.path.exists('CROSS-POL_GAINS.DATA'):
      os.system('mkdir CROSS-POL_GAINS.DATA')
    os.system('rm -f CROSS-POL_GAINS.DATA/*.dat')
    os.system('mv PolConvert.XYGains_IF*.dat CROSS-POL_GAINS.DATA/.')
    os.system('mv POLCAL_OUTPUT_SCAN*.dat CROSS-POL_GAINS.DATA/.')

  # ---------------------------------------
  def CLEAN_FRINGE_BEFORE_POLCONVERT():
    os.system('rm -rf FRINGE.PLOTS/ALL_IFs*.png')
    os.system('rm -rf FRINGE.PLOTS/RL_LR*.png')
    os.system('rm -rf FRINGE.PEAKS')
    if os.path.isdir('FRINGE.PLOTS'):
      os.system('rm -rf FRINGE.PLOTS_NOCALIB')
      os.system('mv FRINGE.PLOTS FRINGE.PLOTS_NOCALIB')

  # ---------------------------------------
  def RUN_STEP1B_AND_MAKE_PLOTS(gains_file, script_name, tag='', use_pcal_for_run=None, save_to_fringe_dir=False):
    CLEAN_FRINGE_BEFORE_POLCONVERT()
    RUN_STEP1B_POLCONVERT(gains_file, script_name, suffix = tag, use_pcal_for_run=use_pcal_for_run)

    with open(gains_file, "rb") as f:
      CALGAINS = pk.load(f)

    PLOT_CROSS_GAINS(CALGAINS, tag=tag)
    PLOT_FRINGE_PEAKS(CALGAINS, tag=tag, scaled=False, save_to_fringe_dir=save_to_fringe_dir)
    PLOT_FRINGE_PEAKS(CALGAINS, tag=tag, scaled=True, save_to_fringe_dir=save_to_fringe_dir)

    return CALGAINS

  # ---------------------------------------
  # 2) Basic run
  TAG0 = '_RUN1' if TWO_ITER_STEP1 else ''
  TAG = TAG0 + '_NO_XPOL_DELAY' if AUTO_XPOL_DELAY else TAG0

  UseRate  = False if AUTO_XPOL_DELAY else USE_RATE
  UseDelay = False if AUTO_XPOL_DELAY else USE_DELAY

  logger.info("Substep 1B: solving POLCAL gains %s" % TAG0)
  
  GAINS_RUN1 = RUN_POLCAL_AND_WRITE_GAINS(
    XPOL_DELAYS,
    tag=TAG,
    fit_antenna=[],
    exclude_antenna=EXCLUDE_ANTENNA,
    use_pcal=USE_PCAL,
    solve_amp=SOLVE_AMP,
    use_rate = UseRate, use_delay = UseDelay,
  )

  if AUTO_XPOL_DELAY:
    logger.info("Substep 1C: estimating XPOL delays %s" % TAG0)
    xpol_res = RUN_XPOL_SCANNER(
      GAINS_RUN1,
      explore_candidates=False,
      out_subdir="XPOL_DELAY_INSPECT_%s%s" % (EXPNAME, TAG0),
    )

    for ant in USE_PCAL.keys():
      if ant in EXCLUDE_ANTENNA:
        continue
      ares = xpol_res.get(ant, {})
      tau_s = ares.get("tau_hat_s", None)
      f0_hz = ares.get("f_zero_hz", None)
      if (tau_s is None) or (f0_hz is None):
        continue
      XPOL_DELAYS[ant] = [float(tau_s), float(f0_hz)]

    logger.info("AUTO_XPOL_DELAY initial estimate produced delays for %d antennas: %s", len(XPOL_DELAYS), sorted(XPOL_DELAYS.keys()))

    logger.info("Substep 1B repeat: solving POLCAL gains with updated XPOL delays %s" % TAG0)

    GAINS_RUN1 = RUN_POLCAL_AND_WRITE_GAINS(
      XPOL_DELAYS,
      tag=TAG0 + '_XPOL_DELAY',
      fit_antenna=[],
      exclude_antenna=EXCLUDE_ANTENNA,
      use_pcal=USE_PCAL,
      solve_amp=SOLVE_AMP,
      use_rate = USE_RATE, use_delay = USE_DELAY,
    )
    kk = RUN_XPOL_SCANNER(
      GAINS_RUN1,
      explore_candidates=False,
      out_subdir="XPOL_DELAY_INSPECT_%s%s"%(EXPNAME,TAG0+'_XPOL_DELAY_corr'),
    )

  with open('XPOL_DELAYS_%s.dat' % EXPNAME, 'wb') as OFF:
    pk.dump(XPOL_DELAYS, OFF, protocol=0)

  save2fringe_dir = False
  if not TWO_ITER_STEP1:
    GAINS_FINAL = 'POLCAL_GAINS_%s.dat' % EXPNAME
    os.system('cp -p "%s" "%s"' % (GAINS_RUN1, GAINS_FINAL))
    save2fringe_dir = True
    logger.info("Final Step 1 gains selected without two-iteration recovery: %s", GAINS_FINAL)

  logger.info("Substep 1D: PolConvert POLCAL scans + diagnostic / fringe plots")

  # RUN1 plots are needed because FRINGE.PEAKS is used below.
  CALGAINS_RUN1 = RUN_STEP1B_AND_MAKE_PLOTS(
    GAINS_RUN1,
    script_name='STEP1B_RUN1',
    tag=TAG0,
    use_pcal_for_run=USE_PCAL,
    save_to_fringe_dir=save2fringe_dir,
  )

  CLEAN_STEP1_OUTPUTS()

  # ---------------------------------------
  # 3) Classify problematic antennas from RUN1 fringe peaks
  if TWO_ITER_STEP1:
    logger.info("Substep 1E: classifying poorly constrained antennas from RUN1 fringe peaks")
    logger.info("Fringe classification thresholds: CROSS_HAND_MAX_FRAC=%s, MAX_BAD_IF_FRAC=%s",
                globals().get("CROSS_HAND_MAX_FRAC", 0.50),
                globals().get("MAX_BAD_IF_FRAC", 0.50))
    if 'CROSS_HAND_MAX_FRAC' not in globals():
      CROSS_HAND_MAX_FRAC = 0.50
    if 'MAX_BAD_IF_FRAC' not in globals():
      MAX_BAD_IF_FRAC = 0.50

    fringe_bad_ifs = {}
    fringe_good_ifs = {}

    ANTS = [ai for ai in sorted(CALGAINS_RUN1['XYadd'].keys())
            if ai != REFANT and ai not in EXCLUDE_ANTENNA]

    print("[STEP 1] Antennas used for fringe classification, excluding REFANT=%s: %s" % (REFANT, ",".join(ANTS)))

    for ant in ANTS:
      if ant in EXCLUDE_ANTENNA:
        continue
      bad_here = []
      good_here = []

      for CURRIF in DOIF:
        if ant in BAD_IF and CURRIF in BAD_IF[ant]:
          continue

        files = sorted(glob.glob('FRINGE.PEAKS/FRINGE.PEAKS_IF*_%s*.dat' % ant))
        if len(files) == 0:
          bad_here.append(CURRIF)
          continue

        ok_if = False
        for iffile in files:
          try:
            with open(iffile) as IFF:
              lines = IFF.readlines()

            currIF = int(float(lines[1].split()[2].replace('#', '').replace('.', '')))
            if currIF != CURRIF:
              continue

            PEAK = float(lines[1].split()[-2])
            RR = float(lines[2].split()[1]) * PEAK
            LL = float(lines[3].split()[1]) * PEAK
            RL = float(lines[4].split()[1]) * PEAK
            LR = float(lines[5].split()[1]) * PEAK

            scale = max(RR, LL)
            if scale <= 0.0:
              continue
            rl_frac = RL / scale
            lr_frac = LR / scale
            if max(rl_frac, lr_frac) < CROSS_HAND_MAX_FRAC:
              ok_if = True
              break

          except Exception:
            pass

        if ok_if:
          good_here.append(CURRIF)
        else:
          bad_here.append(CURRIF)

      fringe_bad_ifs[ant] = bad_here
      fringe_good_ifs[ant] = good_here

    poorly_constrained_ants = []
    for ant in sorted(fringe_bad_ifs.keys()):
      nbad = len(fringe_bad_ifs[ant])
      ngood = len(fringe_good_ifs[ant])
      ntot = nbad + ngood
      if ntot <= 0:
        continue

      bad_frac = float(nbad) / float(ntot)

      print("[STEP 1] Fringe quality %-4s: good IFs=%s bad IFs=%s bad_frac=%.2f" % (ant, str(fringe_good_ifs[ant]), str(fringe_bad_ifs[ant]), bad_frac))

      if bad_frac > MAX_BAD_IF_FRAC:
        poorly_constrained_ants.append(ant)

    exclude_antenna_set1 = sorted(set(EXCLUDE_ANTENNA) | set(poorly_constrained_ants))

    use_pcal_set1 = dict(USE_PCAL)
    for ant in poorly_constrained_ants:
      use_pcal_set1[ant] = False

    fit_ants = sorted([ant for ant in poorly_constrained_ants
                       if ant not in EXCLUDE_ANTENNA])

    use_pcal_set2 = dict(USE_PCAL)
    for ant in fit_ants:
      use_pcal_set2[ant] = bool(REENABLE_BAD_PCAL_USE)

    final_use_pcal = dict(use_pcal_set2)
    USE_PCAL = final_use_pcal

    print("[STEP 1] Poorly constrained antennas from fringe peaks: %s" %
          (",".join(poorly_constrained_ants) or "NONE"))
    print("[STEP 1] Set 1 excludes: %s" %
          (",".join(exclude_antenna_set1) or "NONE"))
    print("[STEP 1] Set 2 fits/recover: %s" %
          (",".join(fit_ants) or "NONE"))
    log_and_print("[STEP 1] Poorly constrained antennas from fringe peaks: %s" %
                  (",".join(poorly_constrained_ants) or "NONE"))
    log_and_print("[STEP 1] Set 1 excludes: %s" %
                  (",".join(exclude_antenna_set1) or "NONE"))
    log_and_print("[STEP 1] Set 2 fits/recover: %s" %
                  (",".join(fit_ants) or "NONE"))

    active_ants = list(ANTS)
    if len(fit_ants) == 0:
      print("[STEP 1] No poorly constrained antennas found; using RUN1 gains as final.")
      log_and_print("[STEP 1] No poorly constrained antennas found; using RUN1 gains as final.")
      GAINS_FINAL = 'POLCAL_GAINS_%s.dat' % EXPNAME
      os.system('cp -p "%s" "%s"' % (GAINS_RUN1, GAINS_FINAL))
      TWO_ITER_STEP1 = False
    elif len(fit_ants) >= len(active_ants):
      print("[STEP 1] All antennas seem to be poorly constrained; using RUN1 gains as final.")
      log_and_print("[STEP 1] All active antennas seem to be poorly constrained; using RUN1 gains as final.")
      GAINS_FINAL = 'POLCAL_GAINS_%s.dat' % EXPNAME
      os.system('cp -p "%s" "%s"' % (GAINS_RUN1, GAINS_FINAL))
      TWO_ITER_STEP1 = False
      
  # ---------------------------------------
  # 3.5) SwConcat again!
  if TWO_ITER_STEP1:
    os.system('rm -rf %s' % CALDIR)
    os.system('mkdir %s' % CALDIR)
    os.system('rm -rf POLCAL_OUTPUT_SCAN-*_IF-*.dat')
    os.system('rm -rf FRINGE.PEAKS FRINGE.PLOTS POLCONVERT.FRINGE')

    swinToCal = ['%s_%s.difx' % (os.path.join(DIFX_DIR, EXPNAME), SI) for SI in POLCAL_SCAN]

    SCRIPT_NAME = 'STEP1_CONCAT'
    keyw = {'SWINs': swinToCal, 'concatName': CALDIR}
    with open('keywords_%s.dat' % SCRIPT_NAME, 'wb') as keys:
      pk.dump(keyw, keys, protocol=0)

    with open('%s.py' % SCRIPT_NAME, 'w') as OFF:
      print(Start % SCRIPT_NAME, file=OFF)
      print('SwConcat.swinConcat(**kww)', file=OFF)

    os.system(PYTHON_CALL % SCRIPT_NAME)

    for SI in POLCAL_SCAN:
      pcals = glob.glob(os.path.join(CALDIR, '%s_%s.difx/PCAL*' % (EXPNAME, SI)))
      for pcal in pcals:
        for CURRIF in DOIF:
          os.system('cp -p "%s" "%s_IF%i"' % (pcal, pcal, CURRIF))

    print("Finished swin concat again!\nWill start the polarization calibration in two sets!\n")

  # ---------------------------------------
  # 4) Two-iteration recovery
  if TWO_ITER_STEP1:
    logger.info("Substep 1F: two-iteration recovery enabled")
    logger.info("Substep 1F.1: SET1 solve using well-constrained antennas only")
    if AUTO_XPOL_DELAY:
        XPOL_DELAYS = {}

    TAG = '_SET1_NO_XPOL_DELAY' if AUTO_XPOL_DELAY else '_SET1'

    GAINS_SET1 = RUN_POLCAL_AND_WRITE_GAINS(
      XPOL_DELAYS,
      tag=TAG,
      fit_antenna=[],
      exclude_antenna=exclude_antenna_set1,
      use_pcal=use_pcal_set1,
      solve_amp=SOLVE_AMP,
      use_rate = UseRate, use_delay = UseDelay,
    )

    if AUTO_XPOL_DELAY:
      logger.info("Substep 1F.2: AUTO_XPOL_DELAY: estimating XPOL delays from SET1 well-constrained antennas.")
      print("\n[STEP 1] AUTO_XPOL_DELAY: single XPOL scan using first set of well-behaved antennas.")
      xpol_res = RUN_XPOL_SCANNER(
        GAINS_SET1,
        explore_candidates=False,
        out_subdir="XPOL_DELAY_INSPECT_%s_set1" % EXPNAME,
      )
      for ant in USE_PCAL.keys():
        if ant in exclude_antenna_set1:
          continue
        ares = xpol_res.get(ant, {})
        tau_s = ares.get("tau_hat_s", None)
        f0_hz = ares.get("f_zero_hz", None)
        if (tau_s is None) or (f0_hz is None):
          continue
        XPOL_DELAYS[ant] = [float(tau_s), float(f0_hz)]

      print("[STEP 1] AUTO_XPOL_DELAY: re-running POLCAL with XPOL_DELAYS for well-behaved antennas.")
      logger.info("Substep 1F.1 repeat: solving POLCAL gains with updated XPOL delays from SET1 well-constrained antennas.")
      GAINS_SET1 = RUN_POLCAL_AND_WRITE_GAINS(
        XPOL_DELAYS,
        tag='_set1_XPOL_DELAY',
        fit_antenna=[],
        exclude_antenna=exclude_antenna_set1,
        use_pcal=use_pcal_set1,
        solve_amp=SOLVE_AMP,
        use_rate = USE_RATE, use_delay = USE_DELAY,
      )

      kk = RUN_XPOL_SCANNER(
        GAINS_SET1,
        explore_candidates=False,
        out_subdir="XPOL_DELAY_INSPECT_%s_set1_XPOL_DELAY_corr" % EXPNAME,
      )

    logger.info("Substep 1F.3: SET2 solve/recovery for poorly constrained antennas: %s", fit_ants)

    TAG = '_SET2_NO_XPOL_DELAY' if AUTO_XPOL_DELAY else '_SET2'

    GAINS_SET2 = RUN_POLCAL_AND_WRITE_GAINS(
      XPOL_DELAYS,
      tag=TAG,
      fit_antenna=fit_ants,
      exclude_antenna=EXCLUDE_ANTENNA,
      use_pcal=use_pcal_set2,
      solve_amp=SET2DOAMP,
   #   use_rate = UseRate, use_delay = UseDelay,
    )

    if AUTO_XPOL_DELAY:
      print("\n[STEP 1] AUTO_XPOL_DELAY: XPOL scan with problematic antennas.")
      logger.info("Substep 1F.4: AUTO_XPOL_DELAY: estimating XPOL delays from SET2 problematic antennas.")
      xpol_res = RUN_XPOL_SCANNER(
        GAINS_SET2,
        explore_candidates=False,
        out_subdir="XPOL_DELAY_INSPECT_%s_set2" % EXPNAME,
      )
      for ant in fit_ants:
        ares = xpol_res.get(ant, {})
        tau_s = ares.get("tau_hat_s", None)
        f0_hz = ares.get("f_zero_hz", None)
        if (tau_s is None) or (f0_hz is None):
          continue
        XPOL_DELAYS[ant] = [float(tau_s), float(f0_hz)]

      print("[STEP 1] AUTO_XPOL_DELAY: re-running POLCAL with XPOL_DELAYS for problematic antennas.")
      logger.info("Substep 1F.3 repeat: solving POLCAL gains with updated XPOL delays from SET2 problematic antennas.")
      GAINS_SET2 = RUN_POLCAL_AND_WRITE_GAINS(
        XPOL_DELAYS,
        tag="_set2_XPOL_DELAY",
        fit_antenna=fit_ants,
        exclude_antenna=EXCLUDE_ANTENNA,
        use_pcal=use_pcal_set2,
        solve_amp=SET2DOAMP,
     #   use_rate = UseRate, use_delay = UseDelay,
      )

      cand_res = RUN_XPOL_SCANNER(
        GAINS_SET2,
        explore_candidates=True,
        out_subdir="XPOL_DELAY_INSPECT_%s_set2_XPOL_DELAY" % EXPNAME,
      )

      rerun_needed = False
      for ant, res in cand_res.items():
        if ant in fit_ants:
          cands = res.get("candidates", [])
          if len(cands) == 1:
            print("[STEP 1] AUTO_XPOL_DELAY: good a priori XPOL_DELAY for problematic antenna %s" % ant)
          else:
            print("[STEP 1] AUTO_XPOL_DELAY: could not get good a priori XPOL_DELAY for problematic antenna %s. Will not provide one..." % ant)
            XPOL_DELAYS.pop(ant, None)
            rerun_needed = True

      if rerun_needed:
        print("[STEP 1] AUTO_XPOL_DELAY: re-running POLCAL without XPOL_DELAYS for problematic antennas.")
        GAINS_SET2 = RUN_POLCAL_AND_WRITE_GAINS(
          XPOL_DELAYS,
          tag="_set2_XPOL_DELAY_v2",
          fit_antenna=fit_ants,
          exclude_antenna=EXCLUDE_ANTENNA,
          use_pcal=use_pcal_set2,
          solve_amp=SET2DOAMP,
       #   use_rate = UseRate, use_delay = UseDelay,
        )

    logger.info("Poorly constrained antennas: %s", poorly_constrained_ants)
    logger.info("SET1 excludes: %s", exclude_antenna_set1)
    logger.info("SET2 fits/recover: %s", fit_ants)
    logger.info("Final USE_PCAL after Step 1: %s", final_use_pcal)

    log_and_print("[STEP 1] Merging SET1 and SET2 POLCAL gains.")
    GAINS_FINAL = 'POLCAL_GAINS_%s.dat' % EXPNAME

    with open(GAINS_SET1, "rb") as f:
      gains_set1 = pk.load(f)

    with open(GAINS_SET2, "rb") as f:
      gains_set2 = pk.load(f)

    for key0 in gains_set1.keys():
      if not isinstance(gains_set1[key0], dict):
        continue
      for ant in fit_ants:
        if ant in gains_set2[key0]:
          gains_set1[key0][ant] = gains_set2[key0][ant]

    with open(GAINS_FINAL, "wb") as f:
      pk.dump(gains_set1, f, protocol=0)

    print("[STEP 1] Final merged gains written to %s" % GAINS_FINAL)
    logger.info("Final merged Step 1 gains written to %s", GAINS_FINAL)

  # ---------------------------------------
  # 5) Final XPOL cleanup
  clean_xpol = {}
  for ant, vals in XPOL_DELAYS.items():
    try:
      tau_s = float(vals[0])
      f0_hz = float(vals[1])
    except Exception:
      continue

    if not np.isfinite(tau_s) or not np.isfinite(f0_hz):
      continue
    if abs(tau_s) < 1e-15:
      continue

    clean_xpol[ant] = [tau_s, f0_hz]

  XPOL_DELAYS = clean_xpol

  with open('XPOL_DELAYS_%s.dat' % EXPNAME, 'wb') as OFF:
    pk.dump(XPOL_DELAYS, OFF)

  with open('XPOL_DELAYS_%s.txt' % EXPNAME, 'w') as O:
    O.write('# ant  tau_s  tau_ns  f_zero_hz\n')
    for a in sorted(XPOL_DELAYS.keys()):
      tau_s = float(XPOL_DELAYS[a][0])
      f0 = float(XPOL_DELAYS[a][1])
      O.write('%s  %.12e  %.6f  %.12e\n' % (a, tau_s, tau_s * 1e9, f0))

  logger.info("Final XPOL_DELAYS for %d antennas: %s", len(XPOL_DELAYS), XPOL_DELAYS)
  logger.info("Wrote XPOL delay files: XPOL_DELAYS_%s.dat/.txt", EXPNAME)

  # ---------------------------------------
  # 6) Final POLCONVERTER run / final plots
  if TWO_ITER_STEP1:
    CALGAINS_FINAL = RUN_STEP1B_AND_MAKE_PLOTS(
      GAINS_FINAL,
      script_name='STEP1B_FINAL',
      tag='',
      use_pcal_for_run=USE_PCAL,
      save_to_fringe_dir=True,
    )
    CLEAN_STEP1_OUTPUTS()
  elif TWO_ITER_STEP1 is False and not save2fringe_dir:
    # We started in two-iteration mode, but RUN1 was already final:
    #   either no antennas were problematic or all were problematic.
    CALGAINS_FINAL = RUN_STEP1B_AND_MAKE_PLOTS(
      GAINS_FINAL,
      script_name='STEP1B_FINAL',
      tag='',
      use_pcal_for_run=USE_PCAL,
      save_to_fringe_dir=True,
    )
    CLEAN_STEP1_OUTPUTS()

  if os.path.exists('POL_CALIBRATE.FAILED'):
    raise Exception('STEP 1 FAILED!')

  # ---------------------------------------
  # 7) Update pipeline auto file
  if os.path.exists(PIPELINE_AUTO_FILE):
    with open(PIPELINE_AUTO_FILE, "rb") as f:
      auto = pk.load(f)
  else:
    auto = {}

  auto["USE_PCAL"] = final_use_pcal
  auto["POORLY_CONSTRAINED_ANTS"] = poorly_constrained_ants
  auto["XPOL_DELAYS"] = XPOL_DELAYS

  with open(PIPELINE_AUTO_FILE, "wb") as f:
    pk.dump(auto, f, protocol=0)

  logger.info("Updated pipeline auto file after Step 1: %s", PIPELINE_AUTO_FILE)
  logger.info("STEP 1 completed successfully. Final gains file: %s", GAINS_FINAL)







# STEP 2: POL-CONVERT THE WHOLE EXPERIMENT.
if 2 in mysteps:

  if len(list(filter(lambda x: 'POLCONVERTER' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

  if not os.path.exists(PCONV_DIR):
    os.system('mkdir %s'%PCONV_DIR)

  NCPU = int(NCPU)
  if NCPU < 1:
    NCPU = multiprocessing.cpu_count() - 1
    
  SCRIPT_NAME = 'STEP2'
  XYG = 'POLCAL_GAINS_%s.dat'%(EXPNAME)

  if not APPLY_AMP:
    IFF = open(XYG,'rb')
    TEMP = pk.load(IFF)
    IFF.close()
    for anti in TEMP['XYratio'].keys():
      for ki in TEMP['XYratio'][anti].keys():
        TEMP['XYratio'][anti][ki][:] = 1.0
    OFF = open('POLCAL_GAINS_NOAMP_%s.dat'%(EXPNAME),'wb')
    pk.dump(TEMP,OFF)
    OFF.close()
    XYG = 'POLCAL_GAINS_NOAMP_%s.dat'%(EXPNAME)

  IFF = open('SOURCES_%s.txt'%EXPNAME)
  lines = IFF.readlines()
  IFF.close()
 
  SCANS = []
  REFANTS = []
  for li,line in enumerate(lines):
    if line.startswith(EXPNAME):
      TEMP = line.split()[0][:-1]
      SCANS.append(TEMP.split('_')[1])
      foundRef = False
      i = li+1
      while not foundRef and i < len(lines):
        if lines[i].startswith(EXPNAME) or 'SNR PASS' in lines[i]:
          REFANTS.append('')
          foundRef = True
        else:
          # skip blank/separator lines safely
          if len(lines[i].split()) == 0:
            i += 1
            continue
          if lines[i].lstrip().startswith('ANT ') and lines[i].split()[-1] == '+':
            foundRef = True
            REFANTS.append(lines[i].split()[1][:-1])
          else:
            i += 1
      if not foundRef:
        REFANTS.append('')

  for sci in range(len(SCANS)):
    if len(REFANTS[sci]) == 0:
      print('WARNING! SCAN %s DOES NOT HAVE ANY VALID ANTENNA!'%SCANS[sci])

  SCRIPT_NAMES = []

  log_header("STEP 2: PolConvert whole experiment")
  logger.info("Using gains file: %s", XYG)
  logger.info("Using XPOL_DELAYS for %d antennas", len(XPOL_DELAYS))
  logger.info("Number of scans to convert: %d", len(SCANS))
  logger.info("Output directory: %s", PCONV_DIR)
  logger.info("USE_PCAL used in Step 2: %s", USE_PCAL)
  logger.info("ZERO_PCALS used in Step 2: %s", ZERO_PCALS)

  for SCAN in SCANS:

    SCRIPT_NAME = 'STEP2_%s'%SCAN

   # os.system('cp -r %s %s'%(os.path.join(DIFX_DIR,'%s_%s*'%(EXPNAME,SCAN)), PCONV_DIR))

    keyw = {'EXPNAME':EXPNAME, 'ORIG_DIR':DIFX_DIR, 'DIFX_DIR':PCONV_DIR, 'XYGAINS':XYG, 
            'SUFFIX': SUFFIX, 'USE_PCAL':USE_PCAL,'SCAN_LIST':[SCAN],'ZERO_PCALS':ZERO_PCALS, 
            'XPOL_DELAYS':XPOL_DELAYS,
            'IF_OFFSET':int(IF_OFFSET), 'AC_WINDOW':int(PCAL_MED_WINDOW), 'XYPCALMODE':XYPCALMODE}
    keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys); keys.close()

    OFF = open('%s.py'%SCRIPT_NAME,'w')
    print(Start%SCRIPT_NAME,file=OFF)
    print('PConv.POLCONVERTER(**kww)',file=OFF)
    OFF.close()
    SCRIPT_NAMES.append(SCRIPT_NAME)

  def DO_PARALLEL(filename):
    print('GOING TO RUN %s'%filename)
    os.system(PYTHON_CALL%filename) 

  if NCPU>1:
    pool = multiprocessing.Pool(processes=NCPU)
    pool.map(DO_PARALLEL,SCRIPT_NAMES)
    pool.close()
    pool.join()
  else:
    for filename in SCRIPT_NAMES:
      DO_PARALLEL(filename)

 # for filename in SCRIPT_NAMES:
 #   os.system('rm -rf %s.py'%filename)
 # os.system('rm -rf keywords_STEP2_*.dat')

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)

  os.system('mv keywords_STEP2_*.dat STEP2*.py STEP_KEYWORDS/.')

  if os.path.exists('POLCONVERTER.FAILED'):
     raise Exception('STEP 2 FAILED!') 

  logger.info("STEP 2 completed successfully")





# STEP 3: PREPARE CF FILE:
if 3 in mysteps:

  log_header("STEP 3: Estimate additive phases and write CF file")
  logger.info("ADDITIVE_PHASE_SCANS: %s", ADDITIVE_PHASE_SCANS)
  logger.info("REFANT: %s", REFANT)
  logger.info("CF output filename: %s", CF_FILENAME)
  logger.info("HOPSNAMES entries: %d", len(HOPSNAMES))
  logger.info("SAMP_DELAYS entries: %d", len(SAMP_DELAYS))
  logger.info("EXCLUDE_BASELINE used in Step 3: %s", EXCLUDE_BASELINE)

  if len(list(filter(lambda x: 'GET_FOURFIT_PHASES' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

  if not os.path.exists('PYPHASES.PLOTS'):
    os.system('mkdir PYPHASES.PLOTS')

  for i,ADDITIVE_PHASE_SCAN in enumerate(ADDITIVE_PHASE_SCANS):
    SCRIPT_NAME = 'STEP3_%i'%i
    SCAN = os.path.join(PCONV_DIR,'%s_%s.difx'%(EXPNAME,ADDITIVE_PHASE_SCAN))

    logger.info("STEP 3 processing additive phase scan %s: %s", ADDITIVE_PHASE_SCAN, SCAN)

    keyw = {'SCAN':SCAN, 'HOPSNAMES': HOPSNAMES, 'IFNAMES': IFNAMES, 'FLAGBAS': EXCLUDE_BASELINE,
           'CALIB_BPASS':CALIB_BPASS, 'PCALDELAYS': PCAL_DELAYS, 'REFANT':REFANT, 
           'FLAG_PCALS':FLAG_PCALS, 'IF_OFFSET':IF_OFFSET, 'SAMP_DELAYS':SAMP_DELAYS,
           "MAX_PCAL_RMS":MAX_PCAL_RMS}

    keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys,protocol=0); keys.close()

    OFF = open('%s.py'%SCRIPT_NAME,'w')
    print(Start%SCRIPT_NAME,file=OFF)
    print('PYF.GET_FOURFIT_PHASES(**kww)',file=OFF)
    OFF.close()
    os.system(PYTHON_CALL%SCRIPT_NAME)
    os.system('mv cf_PyPhases cf_PyPhases_%i'%i)
    os.system('mv PyResults.dat PyResults_%i.dat'%i)

    if not os.path.exists('PYPHASES.PLOTS/SCAN_%i'%i):
      os.system('mkdir PYPHASES.PLOTS/SCAN_%i'%i)
    os.system('mv keywords_STEP3_%i.dat STEP3_%i.py STEP_KEYWORDS/.'%(i,i))
    os.system('mv BandPass_*.png *TEC.png PYPHASES.PLOTS/SCAN_%i/.'%i)

  Nscan = len(ADDITIVE_PHASE_SCANS)

  IFF = open('cf_PyPhases_0','r')
  lines = IFF.readlines()
  for i in range(len(lines)):
    if 'PyPhases. VERSION' in lines[i]:
      HeadLine = i+1
      break
  IFF.close()

  IFF = open('PyResults_0.dat','rb')
  Results = pk.load(IFF)
  IFF.close()

  for i in range(1,Nscan):
    IFF = open('PyResults_%i.dat'%i,'rb')
    ResTemp = pk.load(IFF)
    IFF.close()
    for obgain in ['DEL','PHAS','DEL_OFF','SBD']:
      for ant in ResTemp[obgain].keys():
        if ant not in Results[obgain].keys():
          Results[obgain][ant] = str(ResTemp[obgain][ant])

  OFF = open(CF_FILENAME,'w')
  for i in range(HeadLine):
    print(lines[i][:-1],file=OFF)

  if len(Results['DEL'].keys())>0:
    print('\n\n  *** SAMPLER DELAYS ***\n\n',file=OFF)
    for ant in Results['DEL'].keys():
      print(Results['DEL'][ant],file=OFF)

  if len(Results['SBD'].keys())>0:
    print('\n\n  *** AD-HOC SBDs ***\n\n',file=OFF)
    for ant in Results['SBD'].keys():
      print(Results['SBD'][ant],file=OFF)

  if len(Results['DEL_OFF'].keys())>0:
    print('\n\n  *** OFFSET DELAYS ***\n\n',file=OFF)
    for ant in Results['DEL_OFF'].keys():
      print(Results['DEL_OFF'][ant],file=OFF)
  
  if len(Results['PHAS'].keys())>0:
    print('\n\n  *** ADDITIVE PHASES ***\n\n',file=OFF)
    for ant in Results['PHAS'].keys():
      print(Results['PHAS'][ant],file=OFF)

  OFF.close()

  if os.path.exists('GET_FOURFIT_PHASES.FAILED'):
     raise Exception('STEP 3 FAILED!') 

  logger.info("STEP 3 completed successfully. Wrote CF file: %s", CF_FILENAME)





# STEP 4: CALIBRATE BPASS AND REMOVE IONEX-BASED TEC:
if 4 in mysteps:
  if len(list(filter(lambda x: 'REMOVE_TEC' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

  SCRIPT_NAMES = []

  IFF = open('SOURCES_%s.txt'%EXPNAME)
  lines = IFF.readlines()
  SCANS = []
  REFANTS = []
 
  IFF.close()
 
  for li,line in enumerate(lines):
    if line.startswith(EXPNAME):
      TEMP = line.split()[0][:-1]
      SCANS.append(TEMP.split('_')[1])

  log_header("STEP 4: Calibrate bandpass and remove IONEX TEC")
  logger.info("Input directory: %s", PCONV_DIR)
  logger.info("Output directory: %s", BPCAL_DIR)
  logger.info("APPLY_PHASECAL: %s", APPLY_PHASECAL)
  logger.info("ADDITIVE_PHASE_SCANS: %s", ADDITIVE_PHASE_SCANS)
  logger.info("REFANT: %s", REFANT)
  logger.info("Number of scans: %d", len(SCANS))

  for SCAN in SCANS:

    SCRIPT_NAME = 'STEP4_%s'%SCAN

   # os.system('cp -r %s %s'%(os.path.join(DIFX_DIR,'%s_%s*'%(EXPNAME,SCAN)), PCONV_DIR))

    keyw = {'EXPNAME':EXPNAME, 'SCAN':SCAN, 'ORIG_DIR':PCONV_DIR, 'DEST_DIR':BPCAL_DIR,
            'APPLY_PHASECAL':APPLY_PHASECAL, 'WRITE_DATA':True, 'FLAG_PCALS':FLAG_PCALS,
            'REFSCAN':ADDITIVE_PHASE_SCANS, 'REFANT':REFANT,'FLAGBAS': EXCLUDE_BASELINE, 'IF_OFFSET':IF_OFFSET,
            'SAMP_DELAYS':SAMP_DELAYS, 'HOPS_NAMES':HOPSNAMES}
    keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys); keys.close()

    OFF = open('%s.py'%SCRIPT_NAME,'w')
    print(Start%SCRIPT_NAME,file=OFF)
    print('PYF.removeTEC(**kww)',file=OFF)
    OFF.close()
    SCRIPT_NAMES.append(SCRIPT_NAME)

  def DO_PARALLEL(filename):

    print('GOING TO RUN %s'%filename)
    os.system(PYTHON_CALL%filename) 

  if NCPU>1:
    pool = multiprocessing.Pool(processes=NCPU)
    pool.map(DO_PARALLEL,SCRIPT_NAMES)
    pool.close()
    pool.join()
  else:
    for filename in SCRIPT_NAMES:
      DO_PARALLEL(filename)

  for filename in SCRIPT_NAMES:
    os.system('rm -rf %s.py'%filename)


  os.system('rm -rf keywords_STEP4_*.dat')

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)

  if os.path.exists('REMOVE_TEC.FAILED'):
     raise Exception('STEP 4 FAILED!') 

  logger.info("STEP 4 completed successfully. Output directory: %s", BPCAL_DIR)




## Perform Global Fringe Fitting:
if 5 in mysteps:
  if len(list(filter(lambda x: 'DO_GFF' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

  SCRIPT_NAMES = []

  IFF = open('SOURCES_%s.txt'%EXPNAME)
  lines = IFF.readlines()
  SCANS = []
  REFANTS = []
 
  IFF.close()
 
  for li,line in enumerate(lines):
    if line.startswith(EXPNAME):
      TEMP = line.split()[0][:-1]
      SCANS.append(TEMP.split('_')[1])

  log_header("STEP 5: Broadband global fringe fitting")
  logger.info("Input directory: %s", BPCAL_DIR)
  logger.info("CF file: %s", CF_FILENAME)
  logger.info("GFF_REFANTS: %s", GFF_REFANTS)
  logger.info("ANT_WEIGHTS: %s", ANT_WEIGHTS)
  logger.info("EXCLUDE_BASELINE: %s", EXCLUDE_BASELINE)
  logger.info("Number of scans: %d", len(SCANS))

  for SCAN in sorted(SCANS):

    SCRIPT_NAME = 'STEP5_%s'%SCAN

   # os.system('cp -r %s %s'%(os.path.join(DIFX_DIR,'%s_%s*'%(EXPNAME,SCAN)), PCONV_DIR))

    keyw = {'EXPNAME':EXPNAME, 'SCAN':SCAN, 'DIR':BPCAL_DIR, 'ANT_WEIGHTS':ANT_WEIGHTS,
            'APPLY_PHASECAL':APPLY_PHASECAL, 'FLAG_PCALS':FLAG_PCALS,
            "MAX_PCAL_RMS":MAX_PCAL_RMS,
            'CF_FILE':CF_FILENAME, 'REFANTS':GFF_REFANTS,'FLAGBAS': EXCLUDE_BASELINE, 
            'IF_OFFSET':IF_OFFSET, 'SAMP_DELAYS':SAMP_DELAYS, "HOPS_NAMES":HOPSNAMES}
    keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys); keys.close()

    OFF = open('%s.py'%SCRIPT_NAME,'w')
    print(Start%SCRIPT_NAME,file=OFF)
    print('PYF.DO_GFF(**kww)',file=OFF)
    OFF.close()
    SCRIPT_NAMES.append(SCRIPT_NAME)

  def DO_PARALLEL(filename):
    print('GOING TO RUN %s'%filename)
    os.system(PYTHON_CALL%filename) 

  if NCPU>1:
    pool = multiprocessing.Pool(processes=NCPU)
    pool.map(DO_PARALLEL,SCRIPT_NAMES)
    pool.close()
    pool.join()
  else:
    for filename in SCRIPT_NAMES:
      DO_PARALLEL(filename)

  for filename in SCRIPT_NAMES:
    os.system('rm -rf %s.py'%filename)

  os.system('rm -rf keywords_STEP5_*.dat')

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)

  if os.path.exists('DO_GFF.FAILED'):
     raise Exception('STEP 5 FAILED!') 

  logger.info("STEP 5 completed successfully")




# STEP 6: CALIBRATE THE DATA COMPLETELY:
if 6 in mysteps:
  if len(list(filter(lambda x: 'WRITE_CALIBRATED' not in x, glob.glob('*.FAILED'))))>0:
    raise Exception('ANOTHER TASK FAILED PREVIOUSLY. WILL ABORT UNTIL YOU SOLVE IT!')      

  SCRIPT_NAMES = []

  IFF = open('SOURCES_%s.txt'%EXPNAME)
  lines = IFF.readlines()
  SCANS = []
  REFANTS = []
  IFF.close()
 
  for li,line in enumerate(lines):
    if line.startswith(EXPNAME):
      TEMP = line.split()[0][:-1]
      SCANS.append(TEMP.split('_')[1])

  log_header("STEP 6: Write fully calibrated SWIN files")
  logger.info("Input directory: %s", PCONV_DIR)
  logger.info("Output directory: %s", FINAL_DIR)
  logger.info("CF file: %s", CF_FILENAME)
  logger.info("ADDITIVE_PHASE_SCANS: %s", ADDITIVE_PHASE_SCANS)
  logger.info("REFANT: %s", REFANT)
  logger.info("APPLY_GFF: True")
  logger.info("Number of scans: %d", len(SCANS))

  for SCAN in SCANS:

    SCRIPT_NAME = 'STEP6_%s'%SCAN

   # os.system('cp -r %s %s'%(os.path.join(DIFX_DIR,'%s_%s*'%(EXPNAME,SCAN)), PCONV_DIR))

    keyw = {'EXPNAME':EXPNAME, 'SCAN':SCAN, 'ORIG_DIR':PCONV_DIR, 'DEST_DIR':FINAL_DIR,
            'APPLY_PHASECAL':True, 'WRITE_DATA':True, 'FLAG_PCALS':FLAG_PCALS,
            'REFSCAN':ADDITIVE_PHASE_SCANS, 'REFANT':REFANT,'FLAGBAS': EXCLUDE_BASELINE, 
            'IF_OFFSET':IF_OFFSET, 'SAMP_DELAYS':SAMP_DELAYS, "APPLY_GFF":True, 
            "CF_FILE":CF_FILENAME,"HOPS_NAMES":HOPSNAMES}
    keys = open('keywords_%s.dat'%SCRIPT_NAME,'wb'); pk.dump(keyw, keys); keys.close()

    OFF = open('%s.py'%SCRIPT_NAME,'w')
    print(Start%SCRIPT_NAME,file=OFF)
    print('PYF.removeTEC(**kww)',file=OFF)
    OFF.close()
    SCRIPT_NAMES.append(SCRIPT_NAME)

  def DO_PARALLEL(filename):
    print('GOING TO RUN %s'%filename)
    os.system(PYTHON_CALL%filename) 

  if NCPU>1:
    pool = multiprocessing.Pool(processes=NCPU)
    pool.map(DO_PARALLEL,SCRIPT_NAMES)
    pool.close()
    pool.join()
  else:
    for filename in SCRIPT_NAMES:
      DO_PARALLEL(filename)

  for filename in SCRIPT_NAMES:
    os.system('rm -rf %s.py'%filename)

  os.system('rm -rf keywords_STEP6_*.dat')

  newlogs = glob.glob('*.log')
  for log in newlogs:
    if log not in currlogs:
      os.system('mv %s LOGS/.'%log)

  if os.path.exists('WRITE_CALIBRATED.FAILED'):
     raise Exception('STEP 6 FAILED!') 

  logger.info("STEP 6 completed successfully. Output directory: %s", FINAL_DIR)


logger.info("Pipeline finished for requested mysteps=%s", mysteps)

for handler in logger.handlers:
    handler.flush()

steps_tag = "steps_" + "-".join(str(s) for s in mysteps)
final_master_log = unique_log_copy_name(
    "%s_master_script_%s_final.log_masterscript_run" % (EXPNAME, steps_tag)
)

logger.info("Copying master log to: %s", final_master_log)
for handler in logger.handlers:
    handler.flush()

import shutil
shutil.copy2(MASTER_LOG, final_master_log)
########### +++++++++++++++++++++++++++++++++++++++++++++  ###########
