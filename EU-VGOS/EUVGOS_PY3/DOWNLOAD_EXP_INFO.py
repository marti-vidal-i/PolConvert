# DOWNLOAD_EXP_INFO.py
# Utilities for downloading/parsing antenna coordinates from CDDIS VEX files.
import os, re
from ftplib import FTP_TLS


def _download_cddis_file(exp_code, exp_year, remote_fname, local_path):
    """Download a file from the CDDIS VLBI AUX archive using FTP_TLS."""
    host = "gdc.cddis.eosdis.nasa.gov"
    remote_dir = "vlbi/ivsdata/aux/%04d/%s/" % (int(exp_year), str(exp_code))
    os.makedirs(os.path.dirname(local_path) or ".", exist_ok=True)
    ftps = FTP_TLS(host=host)
    ftps.login(user="anonymous", passwd="anonymous@")
    ftps.prot_p()
    ftps.cwd(remote_dir)
    with open(local_path, "wb") as fh:
        ftps.retrbinary("RETR " + remote_fname, fh.write)
    try:
        ftps.quit()
    except Exception:
        pass
    try:
        if os.path.exists(local_path) and os.path.getsize(local_path) == 0:
            os.remove(local_path)
    except Exception:
        pass


def _parse_ant_coords_from_vex(vex_path):
    """
    Parse VEX blocks for:
      site_ID = Gs;
      site_position = X m : Y m : Z m;
    Returns dict { "GS": (x, y, z), ... }
    """
    re_id  = re.compile(r"site_ID\s*=\s*([A-Za-z0-9]{2})\s*;")
    re_pos = re.compile(
        r"site_position\s*=\s*([-\d.]+)\s*m\s*:\s*([-\d.]+)\s*m\s*:\s*([-\d.]+)\s*m\s*;"
    )
    coords = {}
    in_def = False
    cur_id = None
    cur_pos = None
    with open(vex_path, "r", errors="replace") as f:
        for line in f:
            s = line.strip().lower()
            if s.startswith("def "):
                in_def = True
                cur_id = None
                cur_pos = None
                continue
            if in_def and s.startswith("enddef"):
                if cur_id and cur_pos:
                    coords[cur_id] = cur_pos
                in_def = False
                continue
            if not in_def:
                continue
            m = re_id.search(line)
            if m:
                cur_id = m.group(1).upper()
                continue
            m = re_pos.search(line)
            if m:
                x, y, z = map(float, m.groups())
                cur_pos = (x, y, z)
                continue
    return coords


def ensure_ant_coords_from_vex(
    *,
    exp_code,
    exp_year,
    user_ant_coords=None,
    vex_cache_dir=".",
    verbose=True,
):
    """
      - Always try to obtain VEX coords (download if missing).
      - VEX coordinates form the base.
      - user_ant_coords override any VEX entry.
      - If VEX is unavailable, return only the user coords.
    """
    user_ant_coords = dict(user_ant_coords) if user_ant_coords else {}
    os.makedirs(vex_cache_dir, exist_ok=True)
    vex_local = os.path.join(vex_cache_dir, "%s.vex" % str(exp_code))
    # Ensure VEX exists (download if missing)
    if not os.path.exists(vex_local) or os.path.getsize(vex_local) == 0:
        if verbose:
            print("[VEX] Downloading %s.vex from CDDIS..." % str(exp_code))
        try:
            _download_cddis_file(exp_code, exp_year, "%s.vex" % str(exp_code), vex_local)
        except Exception as e:
            if verbose:
                print("[VEX] WARNING: download failed (%s)." % str(e))
                print("[VEX] Using only user-provided ANT_COORDS.")
            return user_ant_coords
    # Parse VEX
    try:
        vex_coords = _parse_ant_coords_from_vex(vex_local)
    except Exception as e:
        if verbose:
            print("[VEX] WARNING: parsing failed (%s). Using only user coords." % str(e))
        return user_ant_coords
    # Merge: VEX base, user overrides
    #merged = dict(vex_coords)
    #merged.update(user_ant_coords)
    # Merge: user fallback, VEX overrides when available
    merged = dict(user_ant_coords)
    merged.update(vex_coords)
    if verbose:
        print("[VEX] Loaded %d antenna positions (VEX=%d, user_override=%d)" % (len(merged), len(vex_coords), len(user_ant_coords)))
    return merged

def _ensure_corr_file(exp_code, exp_year, corr_cache_dir=".", verbose=True):
    """
    Ensure <exp_code>.corr exists locally (download if missing/empty).
    Returns local file path, or None if cannot be obtained.
    """
    os.makedirs(corr_cache_dir, exist_ok=True)
    corr_local = os.path.join(corr_cache_dir, "%s.corr" % str(exp_code))
    if (not os.path.exists(corr_local)) or os.path.getsize(corr_local) == 0:
        if verbose:
            print("[CORR] Downloading %s.corr from CDDIS..." % str(exp_code))
        try:
            _download_cddis_file(exp_code, exp_year, "%s.corr" % str(exp_code), corr_local)
        except Exception as e:
            if verbose:
                print("[CORR] WARNING: download failed (%s: %s)" % (type(e).__name__, str(e)))
            return None
    return corr_local


def _parse_hopsnames_and_sampler_delays_from_corr(corr_path, verbose=True):
    """
    Parse:
      +STATIONS ... mk4 table  -> HOPSNAMES: { 'GS':'G', ... }
    And sampler delays blocks (mk4-based):
      if station G
        sampler_delay_x -195 130 130 130
    Returns:
      hopsnames: dict[str,str]   # SWIN(2-letter upper) -> mk4(1-letter)
      samp_delays: dict[str,list[float]]  # SWIN(2-letter upper) -> [d1,d2,d3,d4] from sampler_delay_x
    """
    hopsnames = {}
    station_id_by_mk4 = {}   # mk4 -> SWIN two-letter upper
    # --- parse +STATIONS table ---
    in_stations = False
    # Example row:
    # Gs      GGAO12M  G
    # (station_id) (site_name) (mk4)
    row_re = re.compile(r"^\s*([A-Za-z0-9]{2})\s+\S+\s+([A-Za-z0-9])\s*$")
    # --- parse sampler delays blocks ---
    # if station G
    if_station_re = re.compile(r"^\s*if\s+station\s+([A-Za-z0-9])\s*$", re.IGNORECASE)
    sdx_re = re.compile(r"^\s*sampler_delay_x\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$", re.IGNORECASE)
    # We first read whole file lines (simpler state machine)
    with open(corr_path, "r", errors="replace") as f:
        lines = f.readlines()
    # Pass 1: STATIONS
    for line in lines:
        s = line.strip()
        if s.startswith("+STATIONS"):
            in_stations = True
            continue
        if in_stations:
            # End of block (commonly starts with '-' section or another '+SECTION')
            if s.startswith("+") and not s.startswith("+STATIONS"):
                in_stations = False
                continue
            if s.startswith("-") or s.startswith("station") or s.startswith("----") or s == "":
                continue
            m = row_re.match(line)
            if m:
                st2 = m.group(1).upper()
                mk4 = m.group(2).upper()
                hopsnames[st2] = mk4
                station_id_by_mk4[mk4] = st2
    if verbose:
        print("[CORR] Parsed STATIONS: %d entries" % len(hopsnames))
    # Pass 2: sampler_delay_x blocks
    samp_delays = {}
    current_mk4 = None
    for line in lines:
        m = if_station_re.match(line)
        if m:
            current_mk4 = m.group(1).upper()
            continue
        m = sdx_re.match(line)
        if m and current_mk4:
            delays = [float(m.group(i)) for i in range(1, 5)]
            st2 = station_id_by_mk4.get(current_mk4)
            if st2:
                samp_delays[st2] = delays
            # keep current_mk4 until next "if station" (matches your corr style)
    if verbose:
        print("[CORR] Parsed sampler_delay_x: %d stations" % len(samp_delays))
    return hopsnames, samp_delays


_CORR_CACHE = {}  # key: (exp_code, exp_year, corr_cache_dir) -> (hops, delays)
def _load_corr_info(exp_code, exp_year, corr_cache_dir=".", verbose=True):
    key = (str(exp_code), int(exp_year), str(corr_cache_dir))
    if key in _CORR_CACHE:
        return _CORR_CACHE[key]
    corr_path = _ensure_corr_file(exp_code, exp_year, corr_cache_dir=corr_cache_dir, verbose=verbose)
    if corr_path is None:
        _CORR_CACHE[key] = (None, None)
        return None, None
    try:
        hops, delays = _parse_hopsnames_and_sampler_delays_from_corr(corr_path, verbose=verbose)
    except Exception as e:
        if verbose:
            print("[CORR] WARNING: failed parsing (%s: %s)" % (type(e).__name__, str(e)))
        _CORR_CACHE[key] = (None, None)
        return None, None
    _CORR_CACHE[key] = (hops, delays)
    return hops, delays


def ensure_hopsnames_from_corr(
    *,
    exp_code,
    exp_year,
    user_hopsnames=None,
    corr_cache_dir=".",
    verbose=True
):
    """
      - Read HOPSNAMES from <exp_code>.corr (download if needed).
      - CORR values overwrite user_hopsnames (experiment-specific).
      - If CORR unavailable, return user_hopsnames.
    """
    user_hopsnames = dict(user_hopsnames) if user_hopsnames else {}
    corr_hops, _ = _load_corr_info(exp_code, exp_year, corr_cache_dir=corr_cache_dir, verbose=verbose)
    if corr_hops is None:
        return user_hopsnames
    merged = dict(user_hopsnames)
    merged.update(corr_hops)
    return merged


def ensure_samp_delays_from_corr(
    *,
    exp_code,
    exp_year,
    user_samp_delays=None,
    corr_cache_dir=".",
    verbose=True
):
    """
      - Read sampler_delay_x from <exp_code>.corr (download if needed).
      - Returned dict is keyed by SWIN 2-letter station code (uppercase), like 'GS','WF',...
      - CORR values overwrite user_samp_delays (experiment-specific).
      - If CORR unavailable, return user_samp_delays.
    """
    user_samp_delays = dict(user_samp_delays) if user_samp_delays else {}
    _, corr_delays = _load_corr_info(exp_code, exp_year, corr_cache_dir=corr_cache_dir, verbose=verbose)
    if corr_delays is None:
        return user_samp_delays
    merged = dict(user_samp_delays)
    merged.update(corr_delays)
    return merged
