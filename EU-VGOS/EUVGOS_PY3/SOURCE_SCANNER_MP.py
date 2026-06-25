# SOURCE_SCANNER: SIMPLE CHECK OF SNR FROM A SET OF SWIN FILES.
#            Copyright (C) 2022  Ivan Marti-Vidal & Ezequiel Albentosa-Ruiz
#            University of Valencia (Spain)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#
#


import numpy as np
import os, glob, re
import pickle as pk
from multiprocessing import Pool
import math
from collections import defaultdict



# ------------------------------------------------------------------------
# Extra scan-quality metrics for POLCAL selection:
#  - Elevation statistics per antenna (from .im)
#  - Parallactic angle coverage per antenna (from .im + station latitude)
#  - Polarization-product completeness (XX/XY/YX/YY presence)
#  - u/v coverage and baseline diversity (from SWIN rows)
#
# These metrics do NOT change the original SNR-based selection logic.
# They are returned in the result dict under key "metrics" and printed as
# "EXTRA METRICS" in the output text file.
# ------------------------------------------------------------------------

WGS84_A = 6378137.0
WGS84_E2 = 6.69437999014e-3

def ecef_to_geodetic_lat_lon(x, y, z):
    """
    Convert station ECEF coordinates (meters) to geodetic latitude/longitude (radians).
    Only latitude is needed to compute parallactic angle from az/el.
    """
    lon = math.atan2(y, x)
    p = math.sqrt(x*x + y*y)
    # iterative lat
    lat = math.atan2(z, p * (1 - WGS84_E2))
    for _ in range(6):
        sinlat = math.sin(lat)
        N = WGS84_A / math.sqrt(1 - WGS84_E2 * sinlat*sinlat)
        lat = math.atan2(z + WGS84_E2 * N * sinlat, p)
    return lat, lon

def parallactic_angle_from_az_el_lat(az_rad, el_rad, lat_rad):
    """
    Parallactic angle q (radians) from azimuth/elevation and station latitude.
    Used to estimate PA coverage for polarization calibration (D-terms).
    Convention: q = atan2( sin(A), tan(phi)*cos(el) - sin(el)*cos(A) )
    """
    A = az_rad
    e = el_rad
    phi = lat_rad
    num = math.sin(A)
    den = math.tan(phi) * math.cos(e) - math.sin(e) * math.cos(A)
    return math.atan2(num, den)

def parse_station_ecef_from_input(input_path):
    """
    Best-effort parse of antenna ECEF coords from DiFX .input file.
    If the .input doesn't contain positions, dict may be empty.
    Different DiFX builds can format keys differently; this is tolerant:
    - Looks for lines containing TELESCOPE ... X/Y/Z (m) or TELESCOPE ... X (m):
    - Returns dict: antenna_index -> (x,y,z) in meters
    """
    idx_to_xyz = {}
    # several possible patterns:
    # "TELESCOPE 0 X (m):  1234"
    # "TELESCOPE 0 X: 1234"
    # "TELESCOPE X 0: 1234"
    pat1 = re.compile(r"^TELESCOPE\s+(\d+)\s+([XYZ])(?:\s+\(m\))?:\s*([-\d.eE+]+)")
    pat2 = re.compile(r"^TELESCOPE\s+([XYZ])\s+(\d+):\s*([-\d.eE+]+)")
    # Some .input files store position as "TELESCOPE POSITION <i>: x y z"
    pat3 = re.compile(r"^TELESCOPE\s+POSITION\s+(\d+):\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)")

    with open(input_path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            m = pat3.match(line)
            if m:
                i = int(m.group(1))
                x = float(m.group(2)); y = float(m.group(3)); z = float(m.group(4))
                idx_to_xyz[i] = (x, y, z)
                continue

            m = pat1.match(line)
            if m:
                i = int(m.group(1))
                comp = m.group(2)
                val = float(m.group(3))
                if i not in idx_to_xyz:
                    idx_to_xyz[i] = [None, None, None]
                if comp == "X": idx_to_xyz[i][0] = val
                if comp == "Y": idx_to_xyz[i][1] = val
                if comp == "Z": idx_to_xyz[i][2] = val
                continue

            m = pat2.match(line)
            if m:
                comp = m.group(1)
                i = int(m.group(2))
                val = float(m.group(3))
                if i not in idx_to_xyz:
                    idx_to_xyz[i] = [None, None, None]
                if comp == "X": idx_to_xyz[i][0] = val
                if comp == "Y": idx_to_xyz[i][1] = val
                if comp == "Z": idx_to_xyz[i][2] = val
                continue

    # normalize list->tuple if complete
    out = {}
    for i, xyz in idx_to_xyz.items():
        if isinstance(xyz, list):
            if any(v is None for v in xyz):
                continue
            out[i] = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
        else:
            out[i] = xyz
    return out

def parse_im_azel(im_path):
    """
    Parse azimuth and geometric elevation from the DiFX .im file.
    Returns dict: ant_index -> {"az": [deg], "el": [deg]}.
    Note: format can vary; we take the first numeric value on each line.
    """
    az_re = re.compile(r"^SRC\s+\d+\s+ANT\s+(\d+)\s+AZ\b\s*:?\s*(.*)$")
    el_re = re.compile(r"^SRC\s+\d+\s+ANT\s+(\d+)\s+EL\s+GEOM\b\s*:?\s*(.*)$")

    out = defaultdict(lambda: {"az": [], "el": []})

    with open(im_path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            m = az_re.match(line)
            if m:
                ai = int(m.group(1))
                # take first float on RHS
                parts = m.group(2).split()
                if parts:
                    try:
                        out[ai]["az"].append(float(parts[0]))
                    except Exception:
                        pass
                continue
            m = el_re.match(line)
            if m:
                ai = int(m.group(1))
                parts = m.group(2).split()
                if parts:
                    try:
                        out[ai]["el"].append(float(parts[0]))
                    except Exception:
                        pass
                continue
    return out

def parse_output_nchan_by_if_from_lines(lines):
    """
    Return output spectral channels per DiFX frequency index.

    DiFX .input often stores the raw spectral channels as:
        NUM CHANNELS i: N
    and the channels written to SWIN as:
        N / CHANS TO AVG i

    Older/simple jobs may have no CHANS TO AVG lines; then avg=1 and
    this reduces to the old behavior.
    """
    num_channels = {}
    chans_to_avg = {}

    pat_nchan = re.compile(r"^NUM CHANNELS\s+(\d+):\s+(\d+)\s*$")
    pat_avg = re.compile(r"^CHANS TO AVG\s+(\d+):\s+(\d+)\s*$")

    for li in lines:
        s = li.strip()

        m = pat_nchan.match(s)
        if m:
            num_channels[int(m.group(1))] = int(m.group(2))
            continue

        m = pat_avg.match(s)
        if m:
            chans_to_avg[int(m.group(1))] = int(m.group(2))
            continue

    if len(num_channels) == 0:
        raise RuntimeError("No NUM CHANNELS i entries found in .input")

    nchan_by_if = {}
    for ifidx, nraw in num_channels.items():
        avg = int(chans_to_avg.get(ifidx, 1))

        if avg <= 0:
            raise RuntimeError("Bad CHANS TO AVG for IF %d: %d" % (ifidx, avg))

        if int(nraw) % avg != 0:
            raise RuntimeError(
                "NUM CHANNELS not divisible by CHANS TO AVG for IF %d: %d / %d"
                % (ifidx, int(nraw), avg)
            )

        nchan_by_if[int(ifidx)] = int(nraw) // avg

    return nchan_by_if, num_channels, chans_to_avg


def infer_difx_if_base(doif, if_values_in_rows, if_offset=0, requested="auto"):
    """
    Infer whether SWIN row IF/fridx values are 1-based or 0-based relative
    to the user-facing logical IF list.

    Returns 1 or 0.

    With logical DOIF=1..32 and IF_OFFSET=32:
      base=1 tests 1..32 and 33..64
      base=0 tests 0..31 and 32..63

    The score counts how many requested main/shifted physical IFs are
    actually present in the binary rows.
    """
    if str(requested).lower() not in ("auto", "none"):
        return int(requested)

    rows = set(int(x) for x in if_values_in_rows)
    doif_list = [int(x) for x in doif]
    if_offset = int(if_offset)

    def score(base):
        n = 0
        for IFp in doif_list:
            phys = IFp if base == 1 else IFp - 1
            if phys in rows:
                n += 1
            if if_offset != 0 and (phys + if_offset) in rows:
                n += 1
        return n

    score1 = score(1)
    score0 = score(0)

    # Tie-breaker preserves old SOURCE_SCANNER behavior.
    return 0 if score0 > score1 else 1


def physical_if_to_logical_if_with_base(ifc, doif, if_offset=0, difx_if_base=1):
    """
    Map binary SWIN IF/fridx back to logical 1-based IF.

    For base=0, binary IF 0 -> logical IF 1, and binary IF 32 with
    IF_OFFSET=32 -> logical IF 1.
    """
    ifc = int(ifc)

    if int(difx_if_base) == 0:
        return physical_if_to_logical_if(ifc + 1, doif, if_offset)

    return physical_if_to_logical_if(ifc, doif, if_offset)


def if_to_band_range(if_id, nif_per_band=8, nif_total=None):
    """
    Expand one IF to its full band range.
    Example, nif_per_band=8:
      IF 21 -> [17,18,19,20,21,22,23,24]
    """
    if_id = int(if_id)
    nif_per_band = int(nif_per_band)
    b0 = ((if_id - 1) // nif_per_band) * nif_per_band + 1
    b1 = b0 + nif_per_band - 1
    if nif_total is not None and int(nif_total) > 0:
        b1 = min(b1, int(nif_total))
    return list(range(b0, b1 + 1))

def physical_if_to_logical_if(ifc, doif, if_offset=0):
    """
    Map physical DiFX IF number back to the logical IF number used by DOIF.
    Example with IF_OFFSET=32: physical IF 37 -> logical IF 5.
    """
    ifc = int(ifc)
    doif_set = set(int(x) for x in doif)
    if_offset = int(if_offset)

    if ifc in doif_set:
        return ifc

    shifted_back = ifc - if_offset
    if if_offset != 0 and shifted_back in doif_set:
        return shifted_back

    return ifc


def pol_code_for_station_product(product, station_index):
    """
    Return which local pol is required from a station for a visibility product.
    product: "XX", "XY", "YX", "YY"
    station_index: 0 for first antenna in baseline, 1 for second antenna.
    """
    if product == "XX":
        return "X"
    if product == "YY":
        return "Y"
    if product == "XY":
        return "X" if station_index == 0 else "Y"
    if product == "YX":
        return "Y" if station_index == 0 else "X"
    return None


def infer_station_pol_issues(missing_products):
    """
    Infer station-local missing polarization from missing baseline products.
    missing_products example:
      ["XX", "YX"] on baseline A-B implies station B is missing X.
      ["XY", "YY"] on baseline A-B implies station B is missing Y.
      ["XX", "XY"] on baseline A-B implies station A is missing X.
      ["YX", "YY"] on baseline A-B implies station A is missing Y.
    Returns:
      {0: set(["X","Y"]), 1: set(["X","Y"])}
    """
    miss = set(missing_products)
    out = {0: set(), 1: set()}
    for st in (0, 1):
        for pol in ("X", "Y"):
            needed = []
            for prod in ("XX", "XY", "YX", "YY"):
                if pol_code_for_station_product(prod, st) == pol:
                    needed.append(prod)
            if all(prod in miss for prod in needed):
                out[st].add(pol)
    return out


def getFringeSNR(dd, cfg):
    minSNR = float(cfg["SNRCut"])
    DOIF = cfg["DOIF"]
    IF_OFF = int(cfg["IF_OFFSET"])
    MIN_ELEV_DEG = cfg.get("MIN_ELEV_DEG", None)
    REQUIRE_CROSSHANDS = cfg.get("REQUIRE_CROSSHANDS", None)
    MIN_UNIQUE_BASELINES = cfg.get("MIN_UNIQUE_BASELINES", None)
    ANT_COORDS = cfg["ANT_COORDS"]
    MIN_UV_SPAN = cfg.get("MIN_UV_SPAN", None)
    MAX_WEIGHT_ZERO_FRAC = cfg.get("MAX_WEIGHT_ZERO_FRAC", None)
    MIN_WEIGHT_MEDIAN = cfg.get("MIN_WEIGHT_MEDIAN", None)
    DIFX_IF_BASE_REQUESTED = cfg.get("DIFX_IF_BASE", "auto")

    warnings = []
    scan_name = os.path.basename(dd).split(".")[0]
    def _fatal(msg, exc=None):
        if exc is None:
            err = "%s" % msg
        else:
            err = "%s (%s: %s)" % (msg, type(exc).__name__, str(exc))
        return {
            "ok": False,
            "fatal": True,
            "scan": scan_name,
            "result": None,
            "warnings": warnings,
            "error": err,
            "dd": dd,
        }
    
    # Figure out number of antennas, IFs and channels per IF:
    inp = "%s.input" % dd[:-5]
    try:
        with open(inp, "r") as temp:
            lines = temp.readlines()
    except Exception as e:
        return _fatal("Cannot open input file %s" % inp, e)

    Nants = 0
    NIF_scan = 0
    Nchan = 0
    ANT_NAMES = []
    nchan_by_if = {}
    num_channels_raw = {}
    chans_to_avg = {}

    try:
        # Use the actual output channel count written to SWIN.
        nchan_by_if, num_channels_raw, chans_to_avg = parse_output_nchan_by_if_from_lines(lines)

        unique_output_nchan = sorted(set(int(v) for v in nchan_by_if.values()))
        if len(unique_output_nchan) != 1:
            return _fatal(
                "Variable output channel count per IF is not supported by the fast SWIN reader: %s"
                % ",".join(str(x) for x in unique_output_nchan)
            )

        Nchan = int(unique_output_nchan[0])

        for li in lines:
            if li.startswith("TELESCOPE ENTRIES:"):
                Nants = max([Nants, int(li.split()[-1])])
                ANT_NAMES = ["%02i" % i for i in range(1, Nants + 1)]

            if li.startswith("TELESCOPE NAME "):
                tempLine = li.split()
                tid = int(tempLine[2].replace(":", ""))
                NAM = tempLine[-1].replace("\n", "")
                ANT_NAMES[tid] = NAM

            if li.startswith("FREQ ENTRIES:"):
                NIF_scan = int(li.split()[-1])

    except Exception as e:
        return _fatal("Failed parsing .input file", e)

    if Nants <= 1 or Nchan <= 0:
        return _fatal("Bad metadata from input: Nants=%d, Nchan=%d" % (Nants, Nchan))

    if len(chans_to_avg) > 0 and any(int(v) != 1 for v in chans_to_avg.values()):
        warnings.append(
            "%s: using output Nchan=%d from NUM CHANNELS / CHANS TO AVG for SWIN stride"
            % (scan_name, Nchan)
        )

    # try to load station coordinates (ECEF) for parallactic angle
    coords_by_name = {}
    coords_missing = []
    try:
        idx_xyz = parse_station_ecef_from_input(inp)
    except Exception:
        idx_xyz = {}

    # Map to antenna names if possible
    for ai in range(Nants):
        name = ANT_NAMES[ai]
        if name in ANT_COORDS:
            coords_by_name[name] = ANT_COORDS[name]
        elif ai in idx_xyz:
            coords_by_name[name] = idx_xyz[ai]
        else:
            coords_missing.append(name)

    # Read visibilities:
    difx_glob = glob.glob("%s.difx/DIFX*" % dd[:-5])
    if len(difx_glob) == 0:
        return _fatal("Cannot find DIFX file for %s (glob: %s.difx/DIFX*)" % (dd, dd[:-5]))
    DIFXFile = difx_glob[0]
    try:
        frfile = open(DIFXFile, "rb")
    except Exception as e:
        return _fatal("Cannot open DIFX file %s" % DIFXFile, e)
    # consume 8 bytes header
    try:
        _ = frfile.read(8)
    except Exception as e:
        frfile.close()
        return _fatal("Cannot read DIFX header %s" % DIFXFile, e)

    # ID list of baselines:
    SNRBas = []
    for i in range(Nants - 1):
        for j in range(i + 1, Nants):
            SNRBas.append([i, j])

    print("Processing %s" % os.path.basename(DIFXFile))

    dtype2 = np.dtype(
        [
            ("BAS", np.int32),
            ("MJD", np.int32),
            ("SEC", np.float64),
            ("AUX1", np.int32),
            ("AUX2", np.int32),
            ("IF", np.int32),
            ("POL1", np.dtype("b")),
            ("POL2", np.dtype("b")),
            ("AUX3", np.int32),
            ("WEIGHT", np.float64),
            ("U", np.float64),
            ("V", np.float64),
            ("W", np.float64),
            ("VISIB", np.complex64, Nchan + 1),
        ]
    )

    # Stupid numpy error:
    try:
        fringe2 = np.fromfile(frfile, dtype=dtype2)
    except:
        fringe2 = np.fromfile(frfile, dtype=dtype2)

    frfile.close()

    # IF identifiers present in the binary SWIN rows.
    # Use the actual rows to infer whether the binary IF/fridx numbering is 0-based or 1-based.
    if_values_in_rows = set(int(x) for x in np.unique(fringe2["IF"])) if len(fringe2) > 0 else set()
    DIFX_IF_BASE = infer_difx_if_base(
        DOIF,
        if_values_in_rows,
        IF_OFF,
        DIFX_IF_BASE_REQUESTED,
    )
    shifted_ifs_present = []
    for IFp in DOIF:
        IFp = int(IFp)
        IF_main = IFp if DIFX_IF_BASE == 1 else IFp - 1
        IF_shift = IF_main + IF_OFF
        if IF_OFF != 0 and IF_shift in if_values_in_rows:
            shifted_ifs_present.append(IF_shift)
    shifted_ifs_present = sorted(set(shifted_ifs_present))
    warnings.append(
        "%s: SWIN reader uses output Nchan=%d; inferred DIFX_IF_BASE=%d; binary IF range=%s"
        % (scan_name, Nchan, DIFX_IF_BASE, ("%d..%d" % (min(if_values_in_rows), max(if_values_in_rows))) if if_values_in_rows else "EMPTY")
    )
    if shifted_ifs_present:
        warnings.append(
            "%s: binary DIFX rows contain shifted IFs (IF_OFFSET=%d): %s"
            % (scan_name, IF_OFF, ",".join(str(x) for x in shifted_ifs_present))
        )

    SNR_X = [[1.0e18, 0.0] for i in range(Nants)]
    SNR_Y = [[1.0e18, 0.0] for i in range(Nants)]

    any_data = False
    
    # collect extra scan metrics
    pol_present = defaultdict(lambda: {"XX": 0, "XY": 0, "YX": 0, "YY": 0})  # key=(baseline,IF)
    station_if_pol_issues = defaultdict(lambda: defaultdict(set))
    station_if_missing_products = defaultdict(lambda: defaultdict(list))
    empty_if_candidate_rows = []
    ant_pol_present = defaultdict(lambda: {"XX": 0, "XY": 0, "YX": 0, "YY": 0}) # Per-antenna pol presence
    seen_baselines = set()
    seen_baseline_pairs = set()
    n_total_rows = 0

    u_vals = []
    v_vals = []
    w_vals = []
    weight_vals = []

    def _snr_from_pol(vis2d):
        """
        vis2d: complex array of shape (ntime, nchan) or (nrow, nchan).
        Returns SNR float, or None if empty/unusable.
        """
        if vis2d is None or vis2d.size == 0:
            return None

        mat = np.fft.fftshift(np.abs(np.fft.fft2(vis2d)))
        peak = np.unravel_index(np.argmax(mat), np.shape(mat))
        snr = mat[peak[0], peak[1]]
        # zero small patch around peak (avoid biasing noise estimate)
        i0, j0 = peak
        mat[max(i0-1,0):min(i0+2,mat.shape[0]), max(j0-1,0):min(j0+2,mat.shape[1])] = 0.0

        denom = np.sqrt(np.var(mat) + np.average(mat) ** 2.0)
        if denom == 0.0 or not np.isfinite(denom):
            return None

        snr /= denom
        if not np.isfinite(snr):
            return None
        return float(snr)

    for bb in SNRBas:
        BSel = (bb[0] + 1) * 256 + (bb[1] + 1)

        X1M = 0.0
        X2M = 0.0
        Y1M = 0.0
        Y2M = 0.0

        for IFp in DOIF:
            IFp = int(IFp)
            # Logical IFs are user-facing 1-based.
            # Binary SWIN IF/fridx can be either 1-based or 0-based.
            IF_main = IFp if DIFX_IF_BASE == 1 else IFp - 1
            IF_candidates = [int(IF_main)]
            # Only test shifted IF if that physical IF actually exists in the binary rows.
            if IF_OFF != 0:
                IF_shift = int(IF_main + IF_OFF)
                if IF_shift in if_values_in_rows:
                    IF_candidates.append(IF_shift)
            # Remove duplicates in case IF_OFF=0 or something odd.
            IF_candidates = list(dict.fromkeys(IF_candidates))
            any_candidate_has_rows = False
            for IFc in IF_candidates:
                IFl = physical_if_to_logical_if_with_base(IFc, DOIF, IF_OFF, DIFX_IF_BASE)
                MASK = np.logical_and(fringe2["BAS"] == BSel, fringe2["IF"] == IFc)
                if np.sum(MASK) == 0:
                    continue
                any_candidate_has_rows = True

                P1 = np.array([chr(x) for x in fringe2["POL1"][MASK]])
                P2 = np.array([chr(x) for x in fringe2["POL2"][MASK]])

                XXp = np.logical_and(
                    np.logical_or(P1 == "R", P1 == "X"),
                    np.logical_or(P2 == "R", P2 == "X"),
                )
                XYp = np.logical_and(
                    np.logical_or(P1 == "R", P1 == "X"),
                    np.logical_or(P2 == "L", P2 == "Y"),
                )
                YXp = np.logical_and(
                    np.logical_or(P1 == "L", P1 == "Y"),
                    np.logical_or(P2 == "R", P2 == "X"),
                )
                YYp = np.logical_and(
                    np.logical_or(P1 == "L", P1 == "Y"),
                    np.logical_or(P2 == "L", P2 == "Y"),
                )
                missing = []
                if not np.any(XXp): missing.append("XX")
                if not np.any(XYp): missing.append("XY")
                if not np.any(YXp): missing.append("YX")
                if not np.any(YYp): missing.append("YY")
                ant1 = ANT_NAMES[bb[0]]
                ant2 = ANT_NAMES[bb[1]]
                if missing:
                    warnings.append(
                        "%s: baseline %s-%s IF %d missing %s" %
                        (scan_name, ant1, ant2, IFc, ",".join(missing))
                    )
                    station_issue = infer_station_pol_issues(missing)
                    for pol in station_issue[0]:
                        station_if_pol_issues[ant1][int(IFl)].add(pol)
                    for pol in station_issue[1]:
                        station_if_pol_issues[ant2][int(IFl)].add(pol)
                    station_if_missing_products[ant1][int(IFl)].extend(
                        ["%s-%s missing %s" % (ant1, ant2, ",".join(missing))]
                    )
                    station_if_missing_products[ant2][int(IFl)].extend(
                        ["%s-%s missing %s" % (ant1, ant2, ",".join(missing))]
                    )

                # -polarization completeness counts (by baseline+IF)
                key = (int(BSel), int(IFc))
                pol_present[key]["XX"] += int(np.any(XXp))
                pol_present[key]["XY"] += int(np.any(XYp))
                pol_present[key]["YX"] += int(np.any(YXp))
                pol_present[key]["YY"] += int(np.any(YYp))
                # Per-antenna presence:
                # if this baseline has a product at least once,
                # count it for BOTH antennas participating in the baseline.
                a1 = ANT_NAMES[bb[0]]
                a2 = ANT_NAMES[bb[1]]
                if np.any(XXp):
                    ant_pol_present[a1]["XX"] += 1
                    ant_pol_present[a2]["XX"] += 1
                if np.any(XYp):
                    ant_pol_present[a1]["XY"] += 1
                    ant_pol_present[a2]["XY"] += 1
                if np.any(YXp):
                    ant_pol_present[a1]["YX"] += 1
                    ant_pol_present[a2]["YX"] += 1
                if np.any(YYp):
                    ant_pol_present[a1]["YY"] += 1
                    ant_pol_present[a2]["YY"] += 1

                # u/v/w/weight coverage statistics
                try:
                    u_vals.extend(fringe2["U"][MASK].tolist())
                    v_vals.extend(fringe2["V"][MASK].tolist())
                    w_vals.extend(fringe2["W"][MASK].tolist())
                    weight_vals.extend(fringe2["WEIGHT"][MASK].tolist())
                except Exception:
                    pass

                # baseline diversity
                # mark baseline as present for this scan
                seen_baselines.add(int(BSel))
                seen_baseline_pairs.add(tuple(sorted([ant1, ant2])))
                n_total_rows += int(np.sum(MASK))

                SNR_XX = _snr_from_pol(fringe2["VISIB"][MASK, :][XXp, :-1])
                SNR_XY = _snr_from_pol(fringe2["VISIB"][MASK, :][XYp, :-1])
                SNR_YX = _snr_from_pol(fringe2["VISIB"][MASK, :][YXp, :-1])
                SNR_YY = _snr_from_pol(fringe2["VISIB"][MASK, :][YYp, :-1])
                if any(v is not None for v in (SNR_XX, SNR_XY, SNR_YX, SNR_YY)):
                    any_data = True
                cand_X1 = [v for v in (SNR_XX, SNR_XY) if v is not None]
                cand_X2 = [v for v in (SNR_XX, SNR_YX) if v is not None]
                cand_Y1 = [v for v in (SNR_YY, SNR_YX) if v is not None]
                cand_Y2 = [v for v in (SNR_YY, SNR_XY) if v is not None]
                if len(cand_X1) > 0: X1M = max(X1M, max(cand_X1))
                if len(cand_X2) > 0: X2M = max(X2M, max(cand_X2))
                if len(cand_Y1) > 0: Y1M = max(Y1M, max(cand_Y1))
                if len(cand_Y2) > 0: Y2M = max(Y2M, max(cand_Y2))

            if not any_candidate_has_rows:
                ant1 = ANT_NAMES[bb[0]]
                ant2 = ANT_NAMES[bb[1]]
                empty_if_candidate_rows.append((ant1, ant2, int(IFp)))
                warnings.append(
                    "%s: baseline %s-%s IF %d has no visibilities" %
                    (scan_name, ant1, ant2, int(IFp))
                )
        # After scanning all IF candidates for this baseline, update per-antenna
        if X1M > 0.0:
            SNR_X[bb[0]][0] = min(SNR_X[bb[0]][0], X1M)
            SNR_X[bb[0]][1] = max(SNR_X[bb[0]][1], X1M)
        if X2M > 0.0:
            SNR_X[bb[1]][0] = min(SNR_X[bb[1]][0], X2M)
            SNR_X[bb[1]][1] = max(SNR_X[bb[1]][1], X2M)

        if Y1M > 0.0:
            SNR_Y[bb[0]][0] = min(SNR_Y[bb[0]][0], Y1M)
            SNR_Y[bb[0]][1] = max(SNR_Y[bb[0]][1], Y1M)
        if Y2M > 0.0:
            SNR_Y[bb[1]][0] = min(SNR_Y[bb[1]][0], Y2M)
            SNR_Y[bb[1]][1] = max(SNR_Y[bb[1]][1], Y2M)

    del fringe2

    empty_partners = defaultdict(lambda: defaultdict(set))
    for ant1, ant2, IFp in empty_if_candidate_rows:
        empty_partners[ant1][IFp].add(ant2)
        empty_partners[ant2][IFp].add(ant1)

    for ant, ifdict in empty_partners.items():
        for IFp, partners in ifdict.items():
            if len(partners) >= 2:
                station_if_pol_issues[ant][int(IFp)].add("X")
                station_if_pol_issues[ant][int(IFp)].add("Y")
                warnings.append(
                    "%s: station %s IF %d has no rows with %d baselines; flagging both pol. channels" %
                    (scan_name, ant, int(IFp), len(partners))
                )

    # Summarize station/IF polarization issues.
    # If a station is inferred to miss both X and Y for an IF, record it explicitly.
    bad_if_raw = defaultdict(set)
    pol_issue_summary = {}
    for ant, ifdict in station_if_pol_issues.items():
        pol_issue_summary[ant] = {}
        for IFc, pols in ifdict.items():
            pols = set(pols)
            pol_issue_summary[ant][int(IFc)] = sorted(pols)
            if ("X" in pols) or ("Y" in pols):
                bad_if_raw[ant].add(int(IFc))
            if ("X" in pols) and ("Y" in pols):
                warnings.append("%s: station %s IF %d appears to miss BOTH X and Y products" % (scan_name, ant, int(IFc)) )

    bad_if_raw = {ant: sorted(vals) for ant, vals in bad_if_raw.items()}

    # If no usable data at all -> fatal
    if not any_data:
        return _fatal("No usable visibilities found for requested IFs (including IF_OFFSET).")

    # ---- Read sources from .calc ----
    try:
        with open(dd, "r") as IFF:
            clines = IFF.readlines()
    except Exception as e:
        return _fatal("Cannot read calc file %s" % dd, e)
    lsou = None
    Nsou = 0
    for li, line in enumerate(clines):
        if line.startswith("NUM SOURCES:"):
            lsou = li
            Nsou = int(line.split()[-1])
            break

    SCANSOUS = []
    if lsou is not None and Nsou > 0:
        for j in range(Nsou):
            try:
                Snam = clines[lsou + 1 + j * 5].split()[-1]
                SCANSOUS.append(Snam)
            except Exception:
                pass

    # Elevation & parallactic angle metrics
    im_path = "%s.im" % dd[:-5]
    el_stats = {}
    pa_stats = {}
    try:
        if os.path.exists(im_path):
            azel = parse_im_azel(im_path)
        else:
            azel = {}
            warnings.append("%s: .im file not found (%s) -> cannot compute elevation/PA metrics" % (scan_name, im_path))
        for ai in range(Nants):
            name = ANT_NAMES[ai]
            els = azel.get(ai, {}).get("el", [])
            azs = azel.get(ai, {}).get("az", [])
            if els:
                el_stats[name] = {
                    "min": float(np.min(els)),
                    "max": float(np.max(els)),
                    "median": float(np.median(els)),
                    "range": float(np.max(els) - np.min(els)),
                }
            else:
                warnings.append("%s: no EL data for antenna %s in .im" % (scan_name, name))
            # Parallactic angle requires coords (lat) + az+el series
            if name not in coords_by_name:
                # warn once per scan, not per time sample
                continue

            if (azs is None) or (els is None) or (len(azs) == 0) or (len(els) == 0):
                continue

            try:
                x, y, z = coords_by_name[name]
                lat_rad, _lon_rad = ecef_to_geodetic_lat_lon(x, y, z)
                # compute PA for each sample we have (pair az/el by index)
                n = min(len(azs), len(els))
                pa_list = []
                for k in range(n):
                    az_rad = math.radians(float(azs[k]))
                    el_rad = math.radians(float(els[k]))
                    q = parallactic_angle_from_az_el_lat(az_rad, el_rad, lat_rad)
                    pa_list.append(math.degrees(q))
                # unwrap for coverage estimate
                if len(pa_list) >= 1:
                    pa_wrapped = ((np.array(pa_list) + 180.0) % 360.0) - 180.0
                    pa_stats[name] = {
                        "median": float(np.median(pa_wrapped)),
                        "n": int(len(pa_wrapped)),
                    }
            except Exception:
                pass
        # Coordinate warnings
        if coords_missing:
            warnings.append(
                "%s: missing antenna coordinates for PA (%d/%d): %s. "
                "Provide ANT_COORDS or ensure .input has station positions." %
                (scan_name, len(coords_missing), Nants, ",".join(coords_missing[:12]) + ("..." if len(coords_missing) > 12 else ""))
            )
        if len(pa_stats) == 0:
            warnings.append("%s: PA metrics unavailable (missing coords and/or AZ/EL in .im). pa_ok defaults to None." % scan_name)

    except Exception as e:
        warnings.append("%s: failed computing EL/PA metrics (%s: %s)" % (scan_name, type(e).__name__, str(e)))

    IsGood = True
    SNROut = [scan_name, SCANSOUS]
    for j in range(Nants):
        Observes = "+"
        if SNR_X[j][1] < minSNR or SNR_Y[j][1] < minSNR:
            IsGood = False
            if SNR_X[j][1] == 0.00 and SNR_Y[j][1] == 0.0:
                Observes = "-"
        SNROut += [
            ANT_NAMES[j],
            SNR_X[j][0],
            SNR_X[j][1],
            SNR_Y[j][0],
            SNR_Y[j][1],
            Observes,
        ]
    SNROut += ["Y" if IsGood else "N"]

    # summarize pol completeness
    pol_keys = list(pol_present.keys())
    n_keys = len(pol_keys)
    n_full = 0
    n_has_cross = 0
    for k in pol_keys:
        d = pol_present[k]
        full = (d["XX"] > 0 and d["XY"] > 0 and d["YX"] > 0 and d["YY"] > 0)
        has_cross = (d["XY"] > 0 and d["YX"] > 0)
        n_full += int(full)
        n_has_cross += int(has_cross)

    pol_metrics = {
        "baseline_if_pairs": n_keys,
        "full_XX_XY_YX_YY_pairs": n_full,
        "crosshand_pairs": n_has_cross,
        "full_fraction": (float(n_full) / n_keys) if n_keys > 0 else 0.0,
        "crosshand_fraction": (float(n_has_cross) / n_keys) if n_keys > 0 else 0.0,
    }
    # Per-antenna pol summary
    pol_per_ant = {}
    for ai in range(Nants):
        name = ANT_NAMES[ai]
        d = ant_pol_present.get(name, {"XX": 0, "XY": 0, "YX": 0, "YY": 0})
        pol_per_ant[name] = {
            "XX": int(d["XX"]),
            "XY": int(d["XY"]),
            "YX": int(d["YX"]),
            "YY": int(d["YY"]),
            "has_crosshands": bool(d["XY"] > 0 and d["YX"] > 0),
            "has_any_crosshand": bool(d["XY"] > 0 or d["YX"] > 0),
            "has_full_pol": bool(d["XX"] > 0 and d["XY"] > 0 and d["YX"] > 0 and d["YY"] > 0),
        }
    # u/v/weight/baseline metrics
    uv_metrics = {}
    if len(u_vals) > 0 and len(v_vals) > 0:
        umin, umax = float(np.min(u_vals)), float(np.max(u_vals))
        vmin, vmax = float(np.min(v_vals)), float(np.max(v_vals))
        uv_span = float(np.sqrt((umax - umin) ** 2 + (vmax - vmin) ** 2))
        uv_metrics.update({
            "u_min": umin, "u_max": umax,
            "v_min": vmin, "v_max": vmax,
            "uv_span": uv_span,
        })
    wt_metrics = {}
    if len(weight_vals) > 0:
        wv = np.array(weight_vals, dtype=float)
        wt_metrics = {
            "weight_min": float(np.min(wv)),
            "weight_median": float(np.median(wv)),
            "weight_max": float(np.max(wv)),
            "weight_zero_fract": float(np.mean(wv <= 0.0)),
        }
    baseline_metrics = {
        "unique_baselines_seen": int(len(seen_baselines)),
        "total_rows_counted": int(n_total_rows),
        "pairs_seen": sorted(["%s-%s" % p for p in seen_baseline_pairs]),
    }
    # Elevation summary across antennas
    elev_summary = {}
    if el_stats:
        med_els = [el_stats[a]["median"] for a in el_stats if "median" in el_stats[a]]
        min_els = [el_stats[a]["min"] for a in el_stats if "min" in el_stats[a]]
        elev_summary = {
            "median_of_medians": float(np.median(med_els)) if med_els else None,
            "min_of_mins": float(np.min(min_els)) if min_els else None,
        }

    # elevation flag
    elev_thr = cfg.get("MIN_ELEV_DEG", None)
    if elev_thr is None:
        elev_ok = None
    else:
        elev_ok = (elev_summary.get("min_of_mins") is not None) and (elev_summary["min_of_mins"] >= float(elev_thr))
    # crosshands flag
    req_ch = cfg.get("REQUIRE_CROSSHANDS", None)
    if req_ch is None or req_ch is False:
        cross_ok = None
    else:
        cross_ok = (pol_metrics.get("crosshand_fraction", 0.0) > 0.0)
    # baseline diversity flag
    min_bl = cfg.get("MIN_UNIQUE_BASELINES", None)
    if min_bl is None:
        base_ok = None
    else:
        base_ok = (baseline_metrics["unique_baselines_seen"] >= int(min_bl))
    # uv coverage flag
    uv_ok = True
    if MIN_UV_SPAN is not None:
        uv_ok = ("uv_span" in uv_metrics) and (uv_metrics["uv_span"] >= float(MIN_UV_SPAN))
    # weight flag
    wzero_ok = True
    if MAX_WEIGHT_ZERO_FRAC is not None:
        wzero_ok = ("weight_zero_fract" in wt_metrics) and (wt_metrics["weight_zero_fract"] <= float(MAX_WEIGHT_ZERO_FRAC))
    wmed_ok = True
    if MIN_WEIGHT_MEDIAN is not None:
        wmed_ok = ("weight_median" in wt_metrics) and (wt_metrics["weight_median"] >= float(MIN_WEIGHT_MEDIAN))

    # simple pass/fail flags
    flags = {
        "elevation_ok": elev_ok,
        "crosshands_ok": cross_ok,
        "baseline_diversity_ok": base_ok,
        "uv_ok": uv_ok,
        "weights_zero_ok": wzero_ok,
        "weights_median_ok": wmed_ok,
    }
    extra_metrics = {
        "ants_in_input": list(ANT_NAMES),
        "elevation": {"per_ant": el_stats, "summary": elev_summary},
        "parallactic_angle": {"per_ant": pa_stats},
        "pol_completeness": {"summary": pol_metrics, "per_ant": pol_per_ant},
        "station_pol_issues": {
            "per_ant_if": pol_issue_summary,
            "bad_if_raw": bad_if_raw,
        },
        "uv": uv_metrics,
        "weights": wt_metrics,
        "baselines": baseline_metrics,
        "flags": flags,
    }

    return {
        "ok": True,
        "fatal": False,
        "scan": scan_name,
        "result": SNROut,
        "metrics": extra_metrics, # NEW
        "warnings": warnings,
        "error": None,
        "dd": dd,
    }



def getFringeSNR_safe(dd, cfg):
    try:
        return getFringeSNR(dd, cfg)
    except Exception as e:
        # Safety net: anything unexpected is fatal
        return {
            "ok": False,
            "fatal": True,
            "scan": os.path.basename(dd).split(".")[0],
            "result": None,
            "warnings": [],
            "error": "%s: %s" % (type(e).__name__, str(e)),
            "dd": dd,
        }


#######################################
def SOURCE_SCANNER(
    EXPNAME="",
    *,
    DIFX_DIR="",
    SNRCut=10.0,
    SCAN_IF=None,
    NIF=0,
    IF_OFFSET=0,
    MIN_ELEV_DEG=20.0,
    REQUIRE_CROSSHANDS=False,
    MIN_UNIQUE_BASELINES=5,
    MIN_UV_SPAN=1e6,
    MAX_WEIGHT_ZERO_FRAC=None,
    MIN_WEIGHT_MEDIAN=None,
    ANT_COORDS=None,
    MIN_TOTAL_PA_COVERAGE_DEG=20.0,
    PA_SCANS_REQUIRED=0,
    NIF_PER_BAND=8,
    DIFX_IF_BASE="auto",
    NCPU=1,
):
    """
    Reads all the scans in a SWIN directory and creates an ASCII file
    with information about the sources, participating antennas, and SNRs
    of the correlation products.
    
    If minimum SNR is equal or higher than SNRCut, scan is marked as good;
    antennas not reaching SNRCut are marked as bad for that scan.

    Optional scan-quality metrics can be set for further filtering:
    elevation,parallactic-angle, cross-polarization products,
    baseline/uv coverage, and weight statistics.
    These are combined into a FINAL PASS flag to mark scan as good/bad.

    Parameters
    ----------
    EXPNAME : str
        Experiment label used in output filenames (e.g., "h" -> SOURCES_h.txt).
    DIFX_DIR : str
        Directory containing DiFX products (.calc/.input/.im/.difx).
    SNRCut : float
        Minimum SNR threshold for the scan to pass the *core* SNR test.
    SCAN_IF : list[int] or None
        IF indices to test (1-based DiFX IF numbering).
        If empty/None, all IFs in [1..NIF] are used (requires NIF > 0).
    NIF : int
        NUmber of IFs in experiment; only used when SCAN_IF is empty.
    IF_OFFSET : int
        Optional IF shift to also test (useful when cross-hands are stored in a different IF block).
        For each IF in SCAN_IF, the code tests (IF+IF_OFFSET) if in range.
    MIN_ELEV_DEG : float or None
        Minimum allowed geometric elevation (deg).
        If None, elevation is not evaluated.
    REQUIRE_CROSSHANDS : bool
        If True, require that some XY/YX cross-hands are present (otherwise not evaluated).
    MIN_UNIQUE_BASELINES : int or None
        Minimum number of unique baselines seen in the scan.
        If None, not evaluated.
    MIN_UV_SPAN : float or None
        Min. uv-span threshold (same units as stored in DIFX U/V records).
        If None, not evaluated.
    MAX_WEIGHT_ZERO_FRAC : float or None
        Maximum allowed fraction of rows with WEIGHT <= 0.
        If None, not evaluated.
    ANT_COORDS : dict[str, tuple[float,float,float]] or None
        Optional station ECEF coordinates override: {ant_name: (X,Y,Z)} in meters.
        Used for parallactic-angle calculations.
        If None/empty, the code tries to parse station positions from the .input file (best effort).
    MIN_TOTAL_PA_COVERAGE_DEG : float or None
        Target PA range (deg) used to add extra scans for better coverage.
        Requires az/el from .im and station coordinates (via ANT_COORDS or .input).
        If None, PA coverage is not considered.
    MIN_WEIGHT_MEDIAN : float or None
        Minimum allowed median weight. If None, not evaluated.
    PA_SCANS_REQUIRED : int
        Minimum number of PA-usable scans to include per antenna when PA
        enrichment is enabled.
        Values <=1 disable extra PA-scan enrichment.
    NIF_PER_BAND : int
        Number of IFs per frequency band.
        Used to expand detected bad IFs to the whole band in BAD_IF.
    DIFX_IF_BASE : "auto", 0, or 1
        Mapping between logical IFs and binary SWIN IF/fridx values.
        "auto" chooses the base that best matches the SWIN IFs.
    NCPU : int
        Number of worker processes used to scan multiple *.calc scans in parallel.
    """
    global minSNR

    if ANT_COORDS is None:
        ANT_COORDS = {}

    if SCAN_IF is None:
        SCAN_IF = []

    NIF = int(NIF)
    IF_OFF = int(IF_OFFSET)

    if len(SCAN_IF) == 0:
        if NIF <= 0:
            raise ValueError("SOURCE_SCANNER: SCAN_IF is empty, so you must pass a valid NIF (>0).")
        DOIF = list(range(1, NIF + 1))
    else:
        DOIF = [int(a) for a in SCAN_IF]

    # Config dict passed to workers
    cfg = dict(
        MIN_ELEV_DEG=(None if MIN_ELEV_DEG is None else float(MIN_ELEV_DEG)),
        MIN_TOTAL_PA_COVERAGE_DEG=(None if MIN_TOTAL_PA_COVERAGE_DEG is None else float(MIN_TOTAL_PA_COVERAGE_DEG)),
        REQUIRE_CROSSHANDS=(None if REQUIRE_CROSSHANDS is None else bool(REQUIRE_CROSSHANDS)),
        MIN_UNIQUE_BASELINES=(None if MIN_UNIQUE_BASELINES is None else int(MIN_UNIQUE_BASELINES)),
        ANT_COORDS=dict(ANT_COORDS),
        MIN_UV_SPAN=(None if MIN_UV_SPAN is None else float(MIN_UV_SPAN)),
        MAX_WEIGHT_ZERO_FRAC=(None if MAX_WEIGHT_ZERO_FRAC is None else float(MAX_WEIGHT_ZERO_FRAC)),
        MIN_WEIGHT_MEDIAN=(None if MIN_WEIGHT_MEDIAN is None else float(MIN_WEIGHT_MEDIAN)),
        SNRCut=float(SNRCut),
        IF_OFFSET=IF_OFF,
        DOIF=DOIF,
        NIF=NIF,
        NIF_PER_BAND=int(NIF_PER_BAND),
        DIFX_IF_BASE=DIFX_IF_BASE,
    )

    EXP = EXPNAME
    DIRE = DIFX_DIR
    minSNR = float(SNRCut)

    # Read visibilities for SNR estimate:
    OUTPUT = open("SOURCES_%s.txt" % EXP, "w")
    calcs = sorted(
        [
            f[:-1]
            for f in os.popen(
                'find ./%s -maxdepth 1 -name "*.calc" -print' % (DIRE)
            )
        ]
    )

    with Pool(NCPU) as p:
        TORESULT = p.starmap(getFringeSNR_safe, [(calc, cfg) for calc in calcs])

    # Log file for warnings / failures
    LOG = open("SOURCE_SCANNER_FAILED.log", "a")

    # Ensure FAILED folder exists (only used for fatal scans)
    if not os.path.exists("FAILED"):
        os.system("mkdir FAILED")

    # Collect only successful scan result lists for further stats
    OK_RESULTS = []
    RESULT_LINES = []
    SCAN_INFO = []

    NFAILED = 0
    for r in TORESULT:
        # Write warnings (partial issues)
        if r.get("warnings"):
            for w in r["warnings"]:
                print("WARNING " + w, file=LOG)
        # Handle fatal
        if (not r.get("ok", False)) and r.get("fatal", False):
            print("FATAL %s: %s" % (r.get("scan","?"), r.get("error","unknown error")), file=LOG)
            # Move all files related to the scan
            ddpath = r.get("dd", "")
            if ddpath.endswith(".calc"):
                scan_glob = ddpath[:-4] + "*"  # replace ".calc" -> "*"
                os.system("mv -f %s FAILED/ 2>/dev/null" % scan_glob)
            NFAILED += 1
            continue

        # If non-fatal but ok=False (rare), log and skip
        if not r.get("ok", False):
            print("NONFATAL %s: %s" % (r.get("scan","?"), r.get("error","unknown error")), file=LOG)
            continue

        # Normal successful scan
        scan = r["result"]
        OK_RESULTS.append(scan)

        ScanName = scan[0]
        ScanSous = ",".join(scan[1])
        Nants_scan = (len(scan) - 3) // 6

        def _fmt_snr_pair(lo, hi):
            lo = float(lo)
            hi = float(hi)
            # Sentinel from getFringeSNR:
            # [1.0e18, 0.0] means this pol was never updated.
            if hi == 0.0 and lo > 1.0e17:
                return "[NO_DATA]"
            return "[%.2f  %.2f]" % (lo, hi)

        RESULT_LINES.append("%s: %s\n" % (ScanName, ScanSous))
        RESULT_LINES.append("  SNR metrics:")
        for ai in range(Nants_scan):
            ant = scan[2 + ai*6]
            x_txt = _fmt_snr_pair(scan[3 + ai*6], scan[4 + ai*6])
            y_txt = _fmt_snr_pair(scan[5 + ai*6], scan[6 + ai*6])
            obs = scan[7 + ai*6]
            RESULT_LINES.append(
                "    ANT %s: X = %-18s | Y = %-18s  %s  "
                % (ant, x_txt, y_txt, obs)
            )
        RESULT_LINES.append("    SNR PASS: %s\n" % scan[-1])

        m = r.get("metrics") or {}
        flags = m.get("flags", {}) if isinstance(m, dict) else {}
        if m:
            elsum = m.get("elevation", {}).get("summary", {})
            pol_block = m.get("pol_completeness", {})
            pol_sum = pol_block.get("summary", {})
            pol_ant = pol_block.get("per_ant", {})
            basm = m.get("baselines", {})
            uvm = m.get("uv", {})
            wtm = m.get("weights", {})

            # Helpers
            def _thr_str(v, fmt="{}"):
                return "None" if v is None else (fmt.format(v))

            def _eval_str(flag_value):
                # None => not evaluated
                if flag_value is None:
                    return "N/A"
                return "OK" if flag_value else "BAD"

            def _val_str(v, fmt="{:.3g}"):
                if v is None:
                    return "None"
                try:
                    return fmt.format(v)
                except Exception:
                    return str(v)

            def _status_to_bool(status):
                """
                Convert a status token to boolean or None.
                Accepted:
                  - "OK"  -> True
                  - "BAD" -> False
                  - "Y"   -> True
                  - "N"   -> False
                Everything else (None, "unavailable", etc.) -> None
                """
                if status is None:
                    return None
                s = str(status).strip().upper()
                if s == "OK" or s == "Y":
                    return True
                if s == "BAD" or s == "N":
                    return False
                return None

            def final_pass_from_tests(**tests):
                """
                tests: keyword args like snr_pass="Y", elevation="OK", etc
                Only considers tests that map to True/False; ignores None/unknown/unavailable.
                Returns "Y" or "N".
                """
                vals = []
                for _, v in tests.items():
                    b = _status_to_bool(v)
                    if b is not None:
                        vals.append(b)
                return "Y" if all(vals) else "N"

            RESULT_LINES.append("  EXTRA METRICS:")
            # Elevation
            emin = elsum.get("min_of_mins", None)
            emed = elsum.get("median_of_medians", None)
            elev_thr = cfg.get("MIN_ELEV_DEG", None)
            RESULT_LINES.append(
                "    elevation: min=%s deg | med=%s deg | thr=%s | %s" %
                (_val_str(emin, "{:.2f}"),
                 _val_str(emed, "{:.2f}"),
                 _thr_str(elev_thr, "{:.2f}"),
                 _eval_str(flags.get("elevation_ok", None)))
            )
            # Polarization completeness
            ch_frac = 100.0 * float(pol_sum.get("crosshand_fraction", 0.0) or 0.0)
            full_frac = 100.0 * float(pol_sum.get("full_fraction", 0.0) or 0.0)
            req_ch = cfg.get("REQUIRE_CROSSHANDS", None)
            if req_ch is None or req_ch is False:
                # Only show stats, don't label OK/BAD
                RESULT_LINES.append(
                    "    pol prod:  crosshands=%5.1f%% | full=%5.1f%% | (not required)" %
                    (ch_frac, full_frac)
                )
            else:
                RESULT_LINES.append(
                    "    pol prod:  crosshands=%5.1f%% | full=%5.1f%% | %s" %
                    (ch_frac, full_frac, _eval_str(flags.get("crosshands_ok", None)))
                )
            # per-antenna crosshand availability summary (compact)
            if pol_ant:
                missing_xyyx = [a for a, d in pol_ant.items() if not d.get("has_crosshands", False)]
                if len(missing_xyyx) == 0:
                    RESULT_LINES.append("    pol ant:   all ants have XY&YX")
                else:
                    RESULT_LINES.append("    pol ant:   missing XY&YX: %s" % ",".join(sorted(missing_xyyx)))
            # Baselines
            nbl = int(basm.get("unique_baselines_seen", 0) or 0)
            min_bl = cfg.get("MIN_UNIQUE_BASELINES", None)
            RESULT_LINES.append(
                "    baselines: %d unique | thr=%s | %s" %
                (nbl, _thr_str(min_bl, "{}"), _eval_str(flags.get("baseline_diversity_ok", None)))
            )
            # UV coverage
            if "uv_span" in uvm:
                uv_span = uvm.get("uv_span", None)
                uv_thr = cfg.get("MIN_UV_SPAN", None)
                if uv_thr is None:
                    RESULT_LINES.append("    uv span:   %s" % _val_str(uv_span, "{:.3e}"))
                else:
                    RESULT_LINES.append(
                        "    uv span:   %s | thr=%s | %s" %
                        (_val_str(uv_span, "{:.3e}"),
                         _val_str(uv_thr, "{:.3e}"),
                         _eval_str(flags.get("uv_ok", None)))
                    )
            else:
                RESULT_LINES.append("    uv span:   unavailable")
            # Weights
            if wtm:
                med = wtm.get("weight_median", None)
                zero_frac = wtm.get("weight_zero_fract", None)
                zero_pct = None if zero_frac is None else (100.0 * float(zero_frac))
                # median threshold
                parts = ["med=%s" % _val_str(med, "{:.3g}")]
                wmed_thr = cfg.get("MIN_WEIGHT_MEDIAN", None)
                if wmed_thr is not None:
                    parts.append("thr=%s" % _val_str(wmed_thr, "{:.3g}"))
                    parts.append(_eval_str(flags.get("weights_median_ok", None)))
                # zero fraction threshold
                parts.append("N_flag=%s%%" % (_val_str(zero_pct, "{:.1f}") if zero_pct is not None else "None"))
                wzero_thr = cfg.get("MAX_WEIGHT_ZERO_FRAC", None)
                if wzero_thr is not None:
                    parts.append("thr=%s%%" % _val_str(100.0 * float(wzero_thr), "{:.1f}"))
                    parts.append(_eval_str(flags.get("weights_zero_ok", None)))
                RESULT_LINES.append("    weights:   " + " | ".join(parts))
            # FINAL PASS: combine all evaluated tests (ignore N/A)
            snr_pass = scan[-1]  # always present: "Y" or "N"
            elev_status = None
            if flags.get("elevation_ok", None) is not None:
                elev_status = "OK" if flags.get("elevation_ok") else "BAD"
            polprod_status = None
            req_ch = cfg.get("REQUIRE_CROSSHANDS", None)
            if req_ch is not None and req_ch is not False:
                if flags.get("crosshands_ok", None) is not None:
                    polprod_status = "OK" if flags.get("crosshands_ok") else "BAD"
            polant_status = None
            if pol_ant:
                missing_xyyx = [a for a, d in pol_ant.items() if not d.get("has_crosshands", False)]
                polant_status = "OK" if len(missing_xyyx) == 0 else "BAD"
            bl_status = None
            if flags.get("baseline_diversity_ok", None) is not None:
                bl_status = "OK" if flags.get("baseline_diversity_ok") else "BAD"
            spi = m.get("station_pol_issues", {})
            bad_raw = spi.get("bad_if_raw", {})
            if bad_raw:
                txt = []
                for ant in sorted(bad_raw):
                    txt.append("%s:%s" % (ant, ",".join(str(x) for x in bad_raw[ant])))
                RESULT_LINES.append("    bad IF raw: " + " | ".join(txt))
            uv_status = None
            uv_thr = cfg.get("MIN_UV_SPAN", None)
            if uv_thr is not None:
                if "uv_span" in uvm and flags.get("uv_ok", None) is not None:
                    uv_status = "OK" if flags.get("uv_ok") else "BAD"
                else:
                    uv_status = None
            w_status = None
            wmed_thr = cfg.get("MIN_WEIGHT_MEDIAN", None)
            wzero_thr = cfg.get("MAX_WEIGHT_ZERO_FRAC", None)
            w_parts = []
            if wmed_thr is not None and flags.get("weights_median_ok", None) is not None:
                w_parts.append(bool(flags.get("weights_median_ok")))
            if wzero_thr is not None and flags.get("weights_zero_ok", None) is not None:
                w_parts.append(bool(flags.get("weights_zero_ok")))
            if len(w_parts) > 0:
                w_status = "OK" if all(w_parts) else "BAD"
            final_pass = final_pass_from_tests(
                snr_pass=snr_pass,
                elevation=elev_status,
                pol_prod=polprod_status,
                pol_ant=polant_status,
                baselines=bl_status,
                uv_span=uv_status,
                weights=w_status,
            )
            RESULT_LINES.append("\n  FINAL PASS: %s\n" % (final_pass,))
        else:
            # metrics missing -> final_pass falls back to pure SNR pass from getFringeSNR output
            final_pass = scan[-1]   # "Y"/"N" from SNR-only pass
            m = {}                  # normalize

        ant_snr = {}
        ant_present = set()
        for ai in range(Nants_scan):
            ant = scan[2 + ai*6]
            snr_xmax = float(scan[4 + ai*6])
            snr_ymax = float(scan[6 + ai*6])
            obs = scan[7 + ai*6]  # "+" or "-"
            if obs != "-":
                ant_present.add(ant)
            # "worst pol" snr for that antenna in this scan:
            ant_snr[ant] = min(snr_xmax, snr_ymax)

        SCAN_INFO.append({
            "name": ScanName,
            "code": ScanName.split('_')[1] if '_' in ScanName else ScanName,
            "scan_list": scan,           # original list
            "metrics": m,                # your extra metrics dict
            "final_pass": final_pass,    # "Y" or "N"
            "ants": ant_present,         # set of ants observing 
            "ant_snr": ant_snr,          # per-ant worst-pol SNR
        })

    LOG.close()

    # Write the normal text output
    print("\n".join(RESULT_LINES), file=OUTPUT)
    if len(OK_RESULTS) == 0:
        # Nothing to do; consider this a failure of the whole step
        OFF = open("SOURCE_SCANNER.FAILED", "w")
        print("All scans failed fatally. See SOURCE_SCANNER_FAILED.log", file=OFF)
        OFF.close()
        OUTPUT.close()
        return

    ALL_ANTENNAS = np.unique(np.concatenate([result[2:-1:6] for result in OK_RESULTS]))
    
    # Remove FAILED folder if it is empty (and we didn't move anything)
    if NFAILED == 0 and os.path.exists("FAILED"):
        try:
            if len(os.listdir("FAILED")) == 0:
                os.rmdir("FAILED")
        except Exception:
            pass
    if os.path.exists("SOURCE_SCANNER.FAILED"):
        os.system("rm -rf SOURCE_SCANNER.FAILED")
        
    # --------------------------------------------------------------------
    # BEST_SCAN selection:
    #  - guarantee each antenna appears at least once
    #  - prefer FINAL_PASS == "Y" scans
    #  - fallback to FINAL_PASS == "N" only if needed for coverage
    #  - optional PA enrichment: add more scans per antenna if requested
    # --------------------------------------------------------------------
    # Defensive: SCAN_INFO must exist (filled in the main loop above)
    if len(SCAN_INFO) == 0:
        OUTPUT.close()
        raise RuntimeError("SOURCE_SCANNER: SCAN_INFO is empty; cannot build BEST_SCAN selection.")

    def _scan_has_pa(s):
        " Helper: per-scan PA availability check (any antenna has PA stats in this scan) "
        try:
            per_ant = s.get("metrics", {}).get("parallactic_angle", {}).get("per_ant", {})
            return isinstance(per_ant, dict) and (len(per_ant) > 0)
        except Exception:
            return False

    def _scan_has_pa_for_ant(s, ant):
        " Helper: does scan have PA stats for this antenna? "
        try:
            per_ant = s.get("metrics", {}).get("parallactic_angle", {}).get("per_ant", {})
            return (ant in per_ant) and ("median" in per_ant[ant])
        except Exception:
            return False

    def circular_coverage_deg(angles_deg, period=360.0):
        """
            Return circular span (deg) covered by a set of angles
            by subtracting the largest gap on a periodic circle.
        """
        if not angles_deg:
            return 0.0
        a = np.sort(np.mod(np.array(angles_deg, dtype=float), period))
        gaps = np.diff(a)
        wrap_gap = (a[0] + period) - a[-1]
        max_gap = max(np.max(gaps) if len(gaps) else 0.0, wrap_gap)
        return float(period - max_gap)

    def _pa_coverage_for_ant(scan_info, calib_set, ant):
        """
            Compute total PA coverage (deg) for one antenna
            across the currently selected scans using per-scan PA medians.
        """
        angles = []
        for s in scan_info:
            if s.get("code") in calib_set and _scan_has_pa_for_ant(s, ant):
                angles.append(s["metrics"]["parallactic_angle"]["per_ant"][ant]["median"])
        return circular_coverage_deg(angles)

    def _scan_polcal_eligible_for_ant(s, ant):
        try:
            bad_raw = (
                s.get("metrics", {})
                 .get("station_pol_issues", {})
                 .get("bad_if_raw", {})
            )
            return ant not in bad_raw
        except Exception:
            return True

    # Core selection: best scan per antenna by ant_snr[ant], with FINAL_PASS preference
    BEST_SCAN = {}           # ant -> dict(name, code, ant_snr, final_pass, used_fallback)
    FALLBACK_ANTS = []       # antennas that required PASS=N fallback
    for ant in ALL_ANTENNAS:
        # scans where this antenna participates
        cand_all = [s for s in SCAN_INFO if (ant in s.get("ants", set())) and (ant in s.get("ant_snr", {}))]
        if len(cand_all) == 0:
            continue
        cand_pass = [
            s for s in cand_all
            if s.get("final_pass") == "Y"
            and _scan_polcal_eligible_for_ant(s, ant)
        ]
        cand_pass_ineligible = [
            s for s in cand_all
            if s.get("final_pass") == "Y"
            and not _scan_polcal_eligible_for_ant(s, ant)
        ]
        if len(cand_pass) > 0:
            best = max(cand_pass, key=lambda s: float(s["ant_snr"].get(ant, -1e18)))
            used_fallback = False
        elif len(cand_pass_ineligible) > 0:
            best = max(cand_pass_ineligible, key=lambda s: float(s["ant_snr"].get(ant, -1e18)))
            used_fallback = True
            FALLBACK_ANTS.append(ant)
        else:
            best = max(cand_all, key=lambda s: float(s["ant_snr"].get(ant, -1e18)))
            used_fallback = True
            FALLBACK_ANTS.append(ant)
        BEST_SCAN[ant] = {
            "name": best["name"],
            "code": best["code"],
            "ant_snr": float(best["ant_snr"].get(ant, np.nan)),
            "final_pass": best.get("final_pass", "N"),
            "used_fallback": bool(used_fallback),
        }
    # Initial CALIB_SCANS: union of BEST_SCAN scan codes
    CALIB_SCANS = sorted({v["code"] for v in BEST_SCAN.values() if v.get("code")})

    # Optional PA enrichment:
    # - only if user asks PA_SCANS_REQUIRED > 1
    # - only if PA exists at all in the dataset
    # - tries to ensure each antenna has at least K PA-usable scans in CALIB_SCANS
    MIN_TOTAL_PA_COVERAGE_DEG = float(MIN_TOTAL_PA_COVERAGE_DEG) if MIN_TOTAL_PA_COVERAGE_DEG is not None else None
    K = int(PA_SCANS_REQUIRED) if PA_SCANS_REQUIRED is not None else 1
    pa_exists_any = any(_scan_has_pa(s) for s in SCAN_INFO)
    if (K > 1) and pa_exists_any:
        calib_set = set(CALIB_SCANS) # filter for quick checks
        for ant in ALL_ANTENNAS:
            # candidate scans: FINAL_PASS==Y, has PA for ant
            cand = [s for s in SCAN_INFO
                    if (ant in s.get("ants", set()))
                    and s.get("final_pass") == "Y"
                    and _scan_polcal_eligible_for_ant(s, ant)
                    and _scan_has_pa_for_ant(s, ant)
                    and (s.get("code") not in calib_set)]
            if not cand:
                continue
            # enforce min. num. of PA-OK scans if requested
            K = int(PA_SCANS_REQUIRED) if PA_SCANS_REQUIRED is not None else 1
            if K < 1: K = 1
            current_k = sum(1 for s in SCAN_INFO
                            if (s.get("code") in calib_set)
                            and _scan_has_pa_for_ant(s, ant))
            # enforce min. total PA coverage if requested
            current_cov = _pa_coverage_for_ant(SCAN_INFO, calib_set, ant)
            need_cov = None if MIN_TOTAL_PA_COVERAGE_DEG is None else max(0.0, MIN_TOTAL_PA_COVERAGE_DEG - current_cov)

            # keep adding scans while we fail either constraint
            while (current_k < K) or (need_cov is not None and need_cov > 0.0):
                best_s = None
                best_gain = -1.0
                best_snr = -1e18
                base_cov = _pa_coverage_for_ant(SCAN_INFO, calib_set, ant)
                for s in cand:
                    code = s.get("code")
                    new  = _pa_coverage_for_ant(SCAN_INFO, calib_set | {code}, ant)
                    gain = new - base_cov
                    snr = float(s["ant_snr"].get(ant, -1e18))
                    if (gain > best_gain) or (gain == best_gain and snr > best_snr):
                        best_gain = gain
                        best_snr = snr
                        best_s = s

                if best_s is None:
                    break
                if (need_cov is not None and need_cov > 0.0) and best_gain <= 0.0 and current_k >= K:
                    break

                calib_set.add(best_s["code"])
                cand = [s for s in cand if s.get("code") != best_s.get("code")]

                current_k += 1
                current_cov = _pa_coverage_for_ant(SCAN_INFO, calib_set, ant)
                need_cov = None if MIN_TOTAL_PA_COVERAGE_DEG is None else max(0.0, MIN_TOTAL_PA_COVERAGE_DEG - current_cov)
                if not cand:
                    break

        CALIB_SCANS = sorted(calib_set)

    # Optional: print fallback summary
    if len(FALLBACK_ANTS) > 0:
        print("WARNING! These antennas required non-ideal POLCAL fallback for coverage: %s" %
              ",".join(sorted(set(FALLBACK_ANTS))))

    OUTPUT.close()

    # BAD_ANTENNAS: antennas that NEVER reach SNRCut in ANY scan
    BAD_ANTENNAS = sorted(
        ant for ant, info in BEST_SCAN.items()
        if float(info.get("ant_snr", -1e18)) < float(SNRCut)
    )

    # BAD_IF construction with transient-scan handling:
    #   - If station/IF fails in majority of scans where station appears,
    #     expand it to BAD_IF as before.
    #   - If fails only in a minority of scans, do NOT put it in BAD_IF.
    #     Instead, avoid those scans unless needed as fallback.
    #   - Write final text report listing persistent + transient failures.
    BAD_IF_MAJORITY_FRACTION = 0.50  # >50% of station-participating scans => persistent BAD_IF
    # ant -> set(scan codes where antenna participates)
    ant_scan_codes = defaultdict(set)
    # ant -> raw IF -> set(scan codes where that IF had inferred station-pol issues)
    bad_if_scan_hits = defaultdict(lambda: defaultdict(set))
    for s in SCAN_INFO:
        code = s.get("code")
        if not code:
            continue
        metrics = s.get("metrics", {})
        spi = metrics.get("station_pol_issues", {})
        bad_raw = spi.get("bad_if_raw", {})
        denom_ants = set(metrics.get("ants_in_input", []))
        if len(denom_ants) == 0:
            denom_ants = set(s.get("ants", set()))
        # Ensure anything counted as bad is also counted as participating
        denom_ants |= set(bad_raw.keys())
        for ant in denom_ants:
            ant_scan_codes[ant].add(code)
        for ant, ifs in bad_raw.items():
            for IFc in ifs:
                bad_if_scan_hits[ant][int(IFc)].add(code)

    BAD_IF = defaultdict(set)
    TRANSIENT_BAD_IF = defaultdict(dict)
    PERSISTENT_BAD_IF = defaultdict(dict)
    nif_total = cfg.get("NIF", None)
    nif_per_band = int(cfg.get("NIF_PER_BAND", 8))
    for ant in sorted(bad_if_scan_hits):
        n_total = len(ant_scan_codes.get(ant, set()))
        if n_total <= 0:
            continue
        for IFc in sorted(bad_if_scan_hits[ant]):
            scans_bad = sorted(bad_if_scan_hits[ant][IFc])
            n_bad = len(scans_bad)
            frac_bad = float(n_bad) / float(n_total)

            is_persistent = frac_bad > BAD_IF_MAJORITY_FRACTION
            if is_persistent:
                PERSISTENT_BAD_IF[ant][IFc] = {
                    "n_bad": n_bad,
                    "n_total": n_total,
                    "frac_bad": frac_bad,
                }
                can_expand = (
                    nif_total is not None
                    and nif_per_band > 0
                    and 1 <= int(IFc) <= int(nif_total)
                )
                if can_expand:
                    for IFb in if_to_band_range(
                        IFc,
                        nif_per_band=nif_per_band,
                        nif_total=nif_total,
                    ):
                        BAD_IF[ant].add(int(IFb))
                else:
                    BAD_IF[ant].add(int(IFc))
            else:
                TRANSIENT_BAD_IF[ant][IFc] = {
                    "n_bad": n_bad,
                    "n_total": n_total,
                    "frac_bad": frac_bad,
                    "scans": scans_bad,
                }

    BAD_IF = {ant: sorted(vals) for ant, vals in BAD_IF.items()}

    # Final summary report.
    with open("%s_BadIF_Report.txt" % EXPNAME, "w") as f:
        f.write("BAD IF / transient scan report\n")
        f.write("Persistent threshold: > %.1f%% of scans where antenna participates\n\n" % (100.0 * BAD_IF_MAJORITY_FRACTION))
        f.write("PERSISTENT FAILURES -> included in BAD_IF\n")
        if PERSISTENT_BAD_IF:
            for ant in sorted(PERSISTENT_BAD_IF):
                for IFc in sorted(PERSISTENT_BAD_IF[ant]):
                    d = PERSISTENT_BAD_IF[ant][IFc]
                    f.write("%s IF %d: failed %d/%d scans = %.1f%%\n" % (
                            ant, IFc, d["n_bad"], d["n_total"],
                            100.0 * d["frac_bad"], )
                    )
        else:
            f.write("NONE\n")

        f.write("\nTRANSIENT FAILURES -> not included in BAD_IF; avoided as POLCAL candidates unless needed as fallback\n")
        if TRANSIENT_BAD_IF:
            for ant in sorted(TRANSIENT_BAD_IF):
                for IFc in sorted(TRANSIENT_BAD_IF[ant]):
                    d = TRANSIENT_BAD_IF[ant][IFc]
                    f.write("%s IF %d: failed %d/%d scans = %.1f%%; scans=%s\n" % (ant, IFc, d["n_bad"], d["n_total"], 100.0 * d["frac_bad"], ",".join(d["scans"]), ) )
        else:
            f.write("NONE\n")

        f.write("\nFINAL BAD_IF\n")
        if BAD_IF:
            for ant in sorted(BAD_IF):
                f.write("%s: %s\n" % (
                    ant,
                    ",".join(str(x) for x in BAD_IF[ant])
                ))
        else:
            f.write("NONE\n")

    # Write binary ScanStats for the master script:
    #   BEST_SCAN      : dict per antenna
    #   CALIB_SCANS    : list of scan codes used for calibration
    #   FALLBACK_ANTS  : antennas that needed FINAL_PASS="N" fallback
    #   BAD_ANTENNAS   : antennas that never reach SNRCut in any scan
    STATS = open("%s_ScanStats.dat" % EXPNAME, "wb")
    pk.dump([BEST_SCAN, CALIB_SCANS, sorted(set(FALLBACK_ANTS)), BAD_ANTENNAS, BAD_IF, list(ALL_ANTENNAS), SCAN_INFO], STATS)
    STATS.close()

    def write_scanstats_txt(expname, best_scan, calib_scans, fallback_ants, bad_antennas, bad_if, scan_info):
        """
        Write BEST_SCAN and CALIB_SCANS, achieved PA coverage per antenna
        across selected scans, fallback and globally-bad antenna lists.
        """
        by_code = {s["code"]: s for s in scan_info if "code" in s}
        calib_set = set(calib_scans)

        def _pa_entries_for_ant_in_set(ant):
            entries = []
            for code in sorted(calib_set):
                s = by_code.get(code)
                if not s:
                    continue
                try:
                    d = s["metrics"]["parallactic_angle"]["per_ant"][ant]
                    if "median" in d:
                        entries.append((code, float(d["median"])))
                except Exception:
                    pass
            return entries

        def _pa_medians_for_ant_in_set(ant):
            return [pa for _code, pa in _pa_entries_for_ant_in_set(ant)]

        with open("%s_ScanStats.txt" % expname, "w") as f:
            f.write("BEST_SCAN (per antenna)\n")
            for ant in sorted(best_scan):
                b = best_scan[ant]
                code = b.get("code", "")
                snr  = b.get("ant_snr", None)
                fp   = b.get("final_pass", "N")
                fb   = b.get("used_fallback", False)

                pa_entries = _pa_entries_for_ant_in_set(ant)
                angles = [pa for _scan_code, pa in pa_entries]
                if len(angles) == 0:
                    pa_cov_str = "NA"
                else:
                    pa_cov_str = "%.1f deg" % (circular_coverage_deg(angles),)
                pa_scan_codes = [scan_code for scan_code, _pa in pa_entries]
                pa_scans_str = ",".join(pa_scan_codes) if pa_scan_codes else "NONE"
                f.write(
                    "%s  scan=%s  snr=%s  FINAL_PASS=%s  fallback=%s \n"
                    "    PAcov=%s;  NscansPA=%d;  PA_scans=%s\n"
                    % (ant, code, snr, fp, fb, pa_cov_str, len(angles), pa_scans_str)
                )

            f.write("\nCALIB_SCANS (codes): ")
            for c in calib_scans:
                f.write("%s, " % c)
            f.write("\n")

            f.write("\nFALLBACK_ANTS\n")
            f.write(",".join(sorted(set(fallback_ants))) + "\n")

            f.write("\nBAD_ANTENNAS (never reach SNRCut in any scan)\n")
            if bad_antennas:
                f.write(",".join(sorted(set(bad_antennas))) + "\n")
            else:
                f.write("NONE\n")

            f.write("\nBAD_IFS\n")
            if bad_if:
                for ant in sorted(bad_if):
                    f.write("%s: %s\n" % (ant, ",".join(str(x) for x in bad_if[ant])))
            else:
                f.write("NONE\n")

    write_scanstats_txt(EXPNAME, BEST_SCAN, CALIB_SCANS, FALLBACK_ANTS, BAD_ANTENNAS, BAD_IF, SCAN_INFO)


if __name__ == "__main__":
    DIRE = "swin"
    EXP = "h"
    Nchan = 128
    NIF = 32
    minSNR = 10.0
    DOIF = []
    IF_OFF = 32
    SCAN_IF = [4,12,20,28]

    SNRBas = [[0, 1], [0, 2], [1, 2]]
    SNRCut = 50.0
    
    # --- EXTRA METRICS PARAMETERS (NEW) ---
    MIN_ELEV_DEG = 20.0          # elevation threshold for "good"
    MIN_TOTAL_PA_COVERAGE_DEG = 0.1      # minimum parallactic angle range you want (per antenna)
    REQUIRE_CROSSHANDS = True    # require some XY/YX presence
    MIN_UNIQUE_BASELINES = 5     # minimum unique cross baselines present in data

    # Coordinates override: antenna name -> (X,Y,Z) in meters (ECEF)
    ANT_COORDS = {
        "GS": (1130729.877, -4831245.972, 3994228.300),   # GGAO12M
        "HV": (5085502.673, 2668368.171, -2768493.227),   # HARTVGS
        "HB": (-3949991.094, 2522421.259, -4311707.721),  # HOBART12
        "IS": (-3959636.203, 3296825.448, 3747042.571),   # ISHIOKA
        "KE": (-4147354.925, 4581542.266, -1573302.777),  # KATH12M
        "K2": (-5543831.745, -2054585.590, 2387828.974),  # KOKEE12M
        "MG": (-1330788.462, -5328106.593, 3236427.492),  # MACGO12M
        "NN": (1200992.670, 252098.630, 6238038.646),     # NYALE13N
        "OE": (3370889.298, 711571.199, 5349692.048),     # ONSA13NE
        "OW": (3370946.779, 711534.507, 5349660.925),     # ONSA13SW
        "SA": (4618524.302, -2166020.720, 3816270.345),   # RAEGSMAR
        "YJ": (4848831.021, -261629.388, 4122976.576),    # RAEGYEB
        "S6": (-2831646.947, 4675729.541, 3275365.102),   # SESHAN13
        "T1": (-2826801.860, 4679253.933, 3274516.038),   # TIANMA13
        "UM": (228671.597, 4631855.932, 4367130.506),     # URUMQI13
        "WF": (1492206.223, -4458130.552, 4296015.629),   # WESTFORD
        "WN": (4075627.541, 931774.394, 4801552.451),     # WETTZ13N
        "WS": (4075658.807, 931824.883, 4801516.273),     # WETTZ13S
        "YG": (-2388896.500, 5043350.051, -3078590.462),  # YARRA12M
    }
    #### TODO: WHEN IMPLEMENTING IN MASTER SCRIPT, you can find them in the .vex file on https://cddis.nasa.gov/archive/vlbi/ivsdata/aux/

    # Example test run:
    SOURCE_SCANNER(
        EXPNAME=EXP,
        DIFX_DIR=DIRE,
        SNRCut=minSNR,
        SCAN_IF=SCAN_IF,  # empty -> use all IFs from NIF
        NIF=NIF,          # must be >0 if SCAN_IF is empty
        IF_OFFSET=IF_OFF,
        MIN_ELEV_DEG=MIN_ELEV_DEG,
        REQUIRE_CROSSHANDS=REQUIRE_CROSSHANDS,
        MIN_UNIQUE_BASELINES=MIN_UNIQUE_BASELINES,
        MIN_UV_SPAN=1e6,
        MAX_WEIGHT_ZERO_FRAC=None,
        MIN_WEIGHT_MEDIAN=None,
        ANT_COORDS=ANT_COORDS,
        MIN_TOTAL_PA_COVERAGE_DEG=MIN_TOTAL_PA_COVERAGE_DEG,
        PA_SCANS_REQUIRED=2,
        NIF_PER_BAND=8,
        NCPU=4,
    )


