#!/usr/bin/env python3
"""
INSPECT_PCAL_TONES.py  (standalone)

Count IF tones from DiFX PCAL_* files
and visualize what it would do to PCAL tone diff-phases (X-Y or R-L).

"""
import os, re, glob
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

SENTINEL = ("-1", "0", "0", "0")


# ========================================================================
# Helpers for robust stats, phase wrapping/unwrapping, and display shifting.
# ========================================================================
def gap_aware_unwrap(freq_hz, phase_rad, band_gap_hz=1e9):
    """
    Unwrap phase vs frequency,
    restart unwrapping across large frequency gaps (multi-band),
    then align blocks by 2*pi*n offsets.    
    """
    f = np.asarray(freq_hz, float)
    ph = np.asarray(phase_rad, float)
    if f.size == 0:
        return ph
    idx = np.argsort(f)
    f_sorted = f[idx]
    ph_sorted = ph[idx]
    cuts = [0]
    for i in range(1, len(f_sorted)):
        if (f_sorted[i] - f_sorted[i - 1]) > band_gap_hz:
            cuts.append(i)
    cuts.append(len(f_sorted))

    ph_u = ph_sorted.copy()
    for b0, b1 in zip(cuts[:-1], cuts[1:]):
        ph_u[b0:b1] = np.unwrap(ph_u[b0:b1])

    for bi in range(1, len(cuts) - 1):
        prev_end = cuts[bi] - 1
        next_start = cuts[bi]
        jump = ph_u[next_start] - ph_u[prev_end]
        n = np.round(jump / (2.0 * np.pi))
        ph_u[next_start:cuts[bi + 1]] -= n * 2.0 * np.pi

    out = np.empty_like(ph_u)
    out[idx] = ph_u
    return out

def unwrap_relative(phases_rad):
    """
    Simple sequential unwrap (assumes input already ordered)
    unwraps by applying 2*pi jumps: keeps diffs. within [-pi, +pi].
    """
    ph = np.array(phases_rad, dtype=float).copy()
    for i in range(len(ph)-1):
        while True:
            d = ph[i+1] - ph[i]
            if d > np.pi:
                ph[i+1:] -= 2*np.pi
            elif d < -np.pi:
                ph[i+1:] += 2*np.pi
            else:
                break
    return ph

def shift_deg_to_window(deg, lo=-360.0, hi=360.0):
    """ Shift whole IF phases by 360*n so they sit nicely in [lo,hi] limit """
    deg = np.asarray(deg, float)
    m = np.isfinite(deg)
    if m.sum() == 0:
        return deg
    center = 0.5*(lo+hi)
    best = None
    best_score = None
    for n in range(-4, 5):
        d = deg.copy()
        d[m] += 360.0*n
        inside = (d[m] >= lo) & (d[m] <= hi)
        score = (inside.size - inside.sum(), np.nanmean(np.abs(d[m] - center)))
        if best_score is None or score < best_score:
            best_score = score
            best = d
    return best


# ========================================================================
# DiFX .input helpers and IF mapping
# Parse IF from the .input file and map POLCAL IF indices to "unique IFs"
#       by nearest center frequency.
# ========================================================================
def build_unique_ifs_from_input(input_path, round_mhz=1e-6):
    """
    Parse DiFX .input lines (FREQ/BW/SIDEBAND).
    Then, group into a unique-IF list sorted by center frequency.
    Each entry has: IF num (1-based), f_low, f_high, f_center (in MHz)
    """
    freq_re  = re.compile(r"^FREQ \(MHZ\)\s+(\d+):\s+([0-9.+-Ee]+)")
    bw_re    = re.compile(r"^BW \(MHZ\)\s+(\d+):\s+([0-9.+-Ee]+)")
    sb_re    = re.compile(r"^SIDEBAND\s+(\d+):\s+([LU])")

    tmp = {}
    with open(input_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            m = freq_re.match(line)
            if m:
                i = int(m.group(1)); tmp.setdefault(i, {})["f0"] = float(m.group(2)); continue
            m = bw_re.match(line)
            if m:
                i = int(m.group(1)); tmp.setdefault(i, {})["bw"] = float(m.group(2)); continue
            m = sb_re.match(line)
            if m:
                i = int(m.group(1)); tmp.setdefault(i, {})["sb"] = m.group(2); continue

    def r(x):
        return round(x / round_mhz) * round_mhz

    key_to_indices = {}
    for fi, d in tmp.items():
        if not all(k in d for k in ("f0","bw","sb")):
            continue
        key = (r(d["f0"]), r(d["bw"]), d["sb"])
        key_to_indices.setdefault(key, []).append(fi)

    uif_list = []
    for key, idxs in key_to_indices.items():
        f0, bw, sb = key
        if sb == "U":
            lo, hi = f0, f0 + bw
        else:
            lo, hi = f0 - bw, f0
        lo, hi = min(lo,hi), max(lo,hi)
        cen = 0.5*(lo+hi)
        uif_list.append({
            "key": key,
            "indices": sorted(idxs),
            "f_low_mhz": lo,
            "f_high_mhz": hi,
            "f_center_mhz": cen,
        })

    uif_list.sort(key=lambda d: d["f_center_mhz"])
    for i, d in enumerate(uif_list, start=1):
        d["ifnum"] = i
    return uif_list


# ========================================================================
# PCAL file I/O and tone grouping
# Read DiFX PCAL_* files, average complex tones per freq. and pol.,
#       and group tone freqs. into unique-IF bands.
# ========================================================================
def split_blocks_by_sentinel(tokens):
    """
    Given a PCAL line, split the tone records into blocks.
    Returns list of blocks, each a list of (freq_MHz, pol, re, im).
    """
    blocks, cur = [], []
    data = tokens[6:]
    i, n = 0, len(data)
    while i < n:
        if i + 4 <= n and tuple(data[i:i+4]) == SENTINEL:
            if cur:
                blocks.append(cur); cur = []
            i += 4
            continue
        if i + 4 <= n:
            try:
                ff = float(data[i]); pol = data[i+1]
                rr = float(data[i+2]); ii = float(data[i+3])
                cur.append((ff, pol, rr, ii))
                i += 4
            except Exception:
                i += 1
                continue
        else:
            break
    if cur:
        blocks.append(cur)
    return blocks

def read_pcal_file_meanphasor(pcal_path):
    """ Parse PCAL_* file and compute per-(pol,freq) mean tone over all records """
    ant_name = None
    sum_complex = defaultdict(complex)
    sum_amp = defaultdict(float)
    count = defaultdict(int)

    with open(pcal_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                if "Telescope name" in line and "=" in line:
                    ant_name = line.split("=")[-1].strip()
                continue
            tokens = line.split()
            if len(tokens) < 10:
                continue
            if ant_name is None:
                ant_name = tokens[0]
            blocks = split_blocks_by_sentinel(tokens)
            for blk in blocks:
                for (f_mhz, pol, re_, im_) in blk:
                    z = complex(re_, im_)
                    key = (pol, f_mhz)
                    sum_complex[key] += z
                    sum_amp[key] += abs(z)
                    count[key] += 1

    tone_stats = defaultdict(dict)
    for (pol, f_mhz), n in count.items():
        zbar = sum_complex[(pol, f_mhz)] / max(n, 1)
        abar = sum_amp[(pol, f_mhz)] / max(n, 1)
        tone_stats[pol][f_mhz] = {
            "mean_complex": zbar,
            "mean_phase_rad": float(np.angle(zbar)),
            "mean_amp": float(abar),
            "n": int(n),
        }
    return ant_name, tone_stats

def choose_pol_pair(tone_stats):
    pols = sorted(tone_stats.keys())
    polset = set(pols)
    if polset.issuperset({"X","Y"}):
        return "X","Y"
    if polset.issuperset({"R","L"}):
        return "R","L"
    if len(pols) >= 2:
        return pols[0], pols[1]
    if len(pols) == 1:
        return pols[0], None
    return None, None

def group_tones_by_uif(freqs_mhz, uif_list, margin_mhz=1.0):
    freqs = sorted(freqs_mhz)
    groups = []
    assigned = set()
    for d in uif_list:
        lo = d["f_low_mhz"] - margin_mhz
        hi = d["f_high_mhz"] + margin_mhz
        tf = [f for f in freqs if (f >= lo and f <= hi)]
        if tf:
            groups.append((d["ifnum"], sorted(tf)))
            assigned.update(tf)
    unassigned = sorted([f for f in freqs if f not in assigned])
    return groups, unassigned


# ========================================================================
# Plotting
# Diagnostic figures to visualize tone phases and residuals per IF.
# ========================================================================
def plot_pcal_raw_multitone(
    ant, scan_tag, groups, statsA, statsB, polA, polB,
    uif_list, out_png,
    *, band_gap_hz=1e9,
    if_spacing=1.25, intra_if_width=0.75,
    phase_window_deg=180.0,
    use_amp_weight=True,
):
    """
    Multitone-style display of RAW PCAL tone phases.
    """

    def x_positions_for_if(if_idx_1based, freqs_mhz):
        f = np.asarray(freqs_mhz, float)
        if f.size <= 1:
            u = np.array([0.0])
        else:
            u = (f - f.min()) / (f.max() - f.min())
            u = np.clip(u, 0.0, 1.0)
        x0 = if_idx_1based * if_spacing
        return x0 + (u - 0.5) * intra_if_width

    fig = plt.figure(figsize=(11.5, 5.6))
    ax = fig.add_subplot(111)
    fig.suptitle("%s  PCAL %s-%s raw tone phases  [%s]" % (ant, polA, polB, scan_tag), fontsize=14)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.12, top=0.88)

    xticks = []
    xlabels = []

    # groups is [(ifnum, tone_freqs), ...] (skip any non-int labels)
    for ifnum, tone_freqs in groups:
        if not isinstance(ifnum, int):
            continue

        f_mhz = np.asarray(tone_freqs, float)
        if f_mhz.size < 2:
            continue

        phA = np.array([statsA.get(ff, {}).get("mean_phase_rad", np.nan) for ff in f_mhz], float)
        phB = np.array([statsB.get(ff, {}).get("mean_phase_rad", np.nan) for ff in f_mhz], float)
        aA  = np.array([statsA.get(ff, {}).get("mean_amp",       np.nan) for ff in f_mhz], float)
        aB  = np.array([statsB.get(ff, {}).get("mean_amp",       np.nan) for ff in f_mhz], float)

        ok = np.isfinite(phA) & np.isfinite(phB)
        if use_amp_weight:
            ok &= np.isfinite(aA) & np.isfinite(aB) & (aA > 0) & (aB > 0)
        if ok.sum() < 2:
            continue

        f_mhz = f_mhz[ok]
        phA   = phA[ok]
        phB   = phB[ok]
        if use_amp_weight:
            aA = aA[ok]
            aB = aB[ok]

        idx = np.argsort(f_mhz)
        f_mhz = f_mhz[idx]
        phA   = phA[idx]
        phB   = phB[idx]
        if use_amp_weight:
            aA = aA[idx]
            aB = aB[idx]

        f_hz = f_mhz * 1e6

        # raw diff-phase
        d = phA - phB

        # unwrap within IF using your band-gap logic
        d_u = gap_aware_unwrap(f_hz, d, band_gap_hz=band_gap_hz)

        # then unwrap relative in freq-order & constrain to a nice window
        d_u_rel = unwrap_relative(d_u)
        d_deg   = np.degrees(d_u_rel)
        d_deg   = shift_deg_to_window(d_deg, -phase_window_deg, phase_window_deg)

        # x positions in this IF slot
        x = x_positions_for_if(ifnum, f_mhz)

        ax.plot(x, d_deg, marker="o", linestyle="none", markersize=4, alpha=0.75)

        # boundary line after IF slot
        ax.axvline((ifnum * if_spacing) + 0.75*intra_if_width, linestyle=":", linewidth=1, alpha=0.5)

        xticks.append(ifnum * if_spacing)
        xlabels.append(str(ifnum))

    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)
        ax.set_xlim(min(xticks) - 0.75*intra_if_width, max(xticks) + 0.75*intra_if_width)

    ax.set_xlabel("IF number")
    ax.set_ylabel("Phase %s-%s (deg)" % (polA, polB))
    ax.set_ylim(-phase_window_deg, phase_window_deg)
    ax.grid(True, alpha=0.25)

    fig.savefig(out_png, dpi=180)
    plt.close(fig)



# -----------------------------
# To import as a module 
# -----------------------------
def INSPECT_PCAL_TONES(
    input_path,
    scan_difx_dir,
    MAX_NTONES_PER_IF=None,
    make_plots=False,
    plot_outdir=None,
    band_gap_hz=1e9,
    margin_mhz=1.0,
):
    """
    input_path: path to <EXPNAME>.input
    scan_difx_dir: e.g. "<CALDIR>/<EXPNAME>_<scan>.difx"
    MAX_NTONES_PER_IF: tone number threshold
    make_plots: create PNGs or not
    plot_outdir: where to put plots (if make_plots=True)

    Returns:
      results = {
        ant: {
          "pcal_path": str,
          "per_if_counts": {ifnum: n_tones, ...},
          "max_ntones": int,
          "all_below_threshold": bool,
          "plot_path": str or None,
        },
      }
    """
    os.makedirs(plot_outdir, exist_ok=True)
    scan_tag = os.path.basename(os.path.normpath(scan_difx_dir))

    # IF definition from .input
    uif_list = build_unique_ifs_from_input(input_path)

    # Find PCAL_* files
    pcal_files = sorted(
        p for p in glob.glob(os.path.join(scan_difx_dir, "PCAL_*"))
        if "_IF" not in os.path.basename(p)
    )

    results = {}

    for pcal_path in pcal_files:
        ant, tone_stats = read_pcal_file_meanphasor(pcal_path)
        if not ant:
            continue

        polA, polB = choose_pol_pair(tone_stats)
        if not polA or not polB:
            continue

        statsA = tone_stats.get(polA, {})
        statsB = tone_stats.get(polB, {})

        # Only freqs where both pols have tones
        freqs_common = sorted(set(statsA.keys()) & set(statsB.keys()))
        if not freqs_common:
            continue

        # Group by unique-IF using center-frequency overlap
        groups, _ = group_tones_by_uif(freqs_common, uif_list, margin_mhz=margin_mhz)

        per_if_counts = {int(ifnum): len(tone_freqs) for ifnum, tone_freqs in groups}
        max_ntones = max(per_if_counts.values()) if per_if_counts else 0
        all_below = (max_ntones <= MAX_NTONES_PER_IF)

        # Make plot of raw tones vs IF number
        plot_path = None
        if groups and make_plots:
            plot_path = os.path.join(plot_outdir, "pcal_raw_tones_%s_%s.png" % (ant, scan_tag))
            plot_pcal_raw_multitone(
                ant=ant,
                scan_tag=scan_tag,
                groups=groups,
                statsA=statsA,
                statsB=statsB,
                polA=polA,
                polB=polB,
                uif_list=uif_list,
                out_png=plot_path,
                band_gap_hz=band_gap_hz,
            )

        results[ant] = {
            "pcal_path": pcal_path,
            "per_if_counts": per_if_counts,
            "max_ntones": max_ntones,
            "all_below_threshold": all_below,
            "plot_path": plot_path,
        }

    return results




# -----------------------------
# Main 
# -----------------------------
def analyze_pcal_tones_and_plot(
    INPUT_PATH,
    SCAN_DIR,
    OUTDIR,
    MAX_NTONES_PER_IF=5,
    margin_mhz=1.0,
    band_gap_hz=1e9,
):
    """
    For each antenna in SCAN_DIR, using the DiFX INPUT_PATH, counts tones per IF and plots diff-phases (polA-polB) vs IF number

    Returns:
      results = {
        ant: {
          "pcal_path": str,
          "per_if_counts": {ifnum: n_tones, ...},
          "max_ntones": int,
          "all_below_threshold": bool,
          "plot_path": str or None,
        },
      }
    """

    os.makedirs(OUTDIR, exist_ok=True)
    scan_tag = os.path.basename(os.path.normpath(SCAN_DIR))

    # IF definition from .input
    uif_list = build_unique_ifs_from_input(INPUT_PATH)

    # Find PCAL_* files
    pcal_files = sorted(
        p for p in glob.glob(os.path.join(SCAN_DIR, "PCAL_*"))
        if "_IF" not in os.path.basename(p)
    )

    results = {}

    for pcal_path in pcal_files:
        ant, tone_stats = read_pcal_file_meanphasor(pcal_path)
        if not ant:
            continue

        polA, polB = choose_pol_pair(tone_stats)
        if not polA or not polB:
            continue

        statsA = tone_stats.get(polA, {})
        statsB = tone_stats.get(polB, {})

        # Only freqs where both pols have tones
        freqs_common = sorted(set(statsA.keys()) & set(statsB.keys()))
        if not freqs_common:
            continue

        # Group by unique-IF using center-frequency overlap
        groups, _ = group_tones_by_uif(freqs_common, uif_list, margin_mhz=margin_mhz)

        per_if_counts = {int(ifnum): len(tone_freqs) for ifnum, tone_freqs in groups}
        max_ntones = max(per_if_counts.values()) if per_if_counts else 0
        all_below = (max_ntones <= MAX_NTONES_PER_IF)

        # Make plot of raw tones vs IF number
        plot_path = None
        if groups:
            ant_dir = os.path.join(OUTDIR, ant)
            os.makedirs(ant_dir, exist_ok=True)
            plot_path = os.path.join(OUTDIR, "pcal_raw_tones_%s_%s.png" % (ant, scan_tag))

            plot_pcal_raw_multitone(
                ant=ant,
                scan_tag=scan_tag,
                groups=groups,
                statsA=statsA,
                statsB=statsB,
                polA=polA,
                polB=polB,
                uif_list=uif_list,
                out_png=plot_path,
                band_gap_hz=band_gap_hz,
            )

        results[ant] = {
            "pcal_path": pcal_path,
            "per_if_counts": per_if_counts,
            "max_ntones": max_ntones,
            "all_below_threshold": all_below,
            "plot_path": plot_path,
        }

    return results



if __name__ == "__main__":
    # Example usage (adapt paths as needed)
    INPUT_PATH = "swin/h_1186.input"
    SCAN_DIR   = "swin/h_1186.difx"
    OUTDIR     = "PCAL_TONE_PLOTS"
    MAX_NTONES_PER_IF = 5

    res = analyze_pcal_tones_and_plot(
        INPUT_PATH=INPUT_PATH,
        SCAN_DIR=SCAN_DIR,
        OUTDIR=OUTDIR,
        MAX_NTONES_PER_IF=MAX_NTONES_PER_IF,
    )

    print(f"MAX_NTONES_PER_IF = {MAX_NTONES_PER_IF}")
    print("\nAntennas with ALL IFs <= threshold:")
    for ant, d in sorted(res.items()):
        if d["all_below_threshold"]:
            print(f"  {ant:8s}  max_ntones={d['max_ntones']}  per_if={d['per_if_counts']}  plot={d['plot_path']}")

    print("\nAntennas with SOME IFs > threshold:")
    for ant, d in sorted(res.items()):
        if not d["all_below_threshold"]:
            print(f"  {ant:8s}  max_ntones={d['max_ntones']}  per_if={d['per_if_counts']}  plot={d['plot_path']}")

