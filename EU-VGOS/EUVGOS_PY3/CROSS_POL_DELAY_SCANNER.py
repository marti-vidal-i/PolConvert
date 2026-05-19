# CROSS_POL_DELAY_SCANNER: ESTIMATE A PRIORI XPOL_DELAYS FOR POLCONVERT 
#              Copyright (C) 2026  Ezequiel Albentosa-Ruiz
#              University of Valencia (Spain) & 
#              GFZ Helmholtz Centre for Geosciences (Germany)
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
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from typing import Optional

# -----------------------------
# I/O helpers
# -----------------------------
def load_data_file(file_path: str):
    """Load a pickle file and return the object."""
    try:
        with open(file_path, "rb") as f:
            return pickle.load(f)
    except Exception as e:
        print("Error loading file %s: %s" % (file_path, e))
        return None

# -----------------------------
# Phase unwrap + alignment (gap-aware)
# -----------------------------
def gap_aware_unwrap(freq_hz, phase_rad, BAND_GAP_HZ=1e9):
    """
    Unwrap phase in contiguous frequency blocks (split by large gaps),
    then align blocks by adding 2*pi*n offsets to minimize boundary jumps.

    Note: delay estimation uses exp(i*phase) and is tolerant of 2*pi wraps
          this unwrapping mainly stabilizes phi0 and the derived f_zero.
    """
    f = np.asarray(freq_hz, float)
    ph = np.asarray(phase_rad, float)
    if f.size == 0:
        return ph

    idx = np.argsort(f)
    f_sorted = f[idx]
    ph_sorted = ph[idx]
    # Find block boundaries where there is a big gap
    cuts = [0]
    for i in range(1, len(f_sorted)):
        if (f_sorted[i] - f_sorted[i - 1]) > BAND_GAP_HZ:
            cuts.append(i)
    cuts.append(len(f_sorted))

    # Unwrap each block
    ph_u = ph_sorted.copy()
    for b0, b1 in zip(cuts[:-1], cuts[1:]):
        ph_u[b0:b1] = np.unwrap(ph_u[b0:b1])

    for bi in range(1, len(cuts) - 1):
        prev_end = cuts[bi] - 1
        next_start = cuts[bi]
        if prev_end < 0 or next_start >= len(ph_u):
            continue
        jump = ph_u[next_start] - ph_u[prev_end]
        n = np.round(jump / (2.0 * np.pi))
        ph_u[next_start:cuts[bi + 1]] -= n * 2.0 * np.pi

    # Return in original order
    out = np.empty_like(ph_u)
    out[idx] = ph_u
    return out

def band_aware_unwrap(freq_hz, phase_rad, band_id_per_sample):
    """
    Unwrap phase independently within each band (given by band_id_per_sample),
    then align bands by adding 2*pi*n offsets to minimize boundary jumps in
    frequency order (same as gap_aware_unwrap).
    """
    f = np.asarray(freq_hz, float)
    ph = np.asarray(phase_rad, float)
    b = np.asarray(band_id_per_sample)
    if f.size == 0:
        return ph
    
    # sort by frequency
    idx = np.argsort(f)
    ph_sorted = ph[idx]
    b_sorted = b[idx]

    ph_u = ph_sorted.copy()

    # unwrap within each band separately
    for band_id in np.unique(b_sorted):
        sel = (b_sorted == band_id)
        if np.sum(sel) >= 3:
            ph_u[sel] = np.unwrap(ph_u[sel])

    # align bands: walk in frequency order and whenever the band changes,
    # apply a 2*pi*n shift to the entire new band to minimize the jump at the boundary.
    unique_in_order = []
    for band_id in b_sorted:
        if not unique_in_order or band_id != unique_in_order[-1]:
            unique_in_order.append(band_id)

    for i in range(1, len(unique_in_order)):
        b_prev = unique_in_order[i - 1]
        b_next = unique_in_order[i]

        prev_idx = np.where(b_sorted == b_prev)[0]
        next_idx = np.where(b_sorted == b_next)[0]
        if prev_idx.size == 0 or next_idx.size == 0:
            continue

        prev_end = prev_idx[-1]
        next_start = next_idx[0]

        jump = ph_u[next_start] - ph_u[prev_end]
        n = np.round(jump / (2.0 * np.pi))
        ph_u[next_idx] -= n * 2.0 * np.pi

    # return to original order
    out = np.empty_like(ph_u)
    out[idx] = ph_u
    return out

# -----------------------------
# Helpers for flagging bad freq. channels / IFs
# -----------------------------
def mad_sigma(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return 1.4826 * mad

def _quinn_tau(x):
    # x is scalar float
    return 0.25*np.log1p(3.*x*x + 6.*x) - (np.sqrt(6.)/24.)*np.log1p(-2.*np.sqrt(2./3.)/(x+1.+np.sqrt(2./3.)))

def quinn_estimate_complex_1d(Xm1, X0, Xp1):
    """
    Quinn’s second estimator for fractional-bin peak offset using complex FFT bins.
    Returns fractional offset delta in bins (can be +/- ~0.5 typically).
    Based on the same algebra used in your C++ code.

    Xm1 = X[k0-1], X0 = X[k0], Xp1 = X[k0+1]
    """
    denom = (X0.real*X0.real + X0.imag*X0.imag)
    if denom <= 0:
        return 0.0
    Ap = (Xp1.real*X0.real + Xp1.imag*X0.imag) / denom
    Am = (Xm1.real*X0.real + Xm1.imag*X0.imag) / denom

    # Guard against pathological values
    if Ap >= 1.0:
        Ap = 0.999999999
    if Am >= 1.0:
        Am = 0.999999999

    Dp = -Ap/(1.-Ap)
    Dm =  Am/(1.-Am)
    return 0.5*(Dp + Dm) + _quinn_tau(Dp*Dp) - _quinn_tau(Dm*Dm)


def robust_phase_residual_scatter(f_hz, ph_rad):
    f = np.asarray(f_hz, float)
    ph = np.asarray(ph_rad, float)
    ok = np.isfinite(f) & np.isfinite(ph)
    f = f[ok]; ph = ph[ok]
    if f.size < 5:
        return np.nan
    idx = np.argsort(f)
    f = f[idx]; ph = ph[idx]

    ph_u = np.unwrap(ph)  # within IF only

    A = np.vstack([f, np.ones_like(f)]).T
    coef, *_ = np.linalg.lstsq(A, ph_u, rcond=None)     # linear fit
    a, b = coef

    resid = ph_u - (a * f + b)
    return float(mad_sigma(resid))  # radians

def per_if_metrics(f_if, amp_if, ph_if_rad, min_frac=0.3):
    f   = np.asarray(f_if, float)
    amp = np.asarray(amp_if, float)
    ph  = np.asarray(ph_if_rad, float)
    finite = np.isfinite(amp) & np.isfinite(ph) & np.isfinite(f)
    ok = finite & (amp > 0)
    if finite.size == 0:
        return None
    if (np.sum(ok) / np.sum(finite)) < min_frac:
        return None  # insufficient data
    a = amp[ok]
    med_log = float(np.median(np.log(a)))
    scat_lin = float(mad_sigma(a))
    ph_scat = robust_phase_residual_scatter(f[ok], ph[ok])  # radians
    return {"med_log": med_log, "scat_lin": scat_lin, "ph_scat": ph_scat, "n_valid": int(np.sum(ok)),}


def flag_if_edges(nchan: int, edge_frac: float):
    """
    Return a boolean mask of length nchan where True means 'keep'.
    Flags edge_frac of channels at each edge (start and end).
    """
    if nchan <= 0:
        return np.zeros(0, dtype=bool)
    edge_frac = float(edge_frac)
    if edge_frac <= 0:
        return np.ones(nchan, dtype=bool)
    if edge_frac >= 0.5:
        # If user asks to flag >=50% each edge, nothing remains.
        return np.zeros(nchan, dtype=bool)
    n_edge = int(np.ceil(edge_frac * nchan))
    keep = np.ones(nchan, dtype=bool)
    if n_edge > 0:
        keep[:n_edge] = False
        keep[-n_edge:] = False
    return keep


def auto_flag_ifs(if_metrics, if_to_band, k_med=3.5, k_scat=3.5, k_phase=3.5):
    """
    if_metrics: dict {IF: {"med_log":..., "scat_lin":..., "ph_scat":..., "n_valid":...}}
    returns: set of IFs to flag
    """
    flagged = set()
    valid = [(IF, m) for IF, m in if_metrics.items() if m is not None]
    if not valid:
        return flagged

    # Global dispersion flagging
    scats = np.array([m["scat_lin"] for _, m in valid if np.isfinite(m["scat_lin"])])
    if scats.size >= 5:
        exp_center = np.median(scats)
        exp_scale  = mad_sigma(scats)
        if np.isfinite(exp_scale) and exp_scale > 1e-12:
            thr = exp_center + k_scat * exp_scale
            for IF, m in valid:
                if np.isfinite(m["scat_lin"]) and m["scat_lin"] > thr:
                    flagged.add(IF)
    # Global median flagging
    meds = np.array([m["med_log"] for _, m in valid if np.isfinite(m["med_log"])])
    if meds.size >= 5:
        center = np.median(meds)
        scale  = mad_sigma(meds)
        if np.isfinite(scale) and scale > 1e-12:
            lo = center - 1.5 * k_med * scale
            hi = center + 1.5 * k_med * scale
            for IF, m in valid:
                v = m["med_log"]
                if np.isfinite(v) and (v < lo or v > hi):
                    flagged.add(IF)

    # Within-band median flagging
    bands = {}
    for IF, m in valid:
        b = if_to_band.get(IF)
        if b is not None and np.isfinite(m["med_log"]):
            bands.setdefault(b, []).append((IF, m["med_log"]))

    for b, pairs in bands.items():
        vals = np.array([v for _, v in pairs], float)
        c = np.median(vals)
        s = mad_sigma(vals)
        if np.isfinite(s) and s > 1e-12:
            for IF, v in pairs:
                if abs(v - c) > k_med * s:
                    flagged.add(IF)

    # Global phase residual scatter flagging
    phs = np.array([m["ph_scat"] for _, m in valid if np.isfinite(m["ph_scat"])])
    if phs.size >= 5:
        ph_center = np.median(phs)
        ph_scale  = mad_sigma(phs)
        if np.isfinite(ph_scale) and ph_scale > 1e-12:
            thr = ph_center + k_phase * ph_scale
            for IF, m in valid:
                if np.isfinite(m["ph_scat"]) and m["ph_scat"] > thr:
                    flagged.add(IF)

    # Ensure each band keeps at least a fractional minimum of IFs
    band_to_ifs = {}
    for IF, m in valid:
        b = if_to_band.get(IF)
        if b is not None:
            band_to_ifs.setdefault(b, []).append(IF)
    for b, ifs in band_to_ifs.items():
        n_total = len(ifs)
        if n_total == 0:
            continue
        n_keep = sum(IF not in flagged for IF in ifs)
        if n_keep < 2:      # keep at least 2 IFs in each band
            for IF in ifs:  # Drop whole band: flag all IFs in it
                flagged.add(IF)

    return flagged



# -----------------------------
# Band combination helpers
# -----------------------------
def infer_if_to_band_from_data(data, IFs, BAND_GAP_HZ=0.9e9):
    """
    Infer IF->band mapping by clustering IF center-frequencies.
    IFs whose centers are separated by > BAND_GAP_HZ are considered different bands.

    Returns
    -------
    if_to_band : dict {IF: band_id (1..nbands)}
    """
    centers = []
    for IF in IFs:
        f = np.asarray(data["Frequency"][IF], float)
        f = f[np.isfinite(f)]
        if f.size == 0:
            continue
        centers.append((IF, float(np.median(f))))

    # sort by center frequency
    centers.sort(key=lambda x: x[1])

    if_to_band = {}
    band_id = 1
    prev_fc = None
    for IF, fc in centers:
        if prev_fc is not None and (fc - prev_fc) > float(BAND_GAP_HZ):
            band_id += 1
        if_to_band[IF] = band_id
        prev_fc = fc

    return if_to_band

def iter_nonempty_subsets(items):
    """Yield all non-empty subsets of items as tuples (in increasing size order)."""
    items = list(items)
    for r in range(1, len(items) + 1):
        for comb in combinations(items, r):
            yield comb


# -----------------------------
# Helpers for candidate selection
# -----------------------------
def delay_coherence(f_hz, ph_rad, w, tau_s):
    """
    Coherence in [0,1] after removing best-fit constant phase at a given delay.
    """
    f = np.asarray(f_hz, float)
    ph = np.asarray(ph_rad, float)
    w = np.asarray(w, float)
    ok = np.isfinite(f) & np.isfinite(ph) & np.isfinite(w)
    f = f[ok]; ph = ph[ok]; w = w[ok]
    if f.size < 3:
        return np.nan
    r = np.sum(w * np.exp(1j * (ph - 2.0*np.pi*f*tau_s)))  # complex
    denom = np.sum(w)
    if denom <= 0:
        return np.nan
    # remove constant phase optimally => just take magnitude of r
    return float(np.abs(r) / denom)

def refine_tau_local_exact(f_hz, ph_rad, w, tau0_s, dtau_s, ngrid=101):
    """
    Refine a candidate delay by direct exact coherence evaluation around
    the FFT-bin peak. This keeps FFT candidate discovery but avoids relying
    on Quinn for gapped/multi-band spectra.
    """
    f = np.asarray(f_hz, float)
    ph = np.asarray(ph_rad, float)
    w = np.asarray(w, float)
    ok = np.isfinite(f) & np.isfinite(ph) & np.isfinite(w) & (w > 0)
    f = f[ok]
    ph = ph[ok]
    w = w[ok]
    if f.size < 3 or not np.isfinite(tau0_s) or not np.isfinite(dtau_s) or dtau_s <= 0:
        return float(tau0_s)
    grid = tau0_s + np.linspace(-0.75 * dtau_s, 0.75 * dtau_s, int(ngrid))
    vals = np.empty(grid.size, dtype=float)
    wsum = np.sum(w)
    if wsum <= 0:
        return float(tau0_s)
    for i, tau in enumerate(grid):
        vals[i] = abs(np.sum(w * np.exp(1j * (ph - 2.0 * np.pi * f * tau)))) / wsum
    return float(grid[int(np.argmax(vals))])

def spectrum_quality(fringe, cand_sub, ratio_cap=2.0):
    """
    Compute simple robustness indicators for a normalized delay spectrum.
    Combine:
     - snr_like from the whole spectrum (robust noise estimate),
     - uniqueness from candidate peaks (real local maxima, de-duplicated).
    Returns:
      p1      : highest peak
      p2      : second highest peak
      ratio   : p1/p2 (peak uniqueness)
      snr_like: (p1 - median) / MAD_sigma (contrast vs robust noise)
    """
    x = np.asarray(fringe, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan, np.nan, np.nan, np.nan, 0
    # robust "noise floor" / contrast
    med = float(np.median(x))
    sig = float(mad_sigma(x))
    if not np.isfinite(sig) or sig <= 1e-12:
        sig = 1e-12
    # Use peak height from candidates if available, else spectrum max
    if cand_sub:
        p1 = float(cand_sub[0]["fringe"])
        n_cand = len(cand_sub)
        p2 = float(cand_sub[1]["fringe"]) if n_cand >= 2 else 0.0
    else:
        p1 = float(np.max(x))
        p2 = 0.0
        n_cand = 0
    if p2 > 0:
        ratio = p1 / p2
    else:
        ratio = ratio_cap  # treat single-candidate as "very unique" but bounded
    snr_like = float((p1 - med) / sig)
    return float(p1), float(p2), ratio, snr_like, n_cand


def distribute_subset_weight(cand_sub, W_subset, gamma=1.0):
    """
    Distribute a subset-level quality weight across all candidate peaks.
    - Weaker candidates receive a fraction of the weight:
      w_i = W_subset * max(floor_frac, (fringe_i/p1)^gamma)
    - gamma > 1 penalizes weaker peaks strongly.
    """
    if not cand_sub or not np.isfinite(W_subset) or W_subset <= 0:
        return []

    p1 = float(cand_sub[0]["fringe"])
    if not np.isfinite(p1) or p1 <= 0:
        return []

    out = []
    for j, c in enumerate(cand_sub):
        fi = float(c.get("fringe", 0.0))
        if not np.isfinite(fi) or fi <= 0:
            w = 0.0
        else:
            rel = fi / p1
            w = float(W_subset * (rel ** float(gamma)))
        out.append(w)
    return out


def subset_weight_from_quality(p1, ratio, snr_like):
    """
    Reliability weight for a subset spectrum. Higher when:
      - uniqueness: main peak stands clearly above the rest (large p1/p2),
      - spectrum has good contrast relative to noise (snr_like).
    """
    if not np.isfinite(p1) or p1 <= 0:
        return 0.0
    if not np.isfinite(snr_like) or snr_like <= 0:
        return 0.0
    if not np.isfinite(ratio):
        return 0.0
    return float(snr_like * ratio)


def cluster_candidates_by_tau_best(cand_pool, tol_ns):
    """
    Cluster candidates within tol_ns and choose representative tau as the
    best single member (highest per-vote weight), not the weighted mean.

    Returns clusters sorted by total support weight desc.
    """
    if not cand_pool:
        return []

    items = sorted(cand_pool, key=lambda d: d["tau_ns"])
    clusters_raw = []
    cur = [items[0]]

    for d in items[1:]:
        if abs(d["tau_ns"] - cur[-1]["tau_ns"]) <= tol_ns:
            cur.append(d)
        else:
            clusters_raw.append(cur)
            cur = [d]
    clusters_raw.append(cur)

    out = []
    for cl in clusters_raw:
        total_w = float(sum(x["weight"] for x in cl))
        n = len(cl)
        sources = sorted({x["source"] for x in cl})

        # representative tau: choose member with highest weight
        best_member = max(cl, key=lambda x: x["weight"])
        tau_ns_rep = float(best_member["tau_ns"])
        tau_s_rep  = float(best_member.get("tau_s", tau_ns_rep * 1e-9))

        out.append({
            "tau_ns": tau_ns_rep,
            "tau_s": tau_s_rep,
            "weight": total_w,
            "n": n,
            "sources": sources,
            "best_vote_weight": float(best_member["weight"]),
        })

    out.sort(key=lambda d: d["weight"], reverse=True)
    return out



# -----------------------------
# Delay estimation via "delay spectrum" (NUDFT on a delay grid)
# -----------------------------
def find_fringe_peak_candidates(
    taus_s: np.ndarray,
    R: np.ndarray,
    *,
    peak_frac: float = 0.6,
    delay_window_ns: float = 0.5,
    max_candidates: Optional[int] = None,
):
    """
    Find candidate delay peaks in |R(tau)| above peak_frac*max, and
    keep only peaks separated by >= delay_window_ns.

    Returns: list of dicts:
        [{"tau_s":..., "tau_ns":..., "fringe":..., "idx":...}, ...]
    """
    taus_s = np.asarray(taus_s, float)
    fringe = np.abs(np.asarray(R))

    if taus_s.size < 3:
        return []

    m = np.nanmax(fringe)
    if not np.isfinite(m) or m <= 0:
        return []

    thr = float(peak_frac) * float(m)
    if thr <= 0:
        return []

    # local maxima: i where fringe[i-1] < fringe[i] >= fringe[i+1]
    # (use >= on right to keep flat-tops deterministically)
    peaks = []
    for i in range(1, len(fringe) - 1):
        if fringe[i] >= thr and fringe[i] > fringe[i - 1] and fringe[i] >= fringe[i + 1]:
            peaks.append(i)

    if not peaks:
        return []

    # sort by descending peak height
    peaks.sort(key=lambda i: fringe[i], reverse=True)

    # greedy selection with minimum separation
    win_s = float(delay_window_ns) * 1e-9
    chosen = []
    for i in peaks:
        tau_i = float(taus_s[i])
        if all(abs(tau_i - c["tau_s"]) >= win_s for c in chosen):
            chosen.append({
                "tau_s": tau_i,
                "tau_ns": tau_i * 1e9,
                "fringe": float(fringe[i]),
                "idx": int(i),
            })
            if max_candidates is not None and len(chosen) >= int(max_candidates):
                break

    return chosen

'''
def compute_phi0_and_fzero(
    f_sorted_hz: np.ndarray,
    ph_u_sorted_rad: np.ndarray,
    w_sorted: np.ndarray,
    tau_s: float,
):
    """
    Given sorted f and unwrapped phase (sorted), compute phi0 and f_zero for this tau.
    """
    f = np.asarray(f_sorted_hz, float)
    ph_u = np.asarray(ph_u_sorted_rad, float)
    w = np.asarray(w_sorted, float)

    fmin, fmax = float(f[0]), float(f[-1])

    # Estimate phi0 (constant phase) after removing delay term
    resid = ph_u - (2.0 * np.pi * f * tau_s)
    phi0 = float(np.angle(np.sum(w * np.exp(1j * resid))))

    # Frequency where fitted phase is 0 (mod 2*pi):
    #     phi0 + 2*pi tau f = 2*pi m
    # Choose a solution that falls inside [fmin,fmax]    resid = ph_u - (2.0 * np.pi * f * tau_hat)
    if abs(tau_s) < 1e-20:
        return phi0, 0.
    f0_raw = (-phi0) / (2.0 * np.pi * tau_s)  # One branch (m=0)
    step = 1.0 / tau_s  # Solutions differ by 1/tau
    # Choose integer m to land inside band, or closest to band center
    target = 0.5 * (fmin + fmax)
    m = np.round((target - f0_raw) / step)
    f_zero = float(f0_raw + m * step)
    # If still out-of-band due to weird edge cases, clamp to nearest in-band solution
    if f_zero < fmin or f_zero > fmax:
        m2 = np.round((np.clip(f0_raw, fmin, fmax) - f0_raw) / step)
        f_zero = float(f0_raw + m2 * step)

    return phi0, f_zero
'''
def compute_phi0_and_fzero(
    f_sorted_hz: np.ndarray,
    ph_u_sorted_rad: np.ndarray,
    w_sorted: np.ndarray,
    tau_s: float,
):
    """
    Given sorted f and unwrapped phase, compute a stable reference phase
    and f_zero for this tau.

    Returned phi0_rad is the phase at f_ref, not at absolute 0 Hz.
    f_zero_hz is the frequency where:

        phi_ref + 2*pi*tau*(f - f_ref) = 0 mod 2*pi

    choosing the branch closest to the middle of the observed frequency span.
    """
    f = np.asarray(f_sorted_hz, float)
    ph_u = np.asarray(ph_u_sorted_rad, float)
    w = np.asarray(w_sorted, float)

    ok = np.isfinite(f) & np.isfinite(ph_u) & np.isfinite(w) & (w > 0)
    f = f[ok]
    ph_u = ph_u[ok]
    w = w[ok]

    if f.size < 3:
        return np.nan, np.nan

    fmin, fmax = float(np.min(f)), float(np.max(f))
    f_ref = float(np.average(f, weights=w))

    resid = ph_u - 2.0 * np.pi * (f - f_ref) * tau_s
    phi_ref = float(np.angle(np.sum(w * np.exp(1j * resid))))

    if abs(tau_s) < 1e-20:
        return phi_ref, np.nan

    f0_raw = f_ref - phi_ref / (2.0 * np.pi * tau_s)
    step = 1.0 / tau_s

    target = 0.5 * (fmin + fmax)
    m = np.round((target - f0_raw) / step)
    f_zero = float(f0_raw + m * step)

    if f_zero < fmin or f_zero > fmax:
        candidates = f0_raw + step * np.arange(-10, 11)
        j = int(np.argmin(np.abs(candidates - target)))
        f_zero = float(candidates[j])

    return phi_ref, f_zero

# ------------------------------------------------------
def estimate_delay_fft_gridded(
    freq_hz,
    phase_rad,
    amp=None,
    padding=16,
    BAND_GAP_HZ=1e9,
    use_amp_weight=False,
    *,
    band_id_per_sample=None,
    fringe_peak_frac=0.6,
    delay_window_ns=0.5,
    max_candidates=None,
    enable_candidates=True,
    refine_method="quinn",
    apply_freq_window=False,
    window_kind="hann",
):
    """
    Estimate X/Y delay tau from phase vs frequency using an FFT on a
    gridded frequency axis (i.e., a multi-band delay/fringe approach).

    Parameters
    ----------
    freq_hz : array_like
        Sample frequencies in Hz (may include gaps).
    phase_rad : array_like
        Phase samples in radians (can be wrapped).
    amp : array_like or None
        Optional weights per sample (e.g., XY amplitude). If provided and
        use_amp_weight=True, the gridded phasors are weighted by amp.
        Regardless of weighting choice, amp is still used when estimating
        phi0/f_zero and delay-coherence diagnostics.
    padding : int
        Zero-padding factor applied to the gridded freq. array before FFT.
        Higher values produce a finer delay grid.
    BAND_GAP_HZ : float
        Freq gap threshold used ONLY for phase unwrapping/stitching across
        separated bands.
        Affects phi0/f_zero stability, not the FFT delay spectrum itself.
    use_amp_weight : bool
        If True, accumulate amp * exp(i*phi) into frequency bins.
        If False (default), accumulate only exp(i*phi).
    band_id_per_sample : array_like or None
        Optional band identifier for each freq. sample. If given, phase is
        unwrapped independently per band and then stitched by band order.
    fringe_peak_frac : float
        Candidate peak threshold as a fraction of max. fringe amplitude.
    delay_window_ns : float
        Minimum separation, in ns, between returned candidate peaks.
    max_candidates : int or None
        Max. number of candidate peaks to return. None means no cap.
    enable_candidates : bool
        If True, return local fringe-peak candidates.
        If False, only the global maximum is used as fallback.
    refine_method : {"exact", "quinn", "none"}
        Candidate sub-bin refinement method. 
        "quinn" uses Quinn's estimator on neighboring FFT bins. 
        "exact" refines locally by direct coherence evaluation on freqs.
        "none" keeps the FFT-bin delay.
    apply_freq_window : bool
        If True, apply a taper independently to each contiguous occupied
        frequency-grid run before the FFT. 
        This may reduce sidelobes but broadens delay peaks.
    window_kind : {"hann", "blackman"}
        Frequency-domain taper used when apply_freq_window=True.

    Returns
    -------
    tau_hat_s : float
        Best-estimate X/Y delay (seconds). 
        By default this is the strongest non-redundant candidate peak. 
        If no peaks meet candidate threshold,
        falls back to global max. |R(tau)|.
    f_zero_hz : float
        Representative frequency where the fitted phase is 0 (mod 2*pi),
        computed from the selected tau and reference phase.
    taus : ndarray
        Delay axis for the FFT output (seconds), length Npad.
    R : ndarray (complex)
        Complex FFT delay spectrum. Use np.abs(R) for the “fringe amplitude”.
    ph_u : ndarray
        Gap-aware or band-aware unwrapped/stiched phase (radians) in
        sorted-frequency order.
        Used only to stabilize phi0/f_zero estimation.
    phi0_hat : float
        Reference phase offset corresponding to tau_hat_s.
    candidates: list of dict
        List of peak-delay candidates (possible solutions when sidelobes or multi-band ambiguity produce multiple strong peaks). Each entry has:
            - "tau_s"     : refined delay (seconds)
            - "tau_ns"    : refined delay (nanoseconds)
            - "f_zero_hz" : corresponding 0-phase frequency (Hz) for that tau
            - "phi0_rad"  : corresponding constant phase offset (rad)
            - "fringe"    : |R(tau)| amplitude at the *unrefined* peak bin
                           (used for ranking/selection)
        Candidates are selected as local maxima of |R(tau)| above
            `fringe_peak_frac * max(|R|)` and then de-duplicated
            so any two returned candidates differ by at least
            `delay_window_ns` in delay.
    """
    f = np.asarray(freq_hz, float)
    ph = np.asarray(phase_rad, float)

    if amp is None:
        w = np.ones_like(f)
    else:
        w = np.asarray(amp, float)
    # Optional weights: keep them sane (non-negative, finite)
    if w is not None:
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)


    # Keep only finite samples (recommended; otherwise FFT becomes NaN)
    ok = np.isfinite(f) & np.isfinite(ph) & np.isfinite(w)
    f = f[ok]
    ph = ph[ok]
    w = w[ok]
    if f.size < 3:
        raise ValueError("Not enough valid samples")

    # Sort
    idx = np.argsort(f)
    f = f[idx]; ph = ph[idx]; w = w[idx]

    # Unwrap+stitch (for phi0 and f_zero only)
    if band_id_per_sample is None:
        ph_u = gap_aware_unwrap(f, ph, BAND_GAP_HZ=BAND_GAP_HZ)
    else:
        b = np.asarray(band_id_per_sample)[ok]   # apply finite mask
        b = b[idx]                               # apply freq sort
        ph_u = band_aware_unwrap(f, ph, b)

    # Estimate channel spacing (assumes roughly uniform per IF)
    # Robust df estimate:
    #  - ignore non-positive steps
    #  - ignore large gaps (top 20% of steps)
    #  - choose a "typical smallest" spacing
    df_all = np.diff(f)
    df_pos = df_all[np.isfinite(df_all) & (df_all > 0)]
    if df_pos.size == 0:
        raise ValueError("Bad frequency spacing (no positive steps)")
    # Trim big gaps: keep steps <= 80th percentile
    cut = np.percentile(df_pos, 80.0)
    df_small = df_pos[df_pos <= cut]
    if df_small.size < 3:
        df_small = df_pos
    # Use median of lower half to bias toward channel spacing, not gaps
    cut2 = np.percentile(df_small, 50.0)
    df_core = df_small[df_small <= cut2]
    if df_core.size < 3:
        df_core = df_small
    df = float(np.median(df_core))
    if not np.isfinite(df) or df <= 0:
        df = float(np.median(df_pos))
    if not np.isfinite(df) or df <= 0:
        raise ValueError("Bad frequency spacing")

    fmin, fmax = float(f[0]), float(f[-1])
    N = int(np.round((fmax - fmin) / df)) + 1  # Uniform frequency grid

    # Accumulate complex phasors into nearest frequency bins.
    phasor = np.zeros(N, dtype=np.complex128)
    '''
    # If user asked for amp weighting, use it; otherwise, equal weight.
    # If absent, w is all-ones so will keep equal weight.
    if use_amp_weight:
        contrib = w * np.exp(1j * ph)
    else:
        contrib = np.exp(1j * ph)
    # nearest-bin indexing
    k = np.rint((f - fmin) / df).astype(int)
    k = np.clip(k, 0, N - 1)
    np.add.at(phasor, k, contrib)
    '''
    # Use hit normalization so duplicated/near-duplicated samples do not
    # overweight one grid cell.
    hit = np.zeros(N, dtype=np.float64)
    if use_amp_weight:
        contrib = w * np.exp(1j * ph)
        hit_contrib = w
    else:
        contrib = np.exp(1j * ph)
        hit_contrib = np.ones_like(w)
    k = np.rint((f - fmin) / df).astype(int)
    k = np.clip(k, 0, N - 1)
    np.add.at(phasor, k, contrib)
    np.add.at(hit, k, hit_contrib)
    good = hit > 0
    phasor[good] /= hit[good]
    # Optional per-occupied-run taper.
    # For production, keep APPLY_FREQ_WINDOW=False unless sidelobes are a problem.
    if apply_freq_window and np.any(good):
        idx_occ = np.flatnonzero(good)
        splits = np.where(np.diff(idx_occ) > 1)[0] + 1
        runs = np.split(idx_occ, splits)
        for run in runs:
            if run.size < 8:
                continue
            if window_kind == "hann":
                win = np.hanning(run.size)
            elif window_kind == "blackman":
                win = np.blackman(run.size)
            else:
                continue
            phasor[run] *= win

    # Zero-pad (centered)
    Npad = N * int(padding)
    padded = np.zeros(Npad, dtype=np.complex128)
    start = (Npad - N) // 2
    padded[start:start + N] = phasor

    # FFT -> delay spectrum magnitude
    spec = np.fft.fftshift(np.fft.fft(padded))
    fringe = np.abs(spec)

    # --- Normalize fringe to max = 1 for consistent comparison ---
    if fringe.size > 0:
        maxfr = np.nanmax(fringe)
        if np.isfinite(maxfr) and maxfr > 0:
            spec = spec / maxfr
            fringe = fringe / maxfr
    # Delay axis: fftfreq uses cycles per sample; sample spacing is df (Hz)
    # => delay step = 1 / (Npad * df)
    taus = np.fft.fftshift(np.fft.fftfreq(Npad, d=df))  # seconds

    # find peak candidates
    if enable_candidates:
        cand = find_fringe_peak_candidates(
            taus, spec,
            peak_frac=fringe_peak_frac,
            delay_window_ns=delay_window_ns,
            max_candidates=max_candidates,
        )
    else:
        cand = []

    # For each candidate tau, compute phi0 and f_zero
    candidates = []
    for c in cand:
        tau_s = c["tau_s"]
        i = c["idx"]
        tau_ref = tau_s
        if 1 <= i <= (len(spec) - 2):
            dtau = float(abs(taus[1] - taus[0]))
            # Quinn refinement on complex FFT bins (sub-bin)
            if refine_method == "quinn":
                delta = quinn_estimate_complex_1d(spec[i - 1], spec[i], spec[i + 1])
                if not np.isfinite(delta):
                    delta = 0.0
                delta = float(np.clip(delta, -0.6, 0.6))
                tau_ref = float(taus[i] + delta * dtau)
            # Refine locally using exact coherence on the original freqs.
            #       more robust than Quinn for gapped/multi-band spectra.
            elif refine_method == "exact":
                tau_ref = refine_tau_local_exact(f, ph, w, tau_s, dtau, ngrid=101)
            elif refine_method == "none":
                tau_ref = tau_s
            else:
                raise ValueError("Unknown refine_method=%r. Use 'exact', 'quinn', or 'none'." % refine_method)

        phi0, f_zero = compute_phi0_and_fzero(f, ph_u, w, tau_ref)

        candidates.append({
            "tau_s": tau_ref,
            "tau_ns": tau_ref * 1e9,
            "f_zero_hz": float(f_zero),
            "phi0_rad": float(phi0),
            "fringe": float(fringe[i]),
        })

    # Keep the old “best” candidate as the first one (strongest fringe)
    tau_hat = candidates[0]["tau_s"] if candidates else float(taus[np.argmax(np.abs(spec))])
    f_zero_hat = candidates[0]["f_zero_hz"] if candidates else np.nan
    phi0_hat = candidates[0]["phi0_rad"] if candidates else np.nan
    if not candidates:
        phi0_hat, f_zero_hat = compute_phi0_and_fzero(f, ph_u, w, tau_hat)

    return tau_hat, f_zero_hat, taus, spec, ph_u, phi0_hat, candidates



# -----------------------------
# Plotting helpers
# -----------------------------
def plot_amp_phase(
    f_keep, amp_keep, ph_keep_deg,
    f_flagged, amp_flagged, ph_flagged_deg,
    out_png, title_prefix=""
):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))
    # Kept points
    ax1.scatter(f_keep, amp_keep, s=10, label="kept")
    ax2.scatter(f_keep, ph_keep_deg, s=10, label="kept")
    # Flagged points
    if len(f_flagged) > 0:
        ax1.scatter(f_flagged, amp_flagged, s=10, alpha=0.3, c="r", marker="x", linewidths=1.0, label="flagged")
        ax2.scatter(f_flagged, ph_flagged_deg, s=10, alpha=0.3, c="r", marker="x", linewidths=1.0, label="flagged")

    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_ylabel("XY amplitude")
    ax1.set_ylim(0.,np.max(amp_keep)*2.)
    ax1.set_title("%s XY amplitude vs frequency" % (title_prefix,))
    ax1.legend()
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_ylabel("XY phase (deg)")
    ax2.set_title("%s XY phase vs frequency" % (title_prefix,))
    #ax2.legend()

    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def plot_delay_spectrum(taus_s, R, tau_hat_s, out_png, title_prefix="", candidates=None):
    taus_s = np.asarray(taus_s, float)
    fringe = np.abs(np.asarray(R))

    plt.figure(figsize=(10, 4))
    plt.plot(taus_s * 1e9, fringe)
    plt.xlabel(r"Delay $\tau$ (ns)")
    plt.ylabel(r"|R($\tau$)|")
    plt.xlim(-33, 33)

    # Highlight chosen tau_hat (may be refined off-bin)
    if tau_hat_s is not None and np.isfinite(tau_hat_s):
        tau_hat_ns = float(tau_hat_s * 1e9)
        plt.axvline(tau_hat_ns, linestyle="--", color="k", linewidth=1.5, label=("tau_hat %.3f ns" % (tau_hat_ns,)))
    # highlight candidate peaks
    if candidates:
        for k, c in enumerate(candidates, 1):
            tns = float(c.get("tau_ns", np.nan))
            if np.isfinite(tns):
                plt.axvline(tns, linestyle=":", linewidth=1.0, alpha=0.8, label=("cand %d %.3f ns" % (k, tns)))

    plt.title("%s delay spectrum (tau_hat %.3f ns)" % (title_prefix, tau_hat_s*1e9))
    plt.tight_layout()

    plt.savefig(out_png)
    plt.close()



# -----------------------------
# Main
# -----------------------------
def CROSS_POL_DELAY_SCANNER(
    path: str,
    *,
    out_dir=None,
    EXPLORE_CANDIDATES=False,
    FRINGE_PEAK_FRAC=0.9,
    MAX_CANDIDATES=5,     # set max num. of candidates for XPOL_DELAY
    # Flagging options
    FLAG_EDGES=True,
    EDGE_FRAC=0.025,      # flag X% at start AND X% at end of each IF
    FLAG_NOISY_IFS=False,
    BAND_GAP_HZ=0.9e9,    # gap threshold (Hz) to split VGOS bands
    padding=16,           # zero-padding
    # Band-combination exploration
    TRY_BAND_COMBINATIONS=False,
    IF_TO_BAND=None,
    refine_method='quinn',
    apply_freq_window=False,
    fft_window_kind='hann'
):
    """
    Run XPOL delay estimation for all antennas in a POLCAL gain file.

    This pipeline:
      1) Loads POLCAL data from `path` (pickle).
      2) For each IF and antenna, optionally flags edge channels and/or 
         drops entire IFs based on XY amplitude quality metrics.
      3) Concatenates all kept channels across IFs per antenna.
      4) Estimates cross-pol delay using an FFT-gridded delay spectrum,
         with optional candidate refinement and multiple delay candidates.
      5) Writes diagnostic plots and a text summary file.

    Parameters
    ----------
    path : str
        Path to the pickled POLCAL gain file (e.g., "POLCAL_GAINS_h.dat").
    out_dir : str or None, optional
        Output directory for plots and results. If None, defaults to
        "<dirname(path)>/XPOL_DELAY_INSPECT".

    EXPLORE_CANDIDATES: bool, optional
        Set to True to explore candidate fringe peaks

    FRINGE_PEAK_FRAC : float, optional
        Candidate peak threshold as a fraction of maximum fringe amplitude:
            keep peaks with |R(tau)| >= FRINGE_PEAK_FRAC * max(|R|).
        Example: 1/3 keeps peaks above one-third of the maximum.

    MAX_CANDIDATES : int or None, optional
        Max number of candidate peaks returned per antenna.
        None means only one candidate.

    FLAG_EDGES : bool, optional
        If True, flag (drop) EDGE_FRAC of channels at both the low- and high-frequency edges of each IF.
    EDGE_FRAC : float, optional
        Fraction of channels to flag at each IF edge (0 <= EDGE_FRAC < 0.5).
        Example: 0.1 flags 10% at the start and 10% at the end.

    FLAG_NOISY_IFS : bool, optional
        If True, drop entire IFs for an antenna when amplitude quality metrics fail (scatter too large, too few valid channels, too many out-of-range channels).

    BAND_GAP_HZ : float, optional
        Frequency gap threshold (Hz) used for gap-aware phase unwrapping (stitching) inside delay estimation (affects phi0/f_zero stability, not the delay spectrum).
    padding : int, optional
        Zero-padding factor used in the FFT-gridded delay spectrum. Higher values produce a finer delay grid.

    TRY_BAND_COMBINATIONS : bool, optional
        If True, and the full-band spectrum has more than one candidate,
        estimate delays for multi-band subsets and use a
        voting/coherence score to rank the final candidates.
    IF_TO_BAND : dict or None, optional
        Optional mapping {IF: band_id}.
        If None, infers bands from IF center freqs. using BAND_GAP_HZ.
    refine_method : {"exact", "quinn", "none"}, optional
        Candidate sub-bin refinement method passed to estimate_delay_fft_gridded.
        "exact" uses a local direct-coherence scan on the original freqs;
        "quinn" uses Quinn's FFT-bin estimator;
        "none" keeps the FFT-bin delay.
    apply_freq_window : bool, optional
        If True, apply a taper independently to each contiguous occupied
        frequency-grid run before the FFT.
        This may reduce sidelobes but broadens delay peaks.
    fft_window_kind : {"hann", "blackman"}, optional
        Frequency-domain taper used when apply_freq_window=True.

    Returns
    -------
    results : dict
    Dictionary keyed by antenna name. For each antenna:
            results[ant] = {
                "tau_hat_s": float,
                "f_zero_hz": float,
                "candidates": list[dict],
            }
    "candidates" is the list returned by `estimate_delay_fft_gridded()`
        (each entry includes tau_s, tau_ns, f_zero_hz, phi0_rad, fringe).

    Notes
    -----
    - Per-antenna delay estimation uses all kept IFs/channels jointly, so multi-band gaps are handled naturally in the delay spectrum.
    - Candidate peak selection is useful when sidelobes are comparable to the main lobe or when ambiguity peaks appear due to sampling/gaps.
    """
    if  MAX_CANDIDATES is None:
        MAX_CANDIDATES = 1
    data = load_data_file(path)
    if data is None:
        raise SystemExit(1)

    def get_delay_window_ns(f_hz, scale=5.0, fallback_ns=1.0):
        """
        Return delay separation window (in ns) based on bandwidth.
        Uses ~ 1/Bspan (Rayleigh delay resolution).
        Helpful for de-duplicating peaks: candidates closer than  intrinsic resolution are merged.
        """
        f = np.asarray(f_hz, float)
        f = f[np.isfinite(f)]
        if f.size < 2:
            return float(fallback_ns)
        Bspan = float(np.max(f) - np.min(f))  # Hz
        if not np.isfinite(Bspan) or Bspan <= 0:
            return float(fallback_ns)
        return float(2. * scale * 1e9 / Bspan)     # ns

    antennas = list(data["XYadd"].keys())
    IFs_all = sorted(int(k) for k in data["Frequency"].keys())
    # If user didn't provide a mapping, infer from the data
    if IF_TO_BAND is None:
        IF_TO_BAND = infer_if_to_band_from_data(data, IFs_all, BAND_GAP_HZ=BAND_GAP_HZ)
    band_ids_all = sorted({IF_TO_BAND.get(IF) for IF in IFs_all if IF_TO_BAND.get(IF) is not None})


    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(path)), "XPOL_DELAY_INSPECT")
    xpol_dir = out_dir
    os.makedirs(xpol_dir, exist_ok=True)
    combo_plot_dir = None
    if TRY_BAND_COMBINATIONS:
        combo_plot_dir = os.path.join(xpol_dir, "TESTING_BANDS")
        os.makedirs(combo_plot_dir, exist_ok=True)

    # Build per-antenna arrays INCLUDING frequency per sample
    per_ant = {}
    for ant in antennas:
        per_ant[ant] = {
            "f": [], "amp": [], "phase_deg": [],
            "flagged": {"f": [], "amp": [], "phase_deg": []},
            "band": [],
            "bands": {b: {"f": [], "amp": [], "phase_deg": []}
                      for b in band_ids_all},
        }
    # MAIN LOOP — per antenna, then per IF
    for ant in antennas:
        # IFs available for this antenna
        xy_ratio_ant = data.get("XYratio", {}).get(ant, {})
        xy_add_ant   = data.get("XYadd", {}).get(ant, {})
        IFs = set(IFs_all)
        IFs &= set(int(k) for k in xy_ratio_ant.keys())
        IFs &= set(int(k) for k in xy_add_ant.keys())
        IFs = sorted(IFs)
        # get per-IF metrics
        metrics = {}
        for IF in IFs:
            if IF not in data["Frequency"]:
                continue
            f_if = np.asarray(data["Frequency"][IF], float)
            nchan = f_if.size
            edge_keep = flag_if_edges(nchan, EDGE_FRAC) if FLAG_EDGES else np.ones(nchan, bool)

            amp_if = np.asarray(xy_ratio_ant[IF], float)[edge_keep]
            ph_if  = np.asarray(xy_add_ant[IF], float)[edge_keep]

            metrics[IF] = per_if_metrics(f_if[edge_keep], amp_if, ph_if, min_frac=0.3)

        # decide IFs to flag
        flagged_ifs = set()
        if FLAG_NOISY_IFS:
            flagged_ifs = auto_flag_ifs(metrics, IF_TO_BAND, k_med=3., k_scat=3.5, k_phase=5.)
            for IF in flagged_ifs:
                print("Auto-flagging IF %s for %s: %s" % (IF, ant, metrics[IF]))
        # append per-IF data (after flagging)
        for IF in IFs:
            if IF not in data["Frequency"]:
                continue
            f_if = np.asarray(data["Frequency"][IF], float)
            amp_if = np.asarray(xy_ratio_ant[IF], float)
            ph_if  = np.asarray(xy_add_ant[IF], float)
            nchan = f_if.size
            edge_keep = flag_if_edges(nchan, EDGE_FRAC) if FLAG_EDGES else np.ones(nchan, bool)
            # CASE 1 — IF dropped entirely by noise/quality flagging
            if IF in flagged_ifs: # all channels are flagged
                per_ant[ant]["flagged"]["f"].extend(f_if.tolist())
                per_ant[ant]["flagged"]["amp"].extend(amp_if.tolist())
                per_ant[ant]["flagged"]["phase_deg"].extend(ph_if.tolist())
                continue
            # CASE 2 — IF kept, but edges flagged
            f_flag   = f_if[~edge_keep]
            amp_flag = amp_if[~edge_keep]
            ph_flag  = ph_if[~edge_keep]
            if f_flag.size > 0:
                per_ant[ant]["flagged"]["f"].extend(f_flag.tolist())
                per_ant[ant]["flagged"]["amp"].extend(amp_flag.tolist())
                per_ant[ant]["flagged"]["phase_deg"].extend(ph_flag.tolist())
            # CASE 3 — kept channels
            band = IF_TO_BAND.get(IF, None)
            if band is None:
                continue
            f_k   = f_if[edge_keep]
            amp_k = amp_if[edge_keep]
            ph_k  = ph_if[edge_keep]
            per_ant[ant]["f"].extend(f_k.tolist())
            per_ant[ant]["amp"].extend(amp_k.tolist())
            per_ant[ant]["phase_deg"].extend(ph_k.tolist())
            per_ant[ant]["band"].extend([band] * len(f_k))
            if band is not None and band in per_ant[ant]["bands"]:
                bd = per_ant[ant]["bands"][band]
                bd["f"].extend(f_k.tolist())
                bd["amp"].extend(amp_k.tolist())
                bd["phase_deg"].extend(ph_k.tolist())

    # joint plot per antenna (amp/phase vs freq)
    for ant in antennas:
        f_keep = np.asarray(per_ant[ant]["f"], float)
        a_keep = np.asarray(per_ant[ant]["amp"], float)
        p_keep = np.asarray(per_ant[ant]["phase_deg"], float)
        flag = per_ant[ant]["flagged"]
        f_flag = np.asarray(flag["f"], float)
        a_flag = np.asarray(flag["amp"], float)
        p_flag = np.asarray(flag["phase_deg"], float)
        plot_amp_phase(f_keep, a_keep, p_keep, f_flag, a_flag, p_flag, os.path.join(xpol_dir, ("XPOL_%s_amp_phase_allIFs.png" % (ant,))), title_prefix=("%s:" % (ant,)) )

    # Estimate delay per antenna using delay spectrum + Quinn
    results_txt = os.path.join(xpol_dir, "XPOL_DELAY_RESULTS.txt")
    cand_txt    = os.path.join(xpol_dir, "XPOL_DELAY_CANDIDATES_FULLBAND.txt")  # original candidates
    combo_txt = os.path.join(xpol_dir, "XPOL_DELAY_RESULTS_BAND_COMBOS.txt")
    results = {}
    with open(results_txt, "w") as fout, open(cand_txt, "w") as fcand, open(combo_txt, "w") as fcombo:
        fout.write("# Antenna  rank  tau_s  tau_ns  f_zero_hz  score  coh  weight  n_votes  sources\n")
        fcand.write("# Antenna  rank  tau_s  tau_ns  f_zero_hz \n")
        fcombo.write("# Antenna  bands  rank  tau_s  tau_ns  f_zero_hz  nchan\n")

        for ant in antennas:
            fout.write("# ---------------------------------------------\n")
            f = np.asarray(per_ant[ant]["f"], float)
            amp = np.asarray(per_ant[ant]["amp"], float)
            ph_deg = np.asarray(per_ant[ant]["phase_deg"], float)
            ph_rad = np.deg2rad(ph_deg)
            b = np.asarray(per_ant[ant]["band"], int)

            mask = (amp > 1e-12) & np.isfinite(amp) & np.isfinite(ph_rad) & np.isfinite(f)
            f_m = f[mask]
            ph_m = ph_rad[mask]
            amp_m = amp[mask]
            b_m = b[mask]
            if f_m.size < 3:
                fcand.write("%s  1  NaN  NaN  NaN  NaN\n" % (ant,))
                fout.write("%s  1  NaN  NaN  NaN   ---   ---   ---   0   FULL\n" % (ant,))
                continue

            delay_window_ns_main = get_delay_window_ns(f_m)

            tau_hat_s, f_zero_hz, taus, R, ph_u, phi0, candidates = estimate_delay_fft_gridded(
                f_m, ph_m, amp=amp_m,
                padding=padding,
                BAND_GAP_HZ=BAND_GAP_HZ,
                use_amp_weight=False,
                fringe_peak_frac=FRINGE_PEAK_FRAC,
                delay_window_ns=delay_window_ns_main,
                band_id_per_sample=b_m,
                enable_candidates=EXPLORE_CANDIDATES,
                refine_method=refine_method,
                apply_freq_window=apply_freq_window,
                window_kind=fft_window_kind,
            )

            if np.isfinite(tau_hat_s) and np.isfinite(f_zero_hz):
                results[ant] = {
                    "tau_hat_s": float(tau_hat_s),
                    "f_zero_hz": float(f_zero_hz),
                    "candidates": candidates if EXPLORE_CANDIDATES else [],
                }
                tau_ns = tau_hat_s * 1e9
                if EXPLORE_CANDIDATES and candidates:
                    print("%s: %d XPOL_DELAY candidates (sorted by fringe strength)" % (ant, len(candidates)))
                    for j, c in enumerate(candidates, 1):
                        fcand.write("%s  %d  %.9e  %.6f  %.9e  %.6e\n" %
                                    (ant, j, c["tau_s"], c["tau_ns"],
                                     c["f_zero_hz"], c["fringe"]))
                else:
                    fcand.write("%s  1  %.9e  %.6f  %.9e  NaN\n" %
                                (ant, tau_hat_s, tau_ns, f_zero_hz))

                fout.write("%s  1  %.9e  %.6f  %.9e   ---   ---   ---   1   BEST\n" %
                           (ant, tau_hat_s, tau_ns, f_zero_hz))
            else:
                print("Skipping %s: zero/invalid XPOL delay." % ant)
                fcand.write("%s  1  NaN  NaN  NaN  NaN\n" % ant)
                fout.write("%s  1  NaN  NaN  NaN   SKIP_ZERO_DELAY   ---   ---   0   NONE\n" % ant)
                continue

            # Diagnostics
            plot_delay_spectrum(
                taus, R, tau_hat_s,
                os.path.join(xpol_dir, ("XPOL_%s_delay_spectrum.png" % (ant,))),
                title_prefix=("%s:" % (ant,)),
                candidates=candidates
            )

            # explore band combinations ONLY when full-band has ambiguity
            if TRY_BAND_COMBINATIONS and (len(candidates) > 1):
                band_ids = sorted([bb for bb, bd in per_ant[ant]["bands"].items() if len(bd["f"]) >= 3])
                nb = len(band_ids)
                # full-band tolerance for clustering taus (ns)
                tol_ns = 0.5 * get_delay_window_ns(f_m, scale=2.0)

                cand_pool = []
                W_full = 1.0
                for c in candidates:
                    cand_pool.append({
                        "tau_ns": float(c["tau_ns"]),
                        "tau_s": float(c["tau_s"]),
                        "weight": float(W_full * max(0.0, c["fringe"])),
                        "source": "FULL",
                    })

                for band_set in iter_nonempty_subsets(band_ids):
                    if len(band_set) < 2: # explore all subsets size>=2
                        continue
                    if len(band_set) == nb:
                        continue   # skip full band-set (done earlier)
                    bands_label = "+".join(str(bb) for bb in band_set)

                    # Build arrays for this band subset
                    f_sub, amp_sub, ph_deg_sub, b_sub = [], [], [], []
                    for bb in band_set:
                        bd = per_ant[ant]["bands"][bb]
                        f_sub.extend(bd["f"])
                        amp_sub.extend(bd["amp"])
                        ph_deg_sub.extend(bd["phase_deg"])
                        b_sub.extend([bb] * len(bd["f"]))

                    f_sub = np.asarray(f_sub, float)
                    amp_sub    = np.asarray(amp_sub, float)
                    ph_deg_sub = np.asarray(ph_deg_sub, float)
                    b_sub  = np.asarray(b_sub, int)
                    if f_sub.size < 3:
                        continue

                    ph_sub = np.deg2rad(ph_deg_sub)

                    mask_sub = (amp_sub > 1e-12) & np.isfinite(amp_sub) & np.isfinite(ph_sub) & np.isfinite(f_sub)
                    f_m_sub = f_sub[mask_sub]
                    ph_m_sub = ph_sub[mask_sub]
                    amp_m_sub = amp_sub[mask_sub]
                    b_m_sub = np.asarray(b_sub, int)[mask_sub]
                    nchan = int(f_m_sub.size)
                    if nchan < 3:
                        continue

                    # subset-specific delay window
                    delay_window_ns_sub = get_delay_window_ns(f_m_sub)

                    tau_hat_s_sub, f_zero_hz_sub, taus_sub, R_sub, ph_u_sub, phi0_sub, cand_sub = estimate_delay_fft_gridded(
                        f_m_sub, ph_m_sub, amp=amp_m_sub,
                        padding=padding,
                        band_id_per_sample=b_m_sub,
                        BAND_GAP_HZ=BAND_GAP_HZ,
                        use_amp_weight=False,
                        fringe_peak_frac=FRINGE_PEAK_FRAC,
                        delay_window_ns=delay_window_ns_sub,
                        enable_candidates=True,
                        refine_method=refine_method,
                        apply_freq_window=apply_freq_window,
                        window_kind=fft_window_kind,
                    )

                    # Compute subset quality weight from its spectrum
                    fringe_sub = np.abs(R_sub)
                    p1, p2, ratio, snr_like, n_cand = spectrum_quality(fringe_sub, cand_sub)
                    W_subset = subset_weight_from_quality(p1, ratio, snr_like)
                    if W_subset <= 0:
                        continue

                    w_cands = distribute_subset_weight(cand_sub, W_subset, gamma=1.0)

                    for cs, wv in zip(cand_sub, w_cands):
                        if wv <= 0:
                            continue
                        cand_pool.append({
                            "tau_ns": float(cs["tau_ns"]),
                            "tau_s":  float(cs["tau_s"]),
                            "weight": float(wv),
                            "source": bands_label,
                        })
                #        print(
                #            "[%s] subset=%7s  " % (ant, bands_label)
                #            "cand_tau=%9.3f ns  fringe=%.3f  " % (cs["tau_ns"], cs["fringe"])
                #            "w_vote=%.3g  W_subset=%.3g  " % (wv, W_subset)
                #            "p1=%.3f p2=%.3f ratio=%.2f snr=%.2f n_cand=%d" % (p1, p2, ratio, snr_like, n_cand)
                #        )

                    # log + plot
                    tag = bands_label.replace("+", "")
                    plot_delay_spectrum(
                        taus_sub, R_sub, tau_hat_s_sub,
                        os.path.join(combo_plot_dir, ("XPOL_%s_B%s_delay_spectrum.png" % (ant, tag))),
                        title_prefix=("%s B%s:" % (ant, bands_label)),
                        candidates=cand_sub,
                    )

                    # write combo candidates
                    for j, cs in enumerate(cand_sub, 1):
                        fcombo.write("%s  %s  %d  %.9e  %.6f  %.9e\n" % (ant, bands_label, j, cs["tau_s"], cs["tau_ns"], cs["f_zero_hz"]))

                # vote clustering
                clusters = cluster_candidates_by_tau_best(cand_pool, tol_ns=tol_ns)

                final_ranked = []
                for cl in clusters:
                    tau_s_test = cl["tau_s"]
                    coh = delay_coherence(f_m, ph_m, amp_m, tau_s_test)
                    if not np.isfinite(coh):
                        coh = 0.0

                    # normalize away vote count domination
                    alpha = 0.5
                    weight_eff = cl["weight"] / (max(1, cl["n"]) ** alpha)

                    score = float(coh * weight_eff)

                    final_ranked.append({
                        "tau_s": float(tau_s_test),
                        "tau_ns": float(cl["tau_ns"]),
                        "score": score,
                        "coh": float(coh),
                        "weight": float(cl["weight"]),
                        "n_votes": int(cl["n"]),
                        "sources": cl["sources"],
                        "best_vote_weight": float(cl["best_vote_weight"]),
                    })

                final_ranked.sort(key=lambda d: d["score"], reverse=True)
                final_ranked = final_ranked[:MAX_CANDIDATES]

                idxs = np.argsort(f_m)
                f_sorted = f_m[idxs]
                ph_sorted = ph_m[idxs]
                w_sorted = amp_m[idxs]
                b_sorted = b_m[idxs]
                ph_u_sorted = band_aware_unwrap(f_sorted, ph_sorted, b_sorted)

                for d in final_ranked:
                    phi0_hat, f_zero_hat = compute_phi0_and_fzero(f_sorted, ph_u_sorted, w_sorted, d["tau_s"])
                    d["f_zero_hz"] = float(f_zero_hat)

                if final_ranked:
                    results[ant]["tau_hat_s"] = final_ranked[0]["tau_s"]
                    results[ant]["f_zero_hz"] = final_ranked[0]["f_zero_hz"]
                    results[ant]["candidates"] = final_ranked

                if final_ranked:
                    for j, d in enumerate(final_ranked, 1):
                        src = ",".join(sorted(d["sources"]))[:120]  # avoid insane line lengths
                        fout.write("%s  %02d  %.9e  %.6f  %.9e  %.6e  %.4f  %.3f  %d  %s\n" % (ant, j, d["tau_s"], d["tau_ns"], d["f_zero_hz"], d["score"], d["coh"], d["weight"], d["n_votes"], src))
                else:
                    # fallback: at least write something
                    fout.write("%s  01  %.9e  %.6f  %.9e  -1  -1  -1  0  NONE\n" % (ant, tau_hat_s, tau_hat_s*1e9, f_zero_hz))

        #    print("\n>>> Finished antenna %s. Press ENTER to continue to the next antenna..." % (ant,))
        #    input()
    print("\nWrote results to: %s" % (results_txt,))
    print("Plots in: %s" % (xpol_dir,))
    if TRY_BAND_COMBINATIONS:
        print("Wrote band-combo results to: %s" % (combo_txt,))

    return results




def main():
    directory = os.getcwd()
    fname = "POLCAL_GAINS_h.dat"
    path = os.path.join(directory, fname)
    
    EXPLORE_CANDIDATES = True
    FRINGE_PEAK_FRAC = 0.9
    MAX_CANDIDATES = 5

    FLAG_EDGES     = True
    EDGE_FRAC      = 0.025   # flag X% at start AND X% at end of each IF
    FLAG_NOISY_IFS = False
    TRY_BAND_COMBINATIONS = False
    
    FFT_WINDOW=True
    FFT_WINDOW_KIND="hann"
    REFINE_METHOD="exact"
    
    ## If problems identifying IFs corresponding to the different bands
    ##    directly specify them as IF_TO_BAND, e.g.
    # if_to_band = {
    #                        #4 bands, 8 IFs each, VGOS config.
    #                        **{IF: 1 for IF in range(1, 9)},
    #                        **{IF: 2 for IF in range(9, 17)},
    #                        **{IF: 3 for IF in range(17, 25)},
    #                        **{IF: 4 for IF in range(25, 33)},
    #                       }
    if_to_band = None

    CROSS_POL_DELAY_SCANNER(
        path,
        out_dir=os.path.join(directory, "XPOL_DELAY_INSPECT_test"),
        EXPLORE_CANDIDATES=EXPLORE_CANDIDATES,
        FRINGE_PEAK_FRAC=FRINGE_PEAK_FRAC,
        MAX_CANDIDATES=MAX_CANDIDATES,
        BAND_GAP_HZ = 0.9e9,
        FLAG_EDGES=FLAG_EDGES,
        EDGE_FRAC=EDGE_FRAC,
        FLAG_NOISY_IFS=FLAG_NOISY_IFS,
        padding=32,          # zero-padding
        TRY_BAND_COMBINATIONS=TRY_BAND_COMBINATIONS,
        IF_TO_BAND=if_to_band,
        refine_method=REFINE_METHOD,
        apply_freq_window=FFT_WINDOW,
        fft_window_kind=FFT_WINDOW_KIND
    )

if __name__ == "__main__":
    main()

