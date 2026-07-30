#!/usr/bin/env python3
"""
Synthetic model generator for checkerboard validation.

One place that builds a known model and turns it into what two overlapping
subarray inversions would recover from it, so the manuscript figures and the
user-facing checkerboard test are driven by exactly the same code.

The recovery model is deliberately simple but has the two effects that matter
for a merge: away from an array centre, damping and smoothing SUPPRESS anomaly
amplitude and BROADEN anomalies (Barmin et al., 2001; Rawlinson & Spakman,
2016).  Both are produced here by mixing the true field with a blurred copy of
itself, with the mixing weight set by distance from the array centre.
"""
import numpy as np


# ----------------------------------------------------------------------------
# the true model
# ----------------------------------------------------------------------------
def truth_grid(lonv, latv, cell=3.0, amp=0.12, trend=True, pattern="smooth",
               base=3.85):
    """
    A checkerboard of half-wavelength `cell` degrees and amplitude `amp`,
    optionally on a smooth regional trend.

    pattern='smooth' gives a sinusoidal product, which is what a real velocity
    field looks like and what the manuscript uses.  pattern='blocks' gives sharp
    alternating squares, which is the harsher classical checkerboard test: it
    contains high wavenumbers no inversion can recover, so recovered amplitudes
    come out lower.

    Returns a 2-D array indexed [lon, lat].
    """
    X, Y = np.meshgrid(np.asarray(lonv, float), np.asarray(latv, float),
                       indexing="ij")
    board = np.sin(np.pi*(X - X.min())/cell) * np.sin(np.pi*(Y - Y.min())/cell)
    if pattern == "blocks":
        board = np.sign(board)
    elif pattern != "smooth":
        raise ValueError("pattern must be 'smooth' or 'blocks'")
    field = base + amp*board
    if trend:
        span_x = max(X.max() - X.min(), 1e-9)
        span_y = max(Y.max() - Y.min(), 1e-9)
        field = (field
                 + 0.10*np.sin(2*np.pi*(X - X.min())/span_x)
                 - 0.06*np.cos(2*np.pi*(Y - Y.min())/span_y))
    return field


# ----------------------------------------------------------------------------
# what one subarray recovers
# ----------------------------------------------------------------------------
def _gauss_kernel(sigma_nodes):
    r = max(1, int(np.ceil(3*sigma_nodes)))
    x = np.arange(-r, r + 1, dtype=float)
    k = np.exp(-0.5*(x/sigma_nodes)**2)
    return k/k.sum()


def blur2d(F, sigma_x_nodes, sigma_y_nodes):
    """Separable Gaussian blur with edge reflection, numpy only."""
    out = np.asarray(F, float)
    for axis, sig in ((0, sigma_x_nodes), (1, sigma_y_nodes)):
        if sig <= 0:
            continue
        k = _gauss_kernel(sig)
        r = (len(k) - 1)//2
        out = np.apply_along_axis(
            lambda v: np.convolve(np.r_[v[r:0:-1], v, v[-2:-2-r:-1]], k, "valid"),
            axis, out)
    return out


def subarray_recovery(lonv, latv, truth, centre, reach=7.0, blur=1.4,
                      offset=0.0, noise=0.0, seed=0):
    """
    Degrade a true field into one subarray's recovered model.

    reach  : degrees.  Amplitude recovery is exp(-(d/reach)^2), so it is ~1 at
             the array centre and falls off outward.
    blur   : degrees.  The width of the smoothing that replaces the true field
             where recovery is poor, which both broadens and damps anomalies.
    offset : a constant added to the whole model, standing in for a different
             reference level between two independent inversions.
    noise  : standard deviation of independent Gaussian noise.

    Returns a 2-D array the same shape as `truth`.
    """
    lonv = np.asarray(lonv, float); latv = np.asarray(latv, float)
    dlon = np.median(np.diff(lonv)) if len(lonv) > 1 else 1.0
    dlat = np.median(np.diff(latv)) if len(latv) > 1 else 1.0
    X, Y = np.meshgrid(lonv, latv, indexing="ij")

    smoothed = blur2d(truth, blur/abs(dlon), blur/abs(dlat))
    d = np.hypot(X - centre[0], Y - centre[1])
    gain = np.exp(-(d/reach)**2)                    # 1 at the centre, →0 outward
    rec = gain*truth + (1.0 - gain)*smoothed
    if noise:
        rec = rec + noise*np.random.default_rng(seed).standard_normal(rec.shape)
    return rec + offset


def to_xyz(lonv, latv, F):
    """Flatten a 2-D grid to the (lon, lat, value) triple the tool reads."""
    X, Y = np.meshgrid(np.asarray(lonv, float), np.asarray(latv, float),
                       indexing="ij")
    return X.ravel(), Y.ravel(), np.asarray(F, float).ravel()


def save_xyz(path, lonv, latv, F):
    x, y, v = to_xyz(lonv, latv, F)
    np.savetxt(path, np.c_[x, y, v], fmt="%12.6f %12.7f %14.8f")
    return path


def checkerboard_recovery(lonv, latv, truth, rec, cell, at_lon, window=None):
    """
    Fraction of the checkerboard amplitude a recovered field retains near a
    given longitude.  1.0 is perfect recovery; a damped inversion returns less.

    Two things have to be handled or the number is meaningless.  The regional
    trend is removed from both fields first, with a high-pass whose width is the
    checkerboard size - otherwise the trend, which survives smoothing, dominates
    the estimate.  And the slope is measured over a longitude BAND rather than a
    single column, because a column that happens to land on a node of the
    checkerboard contains no anomaly at all, and the ratio there is a division
    of noise by noise.
    """
    lonv = np.asarray(lonv, float); latv = np.asarray(latv, float)
    dlon = np.median(np.diff(lonv)) if len(lonv) > 1 else 1.0
    dlat = np.median(np.diff(latv)) if len(latv) > 1 else 1.0

    def highpass(F):
        F = np.asarray(F, float)
        return F - blur2d(F, cell/abs(dlon), cell/abs(dlat))

    t, r = highpass(truth), highpass(rec)
    w = window if window else max(cell, 1.5)
    m = np.abs(lonv - at_lon) <= w
    if m.sum() < 2:
        m = np.ones(len(lonv), bool)
    tt, rr = t[m, :].ravel(), r[m, :].ravel()
    return float((tt @ rr)/(tt @ tt)) if tt @ tt > 0 else float("nan")
