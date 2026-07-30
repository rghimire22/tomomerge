#!/usr/bin/env python3
"""
merge_tomo.py -- merge two overlapping tomographic (phase-velocity) maps.

Reads two 3-column ASCII files (lon lat vel) from two subarray inversions that
share an overlap zone.  Nodes present in only one map are copied through
untouched.  Nodes present in BOTH maps are replaced by a Gaussian
distance-weighted average that trusts whichever subarray centre is closer:

    w_W = exp(-(d_W / L)^2)          d_W = |lon - centre_W|
    w_E = exp(-(d_E / L)^2)          d_E = |lon - centre_E|
    v   = (w_W*v_W + w_E*v_E) / (w_W + w_E)

This is a pure blend: no bias is removed from either map, so any real
disagreement between the two inversions stays visible in the output.

WHAT L ACTUALLY DOES
--------------------
The normalised weight is an exact logistic in longitude:

    w_W(lon) = 1 / (1 + exp((lon - m) / s))
    m = (cW + cE) / 2            (midpoint between the two centres)
    s = L^2 / (2 * |cE - cW|)    (hand-over scale)

so the hand-over is governed by L^2/(2*sep), NOT by L.  With centres 9 deg
apart, L = 4 gives s = 0.89 deg, not 4 deg.  The 10-90 % transition width is
2*ln(9)*s = 4.39*s.

The blend must reach ~pure-west at the western edge of the overlap and
~pure-east at the eastern edge, otherwise the merged map steps discontinuously
where the overlap ends.  The script computes the largest L that keeps both
overlap edges at >=90 % of the map that continues beyond them, offers it as the
default, and reports the residual seam amplitude in velocity units so the
choice is checkable rather than cosmetic.

Run with no arguments for interactive prompts, or pass flags:

    python merge_tomo.py --west W.xyz --east E.xyz -L 4.0 -o merged.xyz

Requires numpy; matplotlib only for the optional QC plot.
"""

import argparse
import os
import sys

import numpy as np

# Grid coordinates are stored as float32 by most tomography codes, and two
# inversions that nominally share a node can print it as -113.800003 and
# -113.800000.  Nodes are therefore MATCHED on coordinates rounded to
# KEY_DECIMALS (1e-4 deg ~ 11 m: far below any grid spacing, far above float32
# noise) but WRITTEN OUT with their original coordinates, so single-map nodes
# round-trip unchanged.
KEY_DECIMALS = 4

EDGE_TOL = 0.10   # target: overlap edges are >=90 % the map that continues on


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
def read_xyz(path):
    """Read a lon/lat/value file, skipping blank lines and # or > comments."""
    rows = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith(("#", ">")):
                continue
            p = s.split()
            if len(p) < 3:
                raise ValueError(f"{path}:{ln}: expected 3 columns, got {len(p)}")
            try:
                rows.append((float(p[0]), float(p[1]), float(p[2])))
            except ValueError:
                raise ValueError(f"{path}:{ln}: non-numeric value -> {s!r}")
    if not rows:
        raise ValueError(f"{path}: no data rows found")
    a = np.asarray(rows, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{path}: contains NaN/Inf")
    return a[:, 0], a[:, 1], a[:, 2]


MATCH_SNAP = 10 ** -KEY_DECIMALS      # nodes closer than this count as the same


def set_match_tol(deg):
    """Widen or tighten node matching (degrees).  See --match-tol."""
    global MATCH_SNAP
    if not deg > 0:
        raise ValueError("match tolerance must be > 0")
    MATCH_SNAP = float(deg)




# ----------------------------------------------------------------------------
# Blending kernel
# ----------------------------------------------------------------------------
def _centre(c):
    """A centre may be a bare longitude or an (lon, lat) pair."""
    if np.ndim(c) == 0:
        return float(c), 0.0
    return float(c[0]), float(c[1])


def project(lon, lat, cA, cB):
    """
    Signed distance from the midpoint of the two centres, measured along the
    line joining them (positive towards cB), plus the centre separation D.

    This is the only geometry the blend needs.  A difference of two squared
    distances is linear in position:

        d_B^2 - d_A^2 = 2*D*t

    so log(w_A/w_B) = 2*D*t/L^2 and the weight is an exact logistic in t for
    ANY orientation - east-west, north-south or diagonal.  When the two centres
    share a latitude, t is just (lon - midpoint longitude) and everything
    reduces bit-for-bit to the longitude-only formula this started as.
    """
    ax, ay = _centre(cA)
    bx, by = _centre(cB)
    dx, dy = bx - ax, by - ay
    D = float(np.hypot(dx, dy))
    if D == 0.0:
        return np.zeros(np.shape(lon)), 0.0
    ux, uy = dx / D, dy / D
    mx, my = 0.5 * (ax + bx), 0.5 * (ay + by)
    t = ((np.asarray(lon, dtype=float) - mx) * ux
         + (np.asarray(lat, dtype=float) - my) * uy)
    return t, D


def perpendicular(lon, lat, cA, cB):
    """Signed offset from the centre-to-centre line (companion to project)."""
    ax, ay = _centre(cA)
    bx, by = _centre(cB)
    dx, dy = bx - ax, by - ay
    D = float(np.hypot(dx, dy))
    if D == 0.0:
        return np.zeros(np.shape(lon))
    ux, uy = dx / D, dy / D
    mx, my = 0.5 * (ax + bx), 0.5 * (ay + by)
    return (-(np.asarray(lon, dtype=float) - mx) * uy
            + (np.asarray(lat, dtype=float) - my) * ux)


def blend_weight(lon, lat, L, cA, cB):
    """Normalised weight of map A.  Exact logistic, overflow-safe for any L."""
    t, D = project(lon, lat, cA, cB)
    if D == 0.0:
        return np.full(np.shape(t), 0.5)
    z = t / (L * L / (2.0 * D))
    # 1/(1+exp(z)) evaluated as exp(-|z|)/(1+exp(-|z|)) for z>=0 and
    # 1/(1+exp(-|z|)) for z<0: only exp of a non-positive number is ever taken,
    # so this cannot overflow however small L is.
    ez = np.exp(-np.abs(z))
    return np.where(z >= 0, ez / (1.0 + ez), 1.0 / (1.0 + ez))


def west_weight(lon, L, cW, cE):
    """Longitude-only form, kept for the east/west case and old callers."""
    return blend_weight(lon, 0.0, L, cW, cE)


def handover_scale(L, cA, cB):
    _, D = project(0.0, 0.0, cA, cB)
    return float("inf") if D == 0 else L * L / (2.0 * D)


def suggest_L_proj(cA, cB, t_lo, t_hi, p=EDGE_TOL):
    """
    Largest L for which w_A >= 1-p at the near edge of the overlap and <= p at
    the far edge, with the overlap given as a projected extent (see project).

        s <= min(-t_lo, t_hi) / ln((1-p)/p),     L = sqrt(2*D*s)

    None means the midpoint between the centres falls outside the overlap, so
    no L keeps both edges clean and the centres need revisiting.
    """
    _, D = project(0.0, 0.0, cA, cB)
    if D == 0 or t_hi <= t_lo:
        return None
    half = min(-t_lo, t_hi)
    if half <= 0:
        return None
    s = half / np.log((1.0 - p) / p)
    return float(np.sqrt(2.0 * D * s))


def suggest_L(cW, cE, ov_lo, ov_hi, p=EDGE_TOL):
    """Longitude-interval form of suggest_L_proj, for the east/west case."""
    m = 0.5 * (_centre(cW)[0] + _centre(cE)[0])
    return suggest_L_proj(cW, cE, ov_lo - m, ov_hi - m, p)


def _spacing(a):
    """Median node spacing along one axis, used to size plot symbols/bands."""
    u = np.unique(np.round(np.asarray(a, dtype=float) / MATCH_SNAP))
    return float(np.median(np.diff(u)) * MATCH_SNAP) if u.size > 1 else 1.0


def alignment(west, east):
    """
    How well the two grids line up, independent of the match tolerance.

    Returns (n_coincident, worst_gap, pitch): how many nodes of `west` have a
    node of `east` essentially on top of them (within half a node pitch), how
    far apart the worst of those pairs is, and the node pitch itself.

    Comparing n_coincident with the number that actually matched catches the
    worst failure mode of all - a PARTIAL match, where some of the overlap
    blends and the rest quietly does not, leaving scattered holes that look
    like real structure on the finished map.
    """
    pitch = min(_pitch(west[0], west[1]), _pitch(east[0], east[1]))
    pair = match_nodes(west[0], west[1], east[0], east[1], 0.5 * pitch)
    hit = pair >= 0
    if not hit.any():
        return 0, 0.0, pitch
    d = np.hypot(np.asarray(east[0])[pair[hit]] - np.asarray(west[0])[hit],
                 np.asarray(east[1])[pair[hit]] - np.asarray(west[1])[hit])
    # match_nodes accepts a pair exactly ON the radius; "coincident" has to be
    # strictly inside it, or a grid offset by precisely half a pitch would count
    # as lined up.  The browser build uses the same strict test.
    d = d[d < 0.5 * pitch]
    if not len(d):
        return 0, 0.0, pitch
    return int(len(d)), float(d.max()), pitch


def _pitch(lon, lat):
    """Median nearest-neighbour distance inside one map: the real node pitch."""
    lon = np.asarray(lon, dtype=float); lat = np.asarray(lat, dtype=float)
    n = len(lon)
    if n < 2:
        return 1.0
    w, h = lon.max() - lon.min(), lat.max() - lat.min()
    guess = np.sqrt(w * h / n) if w > 0 and h > 0 else max(w, h) / (n - 1)
    if not guess > 0:
        return 1.0
    for attempt in range(5):
        cs = guess * 1.5 * 2 ** attempt
        cell = {}
        for j in range(n):
            cell.setdefault((int(np.floor(lon[j]/cs)), int(np.floor(lat[j]/cs))),
                            []).append(j)
        stride = max(1, n // 2000)
        got = []
        for i in range(0, n, stride):
            cx, cy = int(np.floor(lon[i]/cs)), int(np.floor(lat[i]/cs))
            bd = np.inf
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for j in cell.get((cx+dx, cy+dy), ()):
                        if j == i:
                            continue
                        d = (lon[j]-lon[i])**2 + (lat[j]-lat[i])**2
                        if d < bd:
                            bd = d
            if np.isfinite(bd):
                got.append(np.sqrt(bd))
        if len(got) > 0.5 * len(range(0, n, stride)):
            return float(np.median(got))
    return float(guess)


def _dedupe(lon, lat):
    """Indices of the first occurrence of each exactly repeated coordinate."""
    seen, keep = set(), []
    for i, (x, y) in enumerate(zip(lon, lat)):
        k = (float(x), float(y))
        if k not in seen:
            seen.add(k); keep.append(i)
    return np.asarray(keep, dtype=int)


def match_nodes(lonA, latA, lonB, latB, tol=None):
    """
    For each node of A, the index of the nearest node of B within `tol`
    degrees, or -1.  Ties go to the lower index, so the result is deterministic.

    A hash grid of cell size `tol` with a 3x3 neighbour scan.  Snapping both
    maps onto a common grid and comparing the buckets - the obvious shortcut -
    silently loses any pair that straddles a bucket boundary, which punches
    scattered holes in the overlap that no diagnostic would catch.  Searching a
    true radius cannot do that.
    """
    tol = MATCH_SNAP if tol is None else float(tol)
    cell = {}
    for j, (x, y) in enumerate(zip(lonB, latB)):
        cell.setdefault((int(np.floor(x / tol)), int(np.floor(y / tol))), []).append(j)
    out = np.full(len(lonA), -1, dtype=int)
    t2 = tol * tol
    for i, (x, y) in enumerate(zip(lonA, latA)):
        cx, cy = int(np.floor(x / tol)), int(np.floor(y / tol))
        best, bd = -1, float("inf")
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cell.get((cx + dx, cy + dy), ()):
                    d = (lonB[j] - x) ** 2 + (latB[j] - y) ** 2
                    if d < bd or (d == bd and (best < 0 or j < best)):
                        bd, best = d, j
        out[i] = best if bd <= t2 else -1
    return out


def overlap_projection(west, east, cA, cB, tol=None):
    """
    Projected extent of the nodes the two maps actually share.

    Measuring the overlap from the nodes that match, rather than from where the
    two bounding boxes intersect, is what makes the recommended L correct for
    ragged or non-rectangular coverage.
    """
    pair = match_nodes(west[0], west[1], east[0], east[1], tol)
    hit = pair >= 0
    if not hit.any():
        return None, None, 0
    t, _ = project(np.asarray(west[0])[hit], np.asarray(west[1])[hit], cA, cB)
    return float(t.min()), float(t.max()), int(hit.sum())


# ----------------------------------------------------------------------------
# Levelling
# ----------------------------------------------------------------------------
def level_maps(west, east, mode="constant", to="split", tol=None):
    """
    Remove the systematic part of the disagreement before blending.

    Two overlapping inversions frequently differ by a near-uniform amount, from
    different damping, reference models or data windows.  Blending such a pair
    converts that offset into a smooth ramp right across the hand-over, which is
    indistinguishable from a real velocity gradient.  Estimating the offset on
    the shared nodes and removing it first means the blend acts on structure.

    mode : 'constant' (one number) or 'plane' (a0 + a1*lon + a2*lat)
    to   : 'split' puts half the correction in each map, 'A'/'B' put it all in one

    The RELATIVE offset is measurable; the absolute level is not, so whichever
    choice is made the correction has to be reported, not hidden.  Note that
    levelling necessarily changes nodes only one map covers - the offset has to
    go somewhere, and the alternative is to leave it as a ramp in the overlap.

    Returns (west2, east2, info).
    """
    if mode not in ("none", "constant", "plane"):
        raise ValueError("level mode must be none, constant or plane")
    if mode == "none":
        return west, east, {"mode": "none"}

    pair = match_nodes(west[0], west[1], east[0], east[1], tol)
    hit = pair >= 0
    n = int(hit.sum())
    if n < (3 if mode == "plane" else 1):
        return west, east, {"mode": mode, "n": n,
                            "note": "too few shared nodes to estimate an offset"}

    lo = np.asarray(west[0])[hit]
    la = np.asarray(west[1])[hit]
    d = np.asarray(east[2])[pair[hit]] - np.asarray(west[2])[hit]      # B - A

    if mode == "constant":
        coef = np.array([d.mean()])
        model = lambda x, y: np.full(np.shape(x), coef[0])
    else:
        G = np.column_stack([np.ones_like(lo), lo - lo.mean(), la - la.mean()])
        coef, *_ = np.linalg.lstsq(G, d, rcond=None)
        lo0, la0 = lo.mean(), la.mean()
        model = lambda x, y: coef[0] + coef[1]*(np.asarray(x)-lo0) + coef[2]*(np.asarray(y)-la0)

    fA = {"split": 0.5, "A": 1.0, "B": 0.0}[to]
    corr = model(np.asarray(east[0]), np.asarray(east[1]))
    east2 = (east[0], east[1], np.asarray(east[2]) - (1.0 - fA) * corr)
    corrW = model(np.asarray(west[0]), np.asarray(west[1]))
    west2 = (west[0], west[1], np.asarray(west[2]) + fA * corrW)

    resid = d - model(lo, la)
    info = {"mode": mode, "to": to, "n": n,
            "offset_mean": float(d.mean()),
            "removed_mean": float(model(lo, la).mean()),
            "coef": [float(c) for c in coef],
            "rms_before": float(np.sqrt((d**2).mean())),
            "rms_after": float(np.sqrt((resid**2).mean())),
            "structure_rms": float(resid.std())}
    return west2, east2, info


# ----------------------------------------------------------------------------
# Core merge
# ----------------------------------------------------------------------------
def gaussian_merge(west, east, L, cW, cE):
    """
    west, east : (lon, lat, vel) tuples of 1-D arrays
    L          : smoothing length, degrees
    cW, cE     : longitude of the western / eastern subarray centre

    Returns (lon, lat, vel, wW, source): west nodes in their original file
    order, then east-only nodes.  `source` is 'W', 'E' or 'B' (blended); `wW`
    is the normalised western weight (1 for W-only, 0 for E-only nodes).
    Emitted coordinates are the originals, not the rounded matching keys.
    """
    if L <= 0:
        raise ValueError("smoothing length L must be > 0")

    iW = _dedupe(west[0], west[1])
    iE = _dedupe(east[0], east[1])
    lonW, latW, velW = (np.asarray(a)[iW] for a in west)
    lonE, latE, velE = (np.asarray(a)[iE] for a in east)
    dupW, dupE = len(west[0]) - len(iW), len(east[0]) - len(iE)
    if dupW or dupE:
        print(f"  ! warning: repeated nodes inside a single file "
              f"(west {dupW}, east {dupE}); first occurrence kept")

    pair = match_nodes(lonW, latW, lonE, latE)
    claimed = np.zeros(len(lonE), dtype=bool)
    claimed[pair[pair >= 0]] = True

    out_lon, out_lat, out_vel, out_wW, out_src = [], [], [], [], []

    for i in range(len(lonW)):                           # west file order
        lo, la, j = lonW[i], latW[i], pair[i]
        if j >= 0:
            wWn = float(blend_weight(lo, la, L, cW, cE))
            if not np.isfinite(wWn):                     # pathological L
                ax, ay = _centre(cW); bx, by = _centre(cE)
                wWn = 1.0 if ((lo-ax)**2 + (la-ay)**2) <= ((lo-bx)**2 + (la-by)**2) else 0.0
            v, s = wWn * velW[i] + (1.0 - wWn) * velE[j], "B"
        else:
            v, wWn, s = velW[i], 1.0, "W"
        out_lon.append(lo); out_lat.append(la)
        out_vel.append(v);  out_wW.append(wWn); out_src.append(s)

    for j in range(len(lonE)):                           # then east-only nodes
        if claimed[j]:
            continue
        out_lon.append(lonE[j]); out_lat.append(latE[j])
        out_vel.append(velE[j]); out_wW.append(0.0); out_src.append("E")

    return (np.array(out_lon), np.array(out_lat), np.array(out_vel),
            np.array(out_wW), np.array(out_src))


# ----------------------------------------------------------------------------
# Diagnostics
# ----------------------------------------------------------------------------
def overlap_diff(west, east, lon, lat, src):
    """east-minus-west at every blended node, in output order."""
    ov = src == "B"
    if not ov.any():
        return ov, np.array([])
    pw = match_nodes(lon[ov], lat[ov], west[0], west[1])
    pe = match_nodes(lon[ov], lat[ov], east[0], east[1])
    return ov, np.array([east[2][b] - west[2][a] for a, b in zip(pw, pe)])


def seam_report(lon, wW, ov, diff):
    """
    Amplitude of the step left at each edge of the overlap.

    At the western edge the map continuing west is pure west, while the merged
    value there is w_W*vW + (1-w_W)*vE, so the offset is (1-w_W)*(vE - vW).
    Mirror image at the eastern edge.  This is the seam a plotted map shows.
    """
    if not ov.any():
        return None
    lo_ov, w_ov = lon[ov], wW[ov]
    tol = 10 ** -KEY_DECIMALS
    lo_edge, hi_edge = lo_ov.min(), lo_ov.max()
    wc = np.abs(lo_ov - lo_edge) < tol
    ec = np.abs(lo_ov - hi_edge) < tol
    step_w = (1.0 - w_ov[wc]) * diff[wc]                 # west edge
    step_e = w_ov[ec] * (-diff[ec])                      # east edge
    return dict(lo_edge=lo_edge, hi_edge=hi_edge,
                w_west=float(w_ov[wc].mean()), w_east=float(w_ov[ec].mean()),
                rms_w=float(np.sqrt((step_w ** 2).mean())),
                max_w=float(np.abs(step_w).max()),
                rms_e=float(np.sqrt((step_e ** 2).mean())),
                max_e=float(np.abs(step_e).max()))


# ----------------------------------------------------------------------------
# QC plot
# ----------------------------------------------------------------------------
def qc_plot(west, east, merged, wW, src, diff, ov, L, cW, cE, png):
    # Built through Figure/FigureCanvasAgg rather than pyplot: no global state,
    # no matplotlib.use() call, so this is safe to run on a worker thread and
    # cannot disturb a GUI that already owns a live canvas.
    try:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
    except ImportError:
        print("  ! matplotlib not installed - skipping QC plot")
        return

    lonW, latW, velW = west
    lonE, latE, velE = east
    lonM, latM, velM = merged
    vmin = min(velW.min(), velE.min())
    vmax = max(velW.max(), velE.max())
    asp = 1 / np.cos(np.deg2rad(np.mean(latM)))

    fig = Figure(figsize=(13, 9))
    FigureCanvasAgg(fig)
    gs = fig.add_gridspec(3, 3, height_ratios=[1.35, 1.35, 1.0],
                          hspace=0.45, wspace=0.42)

    def frame(ax):
        ax.set_xlim(lonM.min() - .2, lonM.max() + .2)
        ax.set_ylim(latM.min() - .2, latM.max() + .2)
        ax.set_aspect(asp)
        ax.set_xlabel("lon"); ax.set_ylabel("lat")

    def panel(ax, lo, la, v, title):
        sc = ax.scatter(lo, la, c=v, s=3, marker="s", cmap="jet_r",
                        vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=10); frame(ax)
        return sc

    for col, (lo, la, v, t) in enumerate(((lonW, latW, velW, "West input"),
                                          (lonE, latE, velE, "East input"),
                                          (lonM, latM, velM, "Merged"))):
        ax = fig.add_subplot(gs[0, col])
        fig.colorbar(panel(ax, lo, la, v, t), ax=ax, shrink=.85, label="vel")

    # normalised western weight
    ax = fig.add_subplot(gs[1, 0])
    s4 = ax.scatter(lonM, latM, c=wW, s=3, marker="s", cmap="coolwarm_r",
                    vmin=0, vmax=1)
    ax.set_title("normalised west weight  $w_W$", fontsize=10); frame(ax)
    fig.colorbar(s4, ax=ax, shrink=.85)

    # east-minus-west disagreement inside the overlap
    ax = fig.add_subplot(gs[1, 1])
    if ov.any():
        m = np.abs(diff).max()
        s5 = ax.scatter(lonM[ov], latM[ov], c=diff, s=4, marker="s",
                        cmap="RdBu_r", vmin=-m, vmax=m)
        fig.colorbar(s5, ax=ax, shrink=.85, label="E - W")
        ax.set_title("overlap disagreement (east - west)", fontsize=10)
    else:
        ax.text(.5, .5, "no overlap nodes", ha="center", transform=ax.transAxes)
    frame(ax)

    # weight curve along the blend axis (= longitude when the split is E-W)
    ax = fig.add_subplot(gs[1, 2])
    tM, D = project(lonM, latM, cW, cE)
    ts = np.linspace(tM.min(), tM.max(), 800)
    s_ = handover_scale(L, cW, cE)
    ww = 1.0/(1.0 + np.exp(np.clip(ts/s_, -700, 700)))
    ax.plot(ts, ww, label="$w_A$")
    ax.plot(ts, 1 - ww, label="$w_B$")
    if ov.any():
        ax.axvspan(tM[ov].min(), tM[ov].max(), color="k", alpha=.10,
                   label="overlap")
    ax.axhline(EDGE_TOL, ls="--", lw=.8, c="grey")
    ax.axhline(1 - EDGE_TOL, ls="--", lw=.8, c="grey")
    ax.axvline(-D/2, ls=":", c="C0"); ax.axvline(D/2, ls=":", c="C1")
    ax.set_xlabel("distance along blend axis (deg from midpoint)")
    ax.set_ylabel("normalised weight")
    ax.set_title(f"kernel: L = {L:g}$^\\circ$, scale {s_:.2f}$^\\circ$",
                 fontsize=10)
    ax.legend(fontsize=8); ax.set_ylim(-.05, 1.05); ax.grid(alpha=.3)

    # profile along the blend axis through the midpoint: the seam test
    ax = fig.add_subplot(gs[2, :])
    band = 0.6*max(_spacing(lonM), _spacing(latM))
    for lo, la, v, lab, sty in ((lonW, latW, velW, "map A", "--"),
                                (lonE, latE, velE, "map B", "--"),
                                (lonM, latM, velM, "merged", "-")):
        t = project(lo, la, cW, cE)[0]
        n = perpendicular(lo, la, cW, cE)
        m = np.abs(n) < band
        if not m.any():
            continue
        o = np.argsort(t[m])
        ax.plot(t[m][o], v[m][o], sty, lw=2.2 if lab == "merged" else 1.2,
                label=lab, color="k" if lab == "merged" else None)
    if ov.any():
        ax.axvspan(tM[ov].min(), tM[ov].max(), color="k", alpha=.10)
        for e in (tM[ov].min(), tM[ov].max()):
            ax.axvline(e, lw=.8, c="k", alpha=.5)
    ax.set_xlabel("distance along blend axis (deg from midpoint)")
    ax.set_ylabel("velocity")
    ax.set_title("profile through the centre line   (shaded = overlap; look "
                 "for steps at its edges)", fontsize=10)
    ax.legend(fontsize=9); ax.grid(alpha=.3)

    fig.savefig(png, dpi=140, bbox_inches="tight")
    print(f"  QC plot   -> {png}")


# ----------------------------------------------------------------------------
# Prompt helpers
# ----------------------------------------------------------------------------
def ask(prompt, default=None, cast=str, check=None):
    # Non-interactive run (piped stdin, cron, batch): never block on input.
    if not sys.stdin.isatty():
        if default is None:
            sys.exit(f"error: '{prompt}' is needed but stdin is not a terminal; "
                     f"pass it as a command-line flag")
        return default
    while True:
        d = f" [{default}]" if default is not None else ""
        try:
            s = input(f"{prompt}{d}: ").strip()
        except EOFError:
            if default is None:
                sys.exit("\naborted")
            return default
        if not s:
            if default is None:
                print("  -> required"); continue
            return default
        try:
            v = cast(s)
        except Exception as e:
            print(f"  -> {e}"); continue
        if check:
            msg = check(v)
            if msg:
                print(f"  -> {msg}"); continue
        return v


def existing_file(s):
    p = os.path.expanduser(s.strip().strip('"').strip("'"))
    if not os.path.isfile(p):
        raise ValueError(f"no such file: {p}")
    return p


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Merge two overlapping tomographic .xyz maps with a "
                    "Gaussian distance-weighted blend.")
    ap.add_argument("--west", help="western subarray file (lon lat vel)")
    ap.add_argument("--east", help="eastern subarray file (lon lat vel)")
    ap.add_argument("-L", "--length", type=float,
                    help="smoothing length in degrees")
    ap.add_argument("--cw", type=float, help="first array centre longitude")
    ap.add_argument("--ce", type=float, help="second array centre longitude")
    ap.add_argument("--cwlat", type=float,
                    help="first array centre latitude (default: shared mean "
                         "latitude, i.e. a pure east-west blend)")
    ap.add_argument("--celat", type=float, help="second array centre latitude")
    ap.add_argument("--level", choices=["none", "constant", "plane"], default="none",
                    help="remove the systematic offset between the two maps "
                         "before blending, so the blend acts on structure "
                         "rather than turning a level difference into a ramp")
    ap.add_argument("--level-to", choices=["split", "A", "B"], default="split",
                    help="where the levelling correction is applied")
    ap.add_argument("--match-tol", type=float, default=10 ** -KEY_DECIMALS,
                    help="two nodes count as the same if they agree to this "
                         "many degrees (default 1e-4)")
    ap.add_argument("-o", "--out", help="output .xyz")
    ap.add_argument("--diag", action="store_true",
                    help="also write lon lat vel w_W source to *_diag.txt")
    ap.add_argument("--no-plot", action="store_true", help="skip the QC plot")
    a = ap.parse_args()
    set_match_tol(a.match_tol)

    if a.west is None or a.east is None:
        print("\n" + "=" * 64)
        print("  Merge two overlapping tomographic maps   (lon lat vel .xyz)")
        print("=" * 64)

    fW = existing_file(a.west) if a.west else ask("WEST subarray file",
                                                 cast=existing_file)
    fE = existing_file(a.east) if a.east else ask("EAST subarray file",
                                                 cast=existing_file)

    west, east = read_xyz(fW), read_xyz(fE)
    lonW, lonE = west[0], east[0]

    ov_lo = max(lonW.min(), lonE.min())
    ov_hi = min(lonW.max(), lonE.max())
    ov_w = ov_hi - ov_lo

    print(f"\n  west : {len(lonW):7d} nodes   lon {lonW.min():9.3f} .. "
          f"{lonW.max():9.3f}   lat {west[1].min():7.3f} .. {west[1].max():7.3f}")
    print(f"  east : {len(lonE):7d} nodes   lon {lonE.min():9.3f} .. "
          f"{lonE.max():9.3f}   lat {east[1].min():7.3f} .. {east[1].max():7.3f}")
    if ov_w > 0:
        print(f"  longitude overlap: {ov_lo:.3f} .. {ov_hi:.3f}  "
              f"({ov_w:.3f} deg wide)")
    else:
        print("  ! the maps do not overlap in longitude - nothing to blend")

    # Default centres: midpoint of each map's longitude span.
    print()
    cW = a.cw if a.cw is not None else ask(
        "West array centre lon", default=round((lonW.min() + lonW.max()) / 2, 2),
        cast=float)
    cE = a.ce if a.ce is not None else ask(
        "East array centre lon", default=round((lonE.min() + lonE.max()) / 2, 2),
        cast=float)
    if cW > cE:
        print("  ! the 'west' centre lies east of the 'east' centre - "
              "check --cw/--ce and which file is which")

    # Centre latitudes default to one shared value, which makes the blend axis
    # exactly east-west and the whole thing reduce to the longitude-only case.
    lat0 = float(np.mean(np.concatenate([west[1], east[1]])))
    cA = (cW, a.cwlat if a.cwlat is not None else lat0)
    cB = (cE, a.celat if a.celat is not None else lat0)
    _, sep = project(0.0, 0.0, cA, cB)

    # The overlap is measured from the nodes the two maps actually share, not
    # from where their bounding boxes intersect: ragged coverage otherwise
    # inflates it and the recommended L comes out too smooth.
    t_lo, t_hi, n_match = overlap_projection(west, east, cA, cB)
    if ov_w > 0:
        n_coin, worst, pitch = alignment(west, east)
        if n_coin == 0:
            print(f"\n  ! these maps cover the same ground but sit on DIFFERENT GRIDS "
                  f"(node pitch {pitch:.4g} deg).")
            print(f"    Nothing can be blended and the output will step where they "
                  f"meet.  Regrid one map onto the other; widening --match-tol")
            print(f"    would only pretend that different nodes are the same point.")
        elif n_match < 0.98 * n_coin:
            want = min(worst * 1.05, 0.49 * pitch)
            print(f"\n  ! only {n_match} of about {n_coin} coincident nodes matched: the "
                  f"two files print the same grid to different precision")
            print(f"    (nodes up to {worst:.3g} deg apart, tolerance {a.match_tol:g}). "
                  f"The overlap will come out patchy.")
            print(f"    Re-run with --match-tol {want:.3g} - still far inside the "
                  f"{pitch:.4g} deg node pitch.")
    if n_match:
        bearing = np.degrees(np.arctan2(cB[0]-cA[0], cB[1]-cA[1])) % 360
        axis = ("east-west" if abs(cB[1]-cA[1]) < 1e-9 else
                "north-south" if abs(cB[0]-cA[0]) < 1e-9 else
                f"diagonal, bearing {bearing:.0f} deg")
        print(f"  {n_match} shared nodes; blend axis is {axis}, centres "
              f"{sep:.2f} deg apart")

    Lsug = suggest_L_proj(cA, cB, t_lo, t_hi) if n_match else None

    if a.length is not None:
        L = a.length
    else:
        if n_match and sep > 0:
            print(f"\n  Centres are {sep:.2f} deg apart, so the hand-over scale is")
            print(f"  s = L^2/(2*{sep:.2f}) - much sharper than L itself:")
            for cand in [c for c in (sep / 4, sep / 2, Lsug) if c]:
                w = 1.0/(1.0 + np.exp(t_lo/(cand*cand/(2*sep))))
                print(f"      L = {cand:5.2f}  ->  s = "
                      f"{handover_scale(cand, cA, cB):5.2f} deg,  "
                      f"w_A at the near overlap edge = {w:5.3f}")
            if Lsug:
                print(f"  L = {Lsug:.2f} is the smoothest blend that still keeps "
                      f"both overlap edges >= {1-EDGE_TOL:.0%} pure.")
                # Floor, never round: rounding up would push the default just
                # past its own threshold and trip the warning below.
                default_L = float(np.floor(Lsug * 100) / 100)
            else:
                print("  ! the midpoint of the two centres is outside the "
                      "overlap; no L keeps both edges clean")
                default_L = round(float(np.sqrt(sep)), 2)
        else:
            default_L = 1.0
        L = ask("Smoothing length L (deg)", default=default_L, cast=float,
                check=lambda v: "must be > 0" if v <= 0 else None)

    if sep > 0 and n_match:
        s = handover_scale(L, cA, cB)
        wlo = 1.0/(1.0 + np.exp(t_lo/s))
        whi = 1.0/(1.0 + np.exp(t_hi/s))
        print(f"\n  L = {L:g} deg  ->  hand-over scale {s:.2f} deg, "
              f"10-90 % width {4.394*s:.2f} deg")
        print(f"  weights at the overlap edges:  w_A(near) = {wlo:.3f}   "
              f"w_A(far) = {whi:.3f}")
        # 1e-6 slack so that L sitting exactly on the threshold (the suggested
        # value) does not trip its own warning through rounding.
        if (1 - wlo) > EDGE_TOL + 1e-6 or whi > EDGE_TOL + 1e-6:
            print(f"  ! too smooth for this overlap - the merged map will step at "
                  f"the overlap edges;")
            if Lsug and L > Lsug:
                print(f"    lower L to ~{Lsug:.2f} to keep both edges clean")
            elif not Lsug:
                print(f"    the centres are placed badly relative to the "
                      f"overlap - revisit --cw/--ce")

    out = a.out or ask("Output file", default="merged.xyz")

    if a.level != "none":
        west, east, li = level_maps(west, east, a.level, a.level_to)
        if li.get("note"):
            print(f"\n  ! levelling skipped: {li['note']}")
        else:
            print(f"\n  levelling ({li['mode']}, applied to {li['to']}): removed a "
                  f"mean offset of {li['removed_mean']:+.4f}")
            print(f"    overlap rms {li['rms_before']:.4f} -> {li['rms_after']:.4f}; "
                  f"the remaining {li['structure_rms']:.4f} is structure the blend "
                  f"should handle")
            if li["mode"] == "plane":
                print(f"    plane coefficients [c, d/dlon, d/dlat] = "
                      f"{['%+.5f' % c for c in li['coef']]}")

    lon, lat, vel, wWn, src = gaussian_merge(west, east, L, cA, cB)
    ov, diff = overlap_diff(west, east, lon, lat, src)
    nB = int(ov.sum())

    print(f"\n  merged {len(lon)} nodes:  {int((src=='W').sum())} west-only, "
          f"{int((src=='E').sum())} east-only, {nB} blended")
    if nB:
        print(f"  overlap disagreement (east - west):  mean {diff.mean():+.4f}   "
              f"rms {np.sqrt((diff**2).mean()):.4f}   "
              f"range {diff.min():+.4f} .. {diff.max():+.4f}")
        sr = seam_report(lon, wWn, ov, diff)
        if sr:
            print(f"  residual step at the overlap edges (the visible seam):")
            print(f"      lon {sr['lo_edge']:9.3f}  w_W {sr['w_west']:.3f}   "
                  f"rms {sr['rms_w']:.4f}   max {sr['max_w']:.4f}")
            print(f"      lon {sr['hi_edge']:9.3f}  w_E {1-sr['w_east']:.3f}   "
                  f"rms {sr['rms_e']:.4f}   max {sr['max_e']:.4f}")

    with open(out, "w") as fh:
        for x, y, v in zip(lon, lat, vel):
            fh.write(f"{x:12.6f} {y:12.7f} {v:14.8f}\n")
    print(f"\n  merged    -> {out}")

    if a.diag:
        dpath = os.path.splitext(out)[0] + "_diag.txt"
        with open(dpath, "w") as fh:
            fh.write("# lon lat vel w_W source(W/E/B)\n")
            for x, y, v, w, s_ in zip(lon, lat, vel, wWn, src):
                fh.write(f"{x:12.6f} {y:12.7f} {v:14.8f} {w:9.6f} {s_}\n")
        print(f"  diagnostic-> {dpath}")

    if not a.no_plot:
        qc_plot(west, east, (lon, lat, vel), wWn, src, diff, ov, L, cA, cB,
                os.path.splitext(out)[0] + "_qc.png")
    print()


if __name__ == "__main__":
    main()
