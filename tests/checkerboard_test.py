#!/usr/bin/env python3
"""
Checkerboard validation of the merge, in seven steps you can run and read.

    python tests/checkerboard_test.py                 # all steps
    python tests/checkerboard_test.py --step 3        # just step 3
    python tests/checkerboard_test.py --offset 0      # what if the two agree?
    python tests/checkerboard_test.py --pattern blocks --cell 2

Every step prints what it did, what to look at, and a PASS/FAIL check with the
reason. Files land in tests/checkerboard/.

The point of a checkerboard test is that the answer is known. Two overlapping
"inversions" are built from one true model, merged, and the merged map is
compared with the truth it came from. Anything the merge invents shows up as
error; anything it destroys shows up as lost amplitude.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "checkerboard")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, ROOT)
import merge_tomo as mt                                            # noqa: E402
import synthetic as sy                                             # noqa: E402

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

RESULTS, CHECKS = {}, []


# ============================================================ small utilities
def say(*a):
    print(*a)


def head(n, title):
    say("\n" + "=" * 74)
    say(f"  STEP {n}  {title}")
    say("=" * 74)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    say(f"    [{'PASS' if ok else 'FAIL'}]  {name}" + (f"\n            {detail}" if detail else ""))


def rms(a):
    a = np.asarray(a, float)
    return float(np.sqrt((a**2).mean())) if a.size else float("nan")


CM = None
if HAVE_MPL:
    CM = LinearSegmentedColormap.from_list("seis_r", [
        (i/7, (r/255, g/255, b/255)) for i, (r, g, b) in enumerate(
            [(0, 0, 180), (0, 150, 200), (0, 190, 120), (120, 220, 60),
             (255, 255, 0), (255, 180, 0), (255, 90, 0), (170, 0, 0)])])
    DIV = LinearSegmentedColormap.from_list("rdbu_r", [
        (i/6, (r/255, g/255, b/255)) for i, (r, g, b) in enumerate(
            [(33, 102, 172), (146, 197, 222), (209, 229, 240), (247, 247, 247),
             (244, 165, 130), (214, 96, 77), (178, 24, 43)])])


def imshow(ax, lonv, latv, F, title, vmin, vmax, cm, cblab=""):
    im = ax.imshow(np.asarray(F, float).T, origin="lower", aspect="auto",
                   extent=[lonv[0], lonv[-1], latv[0], latv[-1]],
                   vmin=vmin, vmax=vmax, cmap=cm, interpolation="nearest")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("longitude (°)", fontsize=8)
    ax.set_ylabel("latitude (°)", fontsize=8)
    ax.tick_params(labelsize=7)
    cb = ax.figure.colorbar(im, ax=ax, shrink=.88, pad=.02)
    cb.ax.tick_params(labelsize=7)
    if cblab:
        cb.set_label(cblab, fontsize=8)


def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=140, bbox_inches="tight")
    say(f"    figure -> {os.path.relpath(p, ROOT)}")


# ============================================================ the model setup
def build(args):
    """Grids and the true model. Everything downstream hangs off this."""
    g = args.grid
    latv = np.round(np.arange(args.lat0, args.lat1 + 1e-9, g), 6)
    lonW = np.round(np.arange(args.lon0, args.lon0 + args.width + 1e-9, g), 6)
    lonE = np.round(np.arange(args.lon0 + args.width - args.overlap,
                              args.lon0 + 2*args.width - args.overlap + 1e-9, g), 6)
    lonU = np.round(np.arange(lonW[0], lonE[-1] + 1e-9, g), 6)
    cW = (float(np.mean(lonW)), float(np.mean(latv)))
    cE = (float(np.mean(lonE)), float(np.mean(latv)))
    return dict(g=g, latv=latv, lonW=lonW, lonE=lonE, lonU=lonU, cW=cW, cE=cE)


def truth_on(M, lonv, args):
    """The true model sampled on a given longitude vector, same phase everywhere."""
    full = sy.truth_grid(M["lonU"], M["latv"], cell=args.cell, amp=args.amp,
                         trend=not args.no_trend, pattern=args.pattern)
    idx = [int(np.argmin(np.abs(M["lonU"] - x))) for x in lonv]
    return full[idx, :]


# ============================================================ steps
def step1(M, args):
    head(1, "Build a model whose answer you already know")
    T = truth_on(M, M["lonU"], args)
    RESULTS["truth"] = dict(min=float(T.min()), max=float(T.max()),
                            mean=float(T.mean()),
                            anomaly_pp=float(T.max() - T.min()))
    say(f"    a {args.cell:g}° {args.pattern} checkerboard, amplitude ±{args.amp:g} km/s,"
        f" on a {'sloping' if not args.no_trend else 'flat'} background")
    say(f"    grid {M['g']:g}°, {len(M['lonU'])} x {len(M['latv'])} nodes, "
        f"lon {M['lonU'][0]:.1f} .. {M['lonU'][-1]:.1f}, "
        f"lat {M['latv'][0]:.1f} .. {M['latv'][-1]:.1f}")
    say(f"    values {T.min():.3f} .. {T.max():.3f} km/s")
    sy.save_xyz(os.path.join(OUT, "truth.xyz"), M["lonU"], M["latv"], T)
    say(f"    wrote {os.path.relpath(os.path.join(OUT, 'truth.xyz'), ROOT)}")
    if HAVE_MPL:
        fig = Figure(figsize=(6.2, 3.4)); FigureCanvasAgg(fig)
        ax = fig.add_subplot(111)
        imshow(ax, M["lonU"], M["latv"], T, "true model", T.min(), T.max(), CM,
               "km s$^{-1}$")
        for c, col in ((M["cW"], "#2f6fd0"), (M["cE"], "#d9730d")):
            ax.plot(*c, "o", ms=7, color=col, mec="k", mew=.6)
        save(fig, "step1_truth.png")
    check("the true model has the amplitude you asked for",
          abs((T.max() - T.min())/2 - (args.amp + (0.16 if not args.no_trend else 0))) < 0.12,
          f"peak-to-peak {T.max()-T.min():.3f} km/s")
    say("\n    LOOK AT: step1_truth.png. Squares of alternating fast and slow,")
    say("    plus a gentle regional trend. The two dots are the array centres.")
    return T


def step2(M, args):
    head(2, "Turn it into two overlapping subarray inversions")
    TW = truth_on(M, M["lonW"], args)
    TE = truth_on(M, M["lonE"], args)
    W = sy.subarray_recovery(M["lonW"], M["latv"], TW, M["cW"], args.reach,
                             args.blur, 0.0, args.noise, seed=1)
    E = sy.subarray_recovery(M["lonE"], M["latv"], TE, M["cE"], args.reach,
                             args.blur, args.offset, args.noise, seed=2)
    say(f"    each inversion recovers the truth well near its own centre and")
    say(f"    degrades outward: amplitude fades and anomalies broaden over "
        f"{args.blur:g}°,")
    say(f"    with recovery exp(-(d/{args.reach:g})²).")
    if args.offset:
        say(f"    the eastern model also carries a constant {args.offset:+g} km/s offset,")
        say(f"    standing in for a different reference level between two runs.")
    if args.noise:
        say(f"    plus {args.noise:g} km/s of independent noise in each.")

    # how much checkerboard amplitude each model keeps, at its own centre and at
    # the far side of the overlap
    ov0 = M["lonE"][0]
    ov1 = M["lonW"][-1]
    rec = lambda lonv, T, F, x: sy.checkerboard_recovery(lonv, M["latv"], T, F,
                                                         args.cell, x)
    aW_c = rec(M["lonW"], TW, W, M["cW"][0])
    aW_e = rec(M["lonW"], TW, W, ov1)
    aE_c = rec(M["lonE"], TE, E, M["cE"][0])
    aE_e = rec(M["lonE"], TE, E, ov0)
    RESULTS["amplitude_recovery"] = dict(west_centre=aW_c, west_far_edge=aW_e,
                                        east_centre=aE_c, east_far_edge=aE_e)
    say(f"\n    amplitude recovered, western model:  {aW_c:.2f} at its centre,"
        f"  {aW_e:.2f} at lon {ov1:.1f}")
    say(f"    amplitude recovered, eastern model:  {aE_c:.2f} at its centre,"
        f"  {aE_e:.2f} at lon {ov0:.1f}")
    check("each model recovers more amplitude at its own centre than at its far edge",
          aW_c > aW_e and aE_c > aE_e,
          "if this fails the degradation model is not doing anything")
    check("the overlap is where the two disagree most about amplitude",
          abs(aW_e - aE_c) > 0.05 or abs(aE_e - aW_c) > 0.05,
          "the whole reason for weighting by distance to the array centre")

    sy.save_xyz(os.path.join(OUT, "subarray_west.xyz"), M["lonW"], M["latv"], W)
    sy.save_xyz(os.path.join(OUT, "subarray_east.xyz"), M["lonE"], M["latv"], E)
    say(f"    wrote subarray_west.xyz and subarray_east.xyz in "
        f"{os.path.relpath(OUT, ROOT)}")
    if HAVE_MPL:
        lo, hi = min(W.min(), E.min()), max(W.max(), E.max())
        fig = Figure(figsize=(9.4, 3.2)); FigureCanvasAgg(fig)
        gs = fig.add_gridspec(1, 2, wspace=.28)
        imshow(fig.add_subplot(gs[0]), M["lonW"], M["latv"], W,
               "western subarray recovery", lo, hi, CM, "km s$^{-1}$")
        imshow(fig.add_subplot(gs[1]), M["lonE"], M["latv"], E,
               f"eastern subarray recovery"
               + (f" ({args.offset:+g} offset)" if args.offset else ""),
               lo, hi, CM, "km s$^{-1}$")
        save(fig, "step2_inputs.png")
    say("\n    LOOK AT: step2_inputs.png. The checkerboard is crisp near each")
    say("    centre and washed out towards the far side. That is the thing the")
    say("    merge has to exploit.")
    return W, E


def step3(M, args, W, E):
    head(3, "Merge them, with L chosen by the continuity criterion")
    A = sy.to_xyz(M["lonW"], M["latv"], W)
    B = sy.to_xyz(M["lonE"], M["latv"], E)
    tlo, thi, n = mt.overlap_projection(A, B, M["cW"], M["cE"])
    Lstar = mt.suggest_L_proj(M["cW"], M["cE"], tlo, thi)
    D = mt.project(0.0, 0.0, M["cW"], M["cE"])[1]
    s = mt.handover_scale(Lstar, M["cW"], M["cE"])
    L = args.L if args.L else Lstar
    say(f"    {n} nodes are shared by both models")
    say(f"    array centres {D:.2f}° apart")
    say(f"    criterion gives L* = {Lstar:.3f}°  ->  hand-over scale "
        f"s = L²/(2D) = {s:.3f}°")
    say(f"    10-90% hand-over width = {4.394*s:.2f}°, against an overlap "
        f"{thi-tlo:.2f}° wide")
    if args.L:
        say(f"    (you asked for L = {L:g}° instead of L* = {Lstar:.3f}°)")
    wlo = 1/(1 + np.exp(tlo/mt.handover_scale(L, M["cW"], M["cE"])))
    whi = 1/(1 + np.exp(thi/mt.handover_scale(L, M["cW"], M["cE"])))
    say(f"    weight on the western model: {wlo:.3f} at the near overlap edge, "
        f"{whi:.3f} at the far edge")
    lon, lat, vel, wW, src = mt.gaussian_merge(A, B, L, M["cW"], M["cE"])
    ov, d = mt.overlap_diff(A, B, lon, lat, src)
    sr = mt.seam_report(lon, wW, ov, d)
    RESULTS["merge"] = dict(L=float(L), Lstar=float(Lstar), s=float(s), D=float(D),
                            n_shared=int(n), n_out=int(len(lon)),
                            w_near=float(wlo), w_far=float(whi),
                            seam_near=sr["rms_w"], seam_far=sr["rms_e"])
    np.savetxt(os.path.join(OUT, "merged.xyz"), np.c_[lon, lat, vel],
               fmt="%12.6f %12.7f %14.8f")
    say(f"    merged {len(lon)} nodes -> "
        f"{os.path.relpath(os.path.join(OUT, 'merged.xyz'), ROOT)}")
    # Checked for whatever L is actually in use, not just for L*: forcing a
    # larger L has to be reported as inadmissible, or the test is worthless.
    ok_edges = wlo >= 0.9 - 1e-6 and whi <= 0.1 + 1e-6
    if args.L and not ok_edges:
        detail = (f"you forced L = {L:g}°, above L* = {Lstar:.2f}°. The blend is only "
                  f"{100*wlo:.0f}% western at the near edge and {100*(1-whi):.0f}% "
                  f"eastern at the far edge, so the merged map will step where the "
                  f"overlap ends. This FAIL is the expected result of --L "
                  f"{L:g}, and is what the criterion exists to prevent.")
    else:
        detail = f"w = {wlo:.3f} at the near edge and {whi:.3f} at the far edge"
    check(f"the blend at L = {L:.2f}° is at least 90% pure at both overlap edges",
          ok_edges, detail)
    check("output node count equals the union of the two inputs",
          len(lon) == len(M["lonW"])*len(M["latv"]) + len(M["lonE"])*len(M["latv"]) - n,
          f"{len(lon)} nodes")
    say("\n    NOTE: s is much smaller than L. That is the point of the "
        "criterion —")
    say("    a smoothing length of a few degrees is a hand-over of about one.")
    return A, B, L, Lstar


def _err(M, args, lon, lat, vel, region=None):
    T = truth_on(M, M["lonU"], args)
    tab = {}
    for i, x in enumerate(M["lonU"]):
        for j, y in enumerate(M["latv"]):
            tab[(round(float(x), 6), round(float(y), 6))] = T[i, j]
    keep = np.ones(len(lon), bool) if region is None else region(lon)
    e = np.array([v - tab[(round(float(x), 6), round(float(y), 6))]
                  for x, y, v in zip(lon[keep], lat[keep], vel[keep])])
    return e


def step4(M, args, A, B, L):
    head(4, "Score the merge against the model you started from")
    ov0, ov1 = M["lonE"][0], M["lonW"][-1]
    inside = lambda x: (x >= ov0 - 1e-9) & (x <= ov1 + 1e-9)

    def on_overlap(t):
        """RMS error of a (lon, lat, value) triple over the overlap only."""
        lon, lat, vel = (np.asarray(z) for z in t)
        return rms(_err(M, args, lon, lat, vel, inside))

    lon, lat, vel, wW, src = mt.gaussian_merge(A, B, L, M["cW"], M["cE"])
    e = _err(M, args, lon, lat, vel)
    e_ov = _err(M, args, lon, lat, vel, inside)
    say(f"    RMS error over the whole merged map : {rms(e):.4f} km/s")
    say(f"    RMS error inside the overlap        : {rms(e_ov):.4f} km/s")
    say(f"    mean error (bias)                   : {e.mean():+.4f} km/s")

    # The interesting comparison is not merge-vs-truth on its own but
    # merge-vs-each-input, and how that changes once the level offset is gone.
    say(f"\n    on the overlap, RMS error of each option:")
    say(f"      {'':<18s} {'west only':>10s} {'east only':>10s} {'merge':>10s}")
    A2, B2, li = mt.level_maps(A, B, "constant", "split")
    m2 = mt.gaussian_merge(A2, B2, L, M["cW"], M["cE"])
    rowsA = (on_overlap(A), on_overlap(B), rms(e_ov))
    rowsB = (on_overlap(A2), on_overlap(B2), on_overlap(m2[:3]))
    say(f"      {'as given':<18s} {rowsA[0]:10.4f} {rowsA[1]:10.4f} {rowsA[2]:10.4f}")
    say(f"      {'after levelling':<18s} {rowsB[0]:10.4f} {rowsB[1]:10.4f} {rowsB[2]:10.4f}")
    RESULTS["score"] = dict(rms_all=rms(e), rms_overlap=rms(e_ov),
                            bias=float(e.mean()),
                            as_given=dict(west=rowsA[0], east=rowsA[1], merge=rowsA[2]),
                            levelled=dict(west=rowsB[0], east=rowsB[1], merge=rowsB[2]))

    check("the merge beats the worse of the two inputs",
          rowsA[2] < max(rowsA[0], rowsA[1]),
          "the minimum a distance weighting has to achieve")
    if args.offset:
        check("but with a level offset present it does NOT beat the better input",
              rowsA[2] > min(rowsA[0], rowsA[1]),
              f"merge {rowsA[2]:.4f} against {min(rowsA[0], rowsA[1]):.4f}: blending "
              f"in a model that sits {args.offset:+g} km/s off costs you accuracy. "
              f"This is the case for levelling first, not against merging.")
    check("after levelling, the merge beats BOTH inputs",
          rowsB[2] <= min(rowsB[0], rowsB[1]) + 1e-9,
          f"merge {rowsB[2]:.4f} against {rowsB[0]:.4f} and {rowsB[1]:.4f} — two "
          f"independent estimates of the same field combine to something better "
          f"than either, which is the whole point")
    if args.offset:
        say(f"\n    the residual bias of {e.mean():+.4f} km/s is the absolute level, "
            f"which")
        say(f"    neither input constrains. Levelling removes the RELATIVE offset")
        say(f"    between the two; it cannot know which of them was right.")
    if HAVE_MPL:
        Emap = np.full((len(M["lonU"]), len(M["latv"])), np.nan)
        ix = {round(float(x), 6): i for i, x in enumerate(M["lonU"])}
        iy = {round(float(y), 6): j for j, y in enumerate(M["latv"])}
        for x, y, ee in zip(lon, lat, e):
            Emap[ix[round(float(x), 6)], iy[round(float(y), 6)]] = ee
        a = float(np.nanmax(np.abs(Emap))) or 1.0
        fig = Figure(figsize=(6.4, 3.4)); FigureCanvasAgg(fig)
        ax = fig.add_subplot(111)
        imshow(ax, M["lonU"], M["latv"], Emap, "merged − truth", -a, a, DIV,
               "km s$^{-1}$")
        for xv in (ov0, ov1):
            ax.axvline(xv, color="k", lw=.8, ls=":")
        save(fig, "step4_error.png")
    say("\n    LOOK AT: step4_error.png. Error is largest where BOTH models are")
    say("    far from their centres — the middle of the overlap. If you see a")
    say("    sharp line at either dotted edge, L is too large.")


def step5(M, args, A, B, L, Lstar):
    head(5, "Compare against the alternatives, and against levelling")
    ov0, ov1 = M["lonE"][0], M["lonW"][-1]
    inside = lambda x: (x >= ov0 - 1e-9) & (x <= ov1 + 1e-9)
    mid = 0.5*(ov0 + ov1)

    def step_metric(lon, lat, vel):
        """Largest artificial node-to-node jump across the hand-over."""
        T = truth_on(M, M["lonU"], args)
        tab = {}
        for i, x in enumerate(M["lonU"]):
            for j, y in enumerate(M["latv"]):
                tab[(round(float(x), 6), round(float(y), 6))] = T[i, j]
        band = (ov0 - 1.0, ov1 + 1.0)

        def worst(v):
            out = []
            for c in np.unique(np.round(lat, 6)):
                m = np.abs(lat - c) < 1e-9
                o = np.argsort(lon[m])
                xs, vs = lon[m][o], np.asarray(v)[m][o]
                k = (xs[:-1] >= band[0]) & (xs[:-1] <= band[1])
                if k.any():
                    out.append(np.abs(np.diff(vs))[k].max())
            return float(np.percentile(out, 95)) if out else 0.0
        tv = np.array([tab[(round(float(x), 6), round(float(y), 6))]
                       for x, y in zip(lon, lat)])
        return worst(vel) - worst(tv)

    def score(lon, lat, vel):
        return (rms(_err(M, args, lon, lat, vel)),
                rms(_err(M, args, lon, lat, vel, inside)),
                step_metric(lon, lat, vel))

    rows = []
    for label, (a, b) in (("as given", (A, B)),
                          ("after levelling", mt.level_maps(A, B, "constant", "split")[:2])):
        keepA = a[0] <= mid + 1e-9
        keepB = b[0] > mid + 1e-9
        cut = (np.r_[a[0][keepA], b[0][keepB]], np.r_[a[1][keepA], b[1][keepB]],
               np.r_[a[2][keepA], b[2][keepB]])
        opts = {
            "truncate at the midpoint": cut,
            "50/50 average": mt.gaussian_merge(a, b, 1.0, M["cW"], M["cW"])[:3],
            f"blend, L = {Lstar:.2f}° (L*)": mt.gaussian_merge(a, b, Lstar, M["cW"], M["cE"])[:3],
            "blend, L = 3 L*": mt.gaussian_merge(a, b, 3*Lstar, M["cW"], M["cE"])[:3],
        }
        say(f"\n    {label}:")
        say(f"      {'strategy':<28s} {'RMS all':>9s} {'RMS ovl':>9s} {'false step':>11s}")
        for k, (x, y, v) in opts.items():
            r1, r2, st = score(x, y, v)
            rows.append((label, k, r1, r2, st))
            say(f"      {k:<28s} {r1:9.4f} {r2:9.4f} {st:11.4f}")
    RESULTS["strategies"] = [dict(case=c, strategy=k, rms_all=r1, rms_ovl=r2,
                                  false_step=st) for c, k, r1, r2, st in rows]
    get = lambda case, key, i: [r[i] for r in rows if r[0] == case and r[1].startswith(key)][0]
    check("blending leaves a smaller false step than truncating",
          get("as given", "blend, L = %.2f" % Lstar, 4) < get("as given", "truncate", 4),
          f"{get('as given', 'blend, L = %.2f' % Lstar, 4):.4f} vs "
          f"{get('as given', 'truncate', 4):.4f} km/s")
    if args.offset:
        check("levelling reduces total error far more than any choice of L",
              get("after levelling", "blend, L = %.2f" % Lstar, 2)
              < get("as given", "blend, L = %.2f" % Lstar, 2),
              "removing the offset first is the single biggest improvement")
        _, _, li = mt.level_maps(A, B, "constant", "split")
        say(f"\n    levelling estimated the offset as {li['offset_mean']:+.4f} km/s "
            f"(imposed {args.offset:+g})")
        RESULTS["offset_recovered"] = li["offset_mean"]
        check("the estimated offset is close to the imposed one",
              abs(li["offset_mean"] - args.offset) < 0.35*abs(args.offset) + 0.01,
              "it is contaminated by structure whose mean over the overlap is "
              "not zero, so exact recovery is not expected")
    say("\n    NOTE: compare the two blocks. The gain from levelling is much")
    say("    larger than the spread across merging strategies.")


def step6(M, args, A, B, Lstar):
    head(6, "Scan L, and see why accuracy alone cannot choose it")
    ov0, ov1 = M["lonE"][0], M["lonW"][-1]
    inside = lambda x: (x >= ov0 - 1e-9) & (x <= ov1 + 1e-9)
    A2, B2, _ = mt.level_maps(A, B, "constant", "split")
    Ls = np.linspace(max(0.3, 0.1*Lstar), 3.0*Lstar, 22)
    out = {"L": [], "rms": [], "seam": []}
    for L in Ls:
        lon, lat, vel, wW, src = mt.gaussian_merge(A2, B2, L, M["cW"], M["cE"])
        ov, d = mt.overlap_diff(A2, B2, lon, lat, src)
        sr = mt.seam_report(lon, wW, ov, d)
        out["L"].append(float(L))
        out["rms"].append(rms(_err(M, args, lon, lat, vel)))
        out["seam"].append(float(max(sr["rms_w"], sr["rms_e"])) if sr else 0.0)
    RESULTS["L_scan"] = out
    i_best = int(np.argmin(out["rms"]))
    i_star = int(np.argmin(np.abs(np.array(out["L"]) - Lstar)))
    pen = 100*(out["rms"][i_star] - out["rms"][i_best])/out["rms"][i_best]
    say(f"    {'L (deg)':>9s} {'RMS vs truth':>13s} {'seam':>9s}")
    for k in range(0, len(Ls), 3):
        star = "   <- L*" if k == i_star else ""
        say(f"    {out['L'][k]:9.2f} {out['rms'][k]:13.4f} {out['seam'][k]:9.4f}{star}")
    say(f"\n    lowest RMS error is at L = {out['L'][i_best]:.2f}°, "
        f"not at L* = {Lstar:.2f}°")
    say(f"    but the seam at that L is {out['seam'][i_best]:.4f} km/s versus "
        f"{out['seam'][i_star]:.4f} at L*")
    say(f"    choosing L* costs {pen:.1f}% extra RMS error")
    check("RMS error prefers a sharper hand-over than continuity allows",
          out["L"][i_best] < Lstar,
          "RMS is a point-wise measure and barely sees a discontinuity, so it "
          "cannot be used on its own to set L")
    check("the seam grows with L",
          out["seam"][-1] > out["seam"][0],
          "which is the requirement pulling the other way")
    check("obeying the criterion is cheap in accuracy", pen < 8.0,
          f"{pen:.1f}% extra RMS error")
    if HAVE_MPL:
        fig = Figure(figsize=(8.6, 3.2)); FigureCanvasAgg(fig)
        gs = fig.add_gridspec(1, 2, wspace=.3)
        ax = fig.add_subplot(gs[0])
        ax.plot(out["L"], out["rms"], "-", lw=1.8, c="#067647")
        ax.axvline(Lstar, ls="--", c="#1c3f7c", lw=1)
        ax.plot(out["L"][i_best], out["rms"][i_best], "o", c="#067647")
        ax.set_xlabel("L (°)", fontsize=8)
        ax.set_ylabel("RMS error vs truth (km s$^{-1}$)", fontsize=8)
        ax.set_title("accuracy: wants a sharp hand-over", fontsize=9)
        ax.grid(alpha=.25); ax.tick_params(labelsize=7)
        ax = fig.add_subplot(gs[1])
        ax.plot(out["L"], out["seam"], "-", lw=1.8, c="#b42318")
        ax.axvline(Lstar, ls="--", c="#1c3f7c", lw=1)
        ax.set_xlabel("L (°)", fontsize=8)
        ax.set_ylabel("seam at the overlap edges (km s$^{-1}$)", fontsize=8)
        ax.set_title("continuity: wants a wide one", fontsize=9)
        ax.grid(alpha=.25); ax.tick_params(labelsize=7)
        save(fig, "step6_L_scan.png")
    say("\n    LOOK AT: step6_L_scan.png. The dashed line is L*. The two panels")
    say("    pull in opposite directions; L* is the largest L that keeps the")
    say("    seam small, and it costs almost nothing in accuracy.")


def step7(M, args):
    head(7, "Verdict")
    n = len(CHECKS)
    bad = [c for c, ok in CHECKS if not ok]
    for c, ok in CHECKS:
        say(f"    [{'PASS' if ok else 'FAIL'}]  {c}")
    say("")
    if bad:
        say(f"    {len(bad)} of {n} checks FAILED:")
        for b in bad:
            say(f"      - {b}")
    else:
        say(f"    all {n} checks passed — the merge behaves as the method claims")
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=float)
    say(f"\n    numbers -> {os.path.relpath(os.path.join(OUT, 'results.json'), ROOT)}")
    say("\n    THINGS TO TRY NEXT")
    say("      --offset 0            the two models agree on level; levelling "
        "should stop mattering")
    say("      --offset 0.3          a large level difference; watch the false "
        "gradient grow")
    say("      --cell 1.5            anomalies smaller than the smoothing; "
        "recovery collapses")
    say("      --pattern blocks      sharp squares, the harsher classical test")
    say("      --overlap 1           a narrow overlap forces a much smaller L*")
    say("      --overlap 10          a wide overlap allows a gentler hand-over")
    say("      --noise 0.05          can the merge still be scored through noise?")
    say("      --L 12                deliberately too smooth; step 3 SHOULD fail, "
        "and that failure is the point")
    return 1 if bad else 0


# ============================================================ main
def main():
    p = argparse.ArgumentParser(
        description="Checkerboard validation of the tomographic merge.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--step", type=int, choices=range(1, 8),
                   help="run one step only (later steps re-run what they need)")
    p.add_argument("--cell", type=float, default=3.0, help="checkerboard size, deg")
    p.add_argument("--amp", type=float, default=0.12, help="anomaly amplitude, km/s")
    p.add_argument("--pattern", choices=["smooth", "blocks"], default="smooth")
    p.add_argument("--no-trend", action="store_true", help="drop the regional trend")
    p.add_argument("--offset", type=float, default=0.11,
                   help="level difference imposed on the eastern model, km/s")
    p.add_argument("--noise", type=float, default=0.01, help="noise, km/s")
    p.add_argument("--reach", type=float, default=7.0, help="recovery falloff, deg")
    p.add_argument("--blur", type=float, default=1.4, help="broadening, deg")
    p.add_argument("--L", type=float, default=0.0,
                   help="force a smoothing length instead of using L*")
    p.add_argument("--grid", type=float, default=0.25, help="node spacing, deg")
    p.add_argument("--lon0", type=float, default=-126.0)
    p.add_argument("--lat0", type=float, default=36.0)
    p.add_argument("--lat1", type=float, default=51.0)
    p.add_argument("--width", type=float, default=14.0, help="width of each map, deg")
    p.add_argument("--overlap", type=float, default=5.0, help="overlap width, deg")
    args = p.parse_args()

    if not HAVE_MPL:
        say("  ! matplotlib not found — numbers only, no figures")
    say("\n  CHECKERBOARD VALIDATION OF THE TOMOGRAPHIC MERGE")
    say(f"  outputs in {os.path.relpath(OUT, ROOT)}/")

    M = build(args)
    only = args.step
    step1(M, args)
    if only == 1:
        return 0
    W, E = step2(M, args)
    if only == 2:
        return 0
    A, B, L, Lstar = step3(M, args, W, E)
    if only == 3:
        return 0
    if only in (None, 4):
        step4(M, args, A, B, L)
        if only == 4:
            return 0
    if only in (None, 5):
        step5(M, args, A, B, L, Lstar)
        if only == 5:
            return 0
    if only in (None, 6):
        step6(M, args, A, B, Lstar)
        if only == 6:
            return 0
    return step7(M, args)


if __name__ == "__main__":
    sys.exit(main())
