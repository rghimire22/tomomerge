#!/usr/bin/env python3
"""
Resolution test in the standard tomographic form, and a merge validated on top.

    python tests/resolution_test.py

Unlike tests/checkerboard_test.py, nothing here is an assumed degradation.
Synthetic travel times are computed through a known checkerboard along real
station-pair paths, noise is added, and each subarray's model is recovered by
damped, smoothed least squares.  Amplitude loss, smearing and edge artefacts
emerge from the coverage and the regularisation.

It reports what a tomography paper reports:
  * recovered checkerboards at a range of anomaly sizes, and the amplitude
    recovery curve that follows -- the resolution limit of the experiment
  * the ray-path coverage map
  * the diagonal of the model resolution matrix
  * point-spread functions from the resolution matrix
and then, because that is what this tool is for:
  * whether distance from the array centre -- the quantity the merge weights by
    -- actually predicts the recovered amplitude, against ray density and the
    resolution diagonal as competitors
  * the merge scored against the true model

Outputs land in tests/resolution/.  Takes a couple of minutes; the inversion
operator is factorised once per subarray and reused.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "resolution")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, ROOT)
import merge_tomo as mt                                            # noqa: E402
import tomo_synth as ts                                            # noqa: E402

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

R = {}
CHECKS = []
BASE = 3.85


def say(*a):
    print(*a, flush=True)


def head(t):
    say("\n" + "=" * 76); say("  " + t); say("=" * 76)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    say(f"    [{'PASS' if ok else 'FAIL'}]  {name}" + (f"\n            {detail}" if detail else ""))


if HAVE_MPL:
    CM = LinearSegmentedColormap.from_list("seis_r", [
        (i/7, (r/255, g/255, b/255)) for i, (r, g, b) in enumerate(
            [(0, 0, 180), (0, 150, 200), (0, 190, 120), (120, 220, 60),
             (255, 255, 0), (255, 180, 0), (255, 90, 0), (170, 0, 0)])])
    HOT = LinearSegmentedColormap.from_list("dens", [
        (0, (1, 1, 1)), (.35, (.68, .82, .93)), (.7, (.13, .35, .6)), (1, (.05, .13, .3))])
    DIV = LinearSegmentedColormap.from_list("rdbu_r", [
        (i/6, (r/255, g/255, b/255)) for i, (r, g, b) in enumerate(
            [(33, 102, 172), (146, 197, 222), (209, 229, 240), (247, 247, 247),
             (244, 165, 130), (214, 96, 77), (178, 24, 43)])])


def mapshow(ax, G, m, title, vmin, vmax, cm, cblab="", mask=None, st=None):
    F = G.to_2d(np.asarray(m, float)).astype(float)
    if mask is not None:
        F = F.copy(); F[~G.to_2d(mask).astype(bool)] = np.nan
    im = ax.imshow(F.T, origin="lower", aspect=1/np.cos(np.deg2rad(G.lat_mid)),
                   extent=[G.lonv[0], G.lonv[-1], G.latv[0], G.latv[-1]],
                   vmin=vmin, vmax=vmax, cmap=cm, interpolation="nearest")
    if st is not None:
        ax.plot(st[:, 0], st[:, 1], "k^", ms=1.7, mew=0)
    ax.set_title(title, fontsize=8.5, pad=3)
    ax.set_xlabel("longitude (°)", fontsize=7.5)
    ax.set_ylabel("latitude (°)", fontsize=7.5)
    ax.tick_params(labelsize=7)
    cb = ax.figure.colorbar(im, ax=ax, shrink=.85, pad=.02)
    cb.ax.tick_params(labelsize=6.5)
    if cblab:
        cb.set_label(cblab, fontsize=7)


def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=300, bbox_inches="tight")
    say(f"    figure -> {os.path.relpath(p, ROOT)}")


# ============================================================================
def block_recovery(G, inv, m0, cell, args, fine, NOISE,
                   phases=(0, 0.75, 1.5, 2.25), blk=2.0):
    """
    Amplitude recovery in 2x2 degree blocks, averaged over checkerboard phases.

    Per-cell regression on a single anomaly value is one data point and comes out
    as noise; blocks are what a published figure is actually read at.
    """
    X, Y = np.meshgrid(G.lonv, G.latv, indexing="ij")
    bid = (((X.ravel() - G.lonv[0])//blk).astype(int)*1000
           + ((Y.ravel() - G.latv[0])//blk).astype(int))
    num, den, glob_n, glob_d = {}, {}, 0.0, 0.0
    for sx in phases:
        for sy in phases:
            tv = (BASE + args.amp*np.sign(
                np.sin(np.pi*(X - G.lonv[0] + sx)/cell)
                * np.sin(np.pi*(Y - G.latv[0] + sy)/cell))).ravel()
            rv = 1.0/inv.recover(1.0/tv, m0, noise_vec=NOISE.vector(inv.keys),
                                 fwd_rows=inv.fwd_rows,
                                 fwd_slowness=(1.0/fine(cell, sx, sy)
                                               if inv.fwd_rows is not None else None))
            dt, dr = tv - BASE, rv - rv.mean()
            for b in np.unique(bid):
                m = bid == b
                num[b] = num.get(b, 0.0) + float(dt[m] @ dr[m])
                den[b] = den.get(b, 0.0) + float(dt[m] @ dt[m])
            good = inv.resolution_diag >= args.res_min
            glob_n += float(dt[good] @ dr[good]); glob_d += float(dt[good] @ dt[good])
    per_block = {b: num[b]/den[b] for b in num if den[b] > 0}
    return per_block, bid, (glob_n/glob_d if glob_d > 0 else np.nan)


def build_subarray(G, st, args, tag, Gf=None):
    t0 = time.time()
    pairs = ts.station_pairs(st, args.sep_min, args.sep_max, G.kx)
    rows, _, keys = ts.ray_rows(G, st, pairs)
    inv = ts.Inversion(G, rows, alpha=args.alpha, beta=args.beta)
    inv.resolution_diag = np.diag(inv.resolution())
    inv.keys = keys
    # the same paths traced on the fine grid, so the data can be generated at a
    # resolution the inversion cannot represent
    inv.fwd_rows = ts.ray_rows(Gf, st, pairs)[0] if Gf is not None else None
    say(f"    {tag}: {len(st)} stations, {len(rows)} paths, "
        f"{inv.resolution_diag.max():.2f} peak resolution, "
        f"factorised in {time.time()-t0:.0f} s")
    return inv


# ============================================================================
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--grid", type=float, default=0.5)
    p.add_argument("--fine", type=float, default=4.0,
                   help="forward grid refinement, to avoid an inverse crime")
    p.add_argument("--amp", type=float, default=0.12, help="anomaly amplitude km/s")
    p.add_argument("--noise", type=float, default=1.5, help="travel-time noise, s")
    p.add_argument("--alpha", type=float, default=20.0, help="damping")
    p.add_argument("--beta", type=float, default=150.0, help="smoothing")
    p.add_argument("--stations", type=int, default=170)
    p.add_argument("--sep-min", type=float, default=1.5)
    p.add_argument("--sep-max", type=float, default=10.0)
    p.add_argument("--res-min", type=float, default=0.15,
                   help="publishable threshold on the resolution diagonal")
    p.add_argument("--sizes", type=float, nargs="+",
                   default=[1.0, 1.5, 2.0, 3.0, 4.0, 6.0])
    args = p.parse_args()

    if not HAVE_MPL:
        say("  ! matplotlib missing — numbers only")

    head("SETUP  a real inversion, not an assumed degradation")
    G = ts.Grid(-126, -103, 36, 51, args.grid)
    st = ts.stations(G, n=args.stations)
    W = st[st[:, 0] <= -112.0]
    E = st[st[:, 0] >= -117.0]
    say(f"    grid {G.nx} x {G.ny} = {G.n} cells at {args.grid:g}°")
    say(f"    travel-time noise {args.noise:g} s, damping {args.alpha:g}, "
        f"smoothing {args.beta:g}")
    Gf = ts.Grid(-126, -103, 36, 51, args.grid/args.fine)
    NOISE = ts.MeasurementNoise(args.noise, seed=5)

    def fine(cell, sx=0.0, sy=0.0):
        """The true model on the fine grid, phase-shifted like the coarse one."""
        Xf, Yf = np.meshgrid(Gf.lonv, Gf.latv, indexing="ij")
        return (BASE + args.amp*np.sign(
            np.sin(np.pi*(Xf - Gf.lonv[0] + sx)/cell)
            * np.sin(np.pi*(Yf - Gf.latv[0] + sy)/cell))).ravel()

    say(f"    data are generated on a {args.grid/args.fine:g}° grid "
        f"({Gf.n} cells) and inverted on the {args.grid:g}° grid, so the true "
        f"model is not exactly representable in the model space")
    invW = build_subarray(G, W, args, "western subarray", Gf)
    invE = build_subarray(G, E, args, "eastern subarray", Gf)
    shared = set(invW.keys) & set(invE.keys)
    say(f"    {len(shared)} of {len(invW.keys)} western paths are also eastern "
        f"paths; each carries one noise realisation, not two")
    m0 = np.full(G.n, 1.0/BASE)
    maskW = invW.resolution_diag >= args.res_min
    maskE = invE.resolution_diag >= args.res_min
    say(f"    publishable cells: west {maskW.sum()}, east {maskE.sum()}, "
        f"shared {(maskW & maskE).sum()}")
    R["setup"] = dict(cells=G.n, grid=args.grid, stations=len(st),
                      fwd_cells=Gf.n, fwd_grid=args.grid/args.fine,
                      shared_paths=len(shared),
                      n_west=int(len(W)), n_east=int(len(E)),
                      paths_west=len(invW.rows), paths_east=len(invE.rows),
                      cells_west=int(maskW.sum()), cells_east=int(maskE.sum()),
                      cells_shared=int((maskW & maskE).sum()))
    check("the two subarrays share a substantial region",
          (maskW & maskE).sum() > 100, f"{(maskW & maskE).sum()} cells")

    # ---------------------------------------------------------------- part 1
    head("PART 1  checkerboard recovery against anomaly size")
    say("    the standard resolution test: recover the same experiment at a range")
    say("    of anomaly sizes and see where recovery collapses\n")
    say(f"    {'size (°)':>9s} {'ampl. recovery':>15s} {'verdict':>12s}")
    curve, recovered = [], {}
    for cell in args.sizes:
        _, _, glob = block_recovery(G, invW, m0, cell, args, fine, NOISE,
                                    phases=(0, 1.5))
        curve.append((cell, glob))
        tv = ts.checkerboard(G, cell=cell, amp=args.amp, base=BASE)
        recovered[cell] = (tv, 1.0/invW.recover(
            1.0/tv, m0, noise_vec=NOISE.vector(invW.keys),
            fwd_rows=invW.fwd_rows, fwd_slowness=1.0/fine(cell, 0, 0)))
        verdict = "recovered" if glob > 0.6 else ("marginal" if glob > 0.35 else "lost")
        say(f"    {cell:9.1f} {glob:15.2f} {verdict:>12s}")
    R["resolution_curve"] = [dict(cell=c, recovery=g) for c, g in curve]
    small = [g for c, g in curve if c <= 1.5]
    large = [g for c, g in curve if c >= 4.0]
    check("large anomalies are recovered better than small ones",
          np.mean(large) > np.mean(small) + 0.1,
          f"{np.mean(large):.2f} at >=4° against {np.mean(small):.2f} at <=1.5° — "
          f"this is the resolution limit of the experiment, and it is emergent, "
          f"not imposed")
    lim = next((c for c, g in curve if g > 0.5), None)
    say(f"\n    smallest anomaly recovered at better than 50% amplitude: "
        f"{lim if lim else 'none'}°")
    R["resolution_limit_deg"] = lim

    if HAVE_MPL:
        n = len(args.sizes)
        fig = Figure(figsize=(8.0, 2.2*((n + 1)//2 + 1))); FigureCanvasAgg(fig)
        gs = fig.add_gridspec((n + 1)//2 + 1, 2, hspace=.62, wspace=.42)
        for k, cell in enumerate(args.sizes):
            tv, rv = recovered[cell]
            ax = fig.add_subplot(gs[k//2, k % 2])
            mapshow(ax, G, rv, f"{cell:g}° checkerboard recovered "
                    f"({[g for c, g in curve if c == cell][0]:.2f} amplitude)",
                    BASE - args.amp, BASE + args.amp, CM, "km s$^{-1}$", maskW)
        ax = fig.add_subplot(gs[(n + 1)//2, :])
        ax.plot([c for c, g in curve], [g for c, g in curve], "o-", lw=1.8,
                c="#1c3f7c")
        ax.axhline(0.5, ls="--", lw=1, c="#b42318")
        ax.text(args.sizes[0], 0.52, "50% amplitude", fontsize=7.5, color="#b42318")
        ax.set_xlabel("checkerboard anomaly size (°)", fontsize=8)
        ax.set_ylabel("amplitude recovery", fontsize=8)
        ax.set_title("resolution curve", fontsize=8.5, loc="left")
        ax.grid(alpha=.25); ax.tick_params(labelsize=7); ax.set_ylim(0, 1.05)
        save(fig, "fig_R1_resolution.png")

    # ---------------------------------------------------------------- part 2
    head("PART 2  coverage, resolution matrix and point-spread functions")
    Rmat = invW.resolution()
    diag = invW.resolution_diag
    say(f"    resolution diagonal: median {np.median(diag[maskW]):.2f} over the "
        f"publishable area, peak {diag.max():.2f}")
    psf_pts = [(-119.0, 43.5), (-113.5, 43.5), (-121.0, 48.5)]
    spreads = []
    for pt in psf_pts:
        k = int(G.index(np.array([pt[0]]), np.array([pt[1]]))[0])
        col = np.abs(Rmat[:, k])
        X, Y = np.meshgrid(G.lonv, G.latv, indexing="ij")
        d = np.hypot((X.ravel() - pt[0])*G.kx, Y.ravel() - pt[1])
        # Radius containing half the kernel's mass.  A second moment over the
        # whole grid is dominated by the far-field tails and returns a width of
        # several degrees for a kernel that is plainly one cell across.
        o = np.argsort(d)
        c = np.cumsum(col[o])
        spread = float(d[o][np.searchsorted(c, 0.5*c[-1])]) if c[-1] > 0 else np.nan
        spreads.append(spread)
        say(f"    point-spread at {pt}: half-mass radius {spread:.2f}°, "
            f"peak amplitude {diag[k]:.2f}")
    R["psf"] = [dict(point=list(p), width_deg=s) for p, s in zip(psf_pts, spreads)]
    check("the point-spread function is widest where coverage is thinnest",
          spreads[1] > spreads[0] or spreads[2] > spreads[0],
          f"half-mass radii {['%.2f' % s for s in spreads]}° at the array centre, "
          f"its eastern edge and the north-west corner")
    if HAVE_MPL:
        fig = Figure(figsize=(8.6, 5.0)); FigureCanvasAgg(fig)
        gs = fig.add_gridspec(2, 3, hspace=.62, wspace=.55)
        dens = invW.density/max(invW.density.max(), 1e-9)
        mapshow(fig.add_subplot(gs[0, 0]), G, dens, "(a) ray-path coverage",
                0, 1, HOT, "normalised", None, W)
        mapshow(fig.add_subplot(gs[0, 1]), G, diag,
                "(b) resolution diagonal", 0, float(diag.max()), HOT, "R$_{ii}$")
        tv, rv = recovered[3.0]
        mapshow(fig.add_subplot(gs[0, 2]), G, rv - tv,
                "(c) 3° board: recovered − true",
                -args.amp, args.amp, DIV, "km s$^{-1}$", maskW)
        for k, pt in enumerate(psf_pts):
            kk = int(G.index(np.array([pt[0]]), np.array([pt[1]]))[0])
            col = Rmat[:, kk]
            a = float(np.abs(col).max()) or 1.0
            ax = fig.add_subplot(gs[1, k])
            mapshow(ax, G, col, f"(d{k+1}) kernel at {pt[0]:.0f}°, "
                    f"{pt[1]:.0f}°", -0.22*a, 0.22*a, DIV, "")
            ax.plot(pt[0], pt[1], "k+", ms=8)
            ax.set_xlim(pt[0] - 6.5, pt[0] + 6.5); ax.set_ylim(pt[1] - 5, pt[1] + 5)
        save(fig, "fig_R2_diagnostics.png")

    # ---------------------------------------------------------------- part 3
    head("PART 3  does distance from the array centre predict recovery?")
    say("    this is the assumption the merge weights rest on, so it should be")
    say("    measured rather than asserted\n")
    per_block, bid, _ = block_recovery(G, invW, m0, 3.0, args, fine, NOISE)
    cx, cy = float(np.mean(W[:, 0])), float(np.mean(W[:, 1]))
    X, Y = np.meshgrid(G.lonv, G.latv, indexing="ij")
    rowsB = []
    for b, val in per_block.items():
        m = bid == b
        if diag[m].mean() < args.res_min:
            continue
        lo, la = X.ravel()[m].mean(), Y.ravel()[m].mean()
        rowsB.append((val, np.hypot((lo - cx)*G.kx, la - cy),
                      invW.density[m].mean(), diag[m].mean()))
    Ab = np.array(rowsB)
    cc = lambda col: float(np.corrcoef(Ab[:, col], Ab[:, 0])[0, 1])
    r_dist, r_dens, r_res = cc(1), cc(2), cc(3)
    say(f"    {len(Ab)} publishable 2° blocks, recovery "
        f"{Ab[:, 0].min():.2f} to {Ab[:, 0].max():.2f} (median {np.median(Ab[:, 0]):.2f})")
    say(f"\n    {'predictor':<34s} {'R':>7s} {'R²':>7s}")
    for nm, r in (("distance from array centroid", r_dist), ("ray-path density", r_dens),
                  ("resolution matrix diagonal", r_res)):
        say(f"    {nm:<34s} {r:+7.2f} {r*r:7.2f}")
    R["proxy"] = dict(n_blocks=len(Ab), r_distance=r_dist, r_density=r_dens,
                      r_resolution=r_res)
    check("recovered amplitude does fall off away from the array centre",
          r_dist < -0.4, f"R = {r_dist:+.2f}, so distance is a real predictor")
    check("but the resolution matrix predicts it better than distance does",
          abs(r_res) > abs(r_dist),
          f"R² = {r_res**2:.2f} against {r_dist**2:.2f}. Distance from the centre "
          f"explains only {100*r_dist**2:.0f} per cent of the variance in recovered "
          f"amplitude, so it is a usable but second-best proxy — exactly the "
          f"limitation to state, and the case for weighting by a resolution field "
          f"where one exists")
    if HAVE_MPL:
        fig = Figure(figsize=(7.4, 2.5)); FigureCanvasAgg(fig)
        gs = fig.add_gridspec(1, 3, wspace=.40)
        for k, (col, nm, r) in enumerate(((1, "distance from centroid (°)", r_dist),
                                          (2, "mean ray density (km)", r_dens),
                                          (3, "resolution diagonal", r_res))):
            ax = fig.add_subplot(gs[0, k])
            ax.plot(Ab[:, col], Ab[:, 0], "o", ms=3.2, c="#1c3f7c", alpha=.75)
            ax.set_xlabel(nm, fontsize=8)
            if k == 0:
                ax.set_ylabel("amplitude recovery", fontsize=8)
            ax.set_title(f"R = {r:+.2f},  R² = {r*r:.2f}", fontsize=8.5)
            ax.grid(alpha=.25); ax.tick_params(labelsize=7)
        save(fig, "fig_R3_proxy.png")

    # ---------------------------------------------------------------- part 4
    head("PART 4  merge the two real inversions, scored against the truth")
    cell = 3.0
    tv = ts.checkerboard(G, cell=cell, amp=args.amp, base=BASE)
    rvW = 1.0/invW.recover(1.0/tv, m0, noise_vec=NOISE.vector(invW.keys),
                           fwd_rows=invW.fwd_rows, fwd_slowness=1.0/fine(cell, 0, 0))
    rvE = 1.0/invE.recover(1.0/tv, m0, noise_vec=NOISE.vector(invE.keys),
                           fwd_rows=invE.fwd_rows, fwd_slowness=1.0/fine(cell, 0, 0))
    X, Y = np.meshgrid(G.lonv, G.latv, indexing="ij")
    lon, lat = X.ravel(), Y.ravel()
    A = (lon[maskW], lat[maskW], rvW[maskW])
    B = (lon[maskE], lat[maskE], rvE[maskE])
    cW = (float(np.mean(W[:, 0])), float(np.mean(W[:, 1])))
    cE = (float(np.mean(E[:, 0])), float(np.mean(E[:, 1])))
    tlo, thi, nsh = mt.overlap_projection(A, B, cW, cE)
    Lstar = mt.suggest_L_proj(cW, cE, tlo, thi)
    say(f"    ragged footprints: {maskW.sum()} and {maskE.sum()} nodes, "
        f"{nsh} shared")
    say(f"    array centroids {cW[0]:.2f}° and {cE[0]:.2f}°, "
        f"L* = {Lstar:.2f}°")
    truth = {(round(float(x), 6), round(float(y), 6)): v
             for x, y, v in zip(lon, lat, tv)}
    def err(t):
        x, y, v = (np.asarray(z) for z in t)
        return np.array([vv - truth[(round(float(a), 6), round(float(b), 6))]
                         for a, b, vv in zip(x, y, v)])
    rms = lambda e: float(np.sqrt((e**2).mean()))
    shared = np.zeros(len(lon), bool); shared[maskW & maskE] = True
    ins = lambda x: np.isin(np.round(x, 6), np.round(lon[shared], 6))

    mid = 0.5*(cW[0] + cE[0])
    ka, kb = A[0] <= mid, B[0] > mid
    cut = (np.r_[A[0][ka], B[0][kb]], np.r_[A[1][ka], B[1][kb]],
           np.r_[A[2][ka], B[2][kb]])
    opts = {
        "western model alone": A,
        "eastern model alone": B,
        "truncate at the midpoint": cut,
        "50/50 average": mt.gaussian_merge(A, B, 1.0, cW, cW)[:3],
        f"blend at L* = {Lstar:.2f}°": mt.gaussian_merge(A, B, Lstar, cW, cE)[:3],
    }
    say(f"\n    {'option':<28s} {'RMS vs truth':>13s} {'over the overlap':>18s}")
    scores = {}
    for k, t in opts.items():
        e = err(t)
        x = np.asarray(t[0])
        m = ins(x)
        scores[k] = (rms(e), rms(e[m]) if m.any() else np.nan)
        say(f"    {k:<28s} {scores[k][0]:13.4f} {scores[k][1]:18.4f}")
    R["merge"] = dict(Lstar=float(Lstar), n_shared=int(nsh),
                      scores={k: dict(rms_all=v[0], rms_overlap=v[1])
                              for k, v in scores.items()})
    blend = f"blend at L* = {Lstar:.2f}°"

    # Does distance actually identify the better model? This is the premise of
    # the weighting, and in the overlap - each array's margin - it need not hold.
    lonc = np.unique(np.round(lon[shared], 3))
    agree = tot = 0
    say(f"\n    {'lon':>8s} {'west err':>10s} {'east err':>10s} {'nearer':>8s} "
        f"{'nearer better?':>15s}")
    for L0 in lonc:
        m = shared & (np.abs(lon - L0) < 1e-6)
        if m.sum() < 4:
            continue
        eW = float(np.sqrt(np.mean((rvW[m] - tv[m])**2)))
        eE = float(np.sqrt(np.mean((rvE[m] - tv[m])**2)))
        near = "west" if abs(L0 - cW[0]) < abs(L0 - cE[0]) else "east"
        better = "west" if eW < eE else "east"
        tot += 1; agree += (near == better)
        say(f"    {L0:8.2f} {eW:10.4f} {eE:10.4f} {near:>8s} "
            f"{('yes' if near == better else 'NO'):>15s}")
    frac = agree/max(tot, 1)
    R["merge"]["nearer_is_better_fraction"] = frac
    say(f"\n    the nearer model is the better one in {agree} of {tot} columns "
        f"({100*frac:.0f} per cent)")
    check("the blend beats truncating the models at a line",
          scores[blend][0] < scores["truncate at the midpoint"][0],
          f"{scores[blend][0]:.4f} against "
          f"{scores['truncate at the midpoint'][0]:.4f} km/s over the whole map, and "
          f"{scores[blend][1]:.4f} against {scores['truncate at the midpoint'][1]:.4f} "
          f"over the overlap")
    check("the blend costs nothing measurable against uniform averaging",
          abs(scores[blend][1] - scores["50/50 average"][1])
          < 0.03*scores["50/50 average"][1],
          f"{scores[blend][1]:.4f} against {scores['50/50 average'][1]:.4f} km/s, a "
          f"difference of "
          f"{100*abs(scores[blend][1]-scores['50/50 average'][1])/scores['50/50 average'][1]:.1f}"
          f" per cent. The blend does NOT beat averaging on accuracy, and should not "
          f"be sold as if it did; what it adds is continuity, an explicit criterion "
          f"for the hand-over, and preservation of each model where only it exists.")
    check("distance identifies the better model more often than chance, "
          "but far from always",
          0.5 < frac < 0.95,
          f"{100*frac:.0f} per cent of columns. In the overlap both models sit at "
          f"their own margins and are of similar quality, so a radial proxy has "
          f"little to work with - the same limitation as the R^2 of Part 3, seen "
          f"from the accuracy side. This is why the blend cannot beat averaging here.")

    if HAVE_MPL:
        fig = Figure(figsize=(8.8, 5.2)); FigureCanvasAgg(fig)
        gs = fig.add_gridspec(2, 3, hspace=.50, wspace=.58)
        lo_, hi_ = BASE - args.amp, BASE + args.amp
        mapshow(fig.add_subplot(gs[0, 0]), G, tv, "(a) true model", lo_, hi_, CM,
                "km s$^{-1}$")
        mapshow(fig.add_subplot(gs[0, 1]), G, rvW,
                "(b) western inversion", lo_, hi_, CM, "km s$^{-1}$", maskW, W)
        mapshow(fig.add_subplot(gs[0, 2]), G, rvE,
                "(c) eastern inversion", lo_, hi_, CM, "km s$^{-1}$", maskE, E)
        mg = opts[blend]
        full = np.full(G.n, np.nan); em = np.full(G.n, np.nan)
        ix = {(round(float(x), 6), round(float(y), 6)): i
              for i, (x, y) in enumerate(zip(lon, lat))}
        for x, y, v in zip(*mg):
            full[ix[(round(float(x), 6), round(float(y), 6))]] = v
        e = err(mg)
        for (x, y), ee in zip(zip(mg[0], mg[1]), e):
            em[ix[(round(float(x), 6), round(float(y), 6))]] = ee
        mapshow(fig.add_subplot(gs[1, 0]), G, full,
                f"(d) merged at L* = {Lstar:.2f}°", lo_, hi_, CM, "km s$^{-1}$")
        ec = err(cut); ecm = np.full(G.n, np.nan)
        for (x, y), ee in zip(zip(cut[0], cut[1]), ec):
            ecm[ix[(round(float(x), 6), round(float(y), 6))]] = ee
        a = float(np.nanmax(np.abs(np.r_[em[np.isfinite(em)], ecm[np.isfinite(ecm)]])))
        ax = fig.add_subplot(gs[1, 1])
        mapshow(ax, G, ecm, "(e) truncated − truth", -a, a, DIV, "km s$^{-1}$")
        ax.axvline(mid, c="k", lw=.8, ls="--")
        mapshow(fig.add_subplot(gs[1, 2]), G, em, "(f) merged − truth", -a, a, DIV,
                "km s$^{-1}$")
        save(fig, "fig_R4_merge.png")

    # ---------------------------------------------------------------- part 5
    head("PART 5  choosing L, on inversions that were not built to suit it")
    ov, dif = mt.overlap_diff(A, B, *mt.gaussian_merge(A, B, Lstar, cW, cE)[:2],
                              mt.gaussian_merge(A, B, Lstar, cW, cE)[4])
    say(f"    disagreement between the two inversions over the shared region:")
    say(f"      mean {dif.mean():+.4f}, rms {np.sqrt((dif**2).mean()):.4f} km/s")
    say(f"    this offset is emergent - the two inversions damp towards the same")
    say(f"    starting model but with different coverage, so they disagree on level")
    say(f"    as well as on structure, exactly as the real map pairs do.")
    R["emergent_offset"] = dict(mean=float(dif.mean()),
                                rms=float(np.sqrt((dif**2).mean())))
    Ls = np.linspace(max(0.4, 0.12*Lstar), 3.0*Lstar, 20)
    scan = {"L": [], "rms": [], "seam": []}
    for L in Ls:
        lo2, la2, v2, w2, s2 = mt.gaussian_merge(A, B, L, cW, cE)
        o2, d2 = mt.overlap_diff(A, B, lo2, la2, s2)
        sr = mt.seam_report(lo2, w2, o2, d2)
        scan["L"].append(float(L))
        scan["rms"].append(rms(err((lo2, la2, v2))))
        scan["seam"].append(float(max(sr["rms_w"], sr["rms_e"])) if sr else 0.0)
    R["L_scan"] = scan
    ib = int(np.argmin(scan["rms"]))
    istar = int(np.argmin(np.abs(np.array(scan["L"]) - Lstar)))
    pen = 100*(scan["rms"][istar] - scan["rms"][ib])/scan["rms"][ib]
    say(f"\n    lowest RMS error at L = {scan['L'][ib]:.2f}°, criterion at "
        f"{Lstar:.2f}°, cost of obeying it {pen:.1f}%")
    say(f"    seam at the RMS optimum {scan['seam'][ib]:.4f}, at L* "
        f"{scan['seam'][istar]:.4f} km/s")
    rmsv = 100*(max(scan["rms"]) - min(scan["rms"]))/min(scan["rms"])
    seam_lo = max(min(scan["seam"]), 1e-6)
    seamv = max(scan["seam"])/seam_lo
    R["L_penalty_pct"] = pen
    R["L_min_error"] = scan["L"][ib]
    R["L_span_factor"] = float(max(scan["L"])/min(scan["L"]))
    R["rms_variation_pct"] = float(rmsv)
    R["seam_at_Lstar"] = scan["seam"][istar]
    R["seam_at_3Lstar"] = scan["seam"][-1]
    say(f"    across a factor of {max(scan['L'])/min(scan['L']):.0f} in L, RMS error "
        f"varies by only {rmsv:.1f}% while the seam grows by a factor of {seamv:.0f}")
    check("RMS error is too flat in L to select it", rmsv < 5.0,
          f"only {rmsv:.1f}% variation across a factor of "
          f"{max(scan['L'])/min(scan['L']):.0f} in L. With independent noise in the "
          f"two inversions, wider averaging suppresses variance, so accuracy no "
          f"longer opposes smoothing - it simply stops discriminating.")
    check("the seam is the quantity that does discriminate",
          scan["seam"][-1] > 3*scan["seam"][istar],
          f"{scan['seam'][istar]:.4f} km/s at L* against "
          f"{scan['seam'][-1]:.4f} at 3L*")
    check("the criterion lands at the accuracy optimum, not merely near it",
          pen < 1.0, f"{pen:.2f}% above the best RMS achieved at any L, while "
          f"holding the seam to {scan['seam'][istar]:.4f} km/s")
    if HAVE_MPL:
        fig = Figure(figsize=(7.6, 2.9)); FigureCanvasAgg(fig)
        gs = fig.add_gridspec(1, 2, wspace=.32)
        ax = fig.add_subplot(gs[0, 0])
        ax.plot(scan["L"], scan["rms"], "-", lw=1.9, c="#067647")
        ax.plot(scan["L"][ib], scan["rms"][ib], "o", ms=5, c="#067647")
        ax.axvline(Lstar, ls="--", lw=1, c="#1c3f7c")
        ax.set_xlabel("$L$ (°)", fontsize=8)
        ax.set_ylabel("RMS error vs true model (km s$^{-1}$)", fontsize=8)
        ax.set_title("(a)  accuracy prefers a sharp hand-over", fontsize=8.5, loc="left")
        ax.grid(alpha=.25); ax.tick_params(labelsize=7)
        ax = fig.add_subplot(gs[0, 1])
        ax.plot(scan["L"], scan["seam"], "-", lw=1.9, c="#b42318")
        ax.axvline(Lstar, ls="--", lw=1, c="#1c3f7c")
        ax.text(Lstar, max(scan["seam"])*.96, " $L^*$", fontsize=8, color="#1c3f7c")
        ax.set_xlabel("$L$ (°)", fontsize=8)
        ax.set_ylabel("seam at the overlap edges (km s$^{-1}$)", fontsize=8)
        ax.set_title("(b)  continuity prefers a wide one", fontsize=8.5, loc="left")
        ax.grid(alpha=.25); ax.tick_params(labelsize=7)
        save(fig, "fig_R5_L.png")

    # ---------------------------------------------------------------- verdict
    head("VERDICT")
    bad = [c for c, ok in CHECKS if not ok]
    for c, ok in CHECKS:
        say(f"    [{'PASS' if ok else 'FAIL'}]  {c}")
    say("")
    say(f"    {len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed"
        + ("" if not bad else f"; FAILED: {bad}"))
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=float)
    say(f"    numbers -> {os.path.relpath(os.path.join(OUT, 'results.json'), ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
