#!/usr/bin/env python3
"""
Tests for merge_tomo.py.

    python tests/make_fixtures.py
    python tests/python_test.py
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures")
sys.path.insert(0, ROOT)
import merge_tomo as mt                                          # noqa: E402

F = lambda *p: os.path.join(FIX, *p)
fails = []


def check(name, cond, extra=""):
    print(("  PASS  " if cond else "  FAIL  ") + name + (f"   {extra}" if extra else ""))
    if not cond:
        fails.append(name)


west, east = mt.read_xyz(F("west.xyz")), mt.read_xyz(F("east.xyz"))

# --------------------------------------------------------------------------
print("\n--- the blend kernel ---")
lon = np.linspace(-140, -90, 4001)
for L in (1e-6, 1e-3, 0.5, 4.52, 20, 1e6):
    w = mt.west_weight(lon, L, -119.0, -110.0)
    check(f"L={L:<8g} finite, in [0,1], monotone",
          np.all(np.isfinite(w)) and w.min() >= 0 and w.max() <= 1
          and np.all(np.diff(w) <= 1e-15))

check("matches the raw Gaussian ratio wherever that is computable", (lambda: (
    lambda ra, rb: np.abs(mt.west_weight(lon, 4.52, -119., -110.)[(ra+rb) > 1e-300]
                          - (ra/(ra+rb))[(ra+rb) > 1e-300]).max() < 1e-15)(
        np.exp(-((np.abs(lon+119.)/4.52)**2)), np.exp(-((np.abs(lon+110.)/4.52)**2))))())

X, Y = np.meshgrid(np.linspace(-125, -104, 90), np.linspace(36, 51, 70))
for label, cA, cB in (("east-west", (-119., 43.), (-110., 43.)),
                      ("north-south", (-114.5, 39.), (-114.5, 48.)),
                      ("diagonal", (-119., 40.), (-110., 47.)),
                      ("diagonal reversed", (-110., 47.), (-119., 40.))):
    dA = np.hypot(X - cA[0], Y - cA[1]); dB = np.hypot(X - cB[0], Y - cB[1])
    ea, eb = np.exp(-(dA/4.52)**2), np.exp(-(dB/4.52)**2)
    ok = (ea + eb) > 1e-300
    err = np.abs(mt.blend_weight(X, Y, 4.52, cA, cB)[ok] - (ea/(ea+eb))[ok]).max()
    check(f"{label} split is the exact 2-D Gaussian ratio", err < 1e-14,
          f"max err {err:.1e}")

check("scalar centres == [lon,lat] centres on a shared latitude",
      np.array_equal(mt.west_weight(lon, 4.52, -119., -110.),
                     mt.blend_weight(lon, 43.5, 4.52, (-119., 43.5), (-110., 43.5))))
Ls = mt.suggest_L_proj((-119., 43.), (-110., 43.), -2.5, 2.5)
check("suggested L puts the overlap edges exactly on 90/10",
      abs(float(mt.blend_weight(-117., 43., Ls, (-119., 43.), (-110., 43.))) - 0.9) < 1e-12
      and abs(float(mt.blend_weight(-112., 43., Ls, (-119., 43.), (-110., 43.))) - 0.1) < 1e-12,
      f"L = {Ls:.6f}")
check("no L is offered when the midpoint is outside the overlap",
      mt.suggest_L_proj(-119., -110., 1.0, 5.0) is None)

# --------------------------------------------------------------------------
print("\n--- node matching ---")
rng = np.random.default_rng(3)
lo = np.round(np.arange(-120, -110, 0.05), 4)
la = np.full_like(lo, 41.0)
jit_lo = lo + rng.uniform(-1.6e-5, 1.6e-5, len(lo))
jit_la = la + rng.uniform(-1.6e-5, 1.6e-5, len(la))
pair = mt.match_nodes(lo, la, jit_lo, jit_la, 1e-4)
check("a true radius finds every jittered pair "
      "(coordinate snapping loses the ones straddling a bucket)",
      np.array_equal(pair, np.arange(len(lo))),
      f"{int((pair == np.arange(len(lo))).sum())}/{len(lo)}")
check("nothing is matched beyond the tolerance",
      np.all(mt.match_nodes(lo, la, lo + 0.02, la, 0.01) == -1))
check("ties resolve to the lower index",
      mt.match_nodes(np.array([0.0]), np.array([0.0]),
                     np.array([0.001, -0.001]), np.array([0.0, 0.0]), 0.01)[0] == 0)
n_coin, worst, pitch = mt.alignment(west, east)
# the real files are float32, so the pitch is 0.1 only to about 1e-6
check("alignment finds the pitch and every coincident node of the real pair",
      abs(pitch - 0.1) < 1e-5 and n_coin == 7701, f"pitch {pitch:.8f}, {n_coin} coincident")
base, off = mt.read_xyz(F("base.xyz")), mt.read_xyz(F("offset.xyz"))
check("half-cell-offset grids report zero coincident nodes",
      mt.alignment(base, off)[0] == 0)
jitf = mt.read_xyz(F("jitter.xyz"))
nc, wg, pt = mt.alignment(base, jitf)
check("the same grid printed coarsely is seen as coincident but out of tolerance",
      nc == len(base[0]) and 1e-4 < wg < 1e-3,
      f"{nc} coincident, worst gap {wg:.2e}, pitch {pt:.3f}")

# --------------------------------------------------------------------------
print("\n--- the merge ---")
lonm, latm, velm, wW, src = mt.gaussian_merge(west, east, 4.52, -119.0, -110.0)
ref = np.loadtxt(F("py_default.xyz"))
check("output reproduces the stored reference exactly",
      np.array_equal(np.c_[lonm, latm], ref[:, :2]) and
      np.abs(velm - ref[:, 2]).max() < 5e-9, f"{len(lonm)} nodes")
check("nodes only one map has are passed through untouched", (lambda: (
    lambda kw, ke: all(
        (kw.get((round(x, 4), round(y, 4)), ke.get((round(x, 4), round(y, 4)))) == v)
        for x, y, v, s in zip(lonm, latm, velm, src) if s != "B")
    )({(round(x, 4), round(y, 4)): v for x, y, v in zip(*west)},
      {(round(x, 4), round(y, 4)): v for x, y, v in zip(*east)}))())
check("weights are 1 / 0 / strictly between for A-only / B-only / blended",
      np.all(wW[src == "W"] == 1) and np.all(wW[src == "E"] == 0) and
      np.all((wW[src == "B"] > 0) & (wW[src == "B"] < 1)))
ov, diff = mt.overlap_diff(west, east, lonm, latm, src)
check("every blended node has a disagreement value", len(diff) == int(ov.sum()) == 7701)

# repeated rows inside one file
d = (np.array([-115., -115., -114.9]), np.array([40., 40., 40.]), np.array([3., 9., 4.]))
lo2 = mt.gaussian_merge(d, (np.array([-200.]), np.array([0.]), np.array([1.])),
                        1.0, -115., -114.)[0]
check("a repeated node is kept once, first occurrence winning", len(lo2) == 3,
      f"{len(lo2)} rows out")

# --------------------------------------------------------------------------
print("\n--- levelling ---")
lo_g = np.round(np.arange(-120, -110 + 1e-9, 0.5), 4)
la_g = np.round(np.arange(38, 46 + 1e-9, 0.5), 4)
Xg, Yg = np.meshgrid(lo_g, la_g, indexing="ij")
xs, ys = Xg.ravel(), Yg.ravel()
struct = 0.05*np.sin(xs) * np.cos(ys)
Aw = (xs, ys, 3.8 + struct)
C = 0.137
Be = (xs, ys, 3.8 + struct + C)                       # identical but for a constant
w2, e2, info = mt.level_maps(Aw, Be, "constant", "split")
check("a pure constant offset is recovered", abs(info["offset_mean"] - C) < 1e-12,
      f"{info['offset_mean']:.6f} vs {C}")
check("and removed to machine precision", info["rms_after"] < 1e-12,
      f"rms after {info['rms_after']:.1e}")
check("'split' leaves the pair mean unchanged",
      abs((w2[2].mean() + e2[2].mean()) - (Aw[2].mean() + Be[2].mean())) < 1e-12)
check("node counts are untouched", len(w2[0]) == len(Aw[0]) and len(e2[0]) == len(Be[0]))
wA, eA, _ = mt.level_maps(Aw, Be, "constant", "B")
check("'B' leaves map A alone and moves B onto it",
      np.allclose(wA[2], Aw[2]) and np.allclose(eA[2], Be[2] - C))
wB, eB, _ = mt.level_maps(Aw, Be, "constant", "A")
check("'A' does the reverse",
      np.allclose(eB[2], Be[2]) and np.allclose(wB[2], Aw[2] + C))

Bp = (xs, ys, 3.8 + struct + 0.2 + 0.03*(xs + 115) - 0.02*(ys - 42))
_, _, ip = mt.level_maps(Aw, Bp, "plane", "split")
_, _, ic = mt.level_maps(Aw, Bp, "constant", "split")
check("a linear trend is removed by 'plane' but not by 'constant'",
      ip["rms_after"] < 1e-12 and ic["rms_after"] > 0.01,
      f"plane {ip['rms_after']:.1e}, constant {ic['rms_after']:.4f}")
check("plane coefficients recover the imposed trend",
      abs(ip["coef"][1] - 0.03) < 1e-9 and abs(ip["coef"][2] + 0.02) < 1e-9,
      f"{['%+.5f' % c for c in ip['coef']]}")
check("mode 'none' is a no-op", mt.level_maps(Aw, Be, "none")[2]["mode"] == "none")
try:
    mt.level_maps(Aw, Be, "quadratic"); ok = False
except ValueError:
    ok = True
check("an unknown mode is rejected", ok)
far = (xs - 60, ys, Be[2])
check("no shared nodes -> levelling declines instead of inventing an offset",
      "note" in mt.level_maps(Aw, far, "constant")[2])
# levelling must not change what the blend does to a pair that already agrees
lo_a, la_a, v_a, _, _ = mt.gaussian_merge(Aw, Aw, 3.0, -119., -110.)
w0, e0, _ = mt.level_maps(Aw, Aw, "constant", "split")
lo_b, la_b, v_b, _, _ = mt.gaussian_merge(w0, e0, 3.0, -119., -110.)
check("levelling an already-consistent pair is a no-op", np.allclose(v_a, v_b))

# --------------------------------------------------------------------------
print("\n--- agreement with the original awk ---")
awk = os.path.join(HERE, "awk_reference.d")
if os.path.isfile(awk):
    lo3, la3, v3, _, _ = mt.gaussian_merge(west, east, 4.0, -118.5, -109.5)
    # compare the WRITTEN text, not a double against 8-decimal text, or the awk
    # file's own rounding swamps the thing being measured
    mine = [f"{x:12.6f} {y:12.7f} {v:14.8f}" for x, y, v in zip(lo3, la3, v3)]
    theirs = [ln.rstrip("\n") for ln in open(awk) if ln.strip()]
    check("same number of rows", len(mine) == len(theirs), f"{len(mine)} vs {len(theirs)}")
    bad = [i for i, (p, q) in enumerate(zip(mine, theirs)) if p != q]
    check("agrees with awk on all but a handful of last-digit ties",
          len(bad) < 0.001 * len(mine), f"{len(bad)} of {len(mine)} lines differ")
    if bad:
        dv = np.abs(np.array([float(mine[i].split()[2]) for i in bad])
                    - np.array([float(theirs[i].split()[2]) for i in bad]))
        # one unit in the last printed place; the subtraction itself is not exact
        check("and each is one unit in the last printed place", dv.max() <= 1.5e-8,
              f"max {dv.max():.3e}")
        check("all of them sit exactly at the centre midpoint, where the weight "
              "is 0.5 and the two algebraic forms straddle a rounding tie",
              bool(np.all(np.abs(lo3[bad] - (-114.0)) < 1e-9)))
else:
    print("  skip   tests/awk_reference.d not present")

print("\n" + ("ALL CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
sys.exit(1 if fails else 0)
