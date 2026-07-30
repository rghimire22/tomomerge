#!/usr/bin/env python3
"""
Build the fixtures tests/html_test.js checks against.

Everything here comes out of merge_tomo.py itself - no constants typed by hand -
so if the Python changes, the reference changes with it and the JavaScript is
held to the new behaviour.  Run from the Merging folder:

    python tests/make_fixtures.py
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures")
os.makedirs(FIX, exist_ok=True)
sys.path.insert(0, ROOT)
import merge_tomo as mt                                          # noqa: E402

TOOL = os.path.join(ROOT, "merge_tomo.py")
F = lambda *p: os.path.join(FIX, *p)


def grid(lon, lat, fn):
    X, Y = np.meshgrid(lon, lat, indexing="ij")
    return np.c_[X.ravel(), Y.ravel(), fn(X.ravel(), Y.ravel())]


def save(name, arr):
    np.savetxt(F(name), arr, fmt="%12.6f %12.7f %14.8f")
    return F(name)


def run(args):
    subprocess.run([sys.executable, TOOL, *args, "--no-plot"],
                   check=True, stdin=subprocess.DEVNULL,
                   stdout=subprocess.DEVNULL)


# ---------------------------------------------------------------- real data
src_w = os.path.join(ROOT, "vel02.fg.sa13(west)")
src_e = os.path.join(ROOT, "vel02.fg.sa13(east)")
if os.path.isfile(src_w) and os.path.isfile(src_e):
    for src, dst in ((src_w, "west.xyz"), (src_e, "east.xyz")):
        np.savetxt(F(dst), np.loadtxt(src), fmt="%12.6f %12.7f %14.8f")
else:
    # stand-in with the same geometry as the vel02 sa13 pair
    lat = np.round(np.arange(36, 51.0001, 0.1), 4)
    f = lambda x, y: 3.9 + 0.15*np.sin(x/2.1) + 0.1*np.cos(y/1.7)
    save("west.xyz", grid(np.round(np.arange(-126, -111.9999, 0.1), 4), lat, f))
    save("east.xyz", grid(np.round(np.arange(-117, -102.9999, 0.1), 4), lat,
                          lambda x, y: f(x, y) + 0.11))

# ------------------------------------------------------- synthetic geometry
lon_q = np.round(np.arange(-120, -109.9999, 0.25), 4)
save("south.xyz", grid(lon_q, np.round(np.arange(34, 44.0001, 0.25), 4),
                       lambda x, y: 3.6 + 0.01*np.sin(x) + 0.02*np.cos(y)))
save("north.xyz", grid(lon_q, np.round(np.arange(40, 50.0001, 0.25), 4),
                       lambda x, y: 4.0 + 0.01*np.sin(x) - 0.02*np.cos(y)))
save("diagA.xyz", grid(np.round(np.arange(-124, -111.9999, 0.25), 4),
                       np.round(np.arange(34, 44.0001, 0.25), 4),
                       lambda x, y: 3.7 + 0.02*np.sin(x + y)))
save("diagB.xyz", grid(np.round(np.arange(-118, -105.9999, 0.25), 4),
                       np.round(np.arange(40, 50.0001, 0.25), 4),
                       lambda x, y: 4.1 + 0.02*np.sin(x - y)))
lat_q = np.round(np.arange(36, 46.0001, 0.25), 4)
save("base.xyz", grid(lon_q, lat_q, lambda x, y: np.full(x.shape, 3.5)))
# half a cell out: a genuinely different grid, which must NOT be silently paired
save("offset.xyz", grid(lon_q + 0.125, lat_q, lambda x, y: np.full(x.shape, 4.0)))
# the same grid printed with coarser rounding: 3e-4 out, far inside the 0.25
# node spacing, so widening the tolerance is the right call
rng = np.random.default_rng(7)
jit = grid(lon_q, lat_q, lambda x, y: np.full(x.shape, 4.0))
jit[:, 0] += rng.uniform(-3e-4, 3e-4, len(jit))
jit[:, 1] += rng.uniform(-3e-4, 3e-4, len(jit))
save("jitter.xyz", jit)

# ------------------------------------------------- reference merges, by CLI
run(["--west", F("west.xyz"), "--east", F("east.xyz"),
     "-L", "4.52", "--cw", "-119", "--ce", "-110", "-o", F("py_default.xyz")])
run(["--west", F("west.xyz"), "--east", F("east.xyz"),
     "-L", "4.0", "--cw", "-118.5", "--ce", "-109.5", "-o", F("py_awk.xyz")])
run(["--west", F("south.xyz"), "--east", F("north.xyz"), "-L", "3.3",
     "--cw", "-115", "--cwlat", "39", "--ce", "-115", "--celat", "45",
     "-o", F("py_ns.xyz")])
run(["--west", F("diagA.xyz"), "--east", F("diagB.xyz"), "-L", "3.0",
     "--cw", "-118", "--cwlat", "39", "--ce", "-112", "--celat", "45",
     "-o", F("py_diag.xyz")])

# ------------------------------------------------------- reference numbers
west, east = mt.read_xyz(F("west.xyz")), mt.read_xyz(F("east.xyz"))
lon, lat, vel, wW, src = mt.gaussian_merge(west, east, 4.52, -119.0, -110.0)
ov, diff = mt.overlap_diff(west, east, lon, lat, src)
seam = mt.seam_report(lon, wW, ov, diff)

ref = {
    "suggestL": mt.suggest_L(-119.0, -110.0, -117.0, -112.0),
    "weights": [[x, L, float(mt.west_weight(x, L, -119.0, -110.0))]
                for L in (0.5, 1.0, 4.52, 9.0, 40.0)
                for x in (-130., -119., -117., -114.5, -112., -110., -100., -90.)],
    "weights2d": [[x, y, L, a[0], a[1], b[0], b[1],
                   float(mt.blend_weight(x, y, L, a, b))]
                  for (a, b) in (((-119., 43.5), (-110., 43.5)),      # east-west
                                 ((-115., 39.), (-115., 45.)),        # north-south
                                 ((-118., 39.), (-112., 45.)),        # diagonal
                                 ((-112., 45.), (-118., 39.)))        # reversed
                  for L in (1.0, 3.0, 9.0)
                  for (x, y) in ((-120., 36.), (-115., 42.), (-114.5, 43.5),
                                 (-110., 48.), (-105., 50.))],
    "scales": [[L, mt.handover_scale(L, -119.0, -110.0)] for L in (0.5, 1, 4.52, 9, 40)],
    "diffMean": float(diff.mean()),
    "diffRms": float(np.sqrt((diff**2).mean())),
    "seamRmsW": seam["rms_w"],
    "seamRmsE": seam["rms_e"],
}
json.dump(ref, open(F("ref.json"), "w"), indent=1)

print(f"fixtures written to {FIX}")
print(f"  {len(ref['weights'])} 1-D weight references, "
      f"{len(ref['weights2d'])} 2-D, 4 reference merges")
