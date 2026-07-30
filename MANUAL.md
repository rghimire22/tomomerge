# TomoMerge — short manual

*Developed by Riwaj Ghimire, University of Houston*

---

## 1. What problem it solves

You invert two overlapping subarrays separately and get two maps. In the overlap
each grid node has **two** velocity values that disagree, because the two
inversions had different ray coverage there. You need one map.

Stitching them at a line leaves a visible step. Averaging them 50/50 everywhere
smears the good part of each. TomoMerge does the third thing: in the overlap it
takes a **weighted average that trusts whichever subarray centre is nearer**, so
each map dominates where its resolution is best and hands over smoothly in
between.

```
        map A                        map B
 ┌──────────────────────┐   ┌──────────────────────┐
 │                 ┌────┼───┼────┐                 │
 │   A alone       │  overlap    │      B alone     │
 │   passed        │  blended    │      passed      │
 │   through       └────┼───┼────┘      through     │
 └──────────────────────┘   └──────────────────────┘
        ▲                                  ▲
      centre A          w_A ──►  w_B      centre B
                    100% A    50/50    100% B
```

Nodes only one map has are copied through **untouched**. Only shared nodes change.

---

## 2. How the weighting works

For a shared node, the weight given to map A is

```
w_A = 1 / (1 + exp(t / s))          w_B = 1 − w_A

t = how far the node is from the midpoint between the two array centres,
    measured along the line joining them
s = L² / (2·D)                      D = distance between the two centres
```

Three things to take from that:

- **It is one-dimensional.** Only position along the centre-to-centre line
  matters, so the same formula handles an east–west split, a north–south split
  or a diagonal one.
- **L is not the hand-over width.** The width is set by `s = L²/(2D)`. With
  centres 9° apart, `L = 4` gives `s = 0.89°`, not 4°. The 10–90 % transition
  spans `4.39·s`.
- **The overlap edges are what constrain L.** Where the overlap ends, the blend
  must already be ~pure A on A's side and ~pure B on B's, or the merged map steps
  at that longitude. **Recommend** gives the largest L that keeps both edges at
  least 90 % pure — for the vel02 sa13 pair, `L = 4.52°`.

Bigger L is not "smoother". Past that limit it *creates* a seam at the overlap
boundary. The verdict line goes red when you cross it, and the slider shows the
safe range as a green band.

---

## 3. Components

### The method has five inputs

| Input | What it is | Default |
|---|---|---|
| Map A, Map B | two `lon lat value` text files | — |
| Array centres | one (lon, lat) per map; where each inversion is most trustworthy | each map's centroid |
| **L** | smoothing length, degrees — sets the hand-over width | largest value keeping both edges ≥ 90 % pure |
| Match tolerance | how close two nodes must be to count as the same node | 1e-4° (≈ 11 m) |
| Edge purity target | how pure the overlap edges must stay | 0.10 (90 %) |

### The folder has three front ends over one engine

| File | Role |
|---|---|
| `merge_tomo.html` | **The app.** Double-click, opens in any browser. No Python, no server, no internet. Files are read locally and never uploaded. |
| `merge_tomo.py` | The engine, and a command-line tool. Everything numerical lives here. |
| `merge_tomo_gui.py` | Desktop window (Python + Tk), imports the engine. |
| `README.md` | Full technical reference. |
| `palette_chart.png` | Index of the 21 colour palettes. |
| `overlapping.d` | The original awk script this grew from, kept as the reference. |
| `tests/` | Three test suites. `html_test.js` proves the webpage's output is byte-identical to Python's. |

The webpage reimplements the maths in JavaScript so it can run with nothing
installed. That is a second copy of the numerics, so it is pinned: the test suite
checks the page produces a **byte-identical** `merged.xyz` to `merge_tomo.py` for
east–west, north–south and diagonal splits. Re-run `node tests/html_test.js`
after editing either file.

### What you get out

| Output | Contents |
|---|---|
| `merged.xyz` | `lon lat value` — map A's file order first, then B-only nodes. Coordinates are the originals. |
| `merged_qc.png` | Six panels: both inputs, merged map, the weight field, where A and B disagree, the kernel, and a profile through the seam. |
| `merged_diag.txt` | Adds the per-node weight and a source flag (A / B / blend). |

---

## 4. Using it

1. **Drop both files** on the page. They sort themselves into A and B by
   geometry; `⇄` swaps them if it guessed wrong.
2. **The result is already there.** Centres, L, colour scale and the merge all
   run themselves, and re-run when you touch any control. On a grid large enough
   that live re-merging would lag, it says so and switches to a manual **Merge**
   button.
3. **Read the verdict line.** Green means the overlap edges are clean. Red means
   L is too large for this overlap — press **Recommend**.
4. **Look at the profile panel** in the QC figure for steps at the shaded overlap
   edges. The offset slider moves the profile off the centre line.
5. **Download** `merged.xyz` and the QC figure.

---

## 5. Before you trust the output

Two numbers in the report matter more than the picture.

**Overlap disagreement** — how much the two inversions differ where they see the
same ground. For the vel02 sa13 pair: **+0.114 km/s mean, 0.186 rms**. Blending
cannot remove a systematic offset like that; it only decides where you pay it.
For comparison, `vel02lht` differs by 0.004 mean / 0.111 rms and `vel03` by 0.025
/ 0.114 — on the mean offset, `vel02` is the outlier of your three.

**Residual step at the overlap edges** — the seam that survives, in km/s. At
L = 4.52 on vel02: rms 0.031 at the western edge, 0.008 at the eastern. Lowering
L cleans the edges but concentrates the mismatch into a sharper step mid-overlap.
That trade-off is the one real decision here, which is why it is reported as a
number instead of left to the eye.

**Node matching warnings.** If the two files sit on genuinely different grids the
tool says so and refuses to invent matches — regrid one map first. If they are the
same grid printed to different precision, it widens the tolerance for you and
tells you by how much. Both messages are worth reading rather than dismissing; a
partial match leaves scattered holes in the overlap that look like real structure.

---

## 6. Limits

- L is in degrees in the (lon, lat) plane, not kilometres.
- Weighting uses distance to the array centres, not the actual resolution
  matrices. If you have formal resolution or ray density, that would be a better
  weight than distance.
- Nodes are paired, never interpolated. Incompatible grids are reported, not
  resampled.
