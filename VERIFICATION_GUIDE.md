# Verifying everything yourself

*A step-by-step protocol for checking the method, the code and every number in
the paper on your own machine, before you put your name to any of it.*

Nothing here asks you to take a result on trust. Each step tells you what to
run, what you should see, and — more importantly — **what it would mean if you
saw something else**.

Work through it in order. Steps 1–3 take about a minute each. Step 4 takes a
couple of minutes. Budget an hour to do it properly and read the output.

---

## Step 0 — What you need

```powershell
python --version      # 3.9 or newer
python -c "import numpy, matplotlib; print(numpy.__version__, matplotlib.__version__)"
node --version        # only for the browser-page test
```

Missing anything:

```powershell
pip install numpy matplotlib
```

Node.js is only needed to test the browser version. Get it from
[nodejs.org](https://nodejs.org) if you want that check; everything else runs
without it.

Then, from the `Merging` folder:

```powershell
python tests/make_fixtures.py
```

This writes the reference files the tests compare against. It derives every one
of them **from `merge_tomo.py` itself**, so there are no hand-typed constants
anywhere in the test suite. Expect:

```
fixtures written to .../tests/fixtures
  40 1-D weight references, 60 2-D, 4 reference merges
```

---

## Step 1 — Does the merge engine do what the equations say?

```powershell
python tests/python_test.py
```

**Expect: `ALL CHECKS PASSED`, 41 checks.**

What it is actually proving, and why each matters:

| Check | Why it matters |
|---|---|
| Weights finite and in [0,1] for L from 10⁻⁶ to 10⁶ | The logistic is written to avoid overflow. A NaN here would silently blank part of your map. |
| Matches the raw Gaussian ratio to 10⁻¹⁵ | Confirms equation (4) really is equation (1), not an approximation to it. |
| Exact for east–west, north–south, diagonal, and reversed splits | This is the claim in §2.2 that the reduction to one dimension holds for **any** geometry. |
| Suggested L puts the edges on 90/10 to 10⁻¹² | Equation (6) is solved correctly, not approximately. |
| Radius search finds every jittered pair | §2.5. Coordinate snapping loses pairs; this proves the replacement does not. |
| Levelling recovers an imposed constant exactly | §2.6. A planted 0.137 km s⁻¹ offset is recovered to machine precision. |
| Agrees with the original awk on 34,855 of 34,881 lines | The method reproduces the script it grew from. |

**On those 26 differing lines** — this one deserves your attention rather than a
shrug. They all sit at exactly the midpoint longitude, where the weight is
precisely 0.5, so the merged value is the exact mean of two numbers and lands
exactly halfway between two 8-decimal outputs. The awk computes
`(wA·vA + wB·vB)/(wA + wB)`; this code normalises the weights first. The two are
algebraically identical and differ by 4×10⁻¹⁶ in floating point, which is enough
to fall on opposite sides of a rounding tie. The normalise-first form is kept
because it cannot underflow when both weights are tiny. If you ever see
differences that are *not* at the midpoint, or larger than 10⁻⁸, something is
genuinely wrong.

---

## Step 2 — Is the browser version the same program?

```powershell
node tests/html_test.js
node tests/html_test.js docs/index.html
```

**Expect: `ALL CHECKS PASSED`, 79 checks, both times.**

The browser tool reimplements the maths in JavaScript so it can run with nothing
installed. That is a second copy of the numerics and could drift from the first.
This test extracts the code from the page, runs it, and requires the merged
output to be **byte-identical** to the Python result — same 34,881 lines, same
digits — for east–west, north–south and diagonal geometries.

The second command checks the *published* page at
`https://rghimire22.github.io/tomomerge/`. Run it after any edit to the page. If
it fails, the live tool is not the tested tool.

---

## Step 3 — The resolution test (the one the paper stands on)

```powershell
python tests/resolution_test.py
```

Takes about two minutes; it factorises a 1457×1457 inverse operator twice.
**Expect: `11 of 11 checks passed`.**

### What it does

This is a conventional tomographic resolution test, not a simulation of one:

1. A checkerboard velocity model is specified on a fine 0.125° grid.
2. Synthetic travel times are computed **through that model, along the real
   inter-station paths** of 168 stations (≈4,000 paths per subarray).
3. 1.5 s of Gaussian noise is added — a realistic dispersion-measurement error.
4. Each subarray is inverted for slowness on a **coarser 0.5° grid** by damped,
   smoothed least squares.
5. The recovered models are masked to cells with resolution diagonal > 0.15 —
   what a paper would actually publish.

Nothing about the degradation is imposed. Amplitude loss, smearing, edge
artefacts and the null space all emerge from the coverage and the
regularisation. **This matters:** if the degradation were imposed as a Gaussian
falloff from the array centre, it would have the same functional form as the
merge weights, and the test would be circular. That is exactly why
`tests/checkerboard_test.py` (Step 5) must not be quoted as a resolution test.

### Part 1 — model recovery against anomaly size

```
 size (°)  ampl. recovery      verdict
      1.0            0.24         lost
      1.5            0.48     marginal
      2.0            0.59     recovered
      3.0            0.74     recovered
      4.0            0.82     recovered
      6.0            0.87     recovered
```

**Read it as:** the experiment resolves anomalies down to about 2°; below that,
recovered amplitude is under half the truth and structure should not be
interpreted. This is the standard statement of resolution in a tomography paper,
and `figures/fig_R1_resolution.png` is the standard figure.

*Your numbers should match these to ±0.01.* If recovery is near 1.0 everywhere,
regularisation is too weak or noise too low — the inversion has become trivially
well determined and the test proves nothing.

### Part 2 — coverage, resolution matrix, point-spread functions

```
resolution diagonal: median 0.51 over the publishable area, peak 0.71
point-spread at (-119.0, 43.5): half-mass radius 2.09°
point-spread at (-113.5, 43.5): half-mass radius 2.47°
point-spread at (-121.0, 48.5): half-mass radius 1.24°
```

A resolution diagonal well below 1 means each recovered cell is an **average
over its neighbourhood**, not an independent estimate. The kernel broadens from
2.09° at the array centre to 2.47° at the poorly covered eastern margin — the
degradation the merge is supposed to exploit, appearing on its own.

Look at `fig_R2_diagnostics.png` panels d1–d3: the faint lobes radiating from
each kernel follow the dominant path azimuths. That is anisotropic coverage, and
it is precisely what a radial distance measure cannot see — which is Part 3.

### Part 3 — is distance a fair proxy for resolution? *(the paper's key limitation)*

```
predictor                                R      R²
distance from array centroid         -0.75    0.57
ray-path density                     +0.81    0.65
resolution matrix diagonal           +0.90    0.80
```

**This is the number to understand before you defend the paper.** Distance from
the array centre — the quantity the merge weights by — explains 57 % of the
spatial variation in recovered amplitude. The resolution matrix explains 80 %.
So the proxy is real but coarse, and roughly half the variation is invisible to
it because coverage is ragged and anisotropic while a radial measure is
circular by construction.

A reviewer *will* ask about this. The paper states it in §4.3 and returns to it
in §6, including what would break if you substituted a resolution field: the
reduction to one dimension in equation (4) depends on Gaussian weights about two
points, so the closed-form criterion of equation (6) does not survive and would
have to be imposed numerically.

### Part 4 — how the merge actually performs

```
option                        RMS vs truth   over the overlap
western model alone                 0.0648             0.0716
eastern model alone                 0.0645             0.0645
truncate at the midpoint            0.0642             0.0691
50/50 average                       0.0637             0.0669
blend at L* = 4.61°                 0.0638             0.0672
```

**Do not oversell this.** The blend beats truncation. It does **not** beat a
uniform average, and it does not beat the better single model. The test then
shows why, column by column: the nearer model is the more accurate one in only
**67 %** of longitude columns, because in the overlap both models are at their
own margins and differ in quality by only ~10 %.

The honest claim, and the one the paper makes, is that the method's advantage
over averaging is *control* — continuity at a stated width, exact preservation
of singly covered nodes, a quantified seam — not accuracy.

### Part 5 — why L cannot be chosen by accuracy

```
across a factor of 25 in L, RMS error varies by only 1.0% while the seam
grows by a factor of ~60,000
lowest RMS at L = 8.93°, criterion at L* = 4.61°, cost of obeying it 0.13%
```

Two competing effects nearly cancel: widening the hand-over averages more
independent noise (good) while giving weight to the poorer model (bad). The
error surface is left almost flat, so misfit **cannot** select L. The seam can,
and does. `fig_R5_L.png` shows both curves.

---

## Step 4 — Reproduce every figure and number in the paper

```powershell
python paper/make_paper_figures.py
node paper/build_manuscript.js
```

The first regenerates all eight figures at 300 dpi and writes
`paper/results.json`. The second rebuilds the manuscript, **reading every
quoted number from that JSON**. No figure in the paper is hand-copied, and no
number in the text is typed by hand. Change the data or the method, re-run both,
and the manuscript updates itself.

To satisfy yourself of that, open `paper/results.json` and pick any number in
the text — the 0.0672, the 57 per cent, the 67 per cent — and find it there.

---

## Step 5 — The fast explorer (for intuition, not for the paper)

```powershell
python tests/checkerboard_test.py
python tests/checkerboard_test.py --offset 0
python tests/checkerboard_test.py --cell 1.5
python tests/checkerboard_test.py --overlap 1
python tests/checkerboard_test.py --L 12
```

Ten seconds each, seven annotated steps, 14 checks. Use it to build a feel for
how L, the overlap width and a level offset interact.

**But do not cite it as a resolution test.** Its degradation is *imposed* as
`exp(−(d/reach)²)`, the same functional form as the merge weights, so it cannot
test whether distance weighting is appropriate. `CHECKERBOARD_TEST.md` says so
at the top. Step 3 is the test that counts.

Note that `--L 12` is *designed to fail*: it reports that the blend is only 58 %
pure at the overlap edges, so the map will step there. That failure is the
criterion doing its job.

---

## Step 6 — Run it on your own data

```powershell
python merge_tomo.py --west "vel02.fg.sa13(west)" --east "vel02.fg.sa13(east)" -o merged.xyz --diag
```

Read the report it prints. Four things deserve your attention:

**1. Shared nodes and blend axis.** `7701 shared nodes; blend axis is east-west,
centres 9.00 deg apart`. If the shared count is far below what you expect, the
two grids are not aligned — see the match-tolerance warnings.

**2. The recommended L.** `L = 4.53 is the smoothest blend that still keeps both
overlap edges >= 90% pure`. Larger is not smoother; it puts a step at the
overlap boundary.

**3. The disagreement decomposition.**

```
vel02      offset +0.1139   (37% of the disagreement variance)
vel02lht   offset +0.0044   ( 0%)
vel03      offset +0.0252   ( 5%)
```

`vel02` is the outlier. Its +0.114 km s⁻¹ offset is 15 % of that map's full
dynamic range, and blending converts a constant offset into a **smooth gradient
across the hand-over** that is indistinguishable from real structure. Before
publishing anything derived from a merged `vel02`, either level it —

```powershell
python merge_tomo.py --west "vel02.fg.sa13(west)" --east "vel02.fg.sa13(east)" --level constant -o merged.xyz
```

— or state the offset explicitly. `vel02lht` and `vel03` do not have this
problem, which is itself worth knowing: whatever differs between those
processing runs, it is not the reference level.

**4. The residual seam.** `lon -117.000 rms 0.0309 / lon -112.000 rms 0.0082`.
The western edge is four times worse, because lon −117 is the *eastern* array's
own western boundary, where its coverage is poorest. Your merge carries most of
its risk at one specific place, and the tool identifies it without being told
anything about the deployments.

Then open `merged_qc.png` and check the profile panel for steps at the shaded
overlap edges.

---

## What would make me doubt a result

| Symptom | Likely cause |
|---|---|
| Amplitude recovery ≈ 1.0 at all anomaly sizes | Regularisation too weak or noise too low; the resolution test is not testing anything |
| Resolution diagonal ≈ 1.0 | Same — the inverse problem is over-determined |
| `0 blended` despite overlapping maps | Grids are offset; read the match-tolerance message, do not just widen the tolerance |
| Seam larger than the disagreement rms | L is far too large; press Recommend |
| Python and browser outputs differ | Stop. One of them has been edited without re-running Step 2 |
| Offset > 10 % of the map's dynamic range | Level the models before blending, or report the offset |

---

## Before submitting the paper

Things only you can check, which I could not:

- **Section 5 provenance.** I inferred that the maps are Rayleigh-wave phase
  velocity from western-US subarrays, from the filenames and geometry. Correct
  the period, the method (ambient noise? earthquake? eikonal?), the station
  network, and the inversion parameters.
- **Three citations taken from memory,** not verified against the publisher:
  Shapiro et al. (2005), Deal & Nolet (1996), Smith & Wessel (1990). The other
  nineteen were checked.
- **Abstract length.** Currently ~385 words. GJI allows this; JGR and SRL cap at
  250. Trim the two sentences on the L-selection result first — it is the
  least essential of the three findings.
- **Acknowledgements** is a placeholder.
- **Author list and funding.**
