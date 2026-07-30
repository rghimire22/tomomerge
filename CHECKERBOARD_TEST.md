# Checkerboard explorer — step-by-step

*How to feel out the merge yourself, and what each result means.*

> **This is not a resolution test, and should not be quoted as one.** The
> degradation here is *imposed*: anomaly amplitude is faded as
> `exp(−(d/reach)²)` away from each array centre. That is the same functional
> form as the merge weights, so this experiment cannot test whether distance
> weighting is appropriate — the conclusion would be contained in the premise.
> For the scientific validation, where synthetic travel times are inverted and
> the degradation emerges from coverage and regularisation, run
> `python tests/resolution_test.py`. Use this script to build intuition about how
> L, the overlap width and a level offset interact, which it does far faster.

A checkerboard test works because **you already know the answer**. You build a
true velocity model, degrade it into two overlapping "inversions" the way real
inversions degrade, merge them, and compare the merged map against the model you
started from. Anything the merge invents shows up as error; anything it destroys
shows up as lost amplitude.

## Run it

```
cd Merging
python tests/checkerboard_test.py
```

Needs `numpy`; `matplotlib` adds the figures. Everything is written to
`tests/checkerboard/`. Takes about ten seconds.

Run one step at a time while you read:

```
python tests/checkerboard_test.py --step 1
python tests/checkerboard_test.py --step 2
...
```

Each step prints what it did, tells you which figure to look at, and runs a
`[PASS]`/`[FAIL]` check with the reason spelled out.

---

## What each step does

### Step 1 — build a model whose answer you know

A 3° checkerboard, ±0.12 km/s, on a gentle regional trend. Writes
`truth.xyz` and `step1_truth.png`.

**Look for:** alternating fast and slow squares. The two dots are the array
centres — 119°W and 110°W, 9° apart, the same geometry as your vel02 pair.

### Step 2 — degrade it into two subarray inversions

Each "inversion" recovers the truth well near its own centre and degrades
outward: anomaly amplitude fades and anomalies broaden, mixing in a 1.4°-blurred
copy of the field with weight `1 − exp(−(d/7)²)`. The eastern model also gets a
constant **+0.11 km/s** offset, matching the offset measured in your real vel02
pair. Writes `subarray_west.xyz`, `subarray_east.xyz`, `step2_inputs.png`.

**Look for:** the checkerboard is crisp near each centre and washed out towards
the far side. That contrast is the thing a distance-weighted merge exploits.

**Expected numbers:** amplitude recovery about **0.78** at each centre and
**0.62–0.65** at the far edge of the overlap.

### Step 3 — merge, with L from the criterion

Reports the shared node count, the centre separation *D*, the criterion value
*L\**, the hand-over scale *s = L²/2D*, and the weight at each overlap edge.
Writes `merged.xyz`.

**Expected:** `L* = 4.526°`, `s = 1.138°`, 10–90 % width `5.00°` against a
`5.00°` overlap, weights `0.900` and `0.100` at the edges. Note how much
smaller *s* is than *L* — that is the trap the criterion exists to close.

### Step 4 — score it against the truth

The table here is the most informative thing in the whole test:

```
                    west only  east only      merge
  as given             0.0302     0.1099     0.0675
  after levelling      0.0633     0.0627     0.0616
```

**Read it like this.** With the +0.11 offset still present, the merge is *worse*
than simply keeping the western model — because you are blending in a model that
sits 0.11 km/s off. That is not an argument against merging; it is the argument
for **levelling first**. Once the offset is removed, the merge beats *both*
inputs, which is what combining two independent estimates of the same field
should do.

The residual bias that survives levelling is the **absolute** level, which
neither input constrains. Levelling removes the *relative* offset between the
two; it cannot know which of them was right.

**Look for:** in `step4_error.png`, error is largest in the middle of the
overlap, where both models are far from their centres. A sharp line at either
dotted edge would mean L is too large.

### Step 5 — compare against the alternatives

Truncation, 50/50 averaging, the blend at *L\**, and a deliberately over-smooth
blend at 3*L\**, each scored before and after levelling. The `false step` column
is the largest artificial node-to-node jump introduced across the hand-over,
measured against the same statistic for the true model so real structure is not
counted as an artefact.

**Expected:** truncation leaves a false step of about **0.095 km/s**; the blend
at *L\** cuts that to about **0.017**. And the gain from levelling (0.080 →
0.062 overall) is much larger than any difference between merging strategies.

### Step 6 — scan L

RMS error against truth, and the seam amplitude, as a function of L. Writes
`step6_L_scan.png`.

**Expected, and this is the important one:** RMS error is nearly flat and its
minimum lies *below* *L\**. RMS is a point-wise measure and a discontinuity
occupies almost no points, so **accuracy alone cannot choose L** — it will always
pull you towards the sharp hand-over you are trying to avoid. The seam pulls the
other way. *L\** sits between them and costs under 1 % in RMS.

### Step 7 — verdict

All 14 checks, a `results.json` with every number, and a list of variations to
try.

---

## Things to try, and what should happen

| Command | What it tests | What you should see |
|---|---|---|
| `--offset 0` | the two models agree on level | 11 checks (the offset ones drop out); the merge now beats both inputs immediately |
| `--offset 0.3` | a large level difference | the false gradient grows; levelling matters more |
| `--cell 1.5` | anomalies smaller than the 1.4° smoothing | amplitude recovery collapses — no merge can retrieve what neither inversion resolved |
| `--pattern blocks` | sharp squares, the harsher classical test | lower recovered amplitude; sharp corners are unrecoverable |
| `--overlap 1` | a narrow overlap | *L\** drops sharply — a narrow overlap forces a fast hand-over |
| `--overlap 10` | a wide overlap | *L\** rises; the blend can be gentle |
| `--noise 0.05` | heavy noise | the seam statistic stays usable; RMS rises |
| `--L 12` | forcing an inadmissible L | **step 3 FAILS on purpose**, reporting that the blend is only 58 % pure at the edges and the map will step there. That failure is the criterion doing its job |
| `--reach 3` | very rapid coverage falloff | each model is good only near its own centre; weighting matters more |
| `--no-trend` | checkerboard only | isolates anomaly recovery from the regional field |

## What this test can and cannot tell you

**It can** show that the weighting prefers the better-resolved model, that the
criterion for L is the right way to set it, that levelling matters more than L,
and that the implementation reproduces all of this on data where truth is known.

**It cannot** tell you that the *real* subarray inversions degrade the way this
generator assumes. The degradation model here — amplitude falling off as
`exp(−(d/reach)²)` with distance-dependent broadening — is a reasonable
caricature of what damping does to poorly covered regions, not a simulation of
your inversion. If you have resolution matrices or ray-density maps from the real
inversions, comparing them against `exp(−(d/reach)²)` would be the next honest
step, and would tell you whether distance from the array centre is a fair proxy
in your case.

## Files written

```
tests/checkerboard/
  truth.xyz              the model you started from
  subarray_west.xyz      degraded western "inversion"
  subarray_east.xyz      degraded eastern "inversion"  (carries the offset)
  merged.xyz             the merge, at L*
  results.json           every number the run printed
  step1_truth.png
  step2_inputs.png
  step4_error.png
  step6_L_scan.png
```

`truth.xyz`, `subarray_west.xyz` and `subarray_east.xyz` are ordinary
three-column files, so you can also drop the two subarray files straight into
`merge_tomo.html` and merge them by hand to confirm the browser build agrees.
