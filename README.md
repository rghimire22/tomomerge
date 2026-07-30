# TomoMerge

*Developed by Riwaj Ghimire, University of Houston*


Merges two overlapping maps into one. Nodes seen by only one map pass through
untouched; nodes in the overlap are replaced by a Gaussian distance-weighted
average that trusts whichever array centre is nearer.

Works for an east–west split, a north–south split, or a diagonal one — the
geometry is handled in general, not just for the case it was written for.

## Run it

Three front ends, one method.

```
merge_tomo.html                 # double-click: opens in any browser, no install
python merge_tomo_gui.py        # desktop window: drop files, set L, press MERGE
python merge_tomo.py            # terminal prompts
python merge_tomo.py --west A.xyz --east B.xyz -L 4.52 -o merged.xyz --diag
```

`merge_tomo.html` is self-contained — no Python, no server, no internet, nothing
installed. Your files are read by the browser and never leave the machine, so it
is safe to email to a collaborator. Drop two files and the result appears with no
further clicks: centres, L, colour scale and the merge all run themselves, and
re-run as you move any control. On a grid large enough that live re-merging would
feel sluggish it says so and switches to a manual **Merge** button.

Input files are three columns: `longitude latitude value`, separated by spaces,
commas or semicolons. Blank lines and lines starting with `#`, `>` or `%` are
skipped, and Fortran `D` exponents (`1.5D-3`) are accepted.

Requires `numpy`. `matplotlib` adds the QC figure and the desktop app's live
preview; `pip install tkinterdnd2` turns the desktop app's two panels into real
drop targets (without it they are click-to-browse). The webpage needs none of
this.

## Choosing the smoothing length

The weight given to map A is an exact logistic in one variable — the distance
along the line joining the two array centres:

    w_A(t) = 1 / (1 + exp(t / s))
    t = distance from the midpoint of the centres, along the centre-to-centre line
    s = L² / (2·D)             hand-over scale,  D = separation of the centres

That single formula covers every orientation, because a difference of two squared
distances is linear in position: `d_B² − d_A² = 2·D·t`. When both centres sit on
the same latitude, `t` is simply `lon − midpoint longitude`, which is the plain
east–west form this started as.

**L is not the hand-over width.** With centres 9° apart, L = 4 gives s = 0.89°,
not 4°. The 10–90 % transition spans 4.39·s.

The blend has to reach ~pure A where the overlap ends on A's side, and ~pure B
where it ends on B's — otherwise the merged map steps at those positions.
**Recommend** gives the largest L that keeps both overlap edges at least 90 %
pure; for the vel02 sa13 pair that is L = 4.52°. Larger L is not "smoother", it
is a visible seam at the overlap boundary. The verdict line turns red when you
cross that limit, and the slider shows the safe range as a green band.

The overlap is measured from the nodes the two maps actually *share*, not from
where their bounding boxes cross — otherwise ragged coverage inflates it and the
recommended L comes out too smooth.

## Node matching

Two inversions can print the same node as `-113.800003` and `-113.800000`, so
nodes are paired within a tolerance (default 1e-4° ≈ 11 m) and written out with
their original coordinates.

Pairing is a true radius search — nearest neighbour within the tolerance, via a
hash grid. The obvious shortcut, snapping both maps onto a common grid and
comparing buckets, silently loses any pair that straddles a bucket boundary. With
float32 coordinates that is a few percent of nodes, which would punch scattered
holes through the overlap that no diagnostic would catch and that look like real
structure on the finished map.

Both the page and the CLI check the pairing against how many nodes are
*coincident* (nearest neighbour within half a node pitch) and tell you when they
disagree:

- **overlapping maps, nothing paired** → the two files are on genuinely different
  grids. Nothing can be blended. Regrid one onto the other; widening the
  tolerance would just pretend different nodes are the same point.
- **overlapping maps, only some paired** → the same grid printed to different
  precision. The page widens the tolerance for you, up to a hard ceiling of half
  a node pitch, and says by how much. The CLI prints the `--match-tol` to use.

## Colour palettes

The page ships the 19 standard GMT palettes — `rainbow`, `relief`, `gray`,
`no_green`, `red2green`, `seis`, `jet`, `ocean`, `sealand`, `polar`, `haxby`,
`topo`, `hot`, `cool`, `wysiwyg`, `globe`, `gebco`, `copper`, `split` — plus
`viridis` and `rdbu`. Click the palette chip to open a swatch grid; **reverse**
mirrors any of them (blue-is-fast is just `seis` or `jet` reversed).
`palette_chart.png` in this folder is a printed index of all of them.

**These are reproductions, not the GMT master files.** They are matched by eye and
are close enough to choose a look, but the RGB values are not guaranteed
identical to GMT's `.cpt` masters — `haxby` uses the published 11-colour table
and is exact, while the topographic ones (`relief`, `ocean`, `sealand`, `topo`,
`globe`, `gebco`) are the loosest. When the colours have to match a figure
already in a paper, press **Load .cpt…** and the real GMT file is parsed and used
verbatim. The loader handles both `z r g b z r g b` and `z r/g/b z r/g/b` rows,
skips `B`/`F`/`N` lines and comments, and converts `COLOR_MODEL = +HSV` palettes.

A `✓` on a swatch marks the perceptually uniform maps, where equal steps in
velocity look like equal steps in colour. `rainbow`, `jet` and `seis` are not
uniform — they invent visual edges where the data is smooth, which matters if
someone is going to read structure off the figure. They are here because GMT
figures use them; the choice is yours.

## Reading the output

- `merged.xyz` — `lon lat value`, map A's file order first, then B-only nodes.
  Coordinates are the originals, so single-map nodes round-trip byte-identically.
- `merged_qc.png` (CLI / desktop) or **Download QC figure** (page) — both inputs,
  the merged map, the weight field, where A and B disagree, the blending kernel,
  and a profile along the blend axis. Check the profile for steps at the overlap
  edges; on the page a slider moves the profile off the centre line.
- `merged_diag.txt` (optional) — adds the per-node weight and a source flag.

The run report prints the overlap disagreement and the residual step at each
overlap edge in value units. For the vel02 sa13 pair the two inversions differ by
**+0.114 km/s mean, 0.186 rms** in the overlap. No weighting scheme removes that
— blending only chooses where you pay it. Lowering L cleans the edges but
concentrates the mismatch into a sharper step mid-overlap. That trade-off is the
one real decision here, which is why the seam is reported as a number rather than
left to the eye.

## Tests

```
python tests/make_fixtures.py       # regenerate references from merge_tomo.py
python tests/python_test.py         # the Python implementation
node   tests/html_test.js           # the page, against the Python references
python tests/resolution_test.py     # the scientific validation  (~2 min)
python tests/checkerboard_test.py   # fast interactive explorer   (~10 s)
```

**Which test to trust for what.** `resolution_test.py` is the validation: it
computes synthetic travel times through a checkerboard along real inter-station
paths, adds noise, and recovers each subarray by damped, smoothed least squares,
so amplitude loss and smearing *emerge* from coverage and regularisation. It
reports a resolution curve, ray coverage, the model resolution matrix and
point-spread functions — the standard diagnostics — and only then scores the
merge.

`checkerboard_test.py` is a fast explorer built on an *assumed* degradation. It
is useful for feeling out how L, the overlap width and a level offset interact,
but it is **not** a resolution test and must not be quoted as one: its imposed
error falls off as a Gaussian in distance from the array centre, the same form as
the merge weights, so it cannot test whether that weighting is appropriate. See
`CHECKERBOARD_TEST.md`.

`html_test.js` extracts both `<script>` blocks from the page and runs them
against a stub DOM, so the real merge path *and* the real drawing path execute;
the fake canvas flags any non-finite coordinate, since in a browser those draw
nothing at all. It asserts the page's `merged.xyz` is **byte-identical** to the
Python output for east–west, north–south and diagonal splits. Run it after
editing either file — that byte-identity is the only thing keeping the two
implementations honest.

`make_fixtures.py` derives every reference number from `merge_tomo.py` itself, so
there are no hand-copied constants to go stale.

## Agreement with the original awk

`tests/python_test.py` compares against `tests/awk_reference.d`, produced by the
`overlapping.d` awk. All 34,881 rows match except 26, each differing by one unit
in the last printed decimal (1e-8). Every one of them sits exactly at the
midpoint between the centres, where the weight is precisely 0.5 and the mean
lands exactly on a rounding tie; the awk's `(wA·vA + wB·vB)/(wA + wB)` and this
tool's normalise-first form straddle that tie by 4e-16. The normalise-first form
is kept because it cannot underflow when both weights are tiny.

## Limits

- L is in degrees in the (lon, lat) plane. For a pure east–west split at high
  latitude a degree of longitude is much shorter than a degree of latitude; the
  weighting is self-consistent either way, but L is not a distance in km.
- Blending cannot remove a systematic offset between two inversions, only move
  where it shows.
- Nodes are paired, never interpolated. Two maps on incompatible grids are
  reported, not resampled.
