# shakti-arctic-design

Designing the liquid (non-frozen) region of Shakti artificial spin ice through the way the array is cut.

Ground states of Shakti spin ice are domino tilings, so a finite array can split into frozen regions and a liquid region bounded by an arctic curve. This repository collects the exploratory calculations on how the edges of the array select that curve. It is a follow-up to [shakti-vorticity](https://github.com/cywon-commits/shakti-vorticity) and uses its lattice conventions.

Status: exploratory (November 2026). The notes in `docs/` are working notes in Korean.

## Results so far

1. **Boundary-height dictionary.** In Shakti coordinates the slope polygon is the square `|dh/dX| <= 1, |dh/dY| <= 1`. The boundary height of any cut follows from counting the unit edges whose perpendicular neighbour lies outside the sample (H plaquette +1, V plaquette -1). Only two straight cuts are extremal: axis cuts through the plaquette blocks (a face of the square, Aztec-type edge) and 45-degree cuts at even offset (a vertex, one ice phase attached to the edge).
2. **Stem-island rule.** An edge plaquette whose "stem" island (the island ending at the midpoint of its outer long side) is absent can close into a free vortex and gains `J/2`. The ground-state energy of the open vertex model equals the minimum number of dimers in a boundary-monomer dimer model. This agreed with integer-programming ground states for all 44 samples tested. The edge termination therefore decides the boundary type, and single stem islands act as switches for boundary-height steps of 4.
3. **Box rule and octagons: a movable arctic circle.** Two boundary points whose height difference equals their l1 distance freeze their bounding box. As a corollary, an array with four axis edges and four even 45-degree edges reduces exactly to an Aztec diamond on an inner square; the position and radius of the arctic circle are set independently by the cuts. Each cut continues one adjacent axis edge, chosen by its offset class (d mod 4); four cuts of mixed classes admit no tiling, and convex corner cuts alone never change the circular shape (checked on 8500 samples).
4. **Plus-shaped arrays: a non-circular liquid region.** For arm reach `n` and arm half-width `c`: no tiling for `c < n/3`; a rational curve of class 10 for `n/3 < c < n/sqrt(2)`, with a three-phase frozen pocket at each reflex corner; the full arctic circle for `c >= n/sqrt(2)`. The curve has an explicit one-parameter form, `c/n = (s^2 + 2s + 2)/(6 + 4s - s^2)`, and agrees with exact statistics up to `n = 481` (347,000 plaquettes): no cell farther than 1-2% of `n` from the curve is misclassified (`figs/plus_verify_n241.png`). The width of the boundary grows as `n^0.36`, close to the expected 1/3.
5. **A temperature switch.** At finite temperature the open sample is a dimer model with boundary monomers of weight `exp(J/2T)`, solved exactly by a Pfaffian. A rectangle with two adjacent block-aligned edges and two adjacent through-block edges has a unique ground state and is liquid at high temperature; the liquid fraction falls continuously from about 3/4 to zero on the scale `T ~ J/L` (`figs/switch_corner.png`). The sample freezes completely at `J/(2 T_f) = m ln(1 + sqrt 2)` (`m` freeable plaquettes). The scaling limit is solved by the hodograph method of the companion vacancy paper: one parameter, an explicit Cauchy integral, the edge splits into a frozen, a partially free and an empty zone ending at `X1 = X2^2`; the solution reproduces the freezing constant and the exact free fractions to better than 1% (`figs/corner_theory.png`). Rectangles with two opposite block-aligned edges admit no tiling.

Caveats: the tangent-line family follows from the Kenyon-Okounkov characteristic solution, but the condition at the reflex corners rests on a physical argument, not a proof. It is supported by numerical agreement and by an independent check: applied to the L-shaped sample (one reflex corner, degree 2) the same rules reproduce the curve of Colomo, Pronko and Sportiello (2018) over its whole length (`figs/lshape_vs_cps.png`). Double-precision Kasteleyn inverses are reliable only up to about `n = 61`; larger sizes use the cancellation-free urban-renewal recursion (`shuffle_probs.py`). See `docs/analytic-results.md`.

Relation to prior work (`docs/literature.md`): the Aztec diamond with one cut-off corner is well studied (Colomo-Pronko 2013, third-order transition when the cut reaches the arctic circle; Colomo-Pronko-Sportiello 2018, arctic curve of the L-shaped domain). The threshold `c = n/sqrt(2)` of the plus shape is the same tangency condition. We did not find the four-corner (plus-shaped) curve or its closed form in the literature; the Shakti-specific parts (dictionary, stem-island rule, corner-cut rule) have no counterpart there.

## Layout

```
code/   Python scripts (flat; run them from inside code/)
docs/   boundary-height-dictionary.md   step 1: cuts -> boundary heights
        arctic-checks.md                step 2: ground states of open samples, octagon, plus shape
        analytic-results.md             box rule, derivation and closed form of the curve, precision limits
        literature.md                   prior work on non-convex domino regions
figs/   verification figures
```

## Scripts

| Script | Purpose |
|---|---|
| `ice_phase_slopes.py` | height gradients of the four ice phases |
| `cut_slopes.py` | boundary slope of straight cuts, by direction and offset class |
| `all_paired_counts.py` | tiling counts of small samples by cut class |
| `gen_ilp.py` | open sample of any polygonal shape: vertex-model ground state (ILP) and effective boundary-monomer model (LP) |
| `test_edge_model_axis_diamond.py`, `test_edge_model_shapes.py` | ILP versus effective model, including randomly kept stem islands |
| `eff.py` | ground-state free plaquettes and remaining entropy for clean cuts at larger sizes |
| `peel.py` | forced-domino peeling; scan of balanced octagons |
| `octagon_counts.py` | octagon tiling counts versus Aztec numbers |
| `regimes.py` | plus shape: tiling counts across `c/n` |
| `box_rule.py`, `mixed_octagon.py` | box rule versus peeling: rectangular samples, both offset classes, and a full scan of mixed classes |
| `theory.py`, `scan_theory.py`, `derive_closed_form.py` | class-10 tangent-line family, numerical and closed-form parameters |
| `precision_check.py` | precision of the Kasteleyn inverse (row sums of pairing probabilities) |
| `shuffle_probs.py`, `plus_shuffle.py` | edge probabilities of weighted Aztec diamonds by urban renewal in the log domain; embedding of a plus-shaped sample |
| `verify_shuffle.py`, `width_scaling.py`, `fig_shuffle.py` | curve versus exact statistics at large sizes, boundary-width scaling, figure |
| `lshape_theory.py`, `lshape_cps.py`, `lshape_compare.py`, `fig_lshape.py` | L-shaped sample: degree-2 tangent family versus the published curve (arXiv:1807.07549); needs `mpmath` |
| `bmono.py`, `switch_square.py` | finite temperature: boundary-monomer dimer model by a Pfaffian (checked against exhaustive counts); free-plaquette fraction of the block-aligned square versus `T_on/T` |
| `mixed_scan.py`, `mixed_digest.py` | all 256 rectangles by edge offsets: state at `T = 0` and at `y = 1` (exact rational arithmetic, needs `python-flint`) |
| `zpoly.py` | exact polynomial `Z(y)`: number of states for every number of free plaquettes |
| `corner_theory.py`, `corner_field.py`, `fig_corner_theory.py` | hodograph solution of the corner-switch sample: edge zones, frozen boundary, temperature |
| `profile_width.py` | free-plaquette probability along the switchable edges versus size at fixed scaled field |
| `switch_corner.py`, `fig_switch_corner.py` | corner-switch sample (two adjacent block-aligned edges): free fraction and liquid region versus temperature |
| `maps.py`, `plusmap.py`, `overlay.py` | exact local dimer statistics and frozen-phase maps |
| `verify_big.py`, `fig_verify.py` | quantitative comparison of curve and exact statistics; figure |

`core.py`, `embed.py`, `find.py`, `vort.py` and `shakti_torus.py` are copied unchanged from shakti-vorticity (chessboard-Shakti map, island set, vertex classification).

## Running

```
pip install -r requirements.txt
cd code
mkdir -p data
python cut_slopes.py
python test_edge_model_axis_diamond.py
python verify_big.py 61 21 31 41 && python fig_verify.py 61
python verify_shuffle.py 61,31 121,61 241,121 && python width_scaling.py
```

`pulp` ships the CBC solver used by the integer programs. `derive_closed_form.py` needs `sympy`.

## License

MIT
