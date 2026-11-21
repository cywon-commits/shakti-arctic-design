# shakti-arctic-design

Designing the liquid (non-frozen) region of Shakti artificial spin ice through the way the array is cut.

Ground states of Shakti spin ice are domino tilings, so a finite array can split into frozen regions and a liquid region bounded by an arctic curve. This repository collects the exploratory calculations on how the edges of the array select that curve. It is a follow-up to [shakti-vorticity](https://github.com/cywon-commits/shakti-vorticity) and uses its lattice conventions.

Status: exploratory (November 2026). The notes in `docs/` are working notes in Korean.

## Results so far

1. **Boundary-height dictionary.** In Shakti coordinates the slope polygon is the square `|dh/dX| <= 1, |dh/dY| <= 1`. The boundary height of any cut follows from counting the unit edges whose perpendicular neighbour lies outside the sample (H plaquette +1, V plaquette -1). Only two straight cuts are extremal: axis cuts through the plaquette blocks (a face of the square, Aztec-type edge) and 45-degree cuts at even offset (a vertex, one ice phase attached to the edge).
2. **Stem-island rule.** An edge plaquette whose "stem" island (the island ending at the midpoint of its outer long side) is absent can close into a free vortex and gains `J/2`. The ground-state energy of the open vertex model equals the minimum number of dimers in a boundary-monomer dimer model. This agreed with integer-programming ground states for all 44 samples tested. The edge termination therefore decides the boundary type, and single stem islands act as switches for boundary-height steps of 4.
3. **Octagons: a movable arctic circle.** An array with four axis edges and four even 45-degree edges reduces exactly to an Aztec diamond on an inner square. The position and radius of the arctic circle are set independently by the four 45-degree cuts.
4. **Plus-shaped arrays: a non-circular liquid region.** For arm reach `n` and arm half-width `c`: no tiling for `c < n/3`; a rational curve of class 10 for `n/3 < c < n/sqrt(2)`, with a three-phase frozen pocket at each reflex corner; the full arctic circle for `c >= n/sqrt(2)`. The curve agrees with exact Kasteleyn statistics up to `n = 121` (`figs/plus_verify_n121.png`).

Caveats: the form of the class-10 curve is an ansatz based on the Kenyon-Okounkov structure and is checked numerically, not proven. Prior literature on non-convex domino regions has not been searched yet. See `docs/arctic-checks.md`, sections 4.4 and 6.

## Layout

```
code/   Python scripts (flat; run them from inside code/)
docs/   boundary-height-dictionary.md   step 1: cuts -> boundary heights
        arctic-checks.md                step 2: ground states of open samples, octagon, plus shape
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
| `theory.py`, `scan_theory.py` | class-10 tangent-line family and its parameters |
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
python verify_big.py 61 21 31 41
python verify_big.py 121 45 61 && python fig_verify.py
```

`pulp` ships the CBC solver used by the integer programs. `verify_big.py 121 ...` needs about 3 GB of memory.

## License

MIT
