# Paper figures (draft v1, 2026-11-26)

Story: the cut sets the boundary heights (Fig. 1); convex cuts move the arctic circle (Fig. 2), reflex corners reshape it (Fig. 3);
at finite temperature the edge energies reduce to corner sources, which switch the liquid region on and off (Fig. 4) and make a
through-block square a physical realisation of the boundary-vacancy model (Fig. 5); Fig. 6 checks the boundary model against the full vertex model.
Scripts: `code/paper_fig1.py` ... `paper_fig6.py`, common style `code/figstyle.py` (final width 7.0 in, 8 pt text, Okabe-Ito colours,
liquidity 1 - P_max with the same colour range 0-0.75 in every map).

**Fig. 1.** Cuts as boundary conditions. (a) Shakti islands and plaquettes (orange H, blue V); in a ground state each plaquette shares one unhappy
vertex (dot) with a neighbour, a domino. Red: the stem island of one plaquette. (b) Slope square of the height function: its four corners are the
four ice phases; an axis cut through the blocks fixes one component of the slope (an edge, Aztec type), an even 45-degree cut fixes the whole slope
(a corner, one ice phase). (c) Ground state of a rectangle with block-aligned W and S edges and through-block, fringed N and E edges (L = 12):
plaquettes without their stem island are free (yellow, gain J/2); half-dimer sites (diamonds) cost J/2 when occupied and matter only at T > 0.

**Fig. 2.** Convex cuts move and resize the arctic circle without changing its shape. Exact liquidity of octagonal samples (|X|, |Y| <= 21,
through-block axis edges, even 45-degree cuts; Kasteleyn) with the circle predicted by the box rule (red) inscribed in the remaining square (dashed).
(d) mixes the two offset classes d = 0, 2 (mod 4).

**Fig. 3.** Reflex corners make the liquid region non-circular. (a-c) Plus-shaped samples with arm reach n = 241 and half-width c, exact liquidity
(urban-renewal recursion, about 1.7 x 10^5 plaquettes) with the class-10 curve (red). (d) The family of curves in one quadrant for
1/3 < c/n < 1/sqrt(2) (dotted: sample outline); below c/n = 1/3 there is no tiling, above 1/sqrt(2) the curve is the full circle.

**Fig. 4.** A temperature switch. (a) Ground states of rectangles by edge type: only two adjacent block-aligned edges give a unique ground state
that melts. (b-d) Exact liquidity of that sample (L = 48, boundary model with half-dimer sites) at three temperatures, with the frozen boundaries of
the scaling solution (red). (e) Liquid fraction versus J/(2T)/m (m freeable plaquettes) for L = 32, 48, 64 and the scaling solution; the sample
freezes completely at J/(2T) = m ln(1 + sqrt 2).

**Fig. 5.** Edge energies reduce to corner sources. (a) Sign pattern of the corner terms: alternating on the through-block square, a diagonal pair on
the corner switch. (b, c) Through-block square (Aztec order N = 23) below and above T_on. (d) sigma = 1 - <n_mono>/2N of the Shakti square for
N = 15, 31, 47 against the scaling function of the boundary-monomer Aztec diamond, with J/T = (xN + 2) ln 2; x = 1 is T_on.

**Fig. 6.** The boundary model against the full vertex model (exact density of states, corner switch). (a) Excitation energy for L = 8-20:
vertex model (circles), boundary model with half-dimer sites (lines), freeable sites only (dotted). (b) Remaining difference versus 1/L.

## Checklist (figure guide P1, P3), draft v1

| Item | Status |
|---|---|
| One claim per figure, readable from the figure | yes (titles above) |
| Data as points, predictions as lines | maps are data, red curves are predictions; Fig. 4e, 5d, 6 points vs lines |
| Exact results: no error bars needed | to state in captions (all values exact) |
| Same colour = same object | liquidity map viridis 0-0.75 everywhere; red = scaling prediction everywhere; Okabe-Ito for sizes |
| Final width 7.0 in, 8 pt | yes; Fig. 6 is 4.3 in (could be single column at 3.375 in) |
| Not yet | Fig. 1(c) half-dimer sites are drawn at plaquette centres, not at their boundary vertices; frozen phases are not coloured by ice phase (a phase map would show the three-phase pockets of Fig. 3); compiled-PDF check at 100 % |

## Appendix candidates (existing working figures)
- `figs/lshape_vs_cps.png`: reflex-corner condition reproduces Colomo-Pronko-Sportiello (2018).
- `figs/plus_verify_n241.png`, width scaling (`width_scaling.py`): boundary width ~ n^0.36.
- 256-rectangle table (`mixed_digest.py`).
