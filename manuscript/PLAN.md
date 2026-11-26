# Manuscript plan (writing guide W1-W2), draft v0.1, 2026-11-26

Target: Phys. Rev. B, Regular Article (`revtex4-2`, `prb,reprint`). Companion papers: the vorticity paper (Kim et al., Shakti, T_on)
and "Thermal vacancies melt arctic order" (Won; goes first, cited for the hodograph method).

## Core sentence (W1)

In Shakti artificial spin ice the cut and termination of the array fix the boundary height of the underlying domino tiling,
so the liquid region can be designed: convex cuts move and resize the arctic circle, reflex corners reshape it into an
algebraic curve, and two adjacent block-aligned edges switch it on and off with temperature through corner sources, with
the sample frozen below J/(2T) = m ln(1 + sqrt 2).

Abstract result sentence, last paragraph of the introduction, and first sentence of the conclusion must say this with the
same strength.

## Claims (W1) and evidence grade

Grades: [proof] closed argument; [derived] derivation under stated assumptions; [numerical] exact computation only;
[argument] physical argument with gaps.

| # | Claim | Grade | Where | Evidence |
|---|---|---|---|---|
| C1 | Boundary height of a cut sample follows from counting dangling unit edges (Eq. 2); only through-block axis cuts and even 45-degree cuts are extremal | [proof] for the rule, [numerical] for the table | Sec. II B, Table I, Fig. 1(b) | exact periodic counting; Kasteleyn counts (Aztec 2^{N(N+1)/2}, unique tiling of 45-degree squares) |
| C2 | Stem-island rule: ground states are dimer coverings in which plaquettes without stem island stay unpaired and gain J/2 | [numerical] | Sec. II C, Fig. 1(c) | integer programming of the vertex model, 44 samples, no mismatch |
| C3 | Box rule; convex corner cuts leave a square whose inscribed circle is the liquid region; four cuts of mixed offset class admit no tiling | [proof] + [numerical] | Sec. III, Fig. 2, App. A | 8500 octagons, no mismatch |
| C4 | Plus: tileable iff c >= n/3; circle for c >= n/sqrt 2 | [proof]; [argument] | Sec. IV | Lipschitz bound; restriction of the Aztec solution |
| C5 | Plus, n/3 < c < n/sqrt 2: class-10 curve in closed form, Eq. (8) | [derived] (Kenyon-Okounkov structure + reflex-corner condition) | Sec. IV, Fig. 3, App. B | urban renewal up to n = 481; L-shape reproduces Colomo-Pronko-Sportiello |
| C6 | Among rectangles only two adjacent block-aligned edges give a unique ground state that melts | [numerical] | Sec. V B, Table II, Fig. 4(a) | 256 rectangles |
| C7 | Finite T: boundary model with free plaquettes and half-dimer sites reproduces the vertex model's first excited level | [numerical] | Sec. V A, App. D, Fig. 6 | transfer-matrix density of states, L <= 20 |
| C8 | Edge energies are total derivatives of the height: E_bd = (J/4)[h(NE) - h(SW)] + const (corner switch), alternating corners (through-block square) | [derived] | Sec. V C, Fig. 5(a) | D2 + colour count |
| C9 | Scaling solution of the corner switch; freezing at q = ln(1 + sqrt 2) | [derived] (domain structure assumed) | Sec. V D, Fig. 4, App. E | exact Pfaffian maps L = 32-64; first-level growth rate |
| C10 | Through-block square = boundary-vacancy model with J/T = (xN + 2) ln 2 | [numerical] + [argument] | Sec. VI, Fig. 5, Table IV | first-level counts identical; sigma_N extrapolates to sigma_inf |

Not claimed: proof of the reflex-corner condition; proof that bulk defects vanish in the scaling limit; closed form of the
function cal-K; design gallery (dropped: no separate message).

## Figures (W3)

| Fig. | Claim | Width |
|---|---|---|
| 1 | C1, C2: the cut sets the boundary | 7.0 in |
| 2 | C3: convex cuts move and resize the circle | 7.0 in |
| 3 | C4, C5: reflex corners reshape it | 7.0 in |
| 4 | C6, C9: temperature switch | 7.0 in |
| 5 | C8, C10: corner sources; through-block square = vacancy model | 7.0 in |
| 6 | C7: boundary model vs vertex model | 4.3 in (appendix in v0.1; placement to decide) |

## Outline (W2)

I Introduction. II Shakti ground states as domino tilings (A plaquettes and heights, B boundary heights of cuts, C stem
islands). III Convex cuts move the arctic circle. IV Reflex corners reshape it. V Switching with temperature (A boundary
model, B which rectangles switch, C corner sources, D scaling solution). VI Through-block square and boundary vacancies.
VII Discussion (design rules, validity, experiments). VIII Conclusion. Appendices A box rule and octagon rule; B tangent
family and closed form; C boundary Pfaffian; D vertex-model check; E hodograph solution and fixed-edge variant;
F numerical methods.

## Open decisions for the author

1. Fig. 6 in the main text (Sec. V A) or in Appendix D (v0.1: appendix).
2. Author list (v0.1: corresponding author only).
3. The vorticity paper states that arrays cut through the blocks are free of the T_on mechanism; Sec. VI shows that
   their half-dimer sites give the same temperature scale. One of the two texts must be adjusted.
4. Zenodo DOI for the repository at submission.
