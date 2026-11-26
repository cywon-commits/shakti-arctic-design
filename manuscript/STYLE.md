# Style sheet (completeness guide F1, writing guide W6)

## Terms: one name per object

| Use | Meaning | Do not use |
|---|---|---|
| plaquette (H, V) | rectangular 1x2 loop of the Shakti lattice; H horizontal, V vertical | cell (except for the checkerboard picture) |
| block | 2x2 unit of two parallel plaquettes | |
| unhappy vertex | z = 3 vertex in its excited state; costs J | defect vertex |
| dimer | unhappy vertex shared by two perpendicular plaquettes; equivalently a domino of the checkerboard | tetromino (only once, as background), domino (only in "domino tiling") |
| domino tiling | the dimer covering written on the checkerboard | dimer cover |
| height, h | Thurston height function | height field (allowed for the macroscopic H) |
| slope square | allowed slopes [-1,1]^2 | Newton polygon (once, in parentheses) |
| ice phase | one of the four frozen brick phases, the corners of the slope square | frozen phase, brick phase |
| liquid region | region with nonzero entropy density | arctic region, disordered region |
| frozen region | complement of the liquid region | |
| arctic curve | boundary of the liquid region; arctic circle when circular | frozen boundary (allowed in Sec. V for the scaling solution) |
| liquidity | 1 - P_max of a plaquette | |
| liquid fraction | fraction of plaquettes with P_max < 0.9 | |
| block-aligned edge | axis cut along block boundaries (even offset) | E edge |
| through-block edge | axis cut through the middle of blocks (odd offset) | O edge |
| fringed / clean termination | islands of cut plaquettes kept / removed | tasseled |
| stem island | island ending perpendicularly at the middle of the outer long side of a plaquette | |
| freeable plaquette | edge plaquette whose stem island is missing | |
| free plaquette | freeable plaquette that is unpaired | |
| half-dimer site | boundary z = 3 vertex whose unhappy state covers a single plaquette | H site |
| monomer | unpaired plaquette in the boundary model (free plaquette or occupied half-dimer site) | vacancy (only for the model of Ref. [vacancy]) |
| boundary model | dimer model with monomers on the outer face (free plaquettes and half-dimer sites) | corrected model |
| fixed-edge variant | boundary model without half-dimer sites on N and E edges | original/uncorrected model |
| vertex model | full island-spin model with vertex energies | full model |
| corner switch | square with two adjacent block-aligned edges and two through-block fringed edges | |
| through-block square | square with four through-block fringed edges | Aztec square, OOOO |
| plus | plus-shaped sample, arm reach n, half-width c | cross |
| box rule, octagon rule | Sec. III | |
| corner sources | point sources of the height at corners | corner charges (Ref. [vacancy] term, mention once) |
| source strength q | q = J/(2mT) | phi, (J/2T)/m |
| freezing temperature T_f | J/(2T_f) = m ln(1+sqrt 2) | |
| reflex corner / convex corner | | concave corner |

## Abbreviations
None defined (ASI not used; each candidate appears fewer than three times).

## Symbols

| Symbol | Meaning | Defined |
|---|---|---|
| J | energy of an unhappy vertex | Sec. II A |
| Delta | energy of type-II and charged vertices in units of J | Sec. II A |
| T | temperature, k_B = 1 | Sec. II A |
| (X, Y) | Shakti coordinates | Sec. II A |
| h | height | Sec. II A, Eq. (1) |
| sigma(p) | +1 for H, -1 for V | Eq. (1) |
| d_NE ... | offsets of 45-degree cuts | Sec. III |
| n, c | arm reach and half-width of the plus | Sec. IV |
| tau | parameter of the closed form | Eq. (8) |
| nu, H(nu) | normal angle and support of a tangent line | Eq. (6) |
| L | side of rectangular samples | Sec. V |
| m | number of freeable plaquettes | Sec. V A |
| q | J/(2mT) | Eq. (11) |
| N | order of an Aztec diamond | Sec. VI |
| x | scaling variable of Ref. [vacancy] | Eq. (16) |
| sigma_N | 1 - <n_mono>/(2N) | Sec. VI |
| (u, v) | coordinates scaled to the unit square | Sec. V D |
| omega, zeta | complex slope, hodograph variable | Sec. V D |
| Q, P, P_1, b | auxiliary variables of the solution | Sec. V D |
| cal-K(P) | Cauchy integral, Eq. (13) | Sec. V D |
| X_1 | end of the frozen cap on an edge | Sec. V D |
| a | parameter of the fixed-edge variant | App. E |
| sans-K | Kasteleyn matrix | App. C |

## Conventions
APS cross-references: Fig. 2, Eq. (3), Sec. II, Table I, Appendix B, Ref. [5]. American spelling. "We" for the authors.
Numbers in the body only where they carry a dependence; precise values in tables, captions and appendices.
