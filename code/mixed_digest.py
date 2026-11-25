"""Digest of mixed_scan.py: outcomes grouped by which edges are block-aligned (E) or through-block (O)."""
import sys, pickle, numpy as np
from collections import defaultdict
out=pickle.load(open(f'data/mixed_scan_{sys.argv[1]}.pkl','rb')); par=lambda v:'E' if v%2==0 else 'O'
canon=lambda t: min(t[p[0]]+t[p[1]]+t[p[2]]+t[p[3]] for p in ((0,1,2,3),(1,0,2,3),(0,1,3,2),(1,0,3,2),(2,3,0,1),(3,2,0,1),(2,3,1,0),(3,2,1,0)))
G=defaultdict(list)
for k,r in out.items(): G[canon(''.join(par(v) for v in k))].append((k,r))
print('edges (two X edges, two Y edges) | samples | no state | frozen at T=0 | monomer number fixed | liquid fraction T=0 -> y=1 (range over samples with a state)')
for t,L in sorted(G.items()):
    ok=[(k,r) for k,r in L if r.get('cold')]
    fro=sum(r['S0']<1e-6 for k,r in ok); fix=sum(abs(r['hot']['nmono']-r['cold']['nmono'])<1e-6 for k,r in ok)
    lc=[r['cold']['liqfrac'] for k,r in ok]; lh=[r['hot']['liqfrac'] for k,r in ok]
    print(f'{t[:2]} {t[2:]} | {len(L):3d} | {len(L)-len(ok):3d} | {fro:3d} | {fix:3d} | '+(f'{min(lc):.2f}..{max(lc):.2f} -> {min(lh):.2f}..{max(lh):.2f}' if ok else '-'))
    if len(sys.argv)>2:
        for k,r in sorted(ok):
            c,h=r['cold'],r['hot']; print(f"      {k} box {r['box']} cells {r['cells']} H-V {2*r['H']-r['cells']:+d} freeable H{r['fH']}/V{r['nfree']-r['fH']} free T=0 {round(c['nmono'])} S0/cell {r['S0']/r['cells']:.4f} liq {c['liqfrac']:.2f}->{h['liqfrac']:.2f} <free> y=1 {h['nmono']:.2f}")
