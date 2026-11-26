"""Designed samples: unions of axis rectangles with through-block (odd) edges, fringed termination.  Reports freeable plaquettes,
colour balance and tileability; returns the exact liquidity map (Kasteleyn, all-paired ground state)."""
import math, numpy as np
from gen_ilp import Sample, black
from eff import ntil
from maps import cell_entropy, pcenter
def union(rects):
    return lambda p: any(x0<=p[0]<=x1 and y0<=p[1]<=y1 for x0,x1,y0,y1 in rects)
def build(rects,W=None):
    W=W or max(max(abs(v) for v in r) for r in rects)+4
    S=Sample(union(rects),W=W,mode='A'); nb=sum(black(c) for c in S.cells)
    return S, 2*nb-len(S.cells)
def check(name,rects):
    S,imb=build(rects); lt=ntil(S.cells) if imb==0 else None
    print(f'{name}: cells {len(S.cells)}, freeable {len(S.free)}, H-V {imb}, ln#tilings/ln2 = {None if lt is None else round(lt/math.log(2),1)}')
    return S
if __name__=='__main__':
    s=11
    for g in (3,5,7):
        for w in (3,5,7,9):
            check(f'dumbbell s={s} g={g} w={w}',[(-2*s-g,-g,-s,s),(g,2*s+g,-s,s),(-g,g,-w,w)])
