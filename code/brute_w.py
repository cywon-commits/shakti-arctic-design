"""Exhaustive check of the Pfaffian with per-plaquette monomer weights (corrected_model.poly) on small samples."""
import sys
from functools import lru_cache
from gen_ilp import Sample, halfplanes
from corrected_model import halfsites, poly
def brute(cells,wt):
    """wt[c] = (w0, w1): monomer weight w0 + w1 t.  Returns polynomial in t (dimers carry t)."""
    cl=sorted(cells); idx={c:k for k,c in enumerate(cl)}; n=len(cl)
    nb=[[idx[(c[0]+a,c[1]+b)] for a,b in ((1,0),(-1,0),(0,1),(0,-1)) if (c[0]+a,c[1]+b) in idx] for c in cl]
    W=[wt.get(c,(0,0)) for c in cl]
    def add(p,q,sh=0,f=1):
        r=list(p)+[0]*max(0,len(q)+sh-len(p))
        for i,v in enumerate(q): r[i+sh]+=f*v
        return r
    @lru_cache(None)
    def rec(mask):
        if mask==(1<<n)-1: return (1,)
        k=0
        while mask>>k&1: k+=1
        tot=[0]
        for j in nb[k]:
            if not mask>>j&1: tot=add(tot,rec(mask|1<<k|1<<j),1)
        w0,w1=W[k]
        if w0 or w1:
            r=rec(mask|1<<k)
            if w0: tot=add(tot,r,0,w0)
            if w1: tot=add(tot,r,1,w1)
        return tuple(tot)
    return list(rec(0))
for L in (4,8):
  for var in ('none','both','nonfree'):
    S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); n=halfsites(S,var)
    if len(S.cells)>40: continue
    wt={c:(1 if c in S.free else 0, n.get(c,0)) for c in S.cells}
    pb=brute(S.cells,wt); pp=poly(S,n)
    print(L,var,'brute==pfaffian:',pb==pp[:len(pb)] and all(v==0 for v in pp[len(pb):]), pb[:8], pp[:8])
