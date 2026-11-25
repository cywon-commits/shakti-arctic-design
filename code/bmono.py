"""Boundary-monomer dimer model at finite fugacity: exact partition function and local statistics by a Pfaffian.
Finite-temperature effective model of an open Shakti sample: every freeable edge plaquette may stay unpaired with weight
y = exp(J/2T) (stem-island rule).  Monomers restricted to the outer face keep the model Pfaffian
(Giuliani, Jauslin, Lieb, J. Stat. Phys. 163, 211 (2016)):
    Z(y) = Pf(K + M),   M_ij = G_i G_j y^2  (i < j, freeable cells labelled in cyclic order along the boundary),
where G_i G_j is the Kasteleyn sign of a virtual edge drawn outside the sample between cells i and j."""
import io, contextlib, math, sys, numpy as np
from functools import lru_cache
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, black
def centre(c):
    r=rect(c); return (sum(p[0] for p in r)/2+0.5, sum(p[1] for p in r)/2+0.5)
class BMono:
    def __init__(s,cells,free,c_sign=1):
        cen={c:centre(c) for c in cells}; o=np.mean(list(cen.values()),axis=0)
        F=sorted(free,key=lambda c: math.atan2(cen[c][1]-o[1],cen[c][0]-o[0]))
        s.order=F+sorted(set(cells)-set(free)); s.idx={c:k for k,c in enumerate(s.order)}; s.m=len(F); s.N=len(s.order); s.cen=cen
        K=np.zeros((s.N,s.N))
        for c in cells:
            for d,w in (((1,0),1.0),((0,1),-1.0 if c[0]%2 else 1.0)):
                n=(c[0]+d[0],c[1]+d[1])
                if n in s.idx: K[s.idx[c],s.idx[n]]=w; K[s.idx[n],s.idx[c]]=-w
        s.K=K; s.col=np.array([1 if black(c) else -1 for c in s.order])
        s.g=s._gauge(); s.c=c_sign
    def _gauge(s):
        """G_i = product of the Kasteleyn signs along the boundary path from freeable cell 0 to freeable cell i.
        A virtual outer edge (i,j) then has Kasteleyn sign G_i G_j (the cycle it closes with the boundary encloses no cell)."""
        from collections import deque
        cells=set(s.order); bd={c for c in cells if any((c[0]+a,c[1]+b) not in cells for a in (-1,0,1) for b in (-1,0,1))}
        G=np.ones(s.m)
        for k in range(s.m-1):
            a,b=s.order[k],s.order[k+1]; prev={a:None}; q=deque([a])
            while q and b not in prev:
                c=q.popleft()
                for d in ((1,0),(-1,0),(0,1),(0,-1)):
                    n=(c[0]+d[0],c[1]+d[1])
                    if n in bd and n not in prev: prev[n]=c; q.append(n)
            assert b in prev,'boundary path not found'
            sg=1.0; c=b
            while prev[c] is not None: sg*=np.sign(s.K[s.idx[prev[c]],s.idx[c]]); c=prev[c]
            G[k+1]=G[k]*sg
        return G
    def matrix(s,y):
        m=s.m; A=s.K.copy(); i,j=np.triu_indices(m,1)
        v=s.c*s.g[i]*s.g[j]*y*y; A[i,j]+=v; A[j,i]-=v
        return A
    def lnZ(s,y):
        sg,ld=np.linalg.slogdet(s.matrix(y)); assert sg>0; return ld/2
    def stats(s,y):
        """monomer probability of each freeable cell, pairing probabilities (dict), and per-cell max pairing probability"""
        A=s.matrix(y); Ai=np.linalg.inv(A); m=s.m
        M=A[:m,:m]-s.K[:m,:m]; mono=np.einsum('ij,ji->i',M,Ai[:m,:m])
        P=s.K*Ai.T; 
        return mono,P
def brute(cells,free):
    """Z as a polynomial in y (coefficient list by number of monomers), by exhaustive recursion.  Small samples only."""
    cl=sorted(cells); idx={c:k for k,c in enumerate(cl)}; nb=[[idx[(c[0]+a,c[1]+b)] for a,b in ((1,0),(-1,0),(0,1),(0,-1)) if (c[0]+a,c[1]+b) in idx] for c in cl]
    fr=[c in free for c in cl]; n=len(cl)
    @lru_cache(None)
    def rec(mask):
        if mask==(1<<n)-1: return (1,)
        k=0
        while mask>>k&1: k+=1
        tot=[0]*(n+1)
        for j in nb[k]:
            if not mask>>j&1:
                for a,v in enumerate(rec(mask|1<<k|1<<j)): tot[a]+=v
        if fr[k]:
            for a,v in enumerate(rec(mask|1<<k)): tot[a+1]+=v
        while len(tot)>1 and tot[-1]==0: tot.pop()
        return tuple(tot)
    return rec(0)
if __name__=='__main__':
    from gen_ilp import Sample, halfplanes
    tests=[('block-aligned box 8',halfplanes((1,0,0,8),(0,1,0,8)),'A'),('block-aligned box 8x12',halfplanes((1,0,0,8),(0,1,0,12)),'A'),
           ('mixed box X block-aligned 8, Y through-block 6',halfplanes((1,0,0,8),(0,1,1,7)),'A'),('mixed box 12 x 6',halfplanes((1,0,0,12),(0,1,1,7)),'A'),
           ('through-block box 8, clean',halfplanes((1,0,1,9),(0,1,1,9)),'B'),('45-degree diamond D=4, clean',halfplanes((1,1,0,8),(1,-1,0,8)),'B'),
           ('plus n=6 c=2, clean',lambda p:(abs(p[0])<=6 and abs(p[1])<=2) or (abs(p[0])<=2 and abs(p[1])<=6),'B')]
    for name,ins,mode in tests:
        S=Sample(ins,W=14,mode=mode)
        if len(S.cells)>44: print(name,'too large for brute force',len(S.cells)); continue
        poly=brute(S.cells,S.free)
        print(f'{name}: cells {len(S.cells)}, freeable {len(S.free)}, Z(y) coefficients {poly}')
        if sum(poly)==0: continue
        for c in (1,-1):
            try:
                B=BMono(S.cells,S.free,c); d=[B.lnZ(y)-math.log(sum(v*y**a for a,v in enumerate(poly))) for y in (0.3,1.0,2.5)]
                print(f'   c={c:+d}: ln Z(Pfaffian) - ln Z(brute) at y = 0.3, 1, 2.5: {np.round(d,10)}')
            except AssertionError as e: print(f'   c={c:+d}: failed ({e})')
