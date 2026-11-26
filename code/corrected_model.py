"""Corrected effective model: boundary-monomer dimer model whose monomer sites are the freeable plaquettes (weight 1, energy 0)
and the plaquettes owning a boundary half-dimer site (an unhappy z3 vertex on a boundary point that is not shared by two adjacent plaquettes;
weight t = exp(-J/T) per site).  Dimers carry t.  Exact polynomial in t by integer evaluation of the Pfaffian; compared with the vertex model."""
import sys, math, numpy as np, flint
from collections import defaultdict
from gen_ilp import Sample, halfplanes, ppts
from bmono import BMono
from vertex_dos import dos
from vertex_decomp import fulldeg
def halfsites(S,variant):
    deg=defaultdict(int)
    for a,b in S.isl: deg[a]+=1; deg[b]+=1
    touch=defaultdict(list)
    for c in S.cells:
        for p in ppts(c): touch[p].append(c)
    n=defaultdict(int)
    for p,d in deg.items():
        if d!=3: continue
        t=touch[p]
        if len(t)==1: n[t[0]]+=(2 if fulldeg(p)==3 else 1) if variant=='fit' else 1
        elif len(t)==2 and abs(t[0][0]-t[1][0])+abs(t[0][1]-t[1][1])!=1:
            if variant=='both': n[t[0]]+=1; n[t[1]]+=1
            elif variant=='nonfree':
                for c in t:
                    if c not in S.free: n[c]+=1
            elif variant=='none': pass
            elif variant=='fit':
                for c in t:
                    if c not in S.free: n[c]+=1
    return n
def poly(S,n,k=None,a=1,b=1):
    sites=set(S.free)|set(n); B=BMono(S.cells,sites); N=B.N; m=B.m
    k=k or 3*N+64+N*max(a,b).bit_length(); t=1<<k; nn=N+1 if B.odd else N
    A=[[0]*nn for _ in range(nn)]
    ii,jj=np.nonzero(np.triu(B.K,1))
    for i,j in zip(ii,jj): v=int(B.K[i,j])*t; A[i][j]=v; A[j][i]=-v
    ell=[(a if c in S.free else 0)+b*n.get(c,0)*t for c in B.order[:m]]; G=[int(g) for g in B.g]
    lab=list(range(m)); 
    if B.odd: lab.append(N); ell.append(1); G.append(1)
    for a in range(len(lab)):
        for b in range(a+1,len(lab)):
            v=G[a]*G[b]*ell[a]*ell[b]; A[lab[a]][lab[b]]+=v; A[lab[b]][lab[a]]-=v
    d=int(flint.fmpz_mat(nn,nn,[v for r in A for v in r]).det()); z=math.isqrt(d); assert z*z==d
    c=[]
    while z: c.append(z&(t-1)); z>>=k
    return c
if __name__=='__main__':
    K=int(sys.argv[1]); variant=sys.argv[2]
    for L in [int(v) for v in sys.argv[3:]]:
        S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
        n=halfsites(S,variant); c=poly(S,n); g,_,_=dos(S,Egs+K)
        print(f'L={L} variant={variant}: half-dimer sites {sum(n.values())} on {len(n)} plaquettes')
        for e in range(Egs,Egs+K+1):
            ce=c[e] if e<len(c) else 0
            print(f'   E_gs+{e-Egs}: vertex g/2 = {g[e]/2:.8g}   model = {ce:.8g}   ratio {g[e]/2/ce if ce else float("nan"):.6f}')
