"""Decomposition of the low-lying vertex-model states of the corner-switch sample:
counts by total energy E (units J) and by the number of cost-3 vertices in the bulk (full coordination) and on the boundary.
Same transfer matrix as vertex_dos.py with two extra counters (each kept up to NMAX)."""
import sys, itertools, math, pickle, numpy as np
from collections import defaultdict
from gen_ilp import Sample, halfplanes, COST
from vertex_dos import lab
with __import__('contextlib').redirect_stdout(__import__('io').StringIO()):
    from vort import is_island
NMAX=2
def fulldeg(p):
    x,y=p; return sum(is_island(e) for e in (((x,y),(x+1,y)),((x-1,y),(x,y)),((x,y),(x,y+1)),((x,y-1),(x,y))))
def decomp(S,Emin,K):
    inc=defaultdict(list)
    for k,(a,b) in enumerate(S.isl): inc[a].append((k,'start',a,b)); inc[b].append((k,'end',a,b))
    pts=sorted(inc,key=lambda p:(p[1],p[0]))
    cap=Emin+K; A=np.zeros((cap+1,NMAX+1,NMAX+1)); A[0,0,0]=1.0; axes=[]
    for p in pts:
        L=inc[p]; ks=[k for k,_,_,_ in L]; old=[k for k in ks if k in axes]; new=[k for k in ks if k not in axes]
        keep=[k for k in axes if k not in old]; naxes=keep+new; bulk=len(L)==fulldeg(p)
        B=np.zeros((2,)*len(naxes)+(cap+1,NMAX+1,NMAX+1))
        Am=np.moveaxis(A,[axes.index(k) for k in old],range(len(old))) if old else A
        for bits in itertools.product((1,-1),repeat=len(L)):
            c=COST[lab(L,bits)]
            if c>cap: continue
            val=dict(zip(ks,bits)); src=Am[tuple(0 if val[k]==1 else 1 for k in old)]
            dst=B[(slice(None),)*len(keep)+tuple(0 if val[k]==1 else 1 for k in new)]
            if c==0: dst+=src
            elif c==3:
                if bulk: dst[...,c:,1:,:]+=src[...,:-c,:-1,:]
                else: dst[...,c:,:,1:]+=src[...,:-c,:,:-1]
            else: dst[...,c:,:,:]+=src[...,:-c,:,:]
        A,axes=B,naxes
    return A.reshape(-1,cap+1,NMAX+1,NMAX+1).sum(0)[Emin:]
if __name__=='__main__':
    K=int(sys.argv[1]); out={}
    for L in [int(v) for v in sys.argv[2:]]:
        S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
        G=decomp(S,Egs,K); out[L]=G; g0=G[0].sum()
        print(f'L={L}: m={m}, E_gs={Egs}')
        for k in range(K+1):
            print(f'   E_gs+{k}: total {G[k].sum()/g0:.5g};  by (n3_bulk, n3_edge): '+'  '.join(f'({a},{b}):{G[k,a,b]/g0:.4g}' for a in range(NMAX+1) for b in range(NMAX+1) if G[k,a,b]>0))
        pickle.dump(out,open('data/vertex_decomp.pkl','wb'))
