"""Exact local statistics of a plus-shaped sample at any size, without matrix inversion.
The plus is embedded in the Aztec diamond of the full square; the cells outside the plus are covered by a fixed
perfect matching and all other edges touching them get zero weight (shuffle_probs.py)."""
import sys, math, time, pickle, numpy as np
from collections import defaultdict
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from gen_ilp import Sample, black, ppts
from shuffle_probs import region_probabilities
def region(inside,W): return Sample(inside,W=W,mode='B').cells
def plus(n,c): return lambda p: (abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n)
def perfect_matching(cells):
    B=sorted(c for c in cells if black(c)); W=sorted(c for c in cells if not black(c)); wi={c:k for k,c in enumerate(W)}
    if len(B)!=len(W): return None
    rows=[];cols=[]
    for k,(i,j) in enumerate(B):
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            m=(i+d[0],j+d[1])
            if m in wi: rows.append(k); cols.append(wi[m])
    G=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(B),len(W)))
    mt=maximum_bipartite_matching(G,perm_type='column')
    if (mt<0).any(): return None
    M={}
    for k,b in enumerate(B): M[b]=W[mt[k]]; M[W[mt[k]]]=b
    return M
def plus_probs(n,c):
    S=region(plus(n,c),W=n+6); big=region(lambda p: abs(p[0])<=n and abs(p[1])<=n,W=n+6)
    assert S<=big
    cx=np.mean([a+.5 for a,b in big]); cy=np.mean([b+.5 for a,b in big]); assert abs(cx-round(cx))<1e-9 and abs(cy-round(cy))<1e-9
    cx,cy=int(round(cx)),int(round(cy)); N=int(round(max(abs(a+.5-cx)+abs(b+.5-cy) for a,b in big)))
    assert len(big)==2*N*(N+1),(len(big),N)
    M=perfect_matching(big-S); assert M is not None,'complement of the plus has no perfect matching'
    return S,region_probabilities(S,N,cx,cy,M),N
def cell_stats(S,P):
    mx=defaultdict(float); tot=defaultdict(float)
    for e,p in P.items():
        for q in e: mx[q]=max(mx[q],p); tot[q]+=p
    return mx,tot
if __name__=='__main__':
    for a in sys.argv[1:]:
        n,c=map(int,a.split(',')); t=time.time(); S,P,N=plus_probs(n,c); mx,tot=cell_stats(S,P)
        err=max(abs(tot[q]-1) for q in S)
        print(f'n={n} c={c}: Aztec order {N}, plus cells {len(S)}, max |sum P - 1| = {err:.1e}, {time.time()-t:.1f}s',flush=True)
        if n<=41:      # compare with the Kasteleyn inverse where that is reliable
            from verify_big import maxprob
            mk,ld=maxprob(S); print('   max |P_max(shuffle) - P_max(Kasteleyn)| =',f'{max(abs(mx[q]-mk[q]) for q in S):.1e}')
