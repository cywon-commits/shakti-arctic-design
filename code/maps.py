"""Exact local dimer statistics (Kasteleyn) for the all-paired problem on a Shakti sample; liquid/frozen map in Shakti coordinates."""
from eff import *
import scipy.sparse as sp, scipy.sparse.linalg as spl
def pcenter(c):
    P=ppts(c); return (sum(p[0] for p in P)/len(P), sum(p[1] for p in P)/len(P))
def edge_probs(S):
    B=sorted(c for c in S if black(c)); Wt=sorted(c for c in S if not black(c)); assert len(B)==len(Wt)
    bi={c:k for k,c in enumerate(B)}; wi={c:k for k,c in enumerate(Wt)}
    K=np.zeros((len(B),len(B)),dtype=complex)
    for c in B:
        for d,w in [((1,0),1),((-1,0),1),((0,1),1j),((0,-1),1j)]:
            n=(c[0]+d[0],c[1]+d[1])
            if n in wi: K[bi[c],wi[n]]=w
    Ki=np.linalg.inv(K); P={}
    for c in B:
        for d in [(1,0),(-1,0),(0,1),(0,-1)]:
            n=(c[0]+d[0],c[1]+d[1])
            if n in wi: P[(c,d)]=(K[bi[c],wi[n]]*Ki[wi[n],bi[c]]).real
    s,ld=np.linalg.slogdet(K)
    return P,ld
def cell_entropy(S):
    """per-cell 'liquidity' = 1 - max_d P(cell paired in direction d)"""
    P,ld=edge_probs(S); mx=defaultdict(float)
    for (c,d),p in P.items():
        n=(c[0]+d[0],c[1]+d[1]); mx[c]=max(mx[c],p); mx[n]=max(mx[n],p)
    return {c:1-mx[c] for c in S}, ld
