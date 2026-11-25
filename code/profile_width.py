"""Corner-switch sample: free-plaquette probability along the block-aligned edges at fixed scaled field (J/2T)/m,
for several sizes.  Tests whether the partially free zone narrows relative to the edge (step profile in the scaling limit).
Exact rational arithmetic; only the freeable block of the inverse is computed.  Usage: python profile_width.py r L1 L2 ...  (k = r*m)"""
import sys, math, time, pickle, numpy as np, flint
from gen_ilp import Sample, halfplanes
from bmono import BMono
def mono_exact(B,y):
    m,N=B.m,B.N; assert not B.odd
    A=[[0]*N for _ in range(N)]; ii,jj=np.nonzero(np.triu(B.K,1))
    for i,j in zip(ii,jj): v=int(B.K[i,j]); A[i][j]=v; A[j][i]=-v
    G=[int(g) for g in B.g]
    for a in range(m):
        for b in range(a+1,m): v=G[a]*G[b]*y*y; A[a][b]+=v; A[b][a]-=v
    M=flint.fmpq_mat(flint.fmpz_mat(N,N,[v for r in A for v in r]))
    E=flint.fmpq_mat(N,m,[1 if r==c else 0 for r in range(N) for c in range(m)]); X=M.solve(E)      # X[:, b] = A^{-1} e_b
    p=np.zeros(m)
    for a in range(m):
        s=flint.fmpq(0)
        for b in range(m):
            if a!=b: s+=(G[a]*G[b]*y*y*(1 if a<b else -1))*X[b,a]*(-1)     # M_ab (A^{-1})_ba, (A^{-1})_ba = -X[a,b]... use antisymmetry: (A^{-1})_ba = X[b,a]
        p[a]=-int(s.p)/int(s.q)
    return p
if __name__=='__main__':
    r=float(sys.argv[1]); out={}
    for L in [int(v) for v in sys.argv[2:]]:
        S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); B=BMono(S.cells,S.free); m=B.m; k=round(r*m); t=time.time()
        p=mono_exact(B,1<<k); cen=np.array([B.cen[c] for c in B.order[:m]])
        w=(p*(1-p)).sum()/m                      # width of the partially free zone, as a fraction of the total switchable length
        out[L]=(k,p,cen); print(f'L={L} m={m} k={k} (J/2T)/m={k*math.log(2)/m:.3f}: free fraction {p.sum()/m:.4f}, sum p(1-p)/m = {w:.4f}, sum p(1-p) = {w*m:.3f}   [{time.time()-t:.0f}s]',flush=True)
        print('    ',' '.join(f'{v:.2f}' for v in p),flush=True)
        pickle.dump(out,open(f'data/profile_width_{r}.pkl','wb'))
