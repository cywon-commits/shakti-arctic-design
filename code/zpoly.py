"""Exact polynomial Z(y) = sum_j c_j y^j of the boundary-monomer model: c_j = number of states with j free plaquettes.
Z(2^k) is the exact integer square root of det(K + M); with k large the coefficients are its base-2^k digits."""
import sys, math, pickle, numpy as np, flint
from gen_ilp import Sample, halfplanes
from bmono import BMono
def zpoly(B,k=None):
    m,N=B.m,B.N; assert not B.odd
    k=k or int(0.75*N)+64; y=1<<k
    A=[[0]*N for _ in range(N)]
    ii,jj=np.nonzero(np.triu(B.K,1))
    for i,j in zip(ii,jj): v=int(B.K[i,j]); A[i][j]=v; A[j][i]=-v
    G=[int(g) for g in B.g]
    for a in range(m):
        for b in range(a+1,m): v=G[a]*G[b]*y*y; A[a][b]+=v; A[b][a]-=v
    d=int(flint.fmpz_mat(N,N,[v for r in A for v in r]).det()); z=math.isqrt(d); assert z*z==d
    c=[]
    while z: c.append(z&(y-1)); z>>=k
    assert max(c)<(1<<(k-2)),'k too small'
    return c
if __name__=='__main__':
    out={}
    for L in [int(v) for v in sys.argv[1:]]:
        S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); B=BMono(S.cells,S.free); c=zpoly(B); out[L]=c
        m=B.m; print(f'L={L} cells {B.N} m={m}: top coefficients c_m, c_m-2, ... = {c[m]}, {c[m-2]}, {c[m-4] if m>=4 else ""}, {c[m-6] if m>=6 else ""};  c_0 = {c[0]}',flush=True)
        print(f'      ln c_(m-2) / L = {math.log(c[m-2])/L:.4f}   ln c_0/cells = {math.log(c[0])/B.N:.4f}',flush=True)
        pickle.dump(out,open('data/zpoly_corner.pkl','wb'))
