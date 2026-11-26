"""Local statistics of the corrected boundary model (freeable plaquettes weight 1, half-dimer sites weight n t, dimers t = exp(-J/T))
in floating point, with the row-sum check.  Prints free fraction, half-dimer occupation and liquid fraction versus phi = (J/2T)/m."""
import sys, math, numpy as np, pickle
from gen_ilp import Sample, halfplanes
from bmono import BMono
import corrected_model as cm
def stats(S,n,t):
    sites=set(S.free)|set(n); B=BMono(S.cells,sites); N=B.N; m=B.m; assert not B.odd
    A=B.K*t; ell=np.array([(1.0 if c in S.free else 0.0)+n.get(c,0)*t for c in B.order[:m]]); G=B.g
    i,j=np.triu_indices(m,1); v=G[i]*G[j]*ell[i]*ell[j]; A[i,j]+=v; A[j,i]-=v
    Ai=np.linalg.inv(A); P=(B.K*t)*Ai.T; Mm=A[:m,:m]-B.K[:m,:m]*t; mono=np.einsum('ij,ji->i',Mm,Ai[:m,:m])
    tot=P.sum(1); tot[:m]+=mono; mx=P.max(1); mx[:m]=np.maximum(mx[:m],mono)
    fr=np.array([c in S.free for c in B.order[:m]])
    # split monomer probability of a freeable plaquette that also has half-dimer sites: free part = ell_free/ell share
    share=np.where(fr,1.0/np.maximum(ell,1e-300),0.0)
    return dict(err=abs(tot-1).max(), free=(mono*share).sum()/fr.sum(), half=(mono*(1-share)).sum(), liq=(mx<0.9).mean(), cen=np.array([B.cen[c] for c in B.order]), mx=mx, B=B)
if __name__=='__main__':
    L=int(sys.argv[1]); S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); n=cm.halfsites(S,'fit'); m=len(S.free); out={}
    print(f'L={L}: m={m}, half-dimer weight {sum(n.values())}')
    for phi in [float(v) for v in sys.argv[2:]]:
        t=math.exp(-2*phi*m); st=stats(S,n,t); out[phi]=st
        print(f'  phi={phi:5.3f}: free fraction {st["free"]:.3f}, occupied half-dimer sites {st["half"]:.2f}, liquid fraction {st["liq"]:.3f}, row-sum error {st["err"]:.1e}')
    pickle.dump({k:{kk:vv for kk,vv in v.items() if kk!='B'} for k,v in out.items()},open(f'data/corrected_map_{L}.pkl','wb'))
