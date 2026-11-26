"""Two through-block squares (half side s) overlapping along the diagonal with offset o: liquidity along the diagonal and at the neck."""
import math, pickle, numpy as np
from gen_ilp import Sample
from maps import cell_entropy, pcenter
s=11; out={}
for o in (8,12,14,16,18,20):
    ins=lambda p: (abs(p[0])<=s and abs(p[1])<=s) or (abs(p[0]-o)<=s and abs(p[1]-o)<=s)
    S=Sample(ins,W=o+s+4,mode='A'); liq,ld=cell_entropy(S.cells); P=np.array([pcenter(c) for c in S.cells]); v=np.array([liq[c] for c in S.cells])
    dg=abs(P[:,0]-P[:,1])<1.5; t=(P[dg,0]+P[dg,1])/2; prof=v[dg][np.argsort(t)]; t=np.sort(t)
    neck=v[dg & (abs((P[:,0]+P[:,1])/2-o/2)<1.5)].max()
    out[o]=dict(cen=P,liq=v,t=t,prof=prof,neck=neck); print(f'o={o}: overlap {2*s-o}, max liquidity at the neck {neck:.3f}',flush=True)
pickle.dump(out,open('data/design_valve.pkl','wb'))
