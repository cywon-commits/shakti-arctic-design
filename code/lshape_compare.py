"""Our degree-2 tangent family for the L-shaped sample versus the curve of Colomo-Pronko-Sportiello (2018)."""
import sys, numpy as np
from scipy.spatial import cKDTree
from lshape_theory import solve, curve
from lshape_cps import cps_curve
for a in [float(s) for s in sys.argv[1:]] or [0.4,0.6]:
    q,t,res=solve(1,a,a); th,v,x,y=curve(q,t,m=200001)
    P=cps_curve((1+a)/(1-a)); P=P[(abs(P[:,0])<2)&(abs(P[:,1])<2)]
    best=None
    for sw in (0,1):
        for sx in (1,-1):
            for sy in (1,-1):
                A,B=(x,y) if not sw else (y,x); C=np.c_[(1+sx*A)/2,(1+sy*B)/2]
                d1,_=cKDTree(C).query(P[:,:2]); d2,_=cKDTree(P[:,:2]).query(C)
                if best is None or d1.max()<best[0]: best=(d1.max(),np.median(d1),d2.max(),sw,sx,sy)
    print(f'a={a}: q={q:.5f}, residual {res:.1e}; CPS -> ours: max {best[0]:.2e}, median {best[1]:.2e}; ours -> CPS: max {best[2]:.2e}; map (swap,sx,sy)={best[3:]}',flush=True)
