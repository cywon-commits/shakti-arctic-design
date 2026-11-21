"""Numerical precision of the Kasteleyn inverse: every cell's pairing probabilities must sum to 1."""
import sys, numpy as np
from gen_ilp import Sample, black
def plus(n,c): return lambda p: (abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n)
for a in sys.argv[1:]:
    n,c=map(int,a.split(',')); S=Sample(plus(n,c),W=n+6,mode='B').cells
    B=sorted(q for q in S if black(q)); Wt=sorted(q for q in S if not black(q)); wi={q:k for k,q in enumerate(Wt)}
    K=np.zeros((len(B),len(B)))
    for k,q in enumerate(B):
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            m=(q[0]+d[0],q[1]+d[1])
            if m in wi: K[k,wi[m]]=1.0 if d[1]==0 else (-1.0)**(q[0])
    Ki=np.linalg.inv(K); P=K*Ki.T
    rb=np.abs(P.sum(1)-1).max(); rw=np.abs(P.sum(0)-1).max()
    print(f'n={n} c={c}: max |sum P - 1| black {rb:.2e}, white {rw:.2e}; min P {P.min():.2e}; cond(K) ~ {np.linalg.cond(K,1):.2e}',flush=True)
