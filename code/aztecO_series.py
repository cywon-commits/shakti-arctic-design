"""sigma_N(x) of the through-block Shakti Aztec square (corrected model) for a series of N, and a 1/N extrapolation."""
import sys, math, numpy as np, pickle
from aztecO import shakti, sigma_shakti
xs=[float(v) for v in sys.argv[1].split(',')]; out={}
for L in [int(v) for v in sys.argv[2:]]:
    N=L//2-1; S,c,E0=shakti(L); out[N]=[sigma_shakti(c,E0,N,x) for x in xs]
    print(f'N={N:3d}: '+'  '.join(f'x={x}: {s:.5f}' for x,s in zip(xs,out[N])),flush=True)
pickle.dump((xs,out),open('data/aztecO_series.pkl','wb'))
Ns=sorted(out)
for i,x in enumerate(xs):
    for a,b in zip(Ns[-3:],Ns[-2:]) if len(Ns)>=3 else []:
        pass
    a,b=Ns[-2],Ns[-1]; ext=(b*out[b][i]-a*out[a][i])/(b-a)
    print(f'x={x}: two-point 1/N extrapolation from N={a},{b}: {ext:.4f}')
