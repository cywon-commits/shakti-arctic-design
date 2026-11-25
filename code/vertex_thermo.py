"""Finite-temperature test of the boundary-monomer model against the exact vertex model (vertex_dos.py).
Compares the excitation energy <E - E_gs>/(m J/2) (= 1 - free fraction in the effective model) at fixed scaled field phi = (J/2T)/m.
w: weight (number of spin states) attached to each un-freed edge plaquette, i.e. y -> y/w in the effective model."""
import sys, math, pickle, numpy as np
from gen_ilp import Sample, halfplanes
from bmono import BMono
from zpoly import zpoly
from vertex_dos import dos
K=int(sys.argv[1]); res={}
for L in [int(v) for v in sys.argv[2:]]:
    S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); B=BMono(S.cells,S.free); c=zpoly(B); m=B.m; Egs=(len(S.cells)-m)//2
    g,mo,ni=dos(S,Egs+K); g=g[Egs:]/g[Egs]; res[L]=(g,c)
    print(f'L={L} (m={m}, {ni} islands): g_k/g_0 = '+' '.join(f'{v:.4g}' for v in g[:6]),flush=True)
    print('    phi    J/T   | vertex model: excitation  weight of top level | effective: w=1    w=4')
    for phi in (0.2,0.3,0.433,0.6,0.8,1.0,1.3):
        bJ=2*phi*m; wts=g*np.exp(-bJ*np.arange(K+1)); ev=(wts*np.arange(K+1)).sum()/wts.sum()/(m/2)
        def eff(w):
            y=math.exp(phi*m)/w; z=sum(float(cj)*y**(j-m) for j,cj in enumerate(c)); return 1-sum(j*float(cj)*y**(j-m) for j,cj in enumerate(c))/z/m
        print(f'   {phi:5.3f} {bJ:6.2f}   |      {ev:8.4f}         {wts[-1]/wts.sum():8.1e}        |  {eff(1):7.4f}  {eff(4):7.4f}')
    pickle.dump(res,open('data/vertex_thermo.pkl','wb'))
