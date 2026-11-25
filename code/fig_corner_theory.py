"""Figure: exact liquidity map of the corner-switch sample versus the hodograph solution (frozen boundary and edge zones).
Usage: python fig_corner_theory.py L r   (y = 2^(r m))"""
import sys, math, warnings, pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import brentq
warnings.filterwarnings('ignore')
from gen_ilp import Sample, halfplanes
from bmono import BMono
from mixed_scan import state
from corner_theory import Corner
from corner_field import field
L=int(sys.argv[1]); rs=[float(v) for v in sys.argv[2:]]
S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); B=BMono(S.cells,S.free); m=B.m
fig,ax=plt.subplots(1,len(rs),figsize=(4.2*len(rs),4.4),squeeze=False)
for a,r in zip(ax[0],rs):
    k=round(r*m); phi=k*math.log(2)/m; st=state(B,S,2**k)
    b=brentq(lambda x: field(Corner(x),n=120)-phi,0.05,60,xtol=1e-3); c=Corner(b); f=c.free_fraction()[0]
    a.scatter(st['cen'][:,0],st['cen'][:,1],c=1-st['mx'],s=6*(48/L)**2,marker='s',cmap='viridis',vmin=0,vmax=0.5,lw=0)
    sc=L+1
    for (u,v) in c.curve(n=300):
        a.plot(sc*u,sc*v,color='#ff2d55',lw=2.2,zorder=5); a.plot(sc*v,sc*u,color='#ff2d55',lw=2.2,zorder=5)
    for X,col in ((c.X1(),'k'),(c.X2(),'#f39c12')):
        a.plot([sc*X],[0],'^',color=col,ms=8,clip_on=False); a.plot([0],[sc*X],'>',color=col,ms=8,clip_on=False)
    a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_xlim(-1,sc+1); a.set_ylim(-1,sc+1)
    a.set_title(f'(J/2T)/m = {phi:.3f}, L = {L}\nfree fraction: exact {st["nmono"]/m:.3f}, theory {f:.3f}',fontsize=9)
    print(f'r={r}: field {phi:.4f}, b={b:.4f}, X1={c.X1():.4f}, X2={c.X2():.4f}, free fraction exact {st["nmono"]/m:.4f} theory {f:.4f}',flush=True)
plt.savefig('../figs/corner_theory.png',dpi=170,bbox_inches='tight')
