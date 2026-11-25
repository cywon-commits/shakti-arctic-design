"""Figure: degree-2 tangent family (line) versus Colomo-Pronko-Sportiello (dots) for the L-shaped sample."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from lshape_theory import solve, curve
from lshape_cps import cps_curve
fig,ax=plt.subplots(1,2,figsize=(9,4.6))
for k,a in enumerate((0.4,0.6)):
    q,t,res=solve(1,a,a); th,v,x,y=curve(q,t)
    P=cps_curve((1+a)/(1-a),m=400); P=P[(abs(P[:,0])<2)&(abs(P[:,1])<2)]
    ax[k].plot([-1,1,1,a,a,-1,-1],[-1,-1,a,a,1,1,-1],color='0.4',lw=1)
    ax[k].plot([a,a],[-1,a],color='0.7',lw=.8,ls=':'); ax[k].plot([-1,a],[a,a],color='0.7',lw=.8,ls=':')
    ax[k].plot(x,y,color='#2a6fb0',lw=1.4,label='tangent family, d = 2')
    ax[k].plot(2*P[:,0]-1,1-2*P[:,1],'.',color='#c0392b',ms=3,label='Colomo-Pronko-Sportiello (4.15)')
    ax[k].set_aspect('equal'); ax[k].set_title(f'a/n = {a}'); ax[k].set_xlabel('X/n'); ax[k].set_ylabel('Y/n')
ax[0].legend(fontsize=8,loc='lower left'); plt.tight_layout(); plt.savefig('../figs/lshape_vs_cps.png',dpi=200)
