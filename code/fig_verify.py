"""Figure: exact liquidity map (Kasteleyn) with the class-10 curve. Usage: python fig_verify.py 61   (after verify_big.py 61 21 31 41).
Cells whose probabilities do not sum to 1 within 1e-4 (double-precision failure) are drawn in light grey."""
import pickle, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from theory import curve
n=int(sys.argv[1]); D=pickle.load(open(f'data/verify_n{n}.pkl','rb'))
fig,ax=plt.subplots(1,len(D),figsize=(6.2*len(D),6.6)); ax=np.atleast_1d(ax)
for a,(c,d) in zip(ax,sorted(D.items())):
    q=[k for k in d['cen'] if d['err'][k]<1e-4]; b=[k for k in d['cen'] if d['err'][k]>=1e-4]
    sc=a.scatter([d['cen'][k][0] for k in q],[d['cen'][k][1] for k in q],c=[1-d['mx'][k] for k in q],s=(270/n)**2,marker='s',cmap='viridis',vmin=0,vmax=.6,linewidths=0)
    if b: a.scatter([d['cen'][k][0] for k in b],[d['cen'][k][1] for k in b],c='0.8',s=(270/n)**2,marker='s',linewidths=0)
    th,nu_,cx,cy=curve(*d['sol'],m=20001); a.plot(cx,cy,color='crimson',lw=1.1)
    a.set_aspect('equal'); a.set_xlim(-n-3,n+3); a.set_ylim(-n-3,n+3); a.set_xlabel('X (Shakti units)'); a.set_ylabel('Y')
    a.set_title(f'n = {n}, c = {c} (c/n = {c/n:.2f})'+(f'; {len(b)} cells unreliable' if b else ''),fontsize=10)
fig.colorbar(sc,ax=ax,shrink=.75,label='liquidity 1 - P_max (0 = frozen)')
plt.savefig(f'../figs/plus_verify_n{n}.png',dpi=110,bbox_inches='tight')
