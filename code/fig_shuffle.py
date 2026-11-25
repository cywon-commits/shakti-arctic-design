"""Figure: exact liquidity maps at n = 241 with the class-10 curve, and a zoom on a reflex corner.
Uses data/verify_shuffle.pkl and data/verify_shuffle_shapes.pkl (verify_shuffle.py)."""
import pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from theory import curve
D={}
for f in ('data/verify_shuffle_shapes.pkl','data/verify_shuffle.pkl'): D.update(pickle.load(open(f,'rb')))
keys=[(241,83),(241,121),(241,145)]
fig,ax=plt.subplots(1,4,figsize=(19,5.2))
for a,k in zip(ax,keys):
    d=D[k]; n,c=k; sc=a.scatter(d['cen'][:,0],d['cen'][:,1],c=d['liq'],s=0.35,marker='s',cmap='viridis',vmin=0,vmax=.6,linewidths=0,rasterized=True)
    th,v,x,y=curve(*d['sol'],m=8001); a.plot(x,y,color='crimson',lw=0.8)
    a.set_aspect('equal'); a.set_xlim(-n-6,n+6); a.set_ylim(-n-6,n+6); a.set_title(f'n = {n}, c = {c} (c/n = {c/n:.2f})',fontsize=10); a.set_xlabel('X (Shakti units)')
ax[0].set_ylabel('Y')
a=ax[3]; d=D[(241,83)]; n,c=241,83; m=(abs(d['cen'][:,0]-c)<34)&(abs(d['cen'][:,1]-c)<34)
a.scatter(d['cen'][m,0],d['cen'][m,1],c=d['liq'][m],s=14,marker='s',cmap='viridis',vmin=0,vmax=.6,linewidths=0,rasterized=True)
th,v,x,y=curve(*d['sol'],m=20001); a.plot(x,y,color='crimson',lw=1.0)
a.set_aspect('equal'); a.set_xlim(c-34,c+34); a.set_ylim(c-34,c+34); a.set_title('reflex corner, c/n = 0.34 (zoom)',fontsize=10); a.set_xlabel('X')
fig.colorbar(sc,ax=ax,shrink=.8,label='liquidity 1 - P_max (0 = frozen)')
plt.savefig('../figs/plus_verify_n241.png',dpi=120,bbox_inches='tight')
