import pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from theory import curve
D=pickle.load(open('data/verify_n121.pkl','rb')); n=121
fig,ax=plt.subplots(1,2,figsize=(15,7.6))
for a,(c,d) in zip(ax,sorted(D.items())):
    q=list(d['cen']); x=[d['cen'][k][0] for k in q]; y=[d['cen'][k][1] for k in q]; v=[1-d['mx'][k] for k in q]
    sc=a.scatter(x,y,c=v,s=2.2,marker='s',cmap='viridis',vmin=0,vmax=.6,linewidths=0)
    th,nu_,cx,cy=curve(*d['sol'],m=20001); a.plot(cx,cy,color='crimson',lw=1.1)
    a.set_aspect('equal'); a.set_xlim(-n-4,n+4); a.set_ylim(-n-4,n+4)
    a.set_title(f'plus-shaped Shakti sample, reach n={n}, half-width c={c} (c/n={c/n:.2f})\ncolour: exact 1-P_max (Kasteleyn);  line: class-10 curve',fontsize=10)
    a.set_xlabel('X (Shakti units)'); a.set_ylabel('Y')
fig.colorbar(sc,ax=ax,shrink=.8,label='liquidity 1-P_max (0 = frozen)')
plt.savefig('../figs/plus_verify_n121.png',dpi=110,bbox_inches='tight')
