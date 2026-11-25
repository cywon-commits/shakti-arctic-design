"""Figure: liquidity maps of the corner-switch sample at several temperatures (from data/switch_corner.pkl)."""
import sys, math, pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
out=pickle.load(open('data/switch_corner.pkl','rb')); L=int(sys.argv[1]); ks=[int(v) for v in sys.argv[2:]]
fig,ax=plt.subplots(1,len(ks),figsize=(3.1*len(ks),3.5))
for a,k in zip(ax,ks):
    st=out[(L,k)]; m=round(out[(L,40)]['nmono'])
    sc=a.scatter(st['cen'][:,0],st['cen'][:,1],c=1-st['mx'],s=14*(32/L)**2,marker='s',cmap='viridis',vmin=0,vmax=0.5,lw=0)
    a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_title(f'(J/2T)/m = {k*math.log(2)/m:.2f}\nfree fraction {st["nmono"]/m:.2f}',fontsize=9)
    a.plot([0,L+1,L+1,0,0],[0,0,L+1,L+1,0],color='0.5',lw=.6)
    a.plot([0,0,L+1],[L+1,0,0],color='#c0392b',lw=2)       # block-aligned (switchable) edges
fig.colorbar(sc,ax=ax,fraction=0.02,label='1 - max pairing probability'); plt.savefig('../figs/switch_corner.png',dpi=180,bbox_inches='tight')
