"""Figure: temperature-programmed X (design_switchX.py)."""
import pickle, numpy as np
from figstyle import *
setup()
d=pickle.load(open('data/design_switchX.pkl','rb')); R=d['R']; Ts=[40,16,10,6]
fig,ax=plt.subplots(1,len(Ts),figsize=(W2,2.05),gridspec_kw=dict(wspace=0.04))
lim=((-26,26),(-26,26))
for k,(a,JT) in enumerate(zip(ax,Ts)):
    g=d['maps'][float(JT)]; im=liq_image(a,g['cen'],g['liq'],lim,res=500,rmax=1.1)
    for j,(x0,x1,y0,y1) in enumerate(R):
        col='w' if j<3 else ACCENT
        if j>=3:   # block-aligned outer edges of the switch squares
            if x1-x0>2*11: 
                if x1>0: a.plot([x1,x1],[y0,y1],color=col,lw=1.3); a.plot([x0,x1],[y0,y0],color=col,lw=1.3)
                else: a.plot([x0,x0],[y0,y1],color=col,lw=1.3); a.plot([x0,x1],[y1,y1],color=col,lw=1.3)
    a.set_aspect('equal'); a.axis('off'); a.set_title(f'J/T = {JT}',fontsize=7.5); label(a,'abcd'[k],0.08,0.98)
cb=fig.colorbar(im,ax=list(ax),fraction=0.012,pad=0.01,shrink=0.8); cb.set_label(r'liquidity $1-P_{\max}$')
save(fig,'gallery_switchX')
