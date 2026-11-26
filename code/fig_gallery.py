"""Gallery figure: designed shapes (unions of overlapping through-block squares) and their exact liquid regions."""
import pickle, numpy as np
from figstyle import *
setup()
d=pickle.load(open('data/design_gallery.pkl','rb')); names=['check','chevron','zigzag','X']; titles={'check':'check mark','diag3':'diagonal wire','chevron':'chevron (arrow)','zigzag':'zigzag','X':'cross X'}
fig,ax=plt.subplots(1,len(names),figsize=(W2,2.0),gridspec_kw=dict(wspace=0.05))
for k,(a,nm) in enumerate(zip(ax,names)):
    g=d[nm]; P=g['cen']; s=g['s']; xs=[c[0] for c in g['centres']]; ys=[c[1] for c in g['centres']]
    lim=((min(xs)-s-1,max(xs)+s+1),(min(ys)-s-1,max(ys)+s+1)); im=liq_image(a,P,g['liq'],lim,res=500,rmax=1.1)
    for cx,cy in g['centres']: a.plot([cx-s,cx+s,cx+s,cx-s,cx-s],[cy-s,cy-s,cy+s,cy+s,cy-s],color='w',lw=0.35,alpha=0.6)
    a.set_aspect('equal'); a.axis('off'); a.set_title(titles[nm],fontsize=7.5); label(a,'abcd'[k],0.08,0.98)
cb=fig.colorbar(im,ax=list(ax),fraction=0.012,pad=0.01,shrink=0.8); cb.set_label(r'liquidity $1-P_{\max}$')
save(fig,'gallery_strokes')
