"""Fig. 2. Convex corner cuts move and resize the arctic circle but keep it a circle: exact liquidity maps of octagonal samples
(all-paired, through-block axis edges, even 45-degree cuts) with the circle predicted by the box rule."""
import numpy as np, matplotlib.patches as mp
from figstyle import *
from mixed_octagon import cells, predict
from maps import cell_entropy, pcenter
setup()
a=b=21; cases=[((42,42,42,42),'no cuts'),((32,32,32,32),'four cuts'),((24,24,42,42),'two cuts'),((32,42,38,36),'three cuts, mixed')]
cases=[c for c in cases if predict(a,b,c[0])]
fig,ax=plt.subplots(1,len(cases),figsize=(W2,W2/len(cases)+0.35))
for k,(A,(d,title)) in enumerate(zip(ax,cases)):
    S=cells(a,b,d); liq,ld=cell_entropy(S); P=np.array([pcenter(c) for c in S]); v=np.array([liq[c] for c in S])
    sc=liq_image(A,P,v,((-a-1,a+1),(-b-1,b+1)),res=400,rmax=1.1)
    A.plot([-a,a,a,-a,-a],[-b,-b,b,b,-b],color='0.6',lw=0.5)
    X0,X1,Y0,Y1=predict(a,b,d); cx,cy,r=(X0+X1)/2,(Y0+Y1)/2,(X1-X0)/2
    A.add_patch(mp.Rectangle((X0,Y0),X1-X0,Y1-Y0,fill=False,ls=(0,(2,2)),ec='0.3',lw=0.6))
    th=np.linspace(0,2*np.pi,400); A.plot(cx+r*np.cos(th),cy+r*np.sin(th),color=ACCENT,lw=1.0)
    for di,(sx,sy) in zip(d,((1,1),(-1,1),(-1,-1),(1,-1))):
        if di<a+b:
            # cut line sx X + sy Y = di inside the square
            xs=np.linspace(-a,a,400); ys=(di-sx*xs)/sy; m=(abs(ys)<=b); A.plot(xs[m],ys[m],color='0.2',lw=0.8)
    A.set_aspect('equal'); A.set_xlim(-a-1,a+1); A.set_ylim(-b-1,b+1); A.set_xticks([]); A.set_yticks([]); A.set_title(title,fontsize=7.5)
    for s in A.spines.values(): s.set_visible(False)
    label(A,'abcd'[k],0.06,1.0)
cb=fig.colorbar(sc,ax=ax,fraction=0.012,pad=0.01,shrink=0.75); cb.set_label(r'liquidity $1-P_{\max}$')
save(fig,'fig2')
