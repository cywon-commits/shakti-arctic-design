"""Fig. 1. Shakti islands, plaquettes and dominoes; the slope square; a cut sample whose termination sets the boundary (corner switch, L = 12)."""
import io, contextlib, numpy as np, matplotlib.patches as mp
from figstyle import *
from gen_ilp import Sample, halfplanes, black, R, stem, ppts
from collections import Counter
import corrected_model as cm
from bmono import centre
setup()
def draw_sample(ax,S,dimers=None,free_mark=True,half=None,isl_lw=1.0,shade=True):
    for c in S.cells:
        (x0,y0),(x1,y1)=R(c)[0],R(c)[-1]; X0,Y0=min(x0,x1),min(y0,y1); w=abs(x1-x0)+1; h=abs(y1-y0)+1
        fc=(H_FILL if black(c) else V_FILL) if shade else 'none'
        if free_mark and c in S.free: fc=OI['yellow']
        ax.add_patch(mp.Rectangle((X0,Y0),w,h,fc=fc,ec='none',zorder=0))
    for (a,b) in S.isl:
        d=np.subtract(b,a); p=np.add(a,0.17*d); q=np.subtract(b,0.17*d)
        ax.plot([p[0],q[0]],[p[1],q[1]],color='0.25',lw=isl_lw,solid_capstyle='round',zorder=2)
    if dimers:
        deg=Counter(p for e in S.isl for p in e)
        for a,b in dimers:
            pts=[p for p in ppts(a)&ppts(b) if deg[p]==3]
            for p in pts: ax.plot(*p,'o',ms=3.0,color='k',zorder=4)
    if half:
        for (x,y) in half: ax.plot(x,y,marker='D',ms=3.2,mfc='white',mec=ACCENT,mew=0.9,zorder=4)
    ax.set_aspect('equal'); ax.axis('off')
fig=plt.figure(figsize=(W2,2.45))
# (a) bulk patch with a domino tiling
axa=fig.add_axes([0.0,0.02,0.30,0.88])
S=Sample(halfplanes((1,0,1,9),(0,1,1,9)),W=14,mode='A'); r=S.lp(free=set())
draw_sample(axa,S,dimers=r[2],free_mark=False)
c0=sorted(S.cells,key=lambda c:abs(centre(c)[0]-5)+abs(centre(c)[1]-5))[0]
sa,sb=stem(c0); axa.plot([sa[0],sb[0]],[sa[1],sb[1]],color=ACCENT,lw=2.0,zorder=5)
axa.text(0.5,-0.02,'islands, H/V plaquettes; dots: unhappy vertices = dimers',transform=axa.transAxes,ha='center',va='top',fontsize=7)
label(axa,'a',0.04,0.98)
# (b) slope square
axb=fig.add_axes([0.355,0.14,0.25,0.72])
axb.add_patch(mp.Polygon([(-1,-1),(1,-1),(1,1),(-1,1)],closed=True,fc='0.95',ec='k',lw=0.8))
for (x,y) in [(1,1),(-1,1),(-1,-1),(1,-1)]:
    axb.plot(x,y,'o',ms=5,color='k',zorder=3); axb.text(x*1.0,y*1.17,f'({x:+d},{y:+d})'.replace('+1','1').replace('-1','−1'),ha='center',va='center',fontsize=7)
for (x,y),rot in zip([(0,-1),(1,0),(0,1),(-1,0)],[0,90,0,90]):
    axb.plot(x,y,'s',ms=4,color=OI['orange'],zorder=3)
axb.text(0,-1.45,'edges: through-block axis cut fixes one',ha='center',fontsize=6.5,color=OI['orange'])
axb.text(0,1.38,'corners: ice phases; even 45° cut fixes one',ha='center',fontsize=6.5,color='k')
axb.set_xlim(-1.6,1.6); axb.set_ylim(-1.65,1.6); axb.set_aspect('equal'); axb.axis('off')
axb.text(0,-0.08,r'$(\partial_X h,\ \partial_Y h)$',ha='center',fontsize=7.5)
label(axb,'b',0.06,1.0)
# (c) corner switch L = 12: ground state, free plaquettes, half-dimer sites
axc=fig.add_axes([0.64,0.02,0.36,0.88])
L=12; S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); r=S.lp(); n=cm.halfsites(S,'fit')
def halfpoints(S):
    # boundary z = 3 vertices whose unhappy state covers a single plaquette (half-dimer sites)
    deg=Counter(p for e in S.isl for p in e); touch={}
    for c in S.cells:
        for p in ppts(c): touch.setdefault(p,[]).append(c)
    out=[]
    for p,d in deg.items():
        if d!=3: continue
        t=touch.get(p,[])
        if len(t)==1 and t[0] not in S.free: out.append(p)
        elif len(t)==2 and abs(t[0][0]-t[1][0])+abs(t[0][1]-t[1][1])!=1 and any(c not in S.free for c in t): out.append(p)
    return out
draw_sample(axc,S,dimers=r[2],half=halfpoints(S))
axc.text(0.5,1.01,'W, S: block-aligned  |  N, E: through-block, fringed',transform=axc.transAxes,ha='center',va='bottom',fontsize=7)
axc.plot([],[],'s',color=OI['yellow'],ms=5,label='free plaquette (gains J/2)'); axc.plot([],[],'D',mfc='white',mec=ACCENT,ms=4,label='half-dimer site (costs J/2)')
axc.legend(loc='upper center',bbox_to_anchor=(0.5,-0.07),ncol=2,handletextpad=0.3,columnspacing=0.8,borderaxespad=0)
label(axc,'c',0.0,1.06)
save(fig,'fig1')
