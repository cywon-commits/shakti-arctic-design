"""Fig. 3. Reflex corners make the liquid region non-circular: plus-shaped samples (arm reach n = 241, half-width c), exact pairing
probabilities (urban renewal); frozen cells coloured by ice phase, liquid white; class-10 curve in black.  (d) zoom on a reflex corner:
a pocket where three ice phases meet."""
import pickle, numpy as np, matplotlib.patches as mp
from figstyle import *
from theory import curve
setup()
D=pickle.load(open('data/phase_maps.pkl','rb')); keys=sorted(D,key=lambda k:k[1])
fig,ax=plt.subplots(1,4,figsize=(W2,2.1),gridspec_kw=dict(wspace=0.06))
for k,(a,key) in enumerate(zip(ax,keys)):
    d=D[key]; n,c=key; rgb=phase_rgb(d['prob'])
    rgb_image(a,d['cen'],rgb,((-n-4,n+4),(-n-4,n+4)),res=700)
    th,v,x,y=curve(*d['sol'],m=8001); a.plot(x,y,color='k',lw=0.6)
    a.set_aspect('equal'); a.set_xlim(-n-4,n+4); a.set_ylim(-n-4,n+4); a.axis('off'); a.set_title(f'c/n = {c/n:.2f}',fontsize=7.5)
    label(a,'abc'[k],0.08,0.97)
d=D[keys[0]]; n,c=keys[0]; w=0.30*n
ax[0].add_patch(mp.Rectangle((c-w,c-w),2*w,2*w,fill=False,ec='k',lw=0.6))
a=ax[3]; m=(abs(d['cen'][:,0]-c)<w+2)&(abs(d['cen'][:,1]-c)<w+2)
rgb_image(a,d['cen'][m],phase_rgb(d['prob'][m]),((c-w,c+w),(c-w,c+w)),res=500)
th,v,x,y=curve(*d['sol'],m=20001); a.plot(x,y,color='k',lw=0.8)
a.set_aspect('equal'); a.set_xlim(c-w,c+w); a.set_ylim(c-w,c+w); a.set_xticks([]); a.set_yticks([]); a.set_title('reflex corner (box in a)',fontsize=7.5)
label(a,'d',0.08,0.97)
fig.text(0.33,0.05,'frozen, by ice phase:',fontsize=6.5,va='center',ha='right')
for j,col in enumerate(PHASE_COLS):
    fig.patches.append(mp.Rectangle((0.345+j*0.035,0.035),0.022,0.03,transform=fig.transFigure,fc=col,ec='none'))
fig.text(0.50,0.05,'white: liquid;  black line: class-10 curve',fontsize=6.5,va='center')
save(fig,'fig3')
