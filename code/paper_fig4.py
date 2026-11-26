"""Fig. 4. Temperature switch: (a) which rectangles switch (T = 0 classification by edge type); (b-d) exact liquidity of the corner switch
(corrected boundary model, L = 48) with the scaling solution; (e) liquid fraction versus (J/2T)/m."""
import math, pickle, warnings, numpy as np, matplotlib.patches as mp
warnings.filterwarnings('ignore')
from figstyle import *
from gen_ilp import Sample, halfplanes
import corrected_model as cm
from corrected_map import stats
from corner2_theory import Corner2, liquid_fraction
setup()
fig=plt.figure(figsize=(W2,2.35))
# (a) classification schematic
axa=fig.add_axes([0.0,0.04,0.17,0.86]); axa.axis('off'); axa.set_xlim(0,1); axa.set_ylim(0,1.7); axa.set_aspect('equal',adjustable='datalim',anchor='W')
rows=[('EEEE','Aztec circle'),('EEEO','fixed (no switch)'),('E E\nO O','unique g.s. → liquid'),('EOEO','no tiling'),('OOOO','Aztec circle')]
def sq(x,y,s,edges):   # edges: W,S,E,N : True = block-aligned (thick), False = through-block (thin dashed)
    P=[(x,y+s),(x,y),(x+s,y),(x+s,y+s),(x,y+s)]
    for k,e in enumerate(edges):
        (x0,y0),(x1,y1)=P[k],P[k+1]
        axa.plot([x0,x1],[y0,y1],color='k' if e else '0.5',lw=2.0 if e else 0.8,ls='-' if e else (0,(2,1.5)),solid_capstyle='butt')
cases=[((1,1,1,1),'all block-aligned:\nAztec circle'),((1,1,0,0),'two adjacent:\nswitch'),((1,0,1,0),'two opposite:\nno tiling'),((0,0,0,0),'none:\nAztec circle')]
for k,(e,t) in enumerate(cases):
    y=1.30-k*0.42; sq(0.03,y,0.26,e); axa.text(0.36,y+0.13,t,fontsize=6.5,va='center',color=ACCENT if k==1 else 'k')
axa.text(0.0,1.66,'T = 0; thick edge = block-aligned',fontsize=6.5,va='bottom')
label(axa,'a',0.02,1.07)
# (b-d) maps
L=48; S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); n=cm.halfsites(S,'fit'); m=len(S.free)
bs,tab=pickle.load(open('data/corner2_table.pkl','rb'))
phis=(0.26,0.6,0.8); axs=[fig.add_axes([0.19+k*0.165,0.16,0.16,0.72]) for k in range(3)]
for k,(a,phi) in enumerate(zip(axs,phis)):
    st=stats(S,n,math.exp(-2*phi*m)); b=float(np.exp(np.interp(phi,tab[::-1,0],np.log(bs[::-1])))); c=Corner2(b)
    im=liq_image(a,st['cen'],1-st['mx'],((0,L+1),(0,L+1)),res=300,rmax=1.3)
    for u,v in c.caps():
        for X,Y in ((u,v),(v,u)): a.plot((L+1)*X,(L+1)*Y,color=ACCENT,lw=1.0)
    a.set_aspect('equal'); a.axis('off'); a.set_title(f'(J/2T)/m = {phi:.2f}',fontsize=7.5); label(a,'bcd'[k],0.1,0.98)
    a.plot([0,0,L+1],[L+1,0,0],color='k',lw=1.6,solid_capstyle='butt'); a.plot([0,L+1,L+1],[L+1,L+1,0],color='0.5',lw=0.7,ls=(0,(2,1.5)))
cax=fig.add_axes([0.25,0.06,0.37,0.025]); cb=fig.colorbar(im,cax=cax,orientation='horizontal'); cb.set_label(r'liquidity $1-P_{\max}$',labelpad=1); cb.ax.tick_params(labelsize=6)
# (e) liquid fraction
axe=fig.add_axes([0.765,0.20,0.225,0.66])
axe.plot(tab[:,0],tab[:,1],'k-',lw=1.1,label='scaling solution')
dat={}
for Lk in (32,48,64):
    dd=pickle.load(open(f'data/corrected_map_{Lk}.pkl','rb')); dat[Lk]=[(p,v['liq']) for p,v in sorted(dd.items()) if v['err']<1e-3]
for Lk,mk,col in ((32,'s',OI['sky']),(48,'^',OI['green']),(64,'o',OI['blue'])):
    d=np.array(dat[Lk]); axe.plot(d[:,0],d[:,1],mk,ms=3.2,mfc='none',mec=col,mew=0.8,label=f'L = {Lk}')
axe.axvline(math.log(1+2**.5),color='0.5',ls=':',lw=0.8); axe.text(math.log(1+2**.5)+0.015,0.55,r'$\ln(1+\sqrt{2})$',rotation=90,ha='left',fontsize=6.5)
axe.set_xlabel(r'$J/(2T)\,/\,m$'); axe.set_ylabel('liquid fraction'); axe.set_xlim(0,0.97); axe.set_ylim(0,1.03); axe.legend(loc='lower left',fontsize=6.5,handlelength=1.2)
label(axe,'e',-0.12,1.0)
save(fig,'fig4')
