"""Fig. 5. Edge energies reduce to corner sources: (a) source pattern of the through-block Aztec square and of the corner switch;
(b,c) through-block Aztec square (corrected model, N = 23) above and below T_on; (d) sigma_N(x) versus the scaling function of the
boundary-monomer Aztec diamond, with J/T = (xN + 2) ln 2."""
import sys, math, pickle, warnings, numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/claude/screened-dimer-electrostatics/src')
from figstyle import *
from gen_ilp import Sample, halfplanes
import corrected_model as cm
from corrected_map import stats
setup()
fig=plt.figure(figsize=(W2,2.3))
# (a) schematic
axa=fig.add_axes([0.0,0.06,0.2,0.82]); axa.axis('off'); axa.set_xlim(-0.05,1.05); axa.set_ylim(-0.1,2.3); axa.set_aspect('equal')
def src(ax,x,y,s):
    ax.plot(x,y,'o',ms=8,mfc=ACCENT if s>0 else OI['blue'],mec='none',zorder=3); ax.text(x,y,'+' if s>0 else '−',color='w',ha='center',va='center',fontsize=8,fontweight='bold',zorder=4)
for y0,signs,t in ((1.25,(1,-1,1,-1),'through-block square'),(0.0,(1,0,-1,0),'corner switch')):
    axa.plot([0,0.9,0.9,0,0],[y0,y0,y0+0.9,y0+0.9,y0],color='0.4',lw=0.8)
    for (x,y),s in zip(((0,y0),(0.9,y0),(0.9,y0+0.9),(0,y0+0.9)),signs):
        if s: src(axa,x,y,s)
    axa.text(0.45,y0+0.45,t,ha='center',va='center',fontsize=6.5)
axa.text(0.45,2.33,r'boundary energy $\propto \sum\pm h_{\rm corner}$',ha='center',fontsize=6.5)
label(axa,'a',0.02,1.10)
# (b,c) maps
L=48; N=L//2-1; S=Sample(halfplanes((1,0,1,L+1),(0,1,1,L+1)),W=L+10,mode='A'); n=cm.halfsites(S,'fit')
for k,x in enumerate((1.2,0.4)):
    a=fig.add_axes([0.22+k*0.19,0.12,0.18,0.74]); t=math.exp(-(x*N+2)*math.log(2)); st=stats(S,n,t)
    im=liq_image(a,st['cen'],1-st['mx'],((1,L+1),(1,L+1)),res=300,rmax=1.3)
    a.set_aspect('equal'); a.axis('off'); a.set_title(f'x = {x}' + (r'  ($T<T_{\rm on}$)' if x>1 else r'  ($T>T_{\rm on}$)'),fontsize=7.5); label(a,'bc'[k],0.1,0.98)
cax=fig.add_axes([0.26,0.07,0.30,0.025]); cb=fig.colorbar(im,cax=cax,orientation='horizontal'); cb.set_label(r'liquidity $1-P_{\max}$',labelpad=1); cb.ax.tick_params(labelsize=6)
# (d) sigma
axd=fig.add_axes([0.70,0.20,0.29,0.66])
from exact_sol import solve
Js=np.exp(np.linspace(np.log(1.0003),np.log(3e6),70)); sx=np.array(sorted((solve(J)['x'],solve(J)['sigma']) for J in Js))
axd.plot(np.r_[sx[:,0],1.0,1.3],np.r_[sx[:,1],1.0,1.0],'k-',lw=1.1,label=r'$\sigma_\infty(x)$, vacancy model')
xs=np.linspace(0.15,1.3,24)
for Lk,mk,col in ((32,'s',OI['sky']),(64,'^',OI['green']),(96,'o',OI['blue'])):
    Nk=Lk//2-1; Sk=Sample(halfplanes((1,0,1,Lk+1),(0,1,1,Lk+1)),W=Lk+10,mode='A'); nk=cm.halfsites(Sk,'fit'); pts=[]
    for x in xs:
        st=stats(Sk,nk,math.exp(-(x*Nk+2)*math.log(2)))
        if st['err']<1e-4: pts.append((x,1-st['half']/(2*Nk)))
    pts=np.array(pts); axd.plot(pts[:,0],pts[:,1],mk,ms=3,mfc='none',mec=col,mew=0.8,label=f'Shakti, N = {Nk}')
axd.axvline(1,color='0.5',ls=':',lw=0.8); axd.text(1.0,1.045,r'$T_{\rm on}$',fontsize=7,ha='center')
axd.set_xlabel(r'$x=(J/T-2\ln 2)/(N\ln 2)$'); axd.set_ylabel(r'$\sigma = 1-\langle n_{\rm mono}\rangle/2N$'); axd.set_xlim(0.1,1.3); axd.set_ylim(0.2,1.03)
axd.legend(loc='lower right',fontsize=6,handlelength=1.2)
label(axd,'d',-0.14,1.0)
save(fig,'fig5')
