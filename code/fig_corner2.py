"""Corrected corner switch: exact liquidity maps of the corrected Pfaffian model with the scaling solution (corner2_theory) overlaid,
and the liquid fraction versus (J/2T)/m for several sizes."""
import sys, math, warnings, pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
from scipy.optimize import brentq
from gen_ilp import Sample, halfplanes
import corrected_model as cm
from corrected_map import stats
from corner2_theory import Corner2, field2, liquid_fraction
L=int(sys.argv[1]); phis=[float(v) for v in sys.argv[2:]]
S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); n=cm.halfsites(S,'fit'); m=len(S.free)
bs=np.exp(np.linspace(np.log(0.02),np.log(40),40)); tab=np.array([(field2(Corner2(b),n=250)[0],liquid_fraction(Corner2(b))) for b in bs])
keep=[0]
for i in range(1,len(bs)):
    if 0<tab[i,0]<tab[keep[-1],0]: keep.append(i)      # drop points where the flux quadrature fails (large b)
bs,tab=bs[keep],tab[keep]
pickle.dump((bs,tab),open('data/corner2_table.pkl','wb'))
fig,ax=plt.subplots(1,len(phis)+1,figsize=(3.3*(len(phis)+1),3.6))
for a,phi in zip(ax,phis):
    st=stats(S,n,math.exp(-2*phi*m)); b=float(np.exp(np.interp(phi,tab[::-1,0],np.log(bs[::-1])))); c=Corner2(b)
    a.scatter(st['cen'][:,0],st['cen'][:,1],c=1-st['mx'],s=3.5*(64/L)**2,marker='s',cmap='viridis',vmin=0,vmax=0.5,lw=0)
    sc=L+1
    for u,v in c.caps():
        for X,Y in ((u,v),(v,u)): a.plot(sc*X,sc*Y,color='#ff2d55',lw=1.8)
    a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_title(f'(J/2T)/m = {phi:.2f}\nliquid: exact {st["liq"]:.3f}, theory {liquid_fraction(c):.3f}',fontsize=8.5)
    print(f'phi={phi}: b={b:.3f}, liquid exact {st["liq"]:.4f} (row-sum error {st["err"]:.0e}), theory {liquid_fraction(c):.4f}',flush=True)
a=ax[-1]; a.plot(tab[:,0],tab[:,1],'k-',lw=1.5,label='scaling solution')
dat={32:[(0.1,.996),(0.173,.991),(0.26,.956),(0.347,.926),(0.433,.862),(0.6,.699),(0.7,.562)],
     48:[(0.1,.997),(0.173,.980),(0.26,.951),(0.347,.902),(0.433,.841),(0.52,.759),(0.6,.668),(0.7,.519),(0.8,.315)],
     64:[(0.1,.995),(0.173,.981),(0.26,.945),(0.347,.897),(0.433,.834),(0.52,.749),(0.6,.650),(0.7,.500)]}
for Lk,mk in ((32,'s'),(48,'^'),(64,'o')):
    d=np.array(dat[Lk]); a.plot(d[:,0],d[:,1],mk,ms=4,mfc='none',label=f'exact, L = {Lk}')
a.axvline(math.log(1+2**.5),color='0.6',ls=':'); a.set_xlabel('(J/2T)/m'); a.set_ylabel('liquid fraction'); a.set_xlim(0,0.92); a.set_ylim(0,1.02); a.legend(fontsize=7)
plt.tight_layout(); plt.savefig('../figs/corner2_theory.png',dpi=170)
