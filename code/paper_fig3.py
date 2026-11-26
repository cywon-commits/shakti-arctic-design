"""Fig. 3. Reflex corners make the liquid region non-circular: plus-shaped samples (arm reach n = 241, half-width c), exact liquidity
(urban renewal) with the class-10 curve, and a zoom on a reflex corner."""
import pickle, numpy as np
from figstyle import *
from theory import curve, closed_form, solve
setup()
D={}
for f in ('data/verify_shuffle_shapes.pkl','data/verify_shuffle.pkl'): D.update(pickle.load(open(f,'rb')))
keys=[(241,83),(241,121),(241,165)]
fig,ax=plt.subplots(1,4,figsize=(W2,2.05),gridspec_kw=dict(wspace=0.08))
for k,(a,key) in enumerate(zip(ax,keys)):
    d=D[key]; n,c=key
    sc=liq_image(a,d['cen'],d['liq'],((-n-4,n+4),(-n-4,n+4)),res=700)
    th,v,x,y=curve(*d['sol'],m=8001); a.plot(x,y,color=ACCENT,lw=0.7)
    a.set_aspect('equal'); a.set_xlim(-n-4,n+4); a.set_ylim(-n-4,n+4); a.axis('off'); a.set_title(f'c/n = {c/n:.2f}',fontsize=7.5)
    label(a,'abc'[k],0.08,0.97)
a=ax[3]; cols=[OI['blue'],OI['sky'],OI['green'],OI['orange'],OI['red'],OI['purple']]
for j,rho in enumerate((0.36,0.42,0.50,0.58,0.66,0.70)):
    sol=solve(1.0,rho); sol=min(sol,key=lambda s:abs(s[0])) if isinstance(sol,list) else sol
    th,v,x,y=curve(*sol,m=6001); m_=(x>=-1e-9)&(y>=-1e-9); col=cols[j]; x=np.where(m_,x,np.nan); y=np.where(m_,y,np.nan)
    a.plot(x,y,color=col,lw=0.9,label=f'{rho:.2f}')
    a.plot([1,1,rho,rho,0],[0,rho,rho,1,1],color=col,lw=0.4,ls=':')
a.set_aspect('equal'); a.set_xlim(0,1.02); a.set_ylim(0,1.02); a.set_xlabel('X/n'); a.set_ylabel('Y/n',labelpad=1)
a.set_title(r'curve family, $1/3<c/n<1/\sqrt{2}$',fontsize=7.5); a.legend(title='c/n',fontsize=6,title_fontsize=6.5,loc='lower left',handlelength=1.2,labelspacing=0.2)
label(a,'d',-0.25,1.0)
cb=fig.colorbar(sc,ax=list(ax[:3]),fraction=0.012,pad=0.01,shrink=0.8); cb.set_label(r'liquidity $1-P_{\max}$')
save(fig,'fig3')
