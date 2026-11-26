"""Fig. 6. The boundary model against the full vertex model (exact density of states, corner switch):
(a) excitation energy versus (J/2T)/m for L = 8-20; (b) difference to the corrected model versus 1/L."""
import pickle, numpy as np
from figstyle import *
setup()
d=pickle.load(open('data/corrected_thermo.pkl','rb')); Ls=sorted({k[0] for k in d}); phis=sorted({k[1] for k in d})
fig,(a,b)=plt.subplots(1,2,figsize=(W2*0.62,2.1),gridspec_kw=dict(wspace=0.38))
cols={8:OI['sky'],12:OI['green'],16:OI['orange'],20:OI['blue']}
for L in Ls:
    ev=[d[(L,p)][0] for p in phis]; ec=[d[(L,p)][1] for p in phis]; e0=[d[(L,p)][2] for p in phis]
    a.plot(phis,ev,'o',ms=3,mfc='none',mec=cols[L],mew=0.8); a.plot(phis,ec,'-',color=cols[L],lw=0.9,label=f'L = {L}')
    a.plot(phis,e0,':',color=cols[L],lw=0.8)
a.set_xlabel(r'$q=J/(2mT)$'); a.set_ylabel(r'$\langle E-E_{\rm gs}\rangle/(mJ/2)$'); a.set_yscale('log'); a.set_ylim(5e-4,3)
a.plot([],[],'ko',ms=3,mfc='none',label='vertex model'); a.plot([],[],'k-',label='boundary model'); a.plot([],[],'k:',label='freeable sites only')
a.legend(fontsize=6,ncol=1,loc='lower left',handlelength=1.4)
label(a,'a',-0.12,1.0)
for p,mk in zip(phis,('s','^','o','D','v')):
    if p>0.9: continue
    diff=[d[(L,p)][0]-d[(L,p)][1] for L in Ls]; b.plot([1/L for L in Ls],diff,mk+'-',ms=3,lw=0.7,label=f'{p:.2f}')
b.set_xlabel('1/L'); b.set_ylabel('vertex − boundary model'); b.set_xlim(0,0.135); b.set_ylim(0,None); b.legend(title=r'$q$',fontsize=6,title_fontsize=6)
label(b,'b',-0.18,1.0)
save(fig,'fig6')
