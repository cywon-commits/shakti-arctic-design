"""Excitation energy <E - E_gs>/(m J/2) versus phi = (J/2T)/m: vertex model (exact DOS, data/vertex_thermo*.pkl),
corrected boundary model (half-dimer sites, exact polynomial) and the original model (freeable plaquettes only)."""
import pickle, math, numpy as np
from gen_ilp import Sample, halfplanes
import corrected_model as cm
res={}
for f in ('data/vertex_thermo_8_16.pkl','data/vertex_thermo.pkl'): res.update(pickle.load(open(f,'rb')))
def exc(c,Egs,phi,m):
    bJ=2*phi*m; E=np.arange(len(c)); lw=np.array([math.log(v) if v>0 else -np.inf for v in c],float)-bJ*E
    lw-=lw.max(); w=np.exp(lw); return ((w*(E-Egs)).sum()/w.sum())/(m/2)
print(' L   phi   vertex  corrected  original')
out={}
for L in sorted(res):
    g,c0=res[L]; S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
    cc=cm.poly(S,cm.halfsites(S,'fit'))
    for phi in (0.3,0.433,0.6,0.8,1.0):
        bJ=2*phi*m; K=len(g); wv=g*np.exp(-bJ*np.arange(K)); ev=(wv*np.arange(K)).sum()/wv.sum()/(m/2)
        # original model: c0 indexed by number of free plaquettes j; energy E - E_gs = (m - j)/2
        e0=sum((m-j)/2*float(v)*math.exp(-bJ*(m-j)/2) for j,v in enumerate(c0))/sum(float(v)*math.exp(-bJ*(m-j)/2) for j,v in enumerate(c0))/(m/2)
        ec=exc(cc,Egs,phi,m); out[(L,phi)]=(ev,ec,e0)
        print(f'{L:3d} {phi:5.3f}  {ev:7.4f}  {ec:7.4f}   {e0:7.4f}')
pickle.dump(out,open('data/corrected_thermo.pkl','wb'))
