"""Temperature-programmed symbol: an X of five overlapping squares.  The '/' stroke (centre, NE, SW squares) has through-block edges;
the '\\' stroke squares (NW, SE) have their two outer edges block-aligned (corner-switch type).  Low T: only '/' is liquid; high T: X."""
import math, sys, pickle, numpy as np, warnings
warnings.filterwarnings('ignore')
from gen_ilp import Sample, black
import corrected_model as cm
from corrected_map import stats
o,s=14,11
R=[(-s,s,-s,s),(o-s,o+s,o-s,o+s),(-o-s,-o+s,-o-s,-o+s),(o-s,o+s+1,-o-s-1,-o+s),(-o-s-1,-o+s,o-s,o+s+1)]
ins=lambda p: any(x0<=p[0]<=x1 and y0<=p[1]<=y1 for x0,x1,y0,y1 in R)
S=Sample(ins,W=o+s+5,mode='A'); nb=sum(black(c) for c in S.cells)
print('cells',len(S.cells),'freeable',len(S.free),'H-V',2*nb-len(S.cells))
r=S.lp(); print('T=0 effective model:',None if r is None else r[:2])
n=cm.halfsites(S,'fit'); print('half-dimer weight',sum(n.values()))
out={}
for JT in [float(v) for v in sys.argv[1:]] or [40,25,16,10,6,3]:
    st=stats(S,n,math.exp(-JT)); out[JT]=dict(cen=st['cen'],liq=1-st['mx'],err=st['err'],free=st['free'])
    print(f'J/T={JT}: row-sum error {st["err"]:.1e}, free fraction {st["free"]:.3f}, liquid fraction {st["liq"]:.3f}',flush=True)
pickle.dump(dict(R=R,maps=out),open('data/design_switchX.pkl','wb'))
