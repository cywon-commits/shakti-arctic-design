"""Rectangles with every combination of edge offsets (mod 4), fringed termination: state at T = 0 (y -> infinity) and at y = 1.
Even offset = block-aligned edge (carries freeable plaquettes), odd offset = through-block edge with fringe (none).
Usage: python mixed_scan.py L   (side lengths L + 0..3)"""
import sys, math, pickle, itertools, numpy as np
from gen_ilp import Sample, halfplanes, black
from bmono import BMono
def state(B,S,y):
    r=B.exact(int(y))
    if r is None: return None
    lnZ,mono,pair=r; mx=np.zeros(B.N); tot=np.zeros(B.N); mx[:B.m]=mono; tot[:B.m]=mono
    for (i,j),p in pair.items():
        mx[i]=max(mx[i],p); mx[j]=max(mx[j],p); tot[i]+=p; tot[j]+=p
    liq=mx<0.9; cen=np.array([B.cen[c] for c in B.order])
    return dict(lnZ=lnZ,nmono=mono.sum(),err=abs(tot-1).max(),liqfrac=liq.mean(),liq=liq,cen=cen,mx=mx,mono=mono)
def scan(L,YBIG=2**40):
    out={}
    for x0,y0,dx,dy in itertools.product(range(4),repeat=4):
        x1,y1=x0+L+dx,y0+L+dy; S=Sample(halfplanes((1,0,x0,x1),(0,1,y0,y1)),W=L+10,mode='A')
        key=(x0,x1%4,y0,y1%4); r=dict(box=(x0,x1,y0,y1),cells=len(S.cells),nfree=len(S.free),H=sum(black(c) for c in S.cells),fH=sum(black(c) for c in S.free))
        gs=S.lp(); r['gs_free']=None if gs is None else gs[1]
        if len(S.free)>0:
            B=BMono(S.cells,S.free); hot=state(B,S,1); cold=state(B,S,YBIG)
            r['hot']=hot; r['cold']=cold
            if cold: r['gs_free']=round(cold['nmono']); r['S0']=cold['lnZ']-r['gs_free']*math.log(YBIG)
        else:
            B=BMono(S.cells,set()); c0=state(B,S,1); r['cold']=c0; r['hot']=c0; r['S0']=c0['lnZ'] if c0 else None
        out[key]=r
    return out
if __name__=='__main__':
    L=int(sys.argv[1]); out=scan(L); pickle.dump(out,open(f'data/mixed_scan_{L}.pkl','wb'))
    par=lambda v:'E' if v%2==0 else 'O'
    print('offsets (W,E,S,N mod 4) | type | cells | freeable (H) | free at T=0 | S0/cell | liquid fraction T=0 | at y=1 | <free>/freeable at y=1')
    for k,r in sorted(out.items(),key=lambda kv:(''.join(par(v) for v in kv[0]),kv[0])):
        t=''.join(par(v) for v in k)
        if not r.get('cold'): print(k,t,r['cells'],f"{r['nfree']} ({r['fH']})",'no state in the effective model'); continue
        c=r.get('cold'); h=r.get('hot'); s0=r.get('S0')
        print(k,t,r['cells'],f"{r['nfree']} ({r['fH']})",round(c['nmono']),'-' if s0 is None else f'{s0/r["cells"]:.4f}',
              '-' if not c else f"{c['liqfrac']:.3f}", '-' if not h else f"{h['liqfrac']:.3f}", '-' if not h else f"{h['nmono']/max(r['nfree'],1):.3f}",
              '' if not c else f"err {max(c['err'],h['err']):.0e}")
