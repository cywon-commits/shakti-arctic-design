"""Class-10 curve versus exact statistics at large sizes (urban-renewal recursion, plus_shuffle.py).
Usage: python verify_shuffle.py 61,31 121,61 241,121 481,241"""
import sys, math, time, pickle, io, contextlib, numpy as np
from scipy.spatial import cKDTree
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, black
from shuffle_probs import region_probabilities
from plus_shuffle import perfect_matching, cell_stats
from theory import closed_form, curve
def cells_of(pred,W):
    out=set()
    for i in range(-W,W):
        for j in range(-W,W):
            r=rect((i,j))
            if all(pred((p[0]+a,p[1]+b)) for p in r for a in (0,1) for b in (0,1)): out.add((i,j))
    return out
def centre(c):
    r=rect(c); return (sum(p[0] for p in r)/2+0.5, sum(p[1] for p in r)/2+0.5)
def winding(pts,x,y,chunk=4000):
    w=np.empty(len(pts),int)
    for s in range(0,len(pts),chunk):
        p=pts[s:s+chunk]; a=np.arctan2(y[None,:]-p[:,1:2],x[None,:]-p[:,0:1]); d=np.diff(a,axis=1); d=(d+np.pi)%(2*np.pi)-np.pi
        w[s:s+chunk]=np.rint(d.sum(1)/(2*np.pi)).astype(int)
    return w
def run(n,c):
    t=time.time()
    S=cells_of(lambda p:(abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n),n+6)
    big=cells_of(lambda p: abs(p[0])<=n and abs(p[1])<=n,n+6)
    cx=int(round(np.mean([a+.5 for a,b in big]))); cy=int(round(np.mean([b+.5 for a,b in big]))); N=int(round(max(abs(a+.5-cx)+abs(b+.5-cy) for a,b in big)))
    P=region_probabilities(S,N,cx,cy,perfect_matching(big-S)); mx,tot=cell_stats(S,P)
    q=sorted(S); cen=np.array([centre(k) for k in q]); liq=np.array([1-mx[k] for k in q]); err=max(abs(tot[k]-1) for k in q)
    th,v,x,y=curve(*closed_form(n,c),m=8001)
    wn=winding(cen,x,y); dist,_=cKDTree(np.c_[x,y]).query(cen)
    print(f'n={n} c={c} (c/n={c/n:.3f}): Aztec order {N}, cells {len(q)}, max |sum P - 1| = {err:.1e}  [{time.time()-t:.0f}s]',flush=True)
    frozen=liq<0.1
    for dmin in (0,1,2,3,4,6,8,12):
        far=dist>dmin; bad=int(((wn==0)!=frozen)[far].sum()); print(f'    cells farther than {dmin:2d} from the curve: {int(far.sum()):7d}, misclassified: {bad}')
    # width of the frozen-side tail on the smooth arcs (outside the curve, away from the reflex corners)
    sel=(wn==0)&(np.maximum(abs(cen[:,0]),abs(cen[:,1]))>c+0.1*n)
    bins=np.arange(0,0.2*n,0.5); mid=(bins[:-1]+bins[1:])/2; mean=np.array([liq[sel&(dist>=a)&(dist<b)].mean() if (sel&(dist>=a)&(dist<b)).sum()>8 else np.nan for a,b in zip(bins[:-1],bins[1:])])
    ws=[]
    for level in (1e-1,1e-2,1e-3,1e-4):
        k=np.where(mean<level)[0]; ws.append(float(mid[k[0]]) if len(k) else float('nan'))
    print('    distance at which the mean liquidity outside the curve falls below 1e-1, 1e-2, 1e-3, 1e-4:',ws,flush=True)
    return dict(n=n,c=c,cen=cen,liq=liq,wn=wn,dist=dist,widths=ws,sol=closed_form(n,c))
if __name__=='__main__':
    out={}
    for a in sys.argv[1:]:
        n,c=map(int,a.split(',')); out[(n,c)]=run(n,c); pickle.dump(out,open('data/verify_shuffle.pkl','wb'))
    ns=np.array([k[0] for k in out],float)
    if len(ns)>=3:
        for j,level in enumerate((1e-1,1e-2,1e-3,1e-4)):
            w=np.array([out[k]['widths'][j] for k in out]); ok=~np.isnan(w)
            if ok.sum()>=3: print(f'level {level:g}: widths {w}, fitted exponent {np.polyfit(np.log(ns[ok]),np.log(w[ok]),1)[0]:.2f}')
