"""Finite-size check of the class-10 arctic curve (double precision is reliable only up to about n = 61; see the precision line): exact Kasteleyn statistics on a large plus sample vs winding-number classification by the theoretical curve."""
import sys, math, time, pickle, numpy as np
from collections import defaultdict
from gen_ilp import Sample, black, ppts
from theory import solve, curve
def region(inside,W): return Sample(inside,W=W,mode='B').cells
def plus(n,c): return lambda p: (abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n)
def pcenter(c):
    P=ppts(c); return (sum(p[0] for p in P)/len(P), sum(p[1] for p in P)/len(P))
def maxprob(S):
    B=sorted(c for c in S if black(c)); Wt=sorted(c for c in S if not black(c)); wi={c:k for k,c in enumerate(Wt)}
    K=np.zeros((len(B),len(B)))
    for k,c in enumerate(B):
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            n=(c[0]+d[0],c[1]+d[1])
            if n in wi: K[k,wi[n]]=1.0 if d[1]==0 else (-1.0)**(c[0])
    s,ld=np.linalg.slogdet(K); Ki=np.linalg.inv(K); mx=defaultdict(float)
    P=K*Ki.T; eb=np.abs(P.sum(1)-1); ew=np.abs(P.sum(0)-1)
    maxprob.err={**{c:float(eb[k]) for k,c in enumerate(B)},**{c:float(ew[k]) for k,c in enumerate(Wt)}}   # double precision is not enough for large n: every cell's probabilities must sum to 1
    for k,c in enumerate(B):
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            n=(c[0]+d[0],c[1]+d[1])
            if n in wi:
                p=K[k,wi[n]]*Ki[wi[n],k]; mx[c]=max(mx[c],p); mx[n]=max(mx[n],p)
    return mx,ld
def winding(px,py,x,y):
    a=np.arctan2(y-py,x-px); d=np.diff(a); d=(d+np.pi)%(2*np.pi)-np.pi
    return int(round(d.sum()/(2*np.pi)))
if __name__=='__main__':
    n=int(sys.argv[1]); out={}
    for c in map(int,sys.argv[2:]):
        t=time.time(); S=region(plus(n,c),W=n+6); mx,ld=maxprob(S)
        sol=solve(n,c); assert sol and len(sol)==1, sol
        th,v,x,y=curve(*sol[0],m=20001)
        cen={q:pcenter(q) for q in S}; wn={q:winding(cen[q][0],cen[q][1],x,y) for q in S}
        # distance to curve
        dist={q:float(np.min(np.hypot(x-cen[q][0],y-cen[q][1]))) for q in S}
        tab=defaultdict(lambda:[0,0,0.0])
        for q in S:
            tab[wn[q]][0]+=1; tab[wn[q]][1]+= mx[q]>0.9; tab[wn[q]][2]+=1-mx[q]
        print(f'n={n} c={c} (c/n={c/n:.3f}): cells {len(S)}, lnZ={ld/math.log(2):.3f} ln2, p,a0,a4={tuple(round(float(s),5) for s in sol[0])}  [{time.time()-t:.0f}s]')
        for w in sorted(tab): print(f'    winding {w:+d}: {tab[w][0]:6d} cells, frozen (Pmax>0.9): {tab[w][1]:6d}, mean liquidity 1-Pmax = {tab[w][2]/tab[w][0]:.4f}')
        err=maxprob.err; good={q for q in S if err[q]<1e-4}
        print(f'    precision: max |sum P - 1| = {max(err.values()):.1e}; cells with error < 1e-4: {len(good)} of {len(S)}; unreliable cells: winding 0: {sum(1 for q in S if q not in good and wn[q]==0)}, winding 1: {sum(1 for q in S if q not in good and wn[q]!=0)}')
        for dmin in (0,1,2,3,4,6,8):
            far=[q for q in good if dist[q]>dmin]; bad=sum(1 for q in far if (wn[q]==0)!=(mx[q]>0.9))
            print(f'    reliable cells farther than {dmin} from the curve: {len(far):6d}, misclassified: {bad}')
        for dmin in ():
            far=[q for q in S if dist[q]>dmin]
            bad=sum(1 for q in far if (wn[q]==0)!=(mx[q]>0.9))
            print(f'    cells farther than {dmin} from the curve: {len(far):6d}, misclassified: {bad}')
        out[c]=dict(cen=cen,mx=dict(mx),wn=wn,dist=dist,err=err,sol=sol[0],lnZ=ld)
    pickle.dump(out,open(f'data/verify_n{n}.pkl','wb'))
