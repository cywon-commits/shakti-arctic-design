"""Plus-shaped samples: exact pairing probabilities by direction (urban renewal), to colour frozen cells by ice phase.
Phase of a cell = direction d = (white cell) - (black cell) of its dominant pairing; four directions = four ice phases (domino types)."""
import sys, time, pickle, numpy as np
from collections import defaultdict
from gen_ilp import black
from verify_shuffle import cells_of, centre, perfect_matching, region_probabilities
from theory import closed_form
DIRS=[(1,0),(0,1),(-1,0),(0,-1)]
def run(n,c):
    t=time.time()
    S=cells_of(lambda p:(abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n),n+6)
    big=cells_of(lambda p: abs(p[0])<=n and abs(p[1])<=n,n+6)
    cx=int(round(np.mean([a+.5 for a,b in big]))); cy=int(round(np.mean([b+.5 for a,b in big]))); N=int(round(max(abs(a+.5-cx)+abs(b+.5-cy) for a,b in big)))
    P=region_probabilities(S,N,cx,cy,perfect_matching(big-S))
    pd=defaultdict(lambda: np.zeros(4))
    for (a,b),p in P.items():
        bb,ww=(a,b) if black(a) else (b,a); d=DIRS.index((ww[0]-bb[0],ww[1]-bb[1]))
        pd[a][d]+=p; pd[b][d]+=p
    q=sorted(S); cen=np.array([centre(k) for k in q]); prob=np.array([pd[k] for k in q])
    print(f'n={n} c={c}: cells {len(q)}, max |sum - 1| {abs(prob.sum(1)-1).max():.1e}  [{time.time()-t:.0f}s]',flush=True)
    return dict(n=n,c=c,cen=cen,prob=prob,sol=closed_form(n,c))
if __name__=='__main__':
    out={}
    for a in sys.argv[1:]:
        n,c=map(int,a.split(',')); out[(n,c)]=run(n,c); pickle.dump(out,open('data/phase_maps.pkl','wb'))
