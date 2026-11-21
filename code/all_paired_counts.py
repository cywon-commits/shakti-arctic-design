"""Number of tilings (all plaquettes paired) of small Shakti samples by cut class (Kasteleyn determinant)."""
import io,contextlib,sys,math
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import core
from core import rect, black
def pts(r): return {(p[0]+a,p[1]+b) for p in r for a in (0,1) for b in (0,1)}
def region(pred,W=40):
    return {(i,j) for i in range(-W,W) for j in range(-W,W) if all(pred(p) for p in pts(rect((i,j))))}
def ntilings(S):
    B=[c for c in S if black(c)]; Wt=[c for c in S if not black(c)]
    if len(B)!=len(Wt): return 0
    if not B: return 1
    bi={c:k for k,c in enumerate(B)}; wi={c:k for k,c in enumerate(Wt)}
    K=np.zeros((len(B),len(B)),dtype=complex)
    for c in B:
        for d,w in [((1,0),1),((-1,0),1),((0,1),1j),((0,-1),1j)]:
            n=(c[0]+d[0],c[1]+d[1])
            if n in wi: K[bi[c],wi[n]]=w
    s,ld=np.linalg.slogdet(K)
    return round(math.exp(ld)) if ld<700 else math.exp(ld)
def entropy_per_cell(S):
    t=ntilings(S); return (math.log(t)/len(S) if t else None), t, len(S)
if __name__=='__main__':
    print('axis-aligned squares |x|,|y| <= a  (half-plane class: k odd <-> boundary at odd coordinate)')
    for a0,a1,lab in [(0,8,'block-aligned (even)'),(1,9,'through blocks (odd)'),(0,12,'block-aligned (even)'),(1,13,'through blocks (odd)'),(0,16,'even'),(1,17,'odd')]:
        S=region(lambda p:a0<=p[0]<=a1 and a0<=p[1]<=a1)
        s,t,n=entropy_per_cell(S); print(f'  {lab:22s} side {a1-a0}: cells {n}, tilings {t}, ln Z/cell {s}')
    print('45-degree squares |x+y-c|,|x-y-c2|')
    for D,k in [(8,0),(8,1),(8,2),(8,3),(12,0),(12,1),(12,2),(12,3)]:
        S=region(lambda p,D=D,k=k: 0<=p[0]+p[1]-k<=2*D and -D<=p[0]-p[1]-k<=D)
        s,t,n=entropy_per_cell(S); print(f'  D={D} offset {k}: cells {n}, tilings {t}, ln Z/cell {s}')
    print('block-aligned square with its whole edge layer removed = through-block square (nested check)')
    for a in (8,12,16):
        Se=region(lambda p,a=a:0<=p[0]<=a and 0<=p[1]<=a); So=region(lambda p,a=a:1<=p[0]<=a-1 and 1<=p[1]<=a-1); d=Se-So
        print(f'  side {a}: {len(Se)} -> {len(So)} cells, subset {So<=Se}, removed H {sum(black(c) for c in d)} / V {len(d)-sum(black(c) for c in d)}, tilings {ntilings(So)}')
