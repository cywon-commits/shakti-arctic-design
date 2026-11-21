"""Boundary-height slope of straight cuts of the Shakti lattice (all plaquettes paired): exact count per period for the four offset classes k mod 4."""
import io,contextlib,sys,math
with contextlib.redirect_stdout(io.StringIO()):
    import core
from core import rect, black
cache={}
def R_(c):
    if c not in cache: cache[c]=rect(c)
    return cache[c]
def pts(r): return {(p[0]+a,p[1]+b) for p in r for a in (0,1) for b in (0,1)}
def chess_window(X,Y,k,N):
    # chess cells whose Shakti image lies near the cut line: scan a big chess window
    Rr=int(3*N*math.hypot(X,Y))+12
    return [(i,j) for i in range(-Rr,Rr) for j in range(-Rr,Rr)]
def cut(X,Y,k,N=2):
    """half-plane  -Y x + X y <= k+0.5  (Shakti units, integer form); period 4(X,Y); returns (dh per period, #dangling per period, H plaquettes - V plaquettes per period inside a strip)"""
    per=(4*X,4*Y); L2=X*X+Y*Y
    lo=0; hi=N*16*L2       # along coordinate s = X x + Y y  (integer), window = N periods
    def nrm(p): return -Y*p[0]+X*p[1]
    def alo(p): return X*p[0]+Y*p[1]
    # candidate cells: those whose plaquette lies within normal distance <= 12 of the line and along in [lo-40L2, hi+40L2]
    Rr=int(4*N*math.hypot(X,Y))+16
    S=set(); near=[]
    for i in range(-3*Rr,3*Rr):
        for j in range(-3*Rr,3*Rr):
            r=R_((i,j)); P=pts(r)
            cx=sum(p[0] for p in P)/len(P); cy=sum(p[1] for p in P)/len(P)
            if abs(nrm((cx,cy))-k)>8*math.hypot(X,Y): continue
            a=alo((cx,cy))
            if a<lo-30*L2 or a>hi+30*L2: continue
            near.append((i,j))
            if all(nrm(p)<=k+0.5 for p in P): S.add((i,j))
    dh=0; nd=0
    for c in S:
        for d in [(1,0),(-1,0),(0,1),(0,-1)]:
            n=(c[0]+d[0],c[1]+d[1])
            if n in S: continue
            e=pts(R_(c))&pts(R_(n)); m=(sum(p[0] for p in e)/2,sum(p[1] for p in e)/2)
            if abs(nrm(m)-k)>6*math.hypot(X,Y): continue
            if lo<=alo(m)<hi:
                dh+= 1 if black(c) else -1; nd+=1
    return dh/N, nd/N
if __name__=='__main__':
    print('direction (X,Y) | height rise per primitive step for cut class k=0,1,2,3 | support |X|+|Y|')
    for X,Y in [(1,0),(1,1),(2,1),(3,1),(3,2),(4,1),(4,3),(5,1),(5,2),(5,3),(5,4),(6,1),(7,2)]:
        m=[int(cut(X,Y,k,N=1)[0]/4) for k in range(4)]
        print(f'({X},{Y})   {m}   {X+Y}')
