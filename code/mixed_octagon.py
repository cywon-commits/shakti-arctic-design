"""Rectangle |X|<=a, |Y|<=b with up to four even 45-degree corner cuts of arbitrary offset class: prediction from the box rule versus exact peeling.

Each cut "continues" one of its two adjacent sides and freezes a strip along that side:
    class d = 2 (mod 4): the side before it (counterclockwise order)     NE->right, NW->top, SW->left, SE->bottom
    class d = 0 (mod 4): the side after it                               NE->top,   NW->left, SW->bottom, SE->right
Prediction: the sample is tileable iff no side is continued by two cuts and the remaining rectangle is a square
(side >= 0); the remaining region is then the Aztec diamond on that square.
Samples in which two cuts cross (a side of negative length) are outside the rule and are skipped."""
import io, contextlib, itertools, math, sys, time
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, black
from peel import peel, bbox
from eff import ntil
_pp={}
def ppts(c):
    if c not in _pp: _pp[c]=[(p[0]+i,p[1]+j) for p in rect(c) for i in (0,1) for j in (0,1)]
    return _pp[c]
def cells(a,b,d,W=None):
    W=W or max(a,b)+4; dNE,dNW,dSW,dSE=d
    ok=lambda p: abs(p[0])<=a and abs(p[1])<=b and p[0]+p[1]<=dNE and -p[0]+p[1]<=dNW and -p[0]-p[1]<=dSW and p[0]-p[1]<=dSE
    return {(i,j) for i in range(-W,W) for j in range(-W,W) if all(ok(p) for p in ppts((i,j)))}
SIDE={('NE',2):'R',('NE',0):'T',('NW',2):'T',('NW',0):'L',('SW',2):'L',('SW',0):'B',('SE',2):'B',('SE',0):'R'}
def predict(a,b,d):
    """returns None if predicted untileable, else the remaining box (X0,X1,Y0,Y1)"""
    box={'R':a,'L':-a,'T':b,'B':-b}; used=set()
    for name,di in zip(('NE','NW','SW','SE'),d):
        if di>=a+b: continue                      # no cut at this corner
        s=SIDE[(name,di%4)]
        if s in used: return None
        used.add(s)
        box[s]={'R':di-b,'L':b-di,'T':di-a,'B':a-di}[s]
    w=box['R']-box['L']; h=box['T']-box['B']
    if w!=h or w<0: return None
    return (box['L'],box['R'],box['B'],box['T'])
def actual(a,b,d):
    S=cells(a,b,d); nb=sum(black(c) for c in S)
    if 2*nb!=len(S): return None
    Sr,forced,ok=peel(S)
    if not ok: return None
    if not Sr: return 'empty'
    if ntil(Sr) is None: return None
    return bbox(Sr),len(Sr)
if __name__=='__main__':
    t=time.time(); stats={'both tileable, same box':0,'both untileable':0,'MISMATCH':0}; mixed=0; mixed_tileable=[]; bad=[]; skipped=0
    for a,b in [(13,13),(17,13),(13,9),(17,9)]:
        rng=range(a+b-12,a+b+1,2)
        for d in itertools.product(rng,repeat=4):
            dNE,dNW,dSW,dSE=d
            if min(dNE+dSE-2*a,dNW+dSW-2*a,dNE+dNW-2*b,dSW+dSE-2*b)<0: skipped+=1; continue     # two cuts cross: not an octagon
            pr=predict(a,b,d); ac=actual(a,b,d)
            classes={x%4 for x in d if x<a+b}
            if len(classes)>1: mixed+=1
            if pr is None and ac is None: stats['both untileable']+=1
            elif pr is not None and ac is not None:
                w=pr[1]-pr[0]
                k=w//2-1                                     # Aztec order of a through-block square of side w
                same = (ac=='empty' and k<=0) or (ac!='empty' and ac[0]==pr and ac[1]==2*k*(k+1))
                if same:
                    stats['both tileable, same box']+=1
                    if len(classes)>1: mixed_tileable.append((a,b,d,pr))
                else: stats['MISMATCH']+=1; bad.append((a,b,d,pr,ac))
            else: stats['MISMATCH']+=1; bad.append((a,b,d,pr,ac))
        print(f'a,b=({a},{b}) done, {time.time()-t:.0f}s',stats,flush=True)
    print('skipped (cuts cross):',skipped)
    print('samples with mixed offset classes:',mixed,'; of these tileable:',len(mixed_tileable))
    for x in mixed_tileable[:8]: print('  example',x)
    print('four-cut samples with mixed classes that are tileable:',sum(1 for a,b,d,pr in mixed_tileable if all(x<a+b for x in d)))
    for x in bad[:10]: print('  mismatch',x)
