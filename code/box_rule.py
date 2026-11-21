"""Extremal-pair (box) rule: two boundary points whose height difference equals their l1 distance freeze their bounding box.
Prediction for a rectangle |X|<=a, |Y|<=b with four even 45-degree cuts of one offset class: the unfrozen region is
    d = 2 (mod 4):  X in [b - d_SW, d_NE - b],  Y in [a - d_SE, d_NW - a]
    d = 0 (mod 4):  X in [b - d_NW, d_SE - b],  Y in [a - d_SW, d_NE - a]      (mirror-image pinwheel)
an Aztec diamond (it must be a square). The offset class selects which ice phase a cut carries.
Checked here against forced-domino peeling on samples that are not symmetric octagons."""
from peel import *
def sample(a,b,dNE,dNW,dSW,dSE):
    return region(lambda p: abs(p[0])<=a and abs(p[1])<=b and p[0]+p[1]<=dNE and -p[0]+p[1]<=dNW and -p[0]-p[1]<=dSW and p[0]-p[1]<=dSE, W=max(a,b)+8)
cases=[(21,13,26,34,26,34),(21,13,30,34,22,34),(25,13,26,38,26,38),(21,17,30,38,34,34),(21,17,34,34,30,38),(25,17,38,34,26,42),(21,21,26,42,34,34),(17,13,22,30,26,26),
       (21,21,28,28,36,36),(21,21,28,36,36,28),(21,21,28,32,36,32),(19,19,24,24,32,32),(23,23,28,28,36,36),(21,13,28,36,28,36),(21,21,28,30,32,30)]
ok=0
for a,b,dNE,dNW,dSW,dSE in cases:
    S=sample(a,b,dNE,dNW,dSW,dSE); nb=sum(black(c) for c in S)
    cls={dNE%4,dNW%4,dSW%4,dSE%4}
    if len(cls)>1: print(f'a,b=({a},{b}) d=({dNE},{dNW},{dSW},{dSE}): mixed offset classes, H-V imbalance {2*sum(black(c) for c in S)-len(S)}'); continue
    pred=(b-dSW,dNE-b,a-dSE,dNW-a) if cls=={2} else (b-dNW,dSE-b,a-dSW,dNE-a); square=(pred[1]-pred[0])==(pred[3]-pred[2])
    if 2*nb!=len(S): print(f'a,b=({a},{b}) d=({dNE},{dNW},{dSW},{dSE}): unbalanced (predicted square: {square})'); continue
    Sr,forced,good=peel(S); bb=bbox(Sr) if Sr else None
    R_=region(lambda p,bb=pred: bb[0]<=p[0]<=bb[1] and bb[2]<=p[2-1]<=bb[3],W=max(a,b)+8)
    k=(pred[1]-pred[0])//2-1; ld=ntil(Sr)
    match = bb==pred and Sr==R_
    ok+=match
    print(f'a,b=({a},{b}) d=({dNE},{dNW},{dSW},{dSE}): predicted box {pred}, peeled box {bb}, same cell set {Sr==R_}, lnZ={ld/math.log(2):.3f} ln2 (Aztec order {k}: {k*(k+1)//2})')
print('matches',ok)
