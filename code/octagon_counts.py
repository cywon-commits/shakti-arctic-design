"""Octagon (axis sides + even 45-degree sides) has exactly the tilings of the Aztec diamond of the inner square |X|,|Y|<=d-n; the plus-shaped sample is a different region."""
from eff import *
def region(inside,W=70): return Sample(inside,W=W,mode='B').cells
def plus(n,c): return lambda p: (abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n)
print('octagon |X|,|Y|<=n, |X+-Y|<=d   vs   plus: arms half-width c=d-n, reach n   (all-paired dimer counts)')
for n,d in [(13,18),(13,20),(13,22),(17,24),(17,26),(17,28),(21,30),(21,34),(25,36),(25,40)]:
    octa=region(halfplanes((1,0,-n,n),(0,1,-n,n),(1,1,-d,d),(1,-1,-d,d)))
    pl=region(plus(n,d-n))
    a,b=ntil(octa),ntil(pl)
    print(f' n={n} d={d} c={d-n}: octagon cells {len(octa)}, lnZ={a if a is None else round(a,6)} | plus cells {len(pl)}, lnZ={b if b is None else round(b,6)} | plus subset of octagon: {pl<=octa}; forced cells {len(octa)-len(pl)}; predicted Aztec(order c-1) lnZ={(d-n-1)*(d-n)/2*math.log(2):.4f}')
