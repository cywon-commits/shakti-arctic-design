from eff import *
from peel import peel
def region(inside,W): return Sample(inside,W=W,mode='B').cells
def plus(n,c): return lambda p: (abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n)
L2=math.log(2)
for n in (25,37):
    print(f'n={n}: n/3={n/3:.2f}, n/sqrt2={n/math.sqrt(2):.2f}; Aztec(full square)={(n-1)*n//2} ln2')
    for c in range(3,n+1,2):
        S=region(plus(n,c),W=n+6); ld=ntil(S)
        Sr,forced,ok=peel(S)
        five=5*c*(c-1)//2
        ld=ld if ld is not None else float('nan'); print(f'   c={c:2d} c/n={c/n:.3f}  cells {len(S)} H{sum(black(q) for q in S)}  lnZ={ld/L2:10.4f} ln2   5 Aztec(c): {five:4d}   forced cells by peeling: {2*len(forced)}')
