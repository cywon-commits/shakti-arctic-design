"""Block-aligned square at finite temperature: mean fraction of free edge plaquettes versus x = T_on/T,
with y = exp(J/2T) = 2^(x (n+2)/2) and n the order of the Aztec diamond left at T = 0 (exact, bmono.py)."""
import sys, time, numpy as np
from gen_ilp import Sample, halfplanes
from bmono import BMono
xs=(0,0.25,0.5,0.75,0.9,1.0,1.1,1.25,1.5,2.0)
print('side  cells freeable  n   '+' '.join(f'x={x:<5}' for x in xs))
for a in [int(v) for v in sys.argv[1:]] or [16,24,40,56]:
    S=Sample(halfplanes((1,0,0,a),(0,1,0,a)),W=a+4,mode='A'); n=a//2-2; B=BMono(S.cells,S.free)
    f=[]
    for x in xs:
        mono,P=B.stats(2.0**(x*(n+2)/2)); f.append(mono.sum()/len(S.free))
    print(f'{a:4d} {len(S.cells):6d} {len(S.free):6d} {n:4d}   '+' '.join(f'{v:7.4f}' for v in f),flush=True)
