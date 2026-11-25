"""Corner-switch sample: rectangle with two adjacent block-aligned edges (W, S) and two adjacent through-block fringed edges (E, N).
Exact free-plaquette fraction, entropy and liquid fraction versus y = exp(J/2T) = 2^k.  Usage: python switch_corner.py L1 L2 ..."""
import sys, math, time, pickle, numpy as np
from gen_ilp import Sample, halfplanes
from bmono import BMono
from mixed_scan import state
out={}
for L in [int(v) for v in sys.argv[1:]]:
    S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); B=BMono(S.cells,S.free); m=len(S.free)
    print(f'L={L}: cells {len(S.cells)}, freeable {m}',flush=True)
    for k in (0,1,2,3,4,5,6,8,10,14,20,40):
        t=time.time(); st=state(B,S,2**k); out[(L,k)]=st
        print(f'   y=2^{k:<2d}  (ln y)/m = {k*math.log(2)/m:6.3f}   <free>/freeable = {st["nmono"]/m:.4f}   [ln Z - <free> ln y]/cell = {(st["lnZ"]-st["nmono"]*k*math.log(2))/len(S.cells):.4f}   liquid fraction = {st["liqfrac"]:.3f}   [{time.time()-t:.0f}s]',flush=True)
    pickle.dump(out,open('data/switch_corner.pkl','wb'))
