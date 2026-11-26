"""sigma_N(x) of the through-block Shakti Aztec square, corrected model, in floating point (monomer probabilities, row-sum check)."""
import sys, math, numpy as np, pickle
from gen_ilp import Sample, halfplanes
import corrected_model as cm
from corrected_map import stats
xs=[float(v) for v in sys.argv[1].split(',')]; out={}
for L in [int(v) for v in sys.argv[2:]]:
    N=L//2-1; S=Sample(halfplanes((1,0,1,L+1),(0,1,1,L+1)),W=L+10,mode='A'); n=cm.halfsites(S,'fit'); row=[]
    for x in xs:
        t=math.exp(-(x*N*math.log(2)+2*math.log(2))); st=stats(S,n,t); row.append((1-st['half']/(2*N),st['err']))
    out[N]=row; print(f'N={N:3d}: '+'  '.join(f'x={x}: {s:.5f} (err {e:.0e})' for x,(s,e) in zip(xs,row)),flush=True)
pickle.dump((xs,out),open('data/aztecO_float.pkl','wb'))
