"""Level-1 states: label changes on each edge strip (within 1 of the edge) and in the interior."""
import sys, random
from collections import Counter
from gen_ilp import Sample, halfplanes
from vertex_sample import Sampler
L,k,n=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); rng=random.Random(3)
S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
Sm=Sampler(S,Egs+k); gs=Sm.labels(Sm.sample(Egs,rng))
xs=[p[0] for p in gs]; ys=[p[1] for p in gs]; X0,X1,Y0,Y1=min(xs),max(xs),min(ys),max(ys)
def zone(p):
    x,y=p; z=[]
    if x<=X0+1: z.append('W')
    if y<=Y0+1: z.append('S')
    if x>=X1-1: z.append('E')
    if y>=Y1-1: z.append('N')
    return ''.join(z) or 'in'
tot=Counter(); per=Counter()
for i in range(n):
    labs=Sm.labels(Sm.sample(Egs+k,rng))
    ch=Counter((zone(p),gs[p]+'>'+labs[p]) for p in labs if labs[p]!=gs[p])
    for key,v in ch.items(): tot[key]+=v
    per[tuple(sorted({z for (z,_) in ch}))]+=1
print('zones touched per sample:',per.most_common(8))
print('label changes (zone, gs>excited), total over samples:')
for key,v in sorted(tot.items(),key=lambda kv:-kv[1])[:25]: print('  ',v,key)
