"""Where do sampled excited states differ from the ground state?  Labels of lattice points compared with the (unique) ground state;
differing points grouped into clusters (points within distance 2), and each cluster's contact with the four edges reported."""
import sys, pickle, random
from collections import Counter
from gen_ilp import Sample, halfplanes
from vertex_sample import Sampler
L,k,n=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); rng=random.Random(2)
S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
Sm=Sampler(S,Egs+k); gs=Sm.labels(Sm.sample(Egs,rng))
xs=[p[0] for p in gs]; ys=[p[1] for p in gs]; X0,X1,Y0,Y1=min(xs),max(xs),min(ys),max(ys)
def edges_of(cl):
    e=set()
    for x,y in cl:
        if x<=X0+1: e.add('W')
        if y<=Y0+1: e.add('S')
        if x>=X1-1: e.add('E')
        if y>=Y1-1: e.add('N')
    return ''.join(sorted(e)) or '-'
stat=Counter(); sizes=Counter()
for i in range(n):
    labs=Sm.labels(Sm.sample(Egs+k,rng)); D=[p for p in labs if labs[p]!=gs[p]]
    left=set(D); cls=[]
    while left:
        st=[left.pop()]; cl=set(st)
        while st:
            a=st.pop()
            for b in list(left):
                if abs(a[0]-b[0])<=2 and abs(a[1]-b[1])<=2: left.discard(b); cl.add(b); st.append(b)
        cls.append(cl)
    key=tuple(sorted(edges_of(c) for c in cls)); stat[key]+=1
    for c in cls: sizes[(edges_of(c),len(c)//4*4)]+=1
print(f'L={L}, E_gs+{k}, {n} samples. Clusters of changed points, labelled by the edges they touch (W,S switchable; N,E fixed):')
for key,v in stat.most_common(14): print('  ',v,key)
