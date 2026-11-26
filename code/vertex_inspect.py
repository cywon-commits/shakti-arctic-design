"""Inspect sampled vertex-model states (data/vertex_sample_L_k.pkl): unhappy z3 vertices as dimers between plaquettes
(a z3 vertex sits at a long-side midpoint; it is shared by the plaquettes having it as a long-side midpoint)."""
import sys, pickle
from collections import Counter, defaultdict
from gen_ilp import R
def mids(c):
    r=R(c); xs=sorted({p[0] for p in r}); ys=sorted({p[1] for p in r})
    if len(xs)==2: x=xs[1]; return [(x,ys[0]),(x,ys[0]+1)]          # horizontal plaquette: long sides y = ys[0], ys[0]+1, midpoint x = xs[1]
    y=ys[1]; return [(xs[0],y),(xs[0]+1,y)]
L,k=int(sys.argv[1]),int(sys.argv[2])
cells,free,out=pickle.load(open(f'data/vertex_sample_{L}_{k}.pkl','rb'))
own=defaultdict(list)
for c in cells:
    for p in mids(c): own[p].append(c)
stat=Counter(); ex=[]
for labs,_ in out:
    U=[p for p,l in labs.items() if l=='z3unhappy']; other=Counter(l for p,l in labs.items() if l not in ('z3unhappy','z3happy','z2ok','z1','z4I'))
    nonmid=[p for p in U if p not in own]
    deg=Counter()
    for p in U:
        for c in own.get(p,[]): deg[c]+=1
    cov=Counter(deg[c] for c in cells); half=sum(1 for p in U if len(own.get(p,[]))==1)
    key=(len(nonmid),half,tuple(sorted(cov.items())),sum(deg[c]==0 for c in free),tuple(sorted(other.items())))
    stat[key]+=1
print(f'L={L}, E_gs+{k}: (unhappy z3 not at a midpoint, unhappy at a one-plaquette midpoint, histogram of dimers per plaquette, free plaquettes uncovered, other excited vertices)')
for key,v in stat.most_common(15): print('  ',v,key)
