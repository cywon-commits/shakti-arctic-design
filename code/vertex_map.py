"""Map unhappy z3 vertices to plaquettes.  In the ground state each unhappy vertex should be one dimer of the effective model;
find the geometric rule, then list the boundary points where an unhappy z3 vertex touches a single plaquette."""
import sys, random
from collections import Counter, defaultdict
from gen_ilp import Sample, halfplanes, ppts
from vertex_sample import Sampler
L=int(sys.argv[1]); rng=random.Random(4)
S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
Sm=Sampler(S,Egs+1); gs=Sm.labels(Sm.sample(Egs,rng))
z3=[p for p,l in gs.items() if l in ('z3happy','z3unhappy')]
touch=defaultdict(list)
for c in S.cells:
    for p in ppts(c): touch[p].append(c)
print('z3 points by number of plaquettes touching:',Counter(len(touch[p]) for p in z3))
U=[p for p in z3 if gs[p]=='z3unhappy']; print('unhappy in GS:',len(U),' E_gs =',Egs,'; touching counts',Counter(len(touch[p]) for p in U))
# effective-model dimers
r=S.lp(); dim=r[2]; print('LP dimers',len(dim))
common=Counter()
for a,b in dim:
    cp=ppts(a)&ppts(b); common[tuple(sorted(gs.get(p,'none') for p in cp))]+=1
print('labels on the common points of each GS dimer pair:',common.most_common(5))
