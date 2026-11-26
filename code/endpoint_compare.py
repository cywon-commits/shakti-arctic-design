"""Level-1 states, endpoint by endpoint: vertex model (exact enumeration) versus the Pfaffian model with half-dimer sites.
Vertex endpoints are the two new unhappy boundary z3 points; each is mapped to a site: the freeable plaquette (if the point is shared by a freeable
plaquette and an adjacent one), or the point itself (half-dimer site).  Model counts per site pair come from second differences of the Pfaffian."""
import sys, pickle, itertools, numpy as np
from collections import defaultdict, Counter
from gen_ilp import Sample, halfplanes, ppts
from bmono import BMono
import corrected_model as cm
L=int(sys.argv[1])
S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
cells,free,states,gs=pickle.load(open(f'data/vertex_enum_{L}_1.pkl','rb'))
touch=defaultdict(list)
for c in S.cells:
    for p in ppts(c): touch[p].append(c)
deg=defaultdict(int)
for a,b in S.isl: deg[a]+=1; deg[b]+=1
def site(p):
    t=touch[p]
    if len(t)==2 and abs(t[0][0]-t[1][0])+abs(t[0][1]-t[1][1])==1:
        f=[c for c in t if c in S.free]; return ('F',f[0]) if len(f)==1 else ('D',p)
    return ('H',p)
U0={p for p,l in gs.items() if l=='z3unhappy'}
V=Counter()
for labs in states:
    U={p for p,l in labs.items() if l=='z3unhappy'}; new=[p for p in U-U0 if len(touch[p])<3]
    V[tuple(sorted(site(p) for p in new))]+=1
odd={k:v for k,v in V.items() if len(k)!=2}; print("states with != 2 new boundary points:",odd); V=Counter({k:v for k,v in V.items() if len(k)==2})
for k in V: V[k]//=2
# model: half sites = H points (each owned by the plaquettes it touches; a point touching two non-adjacent plaquettes is tried on each)
H=sorted({s[1] for k in V for s in k if s[0]=='H'})
print('vertex endpoint site types:',Counter(s[0] for k in V for s in k))
own={p:touch[p] for p in H}
print('half-dimer points and their plaquettes:',{p:len(own[p]) for p in H})
# per-pair model counts with one owner choice per point (first owner); level-1 coefficient for activating exactly the two sites
def level1(active):     # active: list of sites; freeable plaquettes always monomer-capable (weight 1), only the active H points add weight t
    n=defaultdict(int)
    for s in active:
        if s[0]=='H':
            for c in own[s[1]]: n[c]+=1      # both owners (variant 'both')
    c=cm.poly(S,n); return c[Egs+1]
base=level1([]); singles={}
rows=[]
for k,v in sorted(V.items(),key=lambda kv:-kv[1]):
    if all(s[0]=='F' for s in k): model=None
    hs=[s for s in k if s[0]=='H']
    if len(hs)==0: tag='F-F'; mc=None
    else:
        key=tuple(hs)
        if len(hs)==1:
            mc=level1(hs)-base
        else:
            mc=level1(hs)-level1(hs[:1])-level1(hs[1:])+base
        tag='-'.join(s[0] for s in k)
    rows.append((k,v,mc,tag))
fF=sum(v for k,v,mc,t in rows if t=='F-F'); print(f'L={L}: F-F (both endpoints un-freed freeable plaquettes): vertex {fF}, model W/S-W/S paths {base}')
for k,v,mc,t in rows[:25]:
    if t!='F-F': print('  ',t,[s[1] for s in k],'vertex',v,' model',mc, ' ratio' ,round(v/mc,3) if mc else None)
