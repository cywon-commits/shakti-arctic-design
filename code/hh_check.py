"""H-H endpoint pairs at level 1: vertex model versus the 'fit' weights (touch-1 point: weight 2; two-plaquette point: weight 1 on each non-freeable owner)."""
import sys, pickle
from collections import defaultdict, Counter
from gen_ilp import Sample, halfplanes, ppts
import corrected_model as cm
L=int(sys.argv[1])
S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
cells,free,states,gs=pickle.load(open(f'data/vertex_enum_{L}_1.pkl','rb'))
touch=defaultdict(list)
for c in S.cells:
    for p in ppts(c): touch[p].append(c)
def site(p):
    t=touch[p]
    if len(t)==2 and abs(t[0][0]-t[1][0])+abs(t[0][1]-t[1][1])==1:
        f=[c for c in t if c in S.free]; return ('F',f[0]) if len(f)==1 else ('D',p)
    return ('H',p)
U0={p for p,l in gs.items() if l=='z3unhappy'}
V=Counter()
for labs in states:
    U={p for p,l in labs.items() if l=='z3unhappy'}; V[tuple(sorted(site(p) for p in U-U0 if len(touch[p])<3))]+=1
for k in V: V[k]/=2
def wts(p):
    t=touch[p]
    return {t[0]:2} if len(t)==1 else {c:1 for c in t if c not in S.free}
def lev(points):
    n=defaultdict(int)
    for p in points:
        for c,w in wts(p).items(): n[c]+=w
    return cm.poly(S,n)[Egs+1]
base=lev([]); single={}
H=sorted({s[1] for k in V for s in k if s[0]=='H'})
for p in H: single[p]=lev([p])
tv=tm=0; rows=[]
for i,p in enumerate(H):
    for q in H[i:]:
        if p==q: v=V.get((('H',p),),0); mc=None
        else:
            v=V.get(tuple(sorted([('H',p),('H',q)])),0); mc=lev([p,q])-single[p]-single[q]+base
        if mc is None: continue
        tv+=v; tm+=mc
        if v!=mc: rows.append((p,q,v,mc))
print(f'L={L}: H-H totals vertex {tv:.0f} model {tm}; mismatching pairs:')
for r in rows[:20]: print('  ',r)
fv=sum(v for k,v in V.items() if len(k)==2 and {k[0][0],k[1][0]}=={'F','H'}); fm=sum(single[p]-base for p in H)
print(f'F-H totals vertex {fv:.0f} model {fm};  single-endpoint states vertex',{k:v for k,v in V.items() if len(k)!=2})
