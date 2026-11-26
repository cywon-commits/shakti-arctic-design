"""For each half-dimer point, find the plaquette it covers and its weight by matching level-1 endpoint counts of the vertex model.
Single-site counts: states with one F endpoint (any freeable plaquette) and this H endpoint, summed over F; and H-H pairs."""
import sys, pickle, numpy as np
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
    U={p for p,l in labs.items() if l=='z3unhappy'}; new=tuple(sorted(site(p) for p in U-U0 if len(touch[p])<3)); V[new]+=1
for k in V: V[k]/=2
H=sorted({s[1] for k in V for s in k if s[0]=='H'})
FH=Counter()
for k,v in V.items():
    if len(k)==2 and k[0][0]=='F' and k[1][0]=='H': FH[k[1][1]]+=v
    if len(k)==2 and k[1][0]=='F' and k[0][0]=='H': FH[k[0][1]]+=v
base=cm.poly(S,{})[Egs+1]
print(f'L={L}: F-F vertex {sum(v for k,v in V.items() if len(k)==2 and k[0][0]=="F" and k[1][0]=="F")}, model {base}')
for p in H:
    res=[]
    for c in touch[p]:
        mc=cm.poly(S,{c:1})[Egs+1]-base
        res.append((c,'free' if c in S.free else '', mc, round(FH[p]/mc,3) if mc else None))
    print(p,'vertex F-H total',FH[p],' per owner (cell, freeable?, model, ratio):',res)
