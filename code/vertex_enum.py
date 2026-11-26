"""Exact enumeration of all vertex-model states at E_gs + k (small samples) by backward traversal of the transfer tables.
Reports, for every state, the unhappy z3 vertices on boundary points (touching fewer than 3 plaquettes) that are not unhappy in the ground state."""
import sys, itertools, pickle
from collections import defaultdict, Counter
from gen_ilp import Sample, halfplanes, COST, ppts
from vertex_dos import lab
from vertex_sample import Sampler
def enum(Sm,E):
    out=[]
    def rec(t,spin,e):
        if t==0:
            if e==0: out.append(dict(spin))
            return
        p=Sm.pts[t-1]; L=Sm.inc[p]; ks=[k for k,_,_,_ in L]; Aprev,axprev=Sm.T[t-1]
        old=[k for k in ks if k in axprev]
        for ob in itertools.product((1,-1),repeat=len(old)):
            val=dict(spin); val.update(zip(old,ob)); c=COST[lab(L,[val[k] for k in ks])]
            if c>e: continue
            if Aprev[tuple(0 if val[k]==1 else 1 for k in axprev)+(e-c,)]>0: rec(t-1,val,e-c)
    rec(len(Sm.pts),{},E); return out
if __name__=='__main__':
    L,k=int(sys.argv[1]),int(sys.argv[2])
    S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
    Sm=Sampler(S,Egs+k); gss=enum(Sm,Egs); sts=enum(Sm,Egs+k); print(f'L={L}: ground states {len(gss)}, states at E_gs+{k}: {len(sts)}')
    touch=defaultdict(list)
    for c in S.cells:
        for p in ppts(c): touch[p].append(c)
    U0={p for p,l in Sm.labels(gss[0]).items() if l=='z3unhappy'}
    stat=Counter()
    for s in sts:
        U={p for p,l in Sm.labels(s).items() if l=='z3unhappy'}
        newb=tuple(sorted((len(touch[p]),p) for p in U-U0 if len(touch[p])<3)); lostb=tuple(sorted((len(touch[p]),p) for p in U0-U if len(touch[p])<3))
        stat[(newb,lostb)]+=1
    print('(new unhappy boundary points, lost unhappy boundary points): number of states')
    for key,v in stat.most_common(40): print('  ',v,key)
    pickle.dump((S.cells,S.free,[Sm.labels(s) for s in sts],Sm.labels(gss[0])),open(f'data/vertex_enum_{L}_{k}.pkl','wb'))
