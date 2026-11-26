"""Exact sampling of vertex-model states of the corner-switch sample at fixed energy E_gs + k (forward transfer matrix, backward sampling),
and a per-plaquette count of unhappy vertices, to identify the excitations missing from the boundary-monomer model."""
import sys, itertools, random, pickle, numpy as np
from collections import defaultdict, Counter
from gen_ilp import Sample, halfplanes, COST, ppts
from vertex_dos import lab
class Sampler:
    def __init__(s,S,cap):
        inc=defaultdict(list)
        for k,(a,b) in enumerate(S.isl): inc[a].append((k,'start',a,b)); inc[b].append((k,'end',a,b))
        s.inc=inc; s.pts=sorted(inc,key=lambda p:(p[1],p[0])); s.cap=cap
        A=np.zeros(cap+1); A[0]=1.0; axes=[]; s.T=[(A,axes)]
        for p in s.pts:
            L=inc[p]; ks=[k for k,_,_,_ in L]; old=[k for k in ks if k in axes]; new=[k for k in ks if k not in axes]
            keep=[k for k in axes if k not in old]; naxes=keep+new
            B=np.zeros((2,)*len(naxes)+(cap+1,))
            Am=np.moveaxis(A,[axes.index(k) for k in old],range(len(old))) if old else A
            for bits in itertools.product((1,-1),repeat=len(L)):
                c=COST[lab(L,bits)]
                if c>cap: continue
                val=dict(zip(ks,bits)); src=Am[tuple(0 if val[k]==1 else 1 for k in old)]
                dst=B[(slice(None),)*len(keep)+tuple(0 if val[k]==1 else 1 for k in new)]
                if c==0: dst+=src
                else: dst[...,c:]+=src[...,:-c]
            A,axes=B,naxes; s.T.append((A,axes))
    def sample(s,E,rng):
        spin={}; e=E
        for t in range(len(s.pts),0,-1):
            p=s.pts[t-1]; L=s.inc[p]; ks=[k for k,_,_,_ in L]; Aprev,axprev=s.T[t-1]
            old=[k for k in ks if k in axprev]; opts=[];w=[]
            for ob in itertools.product((1,-1),repeat=len(old)):
                val=dict(spin); val.update(zip(old,ob)); bits=[val[k] for k in ks]; c=COST[lab(L,bits)]
                if c>e: continue
                idx=tuple(0 if val[k]==1 else 1 for k in axprev)+(e-c,)
                wt=Aprev[idx]
                if wt>0: opts.append((ob,c)); w.append(wt)
            i=rng.choices(range(len(opts)),weights=w)[0]; ob,c=opts[i]; spin.update(zip(old,ob)); e-=c
        assert e==0; return spin
    def labels(s,spin):
        return {p:lab(s.inc[p],[spin[k] for k,_,_,_ in s.inc[p]]) for p in s.pts}
if __name__=='__main__':
    L=int(sys.argv[1]); k=int(sys.argv[2]); n=int(sys.argv[3]); rng=random.Random(1)
    S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
    Sm=Sampler(S,Egs+k); cells=sorted(S.cells)
    pp={c:ppts(c) for c in cells}
    def plaq_counts(labs):
        U={p for p,l in labs.items() if l=='z3unhappy'}
        return {c:len(pp[c]&U) for c in cells},U
    gs=Sm.labels(Sm.sample(Egs,rng)); cg,Ug=plaq_counts(gs)
    print('ground state: unhappy',len(Ug),' plaquette counts',Counter(cg.values()),' free plaquettes with 0:',sum(cg[c]==0 for c in S.free),'of',m)
    stats=Counter(); out=[]
    for i in range(n):
        labs=Sm.labels(Sm.sample(Egs+k,rng)); cc,U=plaq_counts(labs)
        hist=Counter(cc.values()); nfree=sum(cc[c]==0 for c in S.free); n0bulk=sum(cc[c]==0 for c in cells if c not in S.free)
        key=(nfree,tuple(sorted(hist.items())),n0bulk); stats[key]+=1; out.append((labs,cc))
    print(f'E_gs+{k}, {n} samples: (free plaquettes with 0 unhappy, histogram of unhappy per plaquette, non-freeable with 0)')
    for key,v in stats.most_common(12): print('  ',v,key)
    pickle.dump((S.cells,S.free,out),open(f'data/vertex_sample_{L}_{k}.pkl','wb'))
