"""Parity (T-join) model of the low-lying vertex-model states.
Graph: plaquettes; an edge for every z3 point shared by two plaquettes, a half-edge for every z3 point on a single plaquette.
An unhappy z3 vertex = an occupied edge (energy J).  Constraint: every plaquette with its stem island has an odd number of
occupied edges; freeable plaquettes are unconstrained.  Counts are compared with the exact vertex-model density of states."""
import sys, pickle, numpy as np
from collections import defaultdict, Counter
from gen_ilp import Sample, halfplanes, ppts
from vertex_dos import dos
with __import__('contextlib').redirect_stdout(__import__('io').StringIO()):
    from vort import is_island
def graph(S):
    deg=defaultdict(int)
    for a,b in S.isl: deg[a]+=1; deg[b]+=1
    z3=[p for p,d in deg.items() if d==3]
    touch=defaultdict(list)
    for c in S.cells:
        for p in ppts(c): touch[p].append(c)
    E=[]; kinds=Counter()
    for p in z3:
        t=touch[p]
        if len(t)==3:   # bulk: the unhappy vertex pairs the two plaquettes that share an edge through p ... find which pair
            pairs=[(a,b) for i,a in enumerate(t) for b in t[i+1:] if abs(a[0]-b[0])+abs(a[1]-b[1])==1]
            kinds[('bulk',len(pairs))]+=1
            E+= [tuple(pr) for pr in pairs[:1]] if len(pairs)==1 else [('?',p)]
        elif len(t)==2:
            a,b=t; adj=abs(a[0]-b[0])+abs(a[1]-b[1])==1; kinds[('two',adj)]+=1
            E.append((a,b) if adj else ('?',p))
        elif len(t)==1: kinds[('one',)]+=1; E.append((t[0],None))
        else: kinds[('zero',)]+=1
    return E,kinds
def count(S,E,cap):
    cells=sorted(S.cells); idx={c:i for i,c in enumerate(cells)}; odd=[c not in S.free for c in cells]
    last={}
    for k,(a,b) in enumerate(E):
        for c in (a,b):
            if c is not None: last[c]=k
    st={0:np.zeros(cap+1)}; st[0][0]=1.0; active=[]
    for k,(a,b) in enumerate(E):
        bits=[idx[c] for c in (a,b) if c is not None]; mask=sum(1<<i for i in bits)
        new=defaultdict(lambda: np.zeros(cap+1))
        for s,v in st.items():
            new[s]+=v; w=np.zeros(cap+1); w[1:]=v[:-1]; new[s^mask]+=w
        st=dict(new)
        for c in (a,b):                       # retire cells whose last edge is processed: enforce parity
            if c is not None and last[c]==k:
                i=idx[c]; want=1 if odd[i] else None
                nst=defaultdict(lambda: np.zeros(cap+1))
                for s,v in st.items():
                    if want is None or ((s>>i)&1)==want: nst[s&~(1<<i)]+=v
                st=dict(nst)
    for c in cells:
        if c not in last and odd[idx[c]]: return None
    return st[0]
if __name__=='__main__':
    K=int(sys.argv[1])
    for L in [int(v) for v in sys.argv[2:]]:
        S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); m=len(S.free); Egs=(len(S.cells)-m)//2
        E,kinds=graph(S); print(f'L={L}: z3 points by kind {dict(kinds)}')
        if any(a=='?' for a,b in E): print('  ambiguous z3 points; stop'); continue
        tj=count(S,E,Egs+K); g,_,_=dos(S,Egs+K)
        print('   level   vertex g/2        T-join count')
        for k in range(K+1): print(f'   {k:3d}  {g[Egs+k]/2:16.6g}  {tj[Egs+k]:16.6g}   ratio {g[Egs+k]/2/tj[Egs+k] if tj[Egs+k] else float("nan"):.6f}')
        print('   below E_gs: T-join', [tj[e] for e in range(Egs) if tj[e]>0][:5])
