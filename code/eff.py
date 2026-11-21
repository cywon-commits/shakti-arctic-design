"""Effective boundary-monomer model at larger sizes (clean termination, mode B): which freeable plaquettes are free in the ground state, and what region is left."""
from gen_ilp import *
import math
def ntil(S):
    B=[c for c in S if black(c)]; Wt=[c for c in S if not black(c)]
    if len(B)!=len(Wt): return None
    if not B: return 0.0
    bi={c:k for k,c in enumerate(B)}; wi={c:k for k,c in enumerate(Wt)}
    K=np.zeros((len(B),len(B)),dtype=complex)
    for c in B:
        for d,w in [((1,0),1),((-1,0),1),((0,1),1j),((0,-1),1j)]:
            n=(c[0]+d[0],c[1]+d[1])
            if n in wi: K[bi[c],wi[n]]=w
    s,ld=np.linalg.slogdet(K)
    return ld if abs(s)>.5 else None
def analyse(name,inside,W=45,mode='B',verbose=True):
    S=Sample(inside,W=W,mode=mode); F=S.free
    base=S.lp()
    if base is None: print(name,'infeasible'); return S,None
    nd,nf,dim=base
    always=set(); never=set()
    for f in F:
        r=S.lp(free=F-{f})
        if r is None or r[0]>nd: always.add(f)
    matched={c for d in dim for c in d}
    uniq = len(always)==nf
    out=f'{name:30s} plaq {len(S.cells):4d} freeable {len(F):3d} (H{sum(black(c) for c in F)}/V{len(F)-sum(black(c) for c in F)}) free in g.s. {nf:3d}, always-free {len(always):3d}  E_gs={nd}  E_allpaired={(S.lp(free=set()) or ["untileable"])[0]}'
    if uniq:
        Se=S.cells-always; ld=ntil(Se)
        out+=f' | unique free set; ln(#g.s.)={ld:.2f}, per cell {ld/len(Se):.4f}' if ld is not None else ' | S_eff untileable?!'
    else: out+=' | free set NOT unique'
    print(out,flush=True)
    return S,(always if uniq else None)
if __name__=='__main__':
    for X in (16,24,32):
        analyse(f'axis block-aligned {X}',halfplanes((1,0,0,X),(0,1,0,X)))
    for X in (16,18,24,26):
        analyse(f'axis through-block {X}',halfplanes((1,0,1,1+X),(0,1,1,1+X)))
    for a in range(4):
        for b in range(4):
            analyse(f'diamond a={a} b={b} D=12',halfplanes((1,1,a,a+24),(1,-1,b,b+24)))
    for D in (16,20):
        analyse(f'diamond a=0 b=0 D={D}',halfplanes((1,1,0,2*D),(1,-1,0,2*D)))
        analyse(f'diamond a=1 b=1 D={D}',halfplanes((1,1,1,1+2*D),(1,-1,1,1+2*D)))
    for k in range(4):
        analyse(f'tilted (2,1) square off {k}',halfplanes((2,1,k,k+50),(-1,2,k,k+50)))
