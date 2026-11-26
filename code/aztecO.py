"""Through-block (fringed) Aztec square at finite temperature versus the boundary-monomer Aztec diamond of the vacancy paper.
Shakti sample: all four edges through-block; half-dimer sites on every other boundary cell, weight 2, energy J; Aztec order N = L/2 - 1.
Vacancy model: monomers on all 4N boundary cells with fugacity z = 2^(-xN/2).  Matching: 2 exp(-J/2T) = z."""
import sys, math, numpy as np
sys.path.insert(0,'/home/claude/screened-dimer-electrostatics/src')
from gen_ilp import Sample, halfplanes
import corrected_model as cm
from vertex_dos import dos
def shakti(L):
    S=Sample(halfplanes((1,0,1,L+1),(0,1,1,L+1)),W=L+10,mode='A'); n=cm.halfsites(S,'fit'); c=cm.poly(S,n)
    E0=next(i for i,v in enumerate(c) if v); return S,c,E0
def sigma_shakti(c,E0,N,x):
    lt=-(x*N*math.log(2)+2*math.log(2))          # ln t = -J/T with 2 e^{-J/2T} = 2^{-xN/2}
    E=np.arange(len(c)); lw=np.array([math.log(v) if v>0 else -np.inf for v in c])+lt*E; lw-=lw.max(); w=np.exp(lw)
    h=2*((w*(E-E0)).sum()/w.sum()); return 1-h/(2*N)
if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='vertex':
        for L in [int(v) for v in sys.argv[2:]]:
            S,c,E0=shakti(L); g,_,_=dos(S,E0+3)
            print(f'L={L} (N={L//2-1}): '+'  '.join(f'E0+{k}: vertex {g[E0+k]/2:.0f} model {c[E0+k] if E0+k<len(c) else 0}' for k in range(4)))
    else:
        import bpf_exact as bp
        xs=[float(v) for v in sys.argv[2].split(',')]
        for L in [int(v) for v in sys.argv[3:]]:
            N=L//2-1; S,c,E0=shakti(L); R=bp.Rmat(N); Gt=bp.build(N,R)
            for x in xs:
                sv,_=bp.sigma(N,x,Gt,int(6*N)+400); print(f'N={N:3d} x={x:.2f}: sigma Shakti {sigma_shakti(c,E0,N,x):.5f}   vacancy model {sv:.5f}',flush=True)
