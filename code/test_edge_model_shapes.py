from gen_ilp import *
import random, time
def plus(n,c): return lambda p: (abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n)
def octa(n,d): return lambda p: abs(p[0])<=n and abs(p[1])<=n and abs(p[0]+p[1])<=d and abs(p[0]-p[1])<=d
def run(name,inside,mode,W=30,keep=None):
    S=Sample(inside,W=W,mode=mode)
    if keep is not None:      # clean cut, but keep the stem island of a chosen subset of freeable plaquettes
        extra=[stem(c) for c in keep(S)]
        S.isl=sorted(set(S.isl)|set(extra)); S.islset=set(S.isl); S.free={c for c in S.cells if stem(c) not in S.islset}
    t=time.time(); eff=S.lp(); allp=S.lp(free=set()); E,costly,st=S.ilp(tl=150)
    nb=sum(black(c) for c in S.cells)
    print(f'{name:38s} [{mode}] plaq {len(S.cells):3d} (H{nb}/V{len(S.cells)-nb}) freeable {len(S.free):2d} | all-paired {allp[0] if allp else "untileable"} | eff. model E={eff[0] if eff else None} (free {eff[1] if eff else None}) | ILP E={E} {costly} {st} {time.time()-t:.0f}s',flush=True)
# fringed termination (A): balanced tileable shapes must give E = plaq/2
run('plus n=7 c=3',plus(7,3),'A'); run('plus n=5 c=3',plus(5,3),'A'); run('plus n=9 c=1 (c<n/3)',plus(9,1),'A'); run('plus n=11 c=3 (c<n/3)',plus(11,3),'A')
run('octagon n=7 d=10',octa(7,10),'A'); run('octagon n=7 d=12',octa(7,12),'A')
# clean termination (B): effective boundary-monomer model
run('plus n=7 c=3',plus(7,3),'B'); run('octagon n=7 d=10',octa(7,10),'B')
run('tilted (2,1) square',halfplanes((2,1,0,18),(-1,2,0,18)),'B'); run('tilted (2,1) square',halfplanes((2,1,0,18),(-1,2,0,18)),'A')
run('tilted (3,1) square',halfplanes((3,1,1,23),(-1,3,1,23)),'B')
# programmable edge: clean cut with a random subset of stems kept
random.seed(3)
for rep in range(4):
    run(f'axis box 10x10 off0, random stems #{rep}',halfplanes((1,0,0,10),(0,1,0,10)),'B',keep=lambda S: [c for c in sorted(S.free) if random.random()<.5])
for rep in range(4):
    run(f'diamond D=6 (0,0), random stems #{rep}',halfplanes((1,1,0,12),(1,-1,0,12)),'B',keep=lambda S: [c for c in sorted(S.free) if random.random()<.5])
