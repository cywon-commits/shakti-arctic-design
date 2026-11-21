from gen_ilp import *
import time
def report(name,inside,modes=('A','B'),W=30):
    for m in modes:
        S=Sample(inside,W=W,mode=m); t=time.time()
        nb=sum(black(c) for c in S.cells)
        forced=S.lp(free=set()); eff=S.lp()
        E,costly,st=S.ilp()
        print(f'{name:34s} mode {m}: plaq {len(S.cells):3d} (H{nb}/V{len(S.cells)-nb}) freeable {len(S.free):2d} | all-paired dimers: {forced[0] if forced else "untileable"} | eff. model: dimers {eff[0] if eff else None}, free used {eff[1] if eff else None} | ILP E={E} {costly} [{st}, {time.time()-t:.0f}s]',flush=True)
# axis boxes
for a0 in (0,1):
    for X in (8,10):
        report(f'axis box off {a0} side {X}',halfplanes((1,0,a0,a0+X),(0,1,a0,a0+X)))
report('axis box off (0,1) 8x8',halfplanes((1,0,0,8),(0,1,1,9)))
# 45-degree diamonds  lo<=x+y<=lo+2D, lo2<=x-y<=lo2+2D
for a in range(4):
    for b in range(4):
        report(f'diamond a={a} b={b} D=6',halfplanes((1,1,a,a+12),(1,-1,b,b+12)))
