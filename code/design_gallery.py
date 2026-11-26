"""Gallery of designed samples built from overlapping squares with through-block (odd) edges: balance, tileability, liquidity map."""
import math, sys, pickle, numpy as np
from gen_ilp import Sample, black
from eff import ntil
from maps import cell_entropy, pcenter
def squares_inside(centres,s):
    return lambda p: any(abs(p[0]-cx)<=s and abs(p[1]-cy)<=s for cx,cy in centres)
def chain(path,o):
    """path: string of 'u' (up-right) / 'd' (down-right) steps; centres on the diagonal lattice with step o"""
    c=[(0,0)]
    for ch in path:
        x,y=c[-1]; c.append((x+o,y+o) if ch=='u' else (x+o,y-o))
    return c
SHAPES={}
def add(name,centres,s,extra=None):
    SHAPES[name]=(centres,s,extra)
def run(names,compute=True):
    out={}
    for nm in names:
        centres,s,extra=SHAPES[nm]; ins=squares_inside(centres,s) if extra is None else extra
        W=max(max(abs(cx),abs(cy)) for cx,cy in centres)+s+4
        S=Sample(ins,W=W,mode='A'); nb=sum(black(c) for c in S.cells); imb=2*nb-len(S.cells)
        lt=ntil(S.cells) if imb==0 and len(S.free)==0 else None
        print(f'{nm}: cells {len(S.cells)}, freeable {len(S.free)}, H-V {imb}, ln#/ln2 {None if lt is None else round(lt/math.log(2),1)}',flush=True)
        if compute and lt is not None:
            liq,ld=cell_entropy(S.cells); P=np.array([pcenter(c) for c in S.cells]); v=np.array([liq[c] for c in S.cells])
            out[nm]=dict(cen=P,liq=v,centres=centres,s=s,cells=len(S.cells))
    return out
o=14; s=11
add('diag3',chain('uu',o),s); add('V',chain('ddu u'.replace(' ',''),o),s); add('zigzag',chain('udud',o),s)
add('X',[(0,0),(o,o),(-o,-o),(o,-o),(-o,o)],s); add('W',chain('dudu',o),s)
add('chevron',chain('uudd',o),s); add('check',chain('duu',o),s)
if __name__=='__main__':
    names=sys.argv[1].split(',') if len(sys.argv)>1 else list(SHAPES)
    out=run(names,compute='nocompute' not in sys.argv)
    if out: pickle.dump(out,open('data/design_gallery.pkl','wb'))
