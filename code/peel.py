from maps import *
def peel(S):
    """remove forced dominoes (a cell with exactly one neighbour left) until none; returns reduced region, forced dominoes, ok"""
    S=set(S); forced=[]
    nb=lambda c:[n for n in [(c[0]+1,c[1]),(c[0]-1,c[1]),(c[0],c[1]+1),(c[0],c[1]-1)] if n in S]
    Q=[c for c in S if len(nb(c))<=1]
    while Q:
        c=Q.pop()
        if c not in S: continue
        N=nb(c)
        if len(N)==0: return S,forced,False
        if len(N)>1: continue
        m=N[0]; S-={c,m}; forced.append((c,m))
        for x in nb(m)+nb(c):
            if len(nb(x))<=1: Q.append(x)
        for d in [(1,0),(-1,0),(0,1),(0,-1)]:
            for base in (c,m):
                x=(base[0]+d[0],base[1]+d[1])
                if x in S and len(nb(x))<=1: Q.append(x)
    return S,forced,True
def bbox(S):
    pts=[p for c in S for p in ppts(c)]
    return min(p[0] for p in pts),max(p[0] for p in pts),min(p[1] for p in pts),max(p[1] for p in pts)
def region(inside,W=80): return Sample(inside,W=W,mode='B').cells
if __name__=='__main__':
    import itertools
    n=21
    print('general octagons: rectangle |X|<=nx,|Y|<=ny (odd), cuts X+Y<=dNE, -X+Y<=dNW, -X-Y<=dSW, X-Y<=dSE (even)')
    seen=0
    for nx,ny in [(21,21),(21,17),(21,13),(25,17)]:
        for dNE,dNW,dSW,dSE in itertools.product(range(nx+ny-16,nx+ny+1,4),repeat=4):
            ins=lambda p: abs(p[0])<=nx and abs(p[1])<=ny and p[0]+p[1]<=dNE and -p[0]+p[1]<=dNW and -p[0]-p[1]<=dSW and p[0]-p[1]<=dSE
            S=region(ins); nb_=sum(black(c) for c in S)
            if 2*nb_!=len(S): continue
            Sr,forced,ok=peel(S)
            if not ok: continue
            ld=ntil(Sr) if Sr else 0.0
            if ld is None: continue
            if (dNE,dNW,dSW,dSE)==(dNE,)*4 and nx==ny and seen>3: pass
            b=bbox(Sr) if Sr else None
            area=len(Sr)
            rect_cells=None
            if Sr:
                R_=region(lambda p,b=b: b[0]<=p[0]<=b[1] and b[2]<=p[1]<=b[3]); rect_cells=len(R_)
            print(f' nx,ny=({nx},{ny}) d=({dNE},{dNW},{dSW},{dSE}): cells {len(S)} forced {2*len(forced)} reduced {area}  bbox {b}  full-rect cells {rect_cells}  lnZ={ld/math.log(2):.3f} ln2',flush=True)
