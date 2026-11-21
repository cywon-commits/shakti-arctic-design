from maps import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, sys, pickle
def region(inside,W): return Sample(inside,W=W,mode='B').cells
def plus(n,c): return lambda p: (abs(p[0])<=n and abs(p[1])<=c) or (abs(p[0])<=c and abs(p[1])<=n)
def typemap(S):
    P,ld=edge_probs(S); best={}
    for (c,d),p in P.items():
        n=(c[0]+d[0],c[1]+d[1])
        for x in (c,n):
            if p>best.get(x,(0,None))[0]: best[x]=(p,d)
    return best,ld
if __name__=='__main__':
    n=int(sys.argv[1]); cs=[int(a) for a in sys.argv[2:]]
    fig,ax=plt.subplots(1,len(cs),figsize=(9*len(cs),9)); ax=np.atleast_1d(ax)
    col={(1,0):'tab:red',(-1,0):'tab:blue',(0,1):'tab:green',(0,-1):'tab:orange'}
    out={}
    for a,c in zip(ax,cs):
        S=region(plus(n,c),W=n+6); best,ld=typemap(S)
        xs=[pcenter(q)[0] for q in S]; ys=[pcenter(q)[1] for q in S]
        cl=[col[best[q][1]] if best[q][0]>0.97 else '0.85' for q in S]
        a.scatter(xs,ys,c=cl,s=(420/n)**2,marker='s'); a.set_aspect('equal'); a.set_title(f'plus n={n} c={c}: frozen cells coloured by domino type (P>0.97), liquid grey; lnZ/cell={ld/len(S):.4f}')
        out[c]={q:(pcenter(q),best[q]) for q in S}
        print(n,c,len(S),ld,flush=True)
    plt.tight_layout(); plt.savefig(f'../figs/plus_types_n{n}.png',dpi=60)
    pickle.dump(out,open(f'data/plus_n{n}.pkl','wb'))
