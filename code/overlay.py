import pickle, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from theory import *
n=int(sys.argv[1]); D=pickle.load(open(f'data/plus_n{n}.pkl','rb'))
col={(1,0):'tab:red',(-1,0):'tab:blue',(0,1):'tab:green',(0,-1):'tab:orange'}
fig,ax=plt.subplots(1,len(D),figsize=(10*len(D),10)); ax=np.atleast_1d(ax)
for a,(c,M) in zip(ax,sorted(D.items())):
    xs=[v[0][0] for v in M.values()]; ys=[v[0][1] for v in M.values()]
    cl=[col[v[1][1]] if v[1][0]>0.97 else '0.88' for v in M.values()]
    a.scatter(xs,ys,c=cl,s=(460/n)**2,marker='s')
    sol=solve(n,c)
    for s in (sol or []):
        th,v,x,y=curve(*s); a.plot(x,y,'k-',lw=1.2)
        print(n,c,'p,a0,a4=',s)
    a.set_aspect('equal'); a.set_xlim(-n-3,n+3); a.set_ylim(-n-3,n+3); a.set_title(f'n={n}, c={c}')
plt.tight_layout(); plt.savefig(f'../figs/overlay_n{n}.png',dpi=55)
