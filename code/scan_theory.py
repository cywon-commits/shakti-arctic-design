from theory import *
import warnings; warnings.filterwarnings('ignore')
def F(v,n,c):
    p,a0,a4=v
    if abs(p)>=1: return [1e3]*3
    t1=th_at(np.pi/2,p); t2=th_at(np.pi,p)
    return [H(0,p,a0,a4)-n, H(t1,p,a0,a4)-c, H(t2,p,a0,a4)+c]
v=np.array(solve(1,.5)[0]); res={}
for c in np.arange(.5,1.0,.01):
    v,info,ier,msg=fsolve(F,v,args=(1,c),full_output=True)
    if ier!=1 or abs(v[0])>=1: print('up stops at',c); break
    res[round(c,3)]=v.copy()
v=np.array(solve(1,.5)[0])
for c in np.arange(.5,0.0,-.005):
    v,info,ier,msg=fsolve(F,v,args=(1,c),full_output=True)
    if ier!=1 or abs(v[0])>=1: print('down stops at',c); break
    res[round(c,3)]=v.copy()
for c in sorted(res):
    if abs(c*100%5)<1e-6 or c<.34 or c>.93:
        p,a0,a4=res[c]; print(f'c/n={c:.3f}  p={p:+.5f}  a0={a0:+.4f}  a4={a4:+.4f}   (a0-a4)/(1+p)^2 = H(pi/4) = {(a0-a4)/(1+p)**2:+.4f}')
