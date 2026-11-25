"""Width of the frozen-side tail at the class-10 curve versus size.  Outside the curve the liquidity 1 - P_max decays as
exp[-(4/3)(d/w)^(3/2)] (Airy tail), so (-ln liquidity)^(2/3) is linear in the distance d with slope (4/3)^(2/3)/w.
Expected: w ~ n^(1/3).  Uses data/verify_shuffle*.pkl."""
import pickle, glob, numpy as np
rows=[]
for f in sorted(glob.glob('data/verify_shuffle*.pkl')):
    for (n,c),d in pickle.load(open(f,'rb')).items():
        cen,liq,wn,dist=d['cen'],d['liq'],d['wn'],d['dist']
        sel=(wn==0)&(np.maximum(abs(cen[:,0]),abs(cen[:,1]))>c+0.1*n)&(liq>1e-12)&(liq<1e-2)
        y=(-np.log(liq[sel]))**(2/3); s,b=np.polyfit(dist[sel],y,1); w=(4/3)**(2/3)/s
        res=np.std(y-(s*dist[sel]+b))
        rows.append((n,c,int(sel.sum()),w,-b/s,res)); print(f'n={n} c={c}: {int(sel.sum())} cells in the tail, w = {w:.3f}, offset d0 = {-b/s:+.2f}, rms residual of (-ln)^(2/3): {res:.3f}')
rows=sorted(set(rows)); n=np.array([r[0] for r in rows],float); w=np.array([r[3] for r in rows])
print('w / n^(1/3):',np.round(w/n**(1/3),4))
if len(n)>=3: print('fitted exponent d ln w / d ln n = %.3f (all sizes), %.3f (largest three)'%(np.polyfit(np.log(n),np.log(w),1)[0],np.polyfit(np.log(n[-3:]),np.log(w[-3:]),1)[0]))
