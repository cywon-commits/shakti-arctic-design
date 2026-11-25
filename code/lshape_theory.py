"""Tangent-line family for the L-shaped sample: square |X|,|Y| <= n with the corner X > a, Y > b removed.
One reflex corner => degree d = 2 (8 axis tangents: 6 sides + 2 hidden tangents on the extensions X = a, Y = b).
    exp(i nu) = B(u) = e^{i g} u (u - q)/(1 - conj(q) u),   u = e^{i theta}
    H(theta)  = (t0 + t1 cos + t2 sin + t3 cos 2theta + t4 sin 2theta) / |1 - conj(q) u|^2
The rotation of u is fixed by putting nu = 0 at theta = 0 (g is then determined by q).
8 conditions on 7 unknowns (q, t0..t4); they are compatible because the boundary height closes."""
import numpy as np
from scipy.optimize import least_squares, brentq
def nu(th,q):
    u=np.exp(1j*th); B=u*(u-q)/(1-np.conj(q)*u); B0=(1-q)/(1-np.conj(q))
    return np.unwrap(np.angle(B/B0)) if np.ndim(th) else np.angle(B/B0)
def nu_cont(q,m=20001):
    th=np.linspace(0,2*np.pi,m); return th,nu(th,q)
def H(th,q,t):
    u=np.exp(1j*th); return (t[0]+t[1]*np.cos(th)+t[2]*np.sin(th)+t[3]*np.cos(2*th)+t[4]*np.sin(2*th))/abs(1-np.conj(q)*u)**2
def targets(n,a,b): return [n,b,-a,-b,a,n,n,n]      # H at nu = 0, pi/2, ..., 7 pi/2
def thetas(q):
    th,v=nu_cont(q); return np.interp(np.arange(8)*np.pi/2,v,th)
def resid(x,n,a,b):
    r,phi=np.tanh(x[0]),x[1]; q=r*np.exp(1j*phi); return H(thetas(q),q,x[2:])-np.array(targets(n,a,b))
def solve(n,a,b,seed=0):
    rng=np.random.default_rng(seed); best=None
    for k in range(60):
        x0=np.r_[rng.normal(0,1),rng.uniform(0,2*np.pi),rng.normal(0,n,5)]
        s=least_squares(resid,x0,args=(n,a,b),xtol=1e-14,ftol=1e-14,gtol=1e-14)
        if best is None or s.cost<best.cost: best=s
        if best.cost<1e-18: break
    x=best.x; return np.tanh(x[0])*np.exp(1j*x[1]),x[2:],np.sqrt(2*best.cost)
def curve(q,t,m=40001):
    th,v=nu_cont(q,m); h=H(th,q,t); dh=np.gradient(h,th); dv=np.gradient(v,th)
    return th,v,h*np.cos(v)-dh/dv*np.sin(v),h*np.sin(v)+dh/dv*np.cos(v)
if __name__=='__main__':
    for a in (0.2,0.4,0.6,0.7,0.75):
        q,t,res=solve(1,a,a); th,v,x,y=curve(q,t)
        print(f'a=b={a}: q={q:.6f} t={np.round(t,6)} residual={res:.1e}  curve bbox X[{x.min():.4f},{x.max():.4f}] Y[{y.min():.4f},{y.max():.4f}]')
    q,t,res=solve(1,0.3,0.6); print('a=0.3,b=0.6:',q,np.round(t,6),f'{res:.1e}')
