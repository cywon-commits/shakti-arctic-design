"""Degree-5 (class-10) arctic curve of the plus-shaped Shakti sample: arms reach n, half-width c (Shakti units, all edges of face-extremal type).
Tangent-line family  X cos(nu) + Y sin(nu) = H(theta):
    nu(theta) = 5 theta + 2 atan2(p sin 4theta, 1 - p cos 4theta)          (degree-5 Blaschke product, D4-symmetric)
    H(theta)  = (a0 + a4 cos 4theta) / (1 + p^2 - 2 p cos 4theta)
Conditions: H=n at nu=0 (arm end), H=c at nu=pi/2 (arm side), H=-c at nu=pi (ghost extension of the neighbouring arm's side)."""
import numpy as np
from scipy.optimize import brentq, fsolve
def nu(th,p): return 5*th+2*np.arctan2(p*np.sin(4*th),1-p*np.cos(4*th))
def H(th,p,a0,a4): return (a0+a4*np.cos(4*th))/(1+p*p-2*p*np.cos(4*th))
def th_at(target,p): return brentq(lambda t: nu(t,p)-target,0,np.pi/4)
def solve(n,c):
    def F(v):
        p,a0,a4=v; p=np.tanh(p)
        t1=th_at(np.pi/2,p); t2=th_at(np.pi,p)
        return [H(0,p,a0,a4)-n, H(t1,p,a0,a4)-c, H(t2,p,a0,a4)+c]
    best=None
    for p0 in np.linspace(-2.5,2.5,21):
        for s in (1,-1):
            try:
                v,info,ier,msg=fsolve(F,[p0,n*s,0.3*n],full_output=True)
            except Exception: continue
            if ier==1 and max(abs(np.array(F(v))))<1e-9:
                sol=(np.tanh(v[0]),v[1],v[2])
                if best is None or not any(np.allclose(sol,b,atol=1e-6) for b in best): best=(best or [])+[sol]
    return best
def curve(p,a0,a4,m=40001):
    th=np.linspace(0,2*np.pi,m); v=nu(th,p); h=H(th,p,a0,a4)
    dh=np.gradient(h,th); dv=np.gradient(v,th)
    x=h*np.cos(v)-dh/dv*np.sin(v); y=h*np.sin(v)+dh/dv*np.cos(v)
    return th,v,x,y
def closed_form(n,c):
    """Explicit parameters for n/3 < c < n/sqrt(2):  sigma = 2 sin(theta_1) = (2 rho - 1 + sqrt(10 rho^2 - 1))/(1 + rho), rho = c/n."""
    rho=c/n; s=(2*rho-1+np.sqrt(10*rho*rho-1))/(1+rho); q=s*s-4*s-6
    p=(1-s-s*s)/(1+s)
    a0=(s**6+2*s**5-10*s**4-28*s**3-14*s**2+8*s+4)/((s+1)**2*q)
    a4=-2*(s**4+3*s**3+3*s**2+2*s+2)/((s+1)*q)
    return p,n*a0,n*a4
if __name__=='__main__':
    for n,c in [(1,.2),(1,.3443),(1,.4),(1,.5),(1,.6),(1,.8),(1,.95)]:
        print(n,c,solve(n,c))
