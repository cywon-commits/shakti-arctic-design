"""Continuum solution of the corner-switch sample (hodograph method of the vacancy paper, new boundary data).
Unit square, S = (0,0) the corner where the two switchable edges meet, A = (1,0).  Half domain below the diagonal.
zeta-plane image: quarter disc {|zeta|<1, -pi/2 < arg < 0} slit along arg = -pi/4, r_D < |zeta| < 1.
Uniformisation: t = zeta^2, q = (t + 1/t)/2, p = sqrt(q^2 + b^2) (sign of Re q), b = (r_D^-2 - r_D^2)/2, p1 = sqrt(1 + b^2).
    v - i mu (u - 1) = - h(p)/(zeta + 1/zeta),   h = i sqrt(p - p1) K(p),
    K(p) = C + (1/pi) int_{-b}^{b} kappa(p')/(p' - p) dp',  kappa = sqrt( sqrt(1 + b^2 - p'^2) - 1 ) / sqrt(p1 - p'),  K(-p1) = 0.
One parameter b: b -> infinity is zero field (quarter circle), b -> 0 is complete freezing."""
import numpy as np
from scipy.integrate import quad
class Corner:
    def __init__(s,b):
        s.b=b; s.p1=np.sqrt(1+b*b)
        s.C=-s._I(-s.p1)
    def kappa(s,x): return np.sqrt(max(np.sqrt(1+s.b**2-x*x)-1,0.0))/np.sqrt(s.p1-x)
    def _I(s,p):          # (1/pi) int kappa/(p'-p), p real outside (-b,b); substitution p' = b sin a
        f=lambda a: s.kappa(s.b*np.sin(a))*s.b*np.cos(a)/(s.b*np.sin(a)-p)
        return quad(f,-np.pi/2,np.pi/2,limit=200,epsabs=1e-13,epsrel=1e-12)[0]/np.pi
    def K(s,p): return s.C+s._I(p)
    # --- edge v = 0
    def p_of_r(s,r): q=(r*r+1/(r*r))/2; return np.sqrt(q*q+s.b**2)
    def u_ii(s,r): p=s.p_of_r(r); return 1+np.sqrt(p-s.p1)*s.K(p)/(1/r-r)          # partially free zone, zeta = r, density (4/pi) arctan r
    def u_iii(s,t): p=-s.p_of_r(t); return 1+np.sqrt(s.p1-p)*s.K(p)/(t+1/t)        # empty zone, zeta = -i t
    def X1(s): return 1+s.K(s.p1)/np.sqrt(2*s.p1)
    def X2(s): return 1+s.C/np.sqrt(2)
    def free_fraction(s,n=400):
        r=np.cos(np.linspace(0,np.pi/2,n+2)[1:-1]); u=np.array([s.u_ii(x) for x in r]); rho=4/np.pi*np.arctan(r)
        u=np.r_[s.X1(),u,s.X2()]; rho=np.r_[1.0,rho,0.0]; return s.X1()+np.trapezoid(rho,u), u, rho
    # --- frozen boundary: lines (u-1) sin(phi) + v cos(phi) = HA(phi)
    def HA(s,phi):
        q=np.cos(2*phi); p=np.sign(q)*np.sqrt(q*q+s.b**2); return np.sqrt(s.p1-p)*s.K(p)/2
    def curve(s,n=200,eps=1e-3):
        out=[]
        for lo,hi in ((eps,np.pi/4-eps),(np.pi/4+eps,np.pi/2-eps)):
            ph=np.linspace(lo,hi,n); H=np.array([s.HA(x) for x in ph]); dH=np.gradient(H,ph)
            out.append((1+H*np.sin(ph)+dH*np.cos(ph), H*np.cos(ph)-dH*np.sin(ph)))
        return out
if __name__=='__main__':
    print('    b      X1      X2    free fraction   HA(outer, phi=3pi/8) [quarter circle: 1-sin = %.4f]'%(1-np.sin(3*np.pi/8)))
    for b in (0.05,0.2,0.5,1,2,5,20,100):
        c=Corner(b); f,u,rho=c.free_fraction(); print(f'{b:7.2f} {c.X1():7.4f} {c.X2():7.4f} {f:9.4f}      {c.HA(3*np.pi/8):.4f}')
