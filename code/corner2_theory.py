"""Scaling solution of the corrected corner switch (boundary energy ~ h(NE) - h(SW), natural boundary conditions on all four edges).
Unit square, S = (0,0) the SW corner, half-domain below the diagonal.  zeta-plane domain: quarter disc minus the slit arg = -pi/4, r_D < |zeta| < 1
(same as the uncorrected problem); now the imaginary axis is the E edge (u = 1) with natural boundary condition.
    v - i mu (u - 1) = -h/(zeta + 1/zeta),   h = i sqrt(p^2 - p1^2) K(p),
    K(p) = (1/pi) int_{-b}^{b} kappa(p') dp'/(p' - p),   kappa = sqrt(sqrt(p1^2 - p'^2) - 1)/sqrt(p1^2 - p'^2),
p = sqrt(q^2 + b^2) (sign of Re q), q = (zeta^2 + zeta^-2)/2, p1 = sqrt(1 + b^2).  K is odd (point symmetry), no free constant.
Frozen boundaries: lines (u - 1) sin(phi) + v cos(phi) = H(phi) = sqrt(p1^2 - p^2) K(p)/2,  q = cos(2 phi)."""
import numpy as np
from scipy.integrate import quad
class Corner2:
    def __init__(s,b): s.b=b; s.p1=np.sqrt(1+b*b)
    def kappa(s,x): r=np.sqrt(max(s.p1**2-x*x,0)); return np.sqrt(max(r-1,0))/r
    def K(s,p):
        f=lambda a,part: part(s.kappa(s.b*np.sin(a))*s.b*np.cos(a)/(s.b*np.sin(a)-p))
        if np.iscomplexobj(p):
            return (quad(f,-np.pi/2,np.pi/2,args=(np.real,),limit=200)[0]+1j*quad(f,-np.pi/2,np.pi/2,args=(np.imag,),limit=200)[0])/np.pi
        return quad(f,-np.pi/2,np.pi/2,args=(lambda x:x,),limit=200,epsabs=1e-13)[0]/np.pi
    def H(s,phi):
        q=np.cos(2*phi); p=np.sign(q)*np.sqrt(q*q+s.b**2); return np.sqrt(max(s.p1**2-p*p,0))*s.K(p)/2
    def curve(s,lo,hi,n=400):
        ph=lo+(hi-lo)*(1-np.cos(np.linspace(0,np.pi,n)))/2; ph=ph[1:-1]; H=np.array([s.H(x) for x in ph]); dH=np.gradient(H,ph,edge_order=2)
        return 1+H*np.sin(ph)+dH*np.cos(ph), H*np.cos(ph)-dH*np.sin(ph)
    def caps(s,n=400):
        return s.curve(1e-4,np.pi/4-1e-4,n), s.curve(np.pi/4+1e-4,np.pi/2-1e-4,n)
    def edge_S(s,r):     # S edge, zeta = r in (0,1): u and tangential slope
        q=(r*r+1/(r*r))/2; p=np.sqrt(q*q+s.b**2); h=-np.sqrt(p*p-s.p1**2)*s.K(p)      # i sqrt(.) K with K real -> h imaginary: h = i*X; use u = 1 + X/(...)
        # v - i mu (u-1) = -h/(zeta+1/zeta), v = 0, mu real:  mu (u-1) = X/(zeta+1/zeta) with h = i X
        X=np.sqrt(p*p-s.p1**2)*s.K(p); mu=(1-r*r)/(1+r*r); return 1+X/((r+1/r)*mu), 4/np.pi*np.arctan(r)
    def point(s,zeta):
        t=zeta**2; q=(t+1/t)/2; p=np.sqrt(q*q+s.b**2+0j); p=p if p.imag>=0 else -p
        h=1j*np.sqrt(p*p-s.p1**2+0j)*s.K(p)
        if (np.sqrt(p*p-s.p1**2+0j)).imag<0: h=-h      # branch continued through the upper half-plane from p > p1
        mu=(1-zeta**2)/(1+zeta**2); lam=1j*mu; G=-h/(zeta+1/zeta)          # G = v - i mu (u - 1)
        u=1-G.imag/mu.real if abs(mu.imag)<1e-14 else None
        # general: v - lam u = G - i mu  ->  Psi = G - i mu ; u = -Im(Psi)/Im(lam), v = Re(Psi) + Re(lam) u
        Psi=G-1j*mu; u=-Psi.imag/lam.imag; v=Psi.real+lam.real*u
        z=1j*(1+zeta)/(1-zeta); w=(1+z)/(z-1); return u,v,np.log(abs(z)),np.log(abs(w))
def area(xy):
    x,y=xy; return 0.5*abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
if __name__=='__main__':
    for b in (0.2,0.5,1,2,5):
        c=Corner2(b); (u1,v1),(u2,v2)=c.caps()
        print(f'b={b}: SW cap curve from ({u1[0]:.3f},{v1[0]:.3f}) to ({u1[-1]:.3f},{v1[-1]:.3f});  NE cap curve from ({u2[0]:.3f},{v2[0]:.3f}) to ({u2[-1]:.3f},{v2[-1]:.3f})')
        print('   S edge: (u, slope) at r = 0.99, 0.5, 0.01:',[tuple(np.round(c.edge_S(r),3)) for r in (0.99,0.5,0.01)])
def field(c,n=300,d=2e-3):
    """(J/2T)/m: twice the flux of the SW source through a path from the SE corner (zeta = 0) to the centre (slit tip)"""
    rD=np.sqrt(c.p1-c.b) if False else ((np.sqrt(1+c.b**2)-c.b))**0.5
    tip=rD*np.exp(-1j*np.pi/4)
    ts=np.linspace(1e-4,1-d,n); P=np.array([c.point(t*tip) for t in ts]); u,v,lz,lw=P.T; m=lambda a:(a[1:]+a[:-1])/2
    tot=np.sum(m(lz)*np.diff(v))/2-np.sum(m(lw)*np.diff(u))/2
    return 2*tot, (u[0],v[0],u[-1],v[-1])
def liquid_fraction(c,n=400):
    (u1,v1),(u2,v2)=c.caps(n)
    poly=(np.r_[0,u1,0],np.r_[0,v1,0]); return 1-4*area(poly)
def field2(c,n=400,eps=1e-4):
    """same flux, path from the SE corner to a point on the SW half of the diagonal (approached from the S-edge side of the slit)"""
    rD=((np.sqrt(1+c.b**2)-c.b))**0.5; end=0.5*(rD+1)*np.exp(-1j*(np.pi/4-eps))
    ts=np.linspace(1e-5,1,n)**1.0; P=np.array([c.point(t*end) for t in ts]); u,v,lz,lw=P.T; m=lambda a:(a[1:]+a[:-1])/2
    tot=np.sum(m(lz)*np.diff(v))/2-np.sum(m(lw)*np.diff(u))/2
    du=(v[-1]-u[-1])/2; tot+=(lz[-1]*du-lw[-1]*(-du))/2*0     # closing segment is O(eps); ignored
    return 2*tot,(u[-1],v[-1])
def field_closed(c):
    """(J/2T)/m = ln(1+sqrt2) + (1/2) int_0^1 K(sqrt(q^2+b^2)) dq/sqrt(1-q^2)   (flux written as -1/2 int Psi dz/z along the S edge and the SW cap boundary)"""
    f=lambda th: c.K(np.sqrt(np.cos(th)**2+c.b**2))
    return np.log(1+np.sqrt(2))+0.5*quad(f,0,np.pi/2,limit=200)[0]
def X1(c): return 1+c.K(c.p1)
