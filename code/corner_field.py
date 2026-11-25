"""Field (temperature) of the corner-switch solution from flux conservation.
The corner S is a source of the current j = d(ent)/d(grad H) of strength x ln 2 = (J/2T)/m (vacancy-paper normalisation,
j_u = (1/2) ln|z|, j_v = -(1/2) ln|w|).  Half of it crosses any curve from the partially free zone of the edge to the diagonal."""
import sys, numpy as np
from scipy.integrate import quad
from corner_theory import Corner
def Kc(c,p):          # K at complex p (upper half-plane)
    f=lambda a,part: part(c.kappa(c.b*np.sin(a))*c.b*np.cos(a)/(c.b*np.sin(a)-p))
    re=quad(f,-np.pi/2,np.pi/2,args=(np.real,),limit=200)[0]; im=quad(f,-np.pi/2,np.pi/2,args=(np.imag,),limit=200)[0]
    return c.C+(re+1j*im)/np.pi
def point(c,zeta):
    t=zeta**2; q=(t+1/t)/2; p=np.sqrt(q*q+c.b**2+0j); p=p if p.imag>=0 else -p          # upper half-plane
    h=1j*np.sqrt(p-c.p1+0j)*Kc(c,p)       # sqrt branch: positive for p > p1, continued through the upper half-plane
    mu=(1-zeta**2)/(1+zeta**2); lam=1j*mu; Psi=-1j*mu-h/(zeta+1/zeta)
    u=-Psi.imag/lam.imag; v=Psi.real+lam.real*u; z=1j*(1+zeta)/(1-zeta); w=(1+z)/(z-1); return u,v,np.log(abs(z)),np.log(abs(w))
def half_flux(c,rho0,n=200,sign=-1):
    th=np.linspace(0,-np.pi/4,n+2)[1:-1]; P=np.array([point(c,rho0*np.exp(1j*a)) for a in th]); u,v,lz,lw=P.T
    m=lambda a:(a[1:]+a[:-1])/2; return np.sum(m(lz)*np.diff(v))/2, np.sum(m(lw)*np.diff(u))/2, (u[0],v[0],u[-1],v[-1])
if __name__=='__main__':
    for b in [float(v) for v in sys.argv[1:]] or [0.5,2]:
        c=Corner(b); rD=(np.sqrt(1+b*b)-b)**0.5; print('b',b,'r_D',round(rD,4))
        for rho0 in (rD+0.25*(1-rD),rD+0.5*(1-rD),rD+0.75*(1-rD)):
            a,bb,e=half_flux(c,rho0); print(f'   rho0={rho0:.3f}: int ln|z| dv/2 = {a:+.5f}, int ln|w| du/2 = {bb:+.5f}; a-b = {a-bb:+.5f}, a+b = {a+bb:+.5f}; ends (u,v) = ({e[0]:.3f},{e[1]:.3f}) -> ({e[2]:.3f},{e[3]:.3f})')
def field(c,n=300,d=2e-3):
    """(J/2T)/m = x ln 2.  The path stops a small angle d short of the diagonal (the kernel is singular there); the error is O(d)."""
    rD=(c.p1-c.b)**0.5; rho0=(rD+1)/2; th=np.linspace(-1e-9,-np.pi/4+d,n)
    P=np.array([point(c,rho0*np.exp(1j*a)) for a in th]); u,v,lz,lw=P.T; m=lambda a:(a[1:]+a[:-1])/2
    # close the path to the diagonal along the last direction
    tot=np.sum(m(lz)*np.diff(v))/2-np.sum(m(lw)*np.diff(u))/2
    du=(v[-1]-u[-1])/2; tot+=(lz[-1]*(-du)-lw[-1]*du)/2
    return 2*tot
