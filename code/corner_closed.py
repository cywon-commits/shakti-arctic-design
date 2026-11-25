"""Closed form of the corner-switch solution (a = X1 = 1/p1, 0 < a < 1):
    Psi(zeta) = -i sqrt(2 (p - p1)/p1) / (zeta + 1/zeta)
    edge zones: X1 = a, X2 = sqrt(a);  free density rho(u) = (4/pi) arctan r,  (1/r - r)^2 = 4 (u^2 - a^2)/(a^2 - u^4)
    frozen boundary: support function about S, normal (sin phi, cos phi):  h = sin(chi) (inner), cos(chi) (outer),  sin(2 chi) = a sin(2 phi)
    field (J/2T)/m: flux integral along the edge and the outer curve (phi_field)."""
import mpmath as mp
def phi_field(a):
    a=mp.mpf(a); p1=1/a; b=mp.sqrt(p1*p1-1)
    def u3(s):
        q=-(s*s+1/(s*s))/2; p=-mp.sqrt(q*q+b*b); return mp.sqrt(2*(p1-p)/p1)/(s+1/s)
    I3=mp.quad(lambda s: mp.log((1+s)/(1-s))*mp.diff(u3,s),[0,0.5,1])/2
    h=lambda f: mp.cos(mp.asin(a*mp.sin(2*f))/2)
    U=lambda f: h(f)*mp.sin(f)+mp.diff(h,f)*mp.cos(f); V=lambda f: h(f)*mp.cos(f)-mp.diff(h,f)*mp.sin(f)
    lz=lambda f: mp.log(mp.cot(f/2)); lw=lambda f: mp.log(mp.cot(mp.pi/4-f/2))
    Iz=mp.quad(lambda f: lz(f)*mp.diff(V,f),[mp.pi/2,3*mp.pi/8,mp.pi/4])/2
    Iw=mp.quad(lambda f: lw(f)*mp.diff(U,f),[mp.pi/2,3*mp.pi/8,mp.pi/4])/2
    return mp.re(2*(Iz-Iw-I3))
def free_fraction(a):
    a=mp.mpf(a)
    rho=lambda u: 4/mp.pi*mp.atan(mp.exp(-mp.asinh(mp.sqrt((u*u-a*a)/(a*a-u**4)))))
    return a+mp.quad(rho,[a,(a+mp.sqrt(a))/2,mp.sqrt(a)])
if __name__=='__main__':
    mp.mp.dps=20
    for a in (0.05,0.1,0.2,0.4,0.6,0.8,0.95,0.999):
        print(a, mp.nstr(phi_field(a),12), mp.nstr(free_fraction(a),12))
def phi_closed(a):
    """(J/2T)/m = int_0^inf [g(q) - 1] dq/(q^2 - 1),  g = sqrt((1 + sqrt(1 - a^2 + a^2 q^2))/2)
    (the principal value of g/(q^2 - 1); PV int dq/(q^2-1) = 0).  Equals ln(1 + sqrt 2) at a = 1."""
    a=mp.mpf(a); g=lambda q: mp.sqrt((1+mp.sqrt(1-a*a+a*a*q*q))/2)
    f=lambda q: (g(q)-1)/(q*q-1) if abs(q-1)>mp.mpf(10)**(-8) else a*a/8
    return mp.quad(f,[0,0.5,1,2,10,mp.inf])
