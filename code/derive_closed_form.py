"""Symbolic derivation of the closed-form parameters of the class-10 curve (sympy) and comparison with the numerical root-finder."""
import sympy as sp
s=sp.symbols('sigma',positive=True)      # sigma = 2 sin(theta_1)
kap=(2+s)/(1+s)                           # kappa = 2 cos(theta_2)
p=sp.simplify(1-s*kap)
u1=s**2*(4-s**2)/2; u2=kap**2*(4-kap**2)/2
D1=(1-p)**2+2*p*u1; D2=(1-p)**2+2*p*u2
ratio=sp.simplify((D1*u2+D2*u1)/((u2-u1)*(1-p)**2))      # n/c
r=sp.factor(sp.simplify(1/ratio))                         # c/n
print('p =',sp.factor(p)); print('c/n =',r)
print('limits: sigma->0:',sp.limit(r,s,0),'  sigma->sqrt2:',sp.simplify(r.subs(s,sp.sqrt(2))),'  p at sqrt2:',sp.simplify(p.subs(s,sp.sqrt(2))))
# a0,a4 for n=1 : a4 = c (D1+D2)/(u2-u1) ; a0 = c D1 - a4 (1-u1)
c=r; a4=sp.factor(sp.simplify(c*(D1+D2)/(u2-u1))); a0=sp.factor(sp.simplify(c*D1-a4*(1-u1)))
print('a0/n =',a0); print('a4/n =',a4)
print('check H(0)=1:',sp.simplify((a0+a4)/(1-p)**2))
# polynomial relation between rho=c/n and sigma
rho=sp.symbols('rho'); num,den=sp.fraction(sp.together(r)); print('relation:',sp.expand(rho*den-num),'= 0')
# numeric comparison with the root-finder
from theory import solve
for sv in (0.2,0.5,0.8,1.0,1.2,1.35):
    cn=float(r.subs(s,sv)); sol=solve(1,cn)
    print(f'sigma={sv}: c/n={cn:.6f} closed p,a0,a4=({float(p.subs(s,sv)):+.6f},{float(a0.subs(s,sv)):+.6f},{float(a4.subs(s,sv)):+.6f})  solver={tuple(round(float(x),6) for x in sol[0]) if sol else None}')
