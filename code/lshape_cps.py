"""Arctic curve of the L-shaped domain from Colomo, Pronko, Sportiello, J. Stat. Phys. 2018 (arXiv:1807.07549),
eqs. (4.15)-(4.17), symmetric cut (Q = 0), six-vertex coordinates in the unit square.  alpha = 1/2 is the uniform domino measure."""
import numpy as np, mpmath as mp
def cps_curve(R,alpha=0.5,m=4000):
    xx,xy=R/(R+1),1/(R+1); sa=mp.sqrt(alpha)
    u=lambda w:(alpha*w+1-alpha)/w
    M=lambda w:w/((w-1)*(alpha*w+1-alpha))
    def Phi(w,sg):
        U=u(w); rad=xx**2*(U-alpha)*(U-1)+xy**2*(1+sa)**2*U
        return alpha/(2*U)+xx/2+xy*(U-alpha)/(2*U*(U-1))+sg*(U-sa)*mp.sqrt(rad)/(2*U*(U-1))
    pts=[]
    for sg in (1,-1):
        for s in np.linspace(-np.pi/2,np.pi/2,m)[1:-1]:
            w=mp.tan(s)                                  # whole real axis
            try:
                U=u(w); rad=xx**2*(U-alpha)*(U-1)+xy**2*(1+sa)**2*U
                if rad<=0: continue
                dP=mp.diff(lambda v:Phi(v,sg),w); dM=mp.diff(M,w)
                pts.append((float(-M(w)*dP/dM+Phi(w,sg)),float(-dP/dM),sg))
            except ZeroDivisionError: pass
    return np.array(pts)
if __name__=='__main__':
    p=cps_curve(1.5); ok=(abs(p[:,0])<3)&(abs(p[:,1])<3); print(len(p),p[ok][:,0].min(),p[ok][:,0].max(),p[ok][:,1].min(),p[ok][:,1].max())
