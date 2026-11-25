"""Edge probabilities of a weighted Aztec diamond by urban renewal (Propp's generalized domino shuffling), in the log domain.

Every step adds or multiplies positive numbers, so there is no cancellation and double precision is enough at any size.
A sub-region is obtained by embedding: edges that must not be used get log-weight LOG_EPS (effectively zero).

Geometry (chess coordinates).  Aztec diamond of order N centred at the lattice point (cx, cy): cells (a, b) (lower-left
corners) with |a + 1/2 - cx| + |b + 1/2 - cy| <= N.  Faces used by the recursion: lattice points (X, Y) relative to the
centre with |X| + |Y| <= N - 1 and X + Y = N - 1 (mod 2);  i = (X + Y + N - 1)/2,  j = (Y - X + N - 1)/2.
The four edges of face (i, j) sit at  [2i, 2j] bottom, [2i+1, 2j+1] top, [2i, 2j+1] left, [2i+1, 2j] right."""
import numpy as np
LOG_EPS=-600.0
def face_edges(N,cx,cy):
    """dict: (row, col) of the 2N x 2N array -> the two cells of that edge"""
    E={}
    for i in range(N):
        for j in range(N):
            X=i-j+cx; Y=i+j-(N-1)+cy
            sw,se,nw,ne=(X-1,Y-1),(X,Y-1),(X-1,Y),(X,Y)
            E[(2*i,2*j)]=(sw,se); E[(2*i+1,2*j+1)]=(nw,ne); E[(2*i,2*j+1)]=(sw,nw); E[(2*i+1,2*j)]=(se,ne)
    return E
def edge_probabilities(logA):
    """logA: 2N x 2N array of log-weights.  Returns the 2N x 2N array of edge probabilities."""
    N=logA.shape[0]//2; Ws=[logA]
    for k in range(N,1,-1):                       # weight reduction: order k -> k-1
        A=Ws[-1]; w=A[0::2,0::2]; x=A[0::2,1::2]; y=A[1::2,0::2]; z=A[1::2,1::2]
        d=np.logaddexp(w+z,x+y); B=np.empty_like(A)
        B[0::2,0::2]=z-d; B[0::2,1::2]=y-d; B[1::2,0::2]=x-d; B[1::2,1::2]=w-d
        Ws.append(B[1:-1,1:-1])
    P=None
    for A in reversed(Ws):                        # probabilities: order 1 -> N
        k=A.shape[0]//2; w=A[0::2,0::2]; x=A[0::2,1::2]; y=A[1::2,0::2]; z=A[1::2,1::2]
        d=np.logaddexp(w+z,x+y); a=np.exp(w+z-d); b=np.exp(x+y-d)
        Q=np.zeros((2*k,2*k))
        if P is not None: Q[1:-1,1:-1]=P
        p=Q[0::2,0::2]; q=Q[0::2,1::2]; r=Q[1::2,0::2]; s=Q[1::2,1::2]
        c=1-p-q-r-s; Pn=np.empty((2*k,2*k))
        Pn[0::2,0::2]=s+c*a; Pn[0::2,1::2]=r+c*b; Pn[1::2,0::2]=q+c*b; Pn[1::2,1::2]=p+c*a
        P=Pn
    return P
def region_probabilities(cells,N,cx,cy,forced=None):
    """Uniform tilings of `cells` (a subset of the Aztec diamond of order N at (cx,cy)); the other cells are covered by
    the fixed perfect matching `forced` (dict cell -> partner).  Returns dict frozenset({c1,c2}) -> probability."""
    E=face_edges(N,cx,cy); cells=set(cells); forced=forced or {}
    logA=np.empty((2*N,2*N))
    for (r,c),(u,v) in E.items():
        if u in cells and v in cells: logA[r,c]=0.0
        elif forced.get(u)==v: logA[r,c]=0.0
        else: logA[r,c]=LOG_EPS
    P=edge_probabilities(logA)
    return {frozenset(E[k]):P[k] for k in E if E[k][0] in cells and E[k][1] in cells}
if __name__=='__main__':
    # check against a Kasteleyn determinant with random weights
    rng=np.random.default_rng(1)
    for N in (2,3,5,8):
        E=face_edges(N,0,0); logA=rng.normal(size=(2*N,2*N)); P=edge_probabilities(logA)
        cells=sorted({c for e in E.values() for c in e}); B=[c for c in cells if (c[0]+c[1])%2==0]; W=[c for c in cells if (c[0]+c[1])%2]
        bi={c:k for k,c in enumerate(B)}; wi={c:k for k,c in enumerate(W)}; K=np.zeros((len(B),len(W)))
        for (r,c),(u,v) in E.items():
            b,w=(u,v) if u in bi else (v,u); sign=1.0 if u[1]==v[1] else (-1.0)**b[0]
            K[bi[b],wi[w]]=sign*np.exp(logA[r,c])
        Ki=np.linalg.inv(K); err=0
        for (r,c),(u,v) in E.items():
            b,w=(u,v) if u in bi else (v,u); err=max(err,abs(K[bi[b],wi[w]]*Ki[wi[w],bi[b]]-P[r,c]))
        print(f'N={N}: cells {len(cells)} (expected {2*N*(N+1)}), max |P_shuffle - P_Kasteleyn| = {err:.2e}')
