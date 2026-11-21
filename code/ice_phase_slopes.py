import numpy as np
def blk(c): return (c[0]+c[1])%2==0
def dh(dom,x,y,dx,dy):
    if dx==1: L=(x,y); R=(x,y-1)
    elif dx==-1: L=(x-1,y-1); R=(x-1,y)
    elif dy==1: L=(x-1,y); R=(x,y)
    else: L=(x,y-1); R=(x-1,y-1)
    if dom.get(L)==R: return -3 if blk(L) else 3
    if L in dom or R not in dom: return 1 if blk(L) else -1
    return -1 if blk(R) else 1
def tiling(kind,N=60):
    dom={}
    for j in range(-N,N):
        for i in range(-N,N):
            if kind=='Hstag':   # horizontal bricks, row j shifted by j
                if (i-j)%2==0: dom[(i,j)]=(i+1,j); dom[(i+1,j)]=(i,j)
            if kind=='Hstag2':  # shifted the other way
                if (i+j)%2==1: dom[(i,j)]=(i+1,j); dom[(i+1,j)]=(i,j)
            if kind=='Valign':
                if j%2==0: dom[(i,j)]=(i,j+1); dom[(i,j+1)]=(i,j)
            if kind=='Halign':
                if i%2==0: dom[(i,j)]=(i+1,j); dom[(i+1,j)]=(i,j)
            if kind=='Vstag':
                if (j-i)%2==0: dom[(i,j)]=(i,j+1); dom[(i,j+1)]=(i,j)
            if kind=='Vstag2':
                if (j+i)%2==1: dom[(i,j)]=(i,j+1); dom[(i,j+1)]=(i,j)
    return dom
for k in ['Halign','Valign','Hstag','Hstag2','Vstag','Vstag2']:
    d=tiling(k)
    hx=sum(dh(d,x,0,1,0) for x in range(0,40))/40
    hy=sum(dh(d,0,y,0,1) for y in range(0,40))/40
    print(f'{k:7s} grad h = ({hx:+.2f},{hy:+.2f})')
