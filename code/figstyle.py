"""Common style for the paper figures: final widths (3.375 in / 7.0 in), 8 pt text, Okabe-Ito colours, one liquidity colour map."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
W1,W2=3.375,7.0
OI=dict(black='#000000',orange='#E69F00',sky='#56B4E9',green='#009E73',yellow='#F0E442',blue='#0072B2',red='#D55E00',purple='#CC79A7')
H_FILL='#F6D7A7'; V_FILL='#CFE6F5'          # H and V plaquettes (tints of orange / sky)
ACCENT=OI['red']                              # theory curves and the one highlighted object per figure
LIQ_CMAP='viridis'; LIQ_MAX=0.75              # liquidity 1 - P_max, same range in every map
def setup():
    plt.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
        'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,'lines.linewidth':1.1,'lines.markersize':3.5,
        'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans','pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300,'legend.frameon':False})
def label(ax,s,x=-0.02,y=1.02):
    ax.text(x,y,f'({s})',transform=ax.transAxes,fontsize=8,fontweight='bold',ha='right',va='bottom')
def save(fig,name):
    fig.savefig(f'../figs/paper/{name}.pdf',bbox_inches='tight'); fig.savefig(f'../figs/paper/{name}.png',dpi=200,bbox_inches='tight')
def liq_image(ax,cen,val,lim,res=600,rmax=1.6):
    """liquidity map as an image: nearest cell value on a regular grid, blank farther than rmax from any cell centre"""
    import numpy as np
    from scipy.spatial import cKDTree
    (x0,x1),(y0,y1)=lim; xs=np.linspace(x0,x1,res); ys=np.linspace(y0,y1,int(res*(y1-y0)/(x1-x0)))
    X,Y=np.meshgrid(xs,ys); d,i=cKDTree(cen).query(np.c_[X.ravel(),Y.ravel()])
    Z=np.where(d<rmax,np.asarray(val)[i],np.nan).reshape(X.shape)
    return ax.imshow(Z,origin='lower',extent=(x0,x1,y0,y1),cmap=LIQ_CMAP,vmin=0,vmax=LIQ_MAX,interpolation='nearest',rasterized=True)
