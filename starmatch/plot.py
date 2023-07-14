from astropy import visualization as aviz
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path

def show_image(image,figsize=(8, 8),figdpi=300,origin='upper',
               cmap='gray', stretch_mode='linear', clip=True,
               show_colorbar=False, show_ticks=True,mask_rectangle=None,**kwargs):
    """
    Show an image in matplotlib with some basic astronomically-appropriat stretching.

    Inputs:
    image
        The image to show
    percl : number
        The percentile for the lower edge of the stretch (or both edges if ``percu`` is None)
    percu : number or None
        The percentile for the upper edge of the stretch (or None to use ``percl`` for both)
    figsize : 2-tuple
        The size of the matplotlib figure in inches
    """
    percl,percu = 1,99

    fig, ax = plt.subplots(1, 1, figsize=figsize,dpi=figdpi)

    # To preserve details we should *really* downsample correctly and
    # not rely on matplotlib to do it correctly for us (it won't).

    # So, calculate the size of the figure in pixels, block_reduce to
    # roughly that,and display the block reduced image.

    if stretch_mode == 'log':
        stretch = aviz.LogStretch()
    elif stretch_mode == 'sqrt':
        stretch = aviz.SqrtStretch()    
    else:
        stretch = aviz.LinearStretch()

    norm = aviz.ImageNormalize(image,interval=aviz.AsymmetricPercentileInterval(percl, percu),stretch=stretch, clip=clip)
    scale_args = dict(norm=norm)
    im = ax.imshow(image,origin=origin,cmap=cmap, aspect='equal', **scale_args)

    if show_colorbar:
        # I haven't a clue why the fraction and pad arguments below work to make
        # the colorbar the same height as the image, but they do....unless the image
        # is wider than it is tall. Sticking with this for now anyway...
        # Thanks: https://stackoverflow.com/a/26720422/3486425
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        # In case someone in the future wants to improve this:
        # https://joseph-long.com/writing/colorbars/
        # https://stackoverflow.com/a/33505522/3486425
        # https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes

    if not show_ticks:
        ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)     

    if 'mark' in kwargs:
        xy,marker,color,text = kwargs['mark']
        plt.scatter(xy[:,0],xy[:,1],marker=marker, lw=1, facecolors='none', edgecolors=color,s=70)

        for i in range(len(xy)):
            plt.annotate(text[i], (xy[i,0], xy[i,1]),fontsize=10,color='b')

    if mask_rectangle is not None:
        lb_bb,width,height = mask_rectangle
        ax.add_patch(Rectangle(lb_bb, width, height,fill=False,lw=1,color='r',ls='dashdot'))   

    if 'figname' in kwargs: 
        fig_file = kwargs['figname']
        Path(fig_file).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(fig_file,bbox_inches='tight') 
    else:
        plt.show()
            
    plt.close()
