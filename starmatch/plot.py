from astropy import visualization as aviz
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path

def show_image(image, figsize=(8, 8), figdpi=300, origin='upper',
               cmap='gray', stretch_mode='linear', clip=True,
               show_colorbar=False, show_ticks=True, mask_rectangle=None, **kwargs):
    """
    Show an astronomical images.

    Inputs:
        image -> [2d array] The image to display.
        figsize -> [tuple of 2 floats, optional, default=(8, 8)] Size of the matplotlib figure in inches.
        figdpi -> [int, optional, default=300] The resolution of the figure in dots per inch
        origin -> [str, optional, default='upper'] Origin of the image coordinates, such as 'upper' or 'lower'.
        cmap -> [str, optional, default='gray'] Colormap to use for the image.
        stretch_mode -> [str, optional, default='linear'] Stretch mode for the image normalization. Options are 'linear', 'log', 'sqrt'.
        clip -> [bool, optional, default=True] Whether to clip the data to the 0-1 range after stretching.
        show_colorbar -> [bool, optional, default=False] Whether to show the colorbar of the image.
        show_ticks -> [bool, optional, default=True] Whether to show axis ticks.
        mask_rectangle -> 【tuple, optional, default=None] If provided, a rectangle to mask on the image.
        Format is ((x, y), width, height), where (x, y) is the anchor point.
        **kwargs -> [dict, optional] Additional keyword arguments:
            - fig_path (str): If provided, the figure will be saved to this path.
            - mark (tuple): If provided, mark points on the image. Format is (xy, marker, color, text).
                - xy (ndarray): Array of points to mark.
                - marker (str): Marker style.
                - color (str): Marker color.
                - text (list of str): Annotations for the points.
    """
    percl, percu = 1, 99  # Percentiles for the lower and upper edge of the stretch

    # Create a figure and axis with the specified size and DPI
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=figdpi)

    # Determine the stretch function to use based on the stretch_mode
    if stretch_mode == 'log':
        stretch = aviz.LogStretch()
    elif stretch_mode == 'sqrt':
        stretch = aviz.SqrtStretch()
    else:
        stretch = aviz.LinearStretch()

    # Normalize the image using the specified percentiles and stretch function
    norm = aviz.ImageNormalize(image, interval=aviz.AsymmetricPercentileInterval(percl, percu),stretch=stretch, clip=clip)

    # Display the image with the specified normalization and colormap
    im = ax.imshow(image, origin=origin, cmap=cmap, aspect='equal', norm=norm)

    # Optionally add a colorbar to the image
    if show_colorbar:
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Optionally remove the axis ticks
    if not show_ticks:
        ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)

    # Optionally mark points on the image
    if 'mark' in kwargs:
        xy, marker, color, text = kwargs['mark']
        plt.scatter(xy[:, 0], xy[:, 1], marker=marker, lw=1, facecolors='none', edgecolors=color, s=70)
        for i in range(len(xy)):
            plt.annotate(text[i], (xy[i, 0], xy[i, 1]), fontsize=10, color='b')

    # Optionally add a rectangle mask to the image
    if mask_rectangle is not None:
        lb_bb, width, height = mask_rectangle
        ax.add_patch(Rectangle(lb_bb, width, height, fill=False, lw=1, color='r', ls='dashdot'))

    # Save the figure to a file if 'fig_path' is provided, otherwise show the image
    if 'fig_path' in kwargs:
        fig_file = kwargs['fig_path']
        Path(fig_file).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(fig_file, bbox_inches='tight')
    else:
        plt.show()

    # Close the plot to free up memory
    plt.close()

