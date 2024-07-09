import numpy as np
from loess.loess_2d import loess_2d

W = 0.34 # Weights less than W are identified as outliers. Outliers are defined as those with deviations > 4sigma.
def lowess_smooth(xy,UV,frac=0.5):
    """
    Identify outliers in star map matching with the method of LOWESS (Locally Weighted Scatterplot Smoothing).
    Here, LOWESS uses a weighted **linear regression** by default.

    Usage:
        >>> is_outlier = lowess_smooth(xy,UV)
    Inputs:
        xy -> [2D array] Coordinates of sources with shape of (n,2)
        UV -> [2D array] Distortion with shape of (n,2)
        frac -> [float,optional,default=0.5] The fraction of the data used in local regression.
        The value of fraction is between 0 and 1.
    Outputs:
        is_outlier -> [array-like of bool] A boolean array marking the outliers.
    """
    x,y = xy.T
    U,V = UV.T
    _, w_U = loess_2d(x, y, U, frac=frac)
    _, w_V = loess_2d(x, y, V, frac=frac)

    flags = (w_U < W) | (w_V < W)

    return flags

def _iqr_outliers(z, scale=4.0):
    """
    Identifies outliers in a dataset using the Interquartile Range (IQR) method.

    Inputs:
        z -> [array-like] Dataset
        scale -> [float] Scaling factor for the IQR to adjust sensitivity.
        A lower value makes the method more sensitive to outliers.
    Outputs:
        is_outlier -> [array-like,bool] A boolean array marking the outliers.
    """
    # Compute the first and third quartiles (25th and 75th percentiles)
    Q1 = np.nanpercentile(z, 25)
    Q3 = np.nanpercentile(z, 75)

    # Calculate the IQR and determine outlier thresholds
    IQR = Q3 - Q1
    lower_bound = Q1 - scale * IQR
    upper_bound = Q3 + scale * IQR

    # Identify outliers
    flags = (z < lower_bound) | (z > upper_bound)

    return flags

def iqr_outliers(UV, scale=4.0):
    """
    Identify outliers in star map matching with the method of Interquartile Range (IQR).

    Usage:
        is_outlier = iqr_outliers(UV)
    Inputs:
        UV -> [2D array] Distortion with shape of (n,2)
        scale -> [float] Scaling factor for the IQR to adjust sensitivity.
        A lower value makes the method more sensitive to outliers.
    Outputs:
        is_outlier -> [array-like,bool] A boolean array marking the outliers.
    """
    flags_U = _iqr_outliers(UV[:, 0], scale=scale)
    flags_V = _iqr_outliers(UV[:, 1], scale=scale)
    flags = flags_U | flags_V
    return flags

