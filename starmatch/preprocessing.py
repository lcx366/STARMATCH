import numpy as np
from loess.loess_2d import loess_2d
from statsmodels.robust.scale import mad

def lowess_smooth(xy, UV, frac=0.5):
    """
    Identify outliers in star map matching using the method of LOWESS (Locally Weighted Scatterplot Smoothing).
    LOWESS uses a weighted **linear regression** to fit the data locally.

    Usage:
        >>> is_outlier = lowess_smooth(xy, UV, frac=0.5)
    Inputs:
        xy  -> [2D array] Coordinates of sources, with shape (n, 2).
        UV  -> [2D array] Distortion values, with shape (n, 2).
        frac -> [float, optional, default=0.5] The fraction of the data used in each local regression.
                The value of 'frac' should be between 0 and 1. A smaller 'frac' makes the model more sensitive to local changes,
                while a larger 'frac' produces a smoother fit.
    Returns:
        is_outlier -> [array-like of bool] A boolean array of shape (n,), where each element indicates whether the corresponding
                      data point is identified as an outlier (True) or not (False).
    """
    # Extract x, y coordinates from xy array
    x, y = xy.T
    # Extract U, V distortions from UV array
    U, V = UV.T

    # Fit LOWESS models for U and V separately
    U_fit, w_U = loess_2d(x, y, U, frac=frac)
    V_fit, w_V = loess_2d(x, y, V, frac=frac)

    # Compute residuals for U and V
    resi_U = U - U_fit
    resi_V = V - V_fit

    # Calculate the median of the residuals for U and V
    median_U = np.median(resi_U)
    median_V = np.median(resi_V)

    # Calculate the Median Absolute Deviation (MAD) for the residuals of U and V
    sigma_U = mad(resi_U)
    sigma_V = mad(resi_V)

    # Identify outliers: deviations greater than 4 times the MAD from the median
    flags = (np.abs(resi_U - median_U) > 4 * sigma_U) | (np.abs(resi_V - median_V) > 4 * sigma_V)

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

