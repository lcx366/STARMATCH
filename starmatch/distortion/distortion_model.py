
import numpy as np
import warnings

def RadialStandard(coeffs,dc,pixels_xy):
    """
    Calculate the corrected pixel coordinates based on the Standard Radial Distortion Model(SRDM) defined by the distortion center and the distortion coefficients.

    Inputs:
        coeffs -> [list] Coefficients associated with the SRDM, such as [-1e-4,1e-4]
        dc -> [list(2 elements)] Distortion center associated with the SRDM, such as [0.1,0.1]
        pixels_xy -> [list(2 elements),array(n*2)] Distorted pixel coordinates, such as [367,125] or [[2.8671875, 0.9765625], [1.109375, -0.875]]
    Outputs：
        xuyu -> [2d array(n*2)] Distortion-corrected pixel coordinates

    Note:
    The 3rd order polynomial in standard form: ru = rd + k1 * rd**3
    The 5th order polynomial in standard form: ru = rd + k1 * rd**3 + k2 * rd**5
    The 5th order polynomial in all form: ru = rd + k1 * rd**2 + k2 * rd**3 + k3 * rd**4 + k4 * rd**5,
    where rd is the distance between the distorted pixel coordinates and the distortion center, ru is the distance between the distortion-corrected pixel coordinates and the distortion center.

    1. The 3rd order polynomial in standard form only works well for small amounts of distortion.
    2. The 5th order polynomial produce more accurate results, especially for “wave” or “mustache” distortion, which might resemble barrel near the center of the image and pincushion near the corners (or vice-versa).
    3. The 5th order all form use all coefficients up to the maximum instead of alternate coefficients (odd or even-only). We have not observed much advantage to these settings.
    4. Higher order polynomials (7th order or higher) should be used with great care because results can become unstable, especially at the outer parts of the image. The image used to calculate the coefficients should have valid corner points near the edge of the image and there should be sufficient rows or columns. 

    Basic formulas are as follows:
    xu - xc = ru * Cos(theta) = ru * (xd - xc)/rd = ... 
    yu - yc = ru * Sin(theta) = ru * (yd - yc)/rd = ... 
    where (x_c,y_c) is the pixel coordinates of the distortion center.

    For more details, please refer to 
    1. https://en.wikipedia.org/wiki/Distortion_(optics)
    2. https://www.imatest.com/docs/distortion-methods-and-modules/
    3. https://www.imatest.com/support/docs/pre-5-2/geometric-calibration-deprecated/distortion-models/ 
    """

    xdyd = pixels_xy
    xdyd_xcyc = xdyd - dc

    rd2 = np.sum(xdyd_xcyc**2,axis=1) # rd**2
    rd4 = rd2**2 

    n = len(coeffs)
    if n == 1:
        # 3rd order polynomial
        f = coeffs[0]*rd2
    elif n == 2:
        # 5th order polynomial
        f = coeffs[0]*rd2 + coeffs[1]*rd4
    elif n > 2:
        f = coeffs[0]*rd2 + coeffs[1]*rd4
        warnings.warn("Only 5th order standard-form radial distortion is used.")  
    else:
        raise Exception("Radial distortion model requires at least one parameter.")

    xuyu = xdyd + xdyd_xcyc*f[:,None]

    return xuyu

def RadialDivision(coeffs,dc,pixels_xy):
    """
    Calculate the corrected pixel coordinates based on the Division-mode Radial Distortion Model defined by the distortion center and the distortion coefficients.

    Inputs:
        coeffs -> [list] Coefficients associated with the DRDM, such as [-1e-4,1e-4]
        dc -> [list(2 elements)] Distortion center associated with the DRDM, such as [0.1,0.1]
        pixels_xy -> [list(2 elements),array(n*2)] Distorted pixel coordinates, such as [367,125] or [[2.8671875, 0.9765625], [1.109375, -0.875]]
    Outputs：
        xuyu -> [2d array(n*2)] Distortion-corrected pixel coordinates    

    Note:
    The Division-mode Radial Distortion Model typically provides a more accurate approximation than the Standard Radial Distortion Model for the same number of coefficients. 

    The 2nd order polynomial in division form: ru = rd /(1+ k1 * rd**2)
    The 4th order polynomial in division form: ru = rd /(1+ k1 * rd**2 + k2 * rd**4)
    where rd is the distance between the distorted pixel coordinates and the distortion center,
    ru is is the distance between the distortion-corrected pixel coordinates and the distortion center.

    Basic formulas:
    xu - xc = ru * Cos(theta) = ru * (xd - xc)/rd = ... 
    yu - yc = ru * Sin(theta) = ru * (yd - yc)/rd = ... 

    For more details, please refer to 
    1. https://en.wikipedia.org/wiki/Distortion_(optics)
    2. https://www.imatest.com/docs/distortion-methods-and-modules/
    3. https://www.imatest.com/support/docs/pre-5-2/geometric-calibration-deprecated/distortion-models/ 
    """

    xdyd = pixels_xy
    xdyd_xcyc = xdyd - dc

    rd2 = np.sum(xdyd_xcyc**2,axis=1) # rd**2
    rd4 = rd2**2 

    n = len(coeffs)
    if n == 1:
        # 2nd order polynomial
        f = 1 + coeffs[0]*rd2
    elif n == 2:
        # 4th order polynomial
        f = 1 + coeffs[0]*rd2 + coeffs[1]*rd4
    elif n > 2:    
        f = 1 + coeffs[0]*rd2 + coeffs[1]*rd4
        warnings.warn("Only 4th order division-form radial distortion is used.")   
    else:
        raise Exception("Radial distortion model requires at least one parameter.")  

    xuyu = dc + xdyd_xcyc/f[:,None]

    return xuyu   


def Tangential(coeffs,dc,pixels_xy):
    """
    Calculate the corrected pixel coordinates based on the Tangential Distortion Model(also known as the de-centering distortion) defined by the distortion center and the distortion coefficients.

    Inputs:
        coeffs -> [list] Coefficients associated with the TDM, such as [-1e-4,1e-4]
        dc -> [list(2 elements)] Distortion center associated with the TDM, such as [0.1,0.1]
        pixels_xy -> [list(2 elements),array(n*2)] Distorted pixel coordinates, such as [367,125] or [[2.8671875, 0.9765625], [1.109375, -0.875]]
    Outputs：
        xuyu -> [2d array(n*2)] Distortion-corrected pixel coordinates  

    Note:
    Basic formulas:
    xu = xd + (P1 * (rd**2 + 2*(xd-xc)**2) + 2*P2 * (xd-xc)*(yd-yc))*(1 + P3*rd**2 + P4*rd**4 + ...)
    yu = yd + (P2 * (rd**2 + 2*(yd-yc)**2) + 2*P1 * (xd-xc)*(yd-yc))*(1 + P3*rd**2 + P4*rd**4 + ...)

    For more details, please refer to 
    1. https://en.wikipedia.org/wiki/Distortion_(optics)
    2. https://www.imatest.com/docs/distortion-methods-and-modules/
    3. https://www.imatest.com/support/docs/pre-5-2/geometric-calibration-deprecated/distortion-models/ 
    """
    xdyd = pixels_xy
    xd,yd = xdyd.T
    xd_xc,yd_yc = xdyd_xcyc_T = (xdyd - dc).T

    rd2 = np.sum(xdyd_xcyc_T**2,axis=0) # rd**2
    rd4 = rd2**2 

    f = 2*xd_xc*yd_yc
    n = len(coeffs)
    if n == 2:
        # 2nd order tangential distortion coefficient
        xu = xd + coeffs[0]*(rd2+2*xd_xc**2) + coeffs[1]*f
        yu = yd + coeffs[1]*(rd2+2*yd_yc**2) + coeffs[0]*f
    elif n == 3:
        # 3rd order tangential distortion coefficient
        g = 1 + coeffs[2]*rd2
        xu = xd + g * (coeffs[0]*(rd2+2*xd_xc**2) + coeffs[1]*f)
        yu = yd + g * (coeffs[1]*(rd2+2*yd_yc**2) + coeffs[0]*f)
    elif n == 4:
        # 4th order tangential distortion coefficient
        g = 1 + coeffs[2]*rd2 + coeffs[3]*rd4
        xu = xd + g * (coeffs[0]*(rd2+2*xd_xc**2) + coeffs[1]*f)
        yu = yd + g * (coeffs[1]*(rd2+2*yd_yc**2) + coeffs[0]*f) 
    elif n > 4:
        g = 1 + coeffs[2]*rd2 + coeffs[3]*rd4
        xu = xd + g * (coeffs[0]*(rd2+2*xd_xc**2) + coeffs[1]*f)
        yu = yd + g * (coeffs[1]*(rd2+2*yd_yc**2) + coeffs[0]*f) 
        warnings.warn("Only 4th order tangential distortion is used.")            
    else:
        raise Exception("Tangential distortion model requires at least two parameters.")  

    return np.stack([xu,yu]).T      

def Brown_Conrady(coeffs,dc,pixels_xy):
    """
    Calculate the corrected pixel coordinates based on the Brown-Conrady Distortion Model(BCDM) defined by the distortion center and the distortion coefficients.

    Inputs:
        coeffs -> [list] Coefficients associated with the BCDM, such as [[-1e-4,1e-4],[1e-3,1e-3,1e-4,1e-5]]
        dc -> [list(2 elements)] Distortion center associated with the BCDM, such as [0.1,0.1]
        pixels_xy -> [list(2 elements),array(n*2)] Distorted pixel coordinates, such as [367,125] or [[2.8671875, 0.9765625], [1.109375, -0.875]]
    Outputs：
        xuyu -> [2d array(n*2)] Distortion-corrected pixel coordinates      

    Note:
    The Brown–Conrady model corrects both the radial distortion and the tangential distortion caused by physical elements in a lens not being perfectly aligned.

    Basic formulas:
    xu = xd + (xd - xc) * (K1*rd**2 + K2*rd**4 + ...) + (P1 * (rd**2 + 2*(xd-xc)**2) + 2*P2 * (xd-xc)*(yd-yc))*(1 + P3*rd**2 + P4*rd**4 + ...)
    yu = yd + (yd - xc) * (K1*rd**2 + K2*rd**4 + ...) + (P2 * (rd**2 + 2*(yd-yc)**2) + 2*P1 * (xd-xc)*(yd-yc))*(1 + P3*rd**2 + P4*rd**4 + ...)

    For more details, please refer to 
    1. https://en.wikipedia.org/wiki/Distortion_(optics)
    2. https://www.imatest.com/docs/distortion-methods-and-modules/
    3. https://www.imatest.com/support/docs/pre-5-2/geometric-calibration-deprecated/distortion-models/ 
    """
    xc,yc = dc
    xdyd = pixels_xy
    xdyd_xcyc = xdyd - dc
    xd_xc,yd_yc = xdyd_xcyc.T
    try:
        coeffs_r,coeffs_t = coeffs
    except:
        raise Exception("'coeffs' should be in form of [[coeffs_r],[coeffs_t]] for Brown–Conrady Distortion model.")    

    rd2 = np.sum(xdyd_xcyc**2,axis=1) # rd**2
    rd4 = rd2**2

    n_r,n_t = len(coeffs_r),len(coeffs_t)

    # for radial distortion
    if n_r == 1:
        # 3rd order polynomial
        f_r = coeffs_r[0]*rd2
    elif n_r == 2:
        # 5th order polynomial
        f_r = coeffs_r[0]*rd2 + coeffs_r[1]*rd4
    elif n_r > 2:
        f_r = coeffs_r[0]*rd2 + coeffs_r[1]*rd4
        warnings.warn("Only 5th order standard-form radial distortion is used.")  
    else:
        raise Exception("Radial distortion model requires at least one parameter.")  

    xuyu_r = xdyd + xdyd_xcyc*f_r[:,None]

    # for tangential distortion
    f_t = 2*xd_xc*yd_yc
    if n_t == 2:
        # 2nd order tangential distortion coefficient
        xu_t = coeffs_t[0]*(rd2+2*xd_xc**2) + coeffs_t[1]*f_t
        yu_t = coeffs_t[1]*(rd2+2*yd_yc**2) + coeffs_t[0]*f_t
    elif n_t == 3:
        # 3rd order tangential distortion coefficient
        g = 1 + coeffs_t[2]*rd2
        xu_t = g * (coeffs_t[0]*(rd2+2*xd_xc**2) + coeffs_t[1]*f_t)
        yu_t = g * (coeffs_t[1]*(rd2+2*yd_yc**2) + coeffs_t[0]*f_t)
    elif n_t == 4:
        # 4th order tangential distortion coefficient
        g = 1 + coeffs_t[2]*rd2 + coeffs_t[3]*rd4
        xu_t = g * (coeffs_t[0]*(rd2+2*xd_xc**2) + coeffs_t[1]*f_t)
        yu_t = g * (coeffs_t[1]*(rd2+2*yd_yc**2) + coeffs_t[0]*f_t)    
    elif n_t > 4:
        g = 1 + coeffs_t[2]*rd2 + coeffs_t[3]*rd4
        xu_t = g * (coeffs_t[0]*(rd2+2*xd_xc**2) + coeffs_t[1]*f_t)
        yu_t = g * (coeffs_t[1]*(rd2+2*yd_yc**2) + coeffs_t[0]*f_t) 
        warnings.warn("Only 4th order tangential distortion is used.")         
    else:
        raise Exception("Tangential distortion model requires at least two parameters.")  
    xuyu_t = np.stack([xu_t,yu_t]).T      
    xuyu = xuyu_r + xuyu_t

    return xuyu