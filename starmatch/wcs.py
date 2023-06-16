import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

def wcs_trans(pixel_width,fp_radec):
    """
    Define the WCS transformation between the celestial coordinates and the pixel coordinates.

    Usage:
        >>> wcs = wcs_trans(0.01,[20,30]):
    
    Inputs: 
        pixel_width -> [list or tuple of int] Pixel width in [deg], such as [0.001,0.001]
        fp_radec -> [list or tuple of int] Fiducial point(center pointing) of the camera in form of [Ra,Dec] in [deg]

    Outputs:
        wcs -> object of class WCS   
    """
    if np.isscalar(pixel_width):
        pixel_width_axis1 = pixel_width_axis2 = pixel_width
    else:
        pixel_width_axis1, pixel_width_axis2 = pixel_width
    
    fp_ra,fp_dec = fp_radec
    
    wcs_input_dict = {
        'CTYPE1': 'RA---TAN',
        'CUNIT1': 'deg',
        'CDELT1': pixel_width_axis1,
        'CRPIX1': 1,
        'CRVAL1': fp_ra,
        'CTYPE2': 'DEC--TAN',
        'CUNIT2': 'deg',
        'CDELT2': pixel_width_axis2,
        'CRPIX2': 1,
        'CRVAL2': fp_dec  
    }
    wcs = WCS(wcs_input_dict)
    return wcs    

def xy_catalog(fp_radec,radec,pixel_width):
    """
    Calculate the pixel coordinates of stars in a sky area.

    Usage:
        >>> x,y = xy_catalog([10,20],[[11,15],[22,-4]],0.01)

    Inputs:
        fp_radec -> [list or tuple of int] Fiducial point(center pointing) of the camera in form of [Ra,Dec] in [deg]
        radec -> [2d array of float] Celestial coordinates of stars in format of [[Ra0,Dec0],..,[Ran,Decn]] in [deg]
        pixel_width -> [float] Pixel width in [deg]

    Outputs:
        x -> [array of float] x components of the star pixel coordinates       
        y -> [array of float] y components of the star pixel coordinates   
        wcs -> Instance of class WCS 
    """
    wcs = wcs_trans(pixel_width,fp_radec)
    radec = SkyCoord(radec, unit='deg')
    x,y = np.array(wcs.world_to_pixel(radec))

    return x,y,wcs 