from starcatalogquery import StarCatalog
from scipy.spatial import KDTree
from colorama import Fore

from .astroalign import find_transform_tree,matrix_transform
from .wcs import wcs_trans
      
def get_orientation(camera_xy,camera_asterisms,camera_invariant_tree,pixel_width,catalogfile_h5):
    """
    Obtain the center pointing and pixel width of the camera through blind matching of star maps.

    Inputs:
        camera_xy -> [2d array, n*2] Pixel coordinates of sources in camera.
        camera_asterisms -> [2d array, m*3] Array of the indices of sources that correspond to each invariant triangle.
        camera_invariant_tree -> [Instance of class scipy.spatial.KDTree] 2d-tree constructed from array of (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of sources.
        pixel_width -> [float] Pixel size in deg, such as 0.01
        catalogfile_h5 -> Star catalog index file in h5 format, which records the center pointing of each sky area, the pixel coordinates, the triangle invariants and the asterism indices of the stars.
    Outputs:
        fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in [deg]
        pixel_width_estimate -> [float] Pixel width of camera in [deg]
    """

    # Read the star catalog index file
    fp_radecs,stars_xy,stars_invariants,stars_asterisms = StarCatalog.read_h5_indices(catalogfile_h5)

    camera = (camera_xy,camera_asterisms,camera_invariant_tree)

    n = len(fp_radecs)
    for i in range(n):
        desc = 'Iterate over sky area index {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,i+1,Fore.RESET,n)
        print(desc,end='\r')
        ra_c,dec_c = fp_radecs[i]
        catalog_xy = stars_xy[i]
        catalog_asterisms = stars_asterisms[i]
        catalog_invariant_tree = KDTree(stars_invariants[i])
        stars = (catalog_xy,catalog_asterisms,catalog_invariant_tree)
        # Align the sources in camera and the stars in catalog, and establish the mapping relationship.
        try:
            transf, (pixels_camera_match, pixels_catalog_match),_s,_d = find_transform_tree(camera, stars)
            # Roughly calibrate the center pointing of the camera
            pixels_cc_affine = matrix_transform([0,0],transf.params)  
            pixels_cc_affine_x,pixels_cc_affine_y = pixels_cc_affine[:,0],pixels_cc_affine[:,1]
            wcs = wcs_trans(pixel_width,fp_radecs[i])
            cc_radec_estimate = wcs.pixel_to_world(pixels_cc_affine_x,pixels_cc_affine_y)
            fp_radec = fp_ra,fp_dec = cc_radec_estimate.ra.deg[0],cc_radec_estimate.dec.deg[0]
            pixel_width_estimate = pixel_width * transf.scale
            break
        except:
            fp_radec,pixel_width_estimate = None,None
            continue
    print('')        

    return fp_radec,pixel_width_estimate 

def get_orientation_mp(camera_xy,camera_asterisms,camera_invariant_tree,pixel_width,catalogfile_h5):
    """
    Obtain the center pointing and pixel width of the camera through blind matching of star maps with the multi-core parallel computing.

    Inputs:
        camera_xy -> [2d array, n*2] Pixel coordinates of sources in camera.
        camera_asterisms -> [2d array, m*3] Array of the indices of sources that correspond to each invariant triangle.
        camera_invariant_tree -> [Instance of class scipy.spatial.KDTree] 2d-tree constructed from array of (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of sources.
        pixel_width -> [float] Pixel size in deg, such as 0.01
        catalogfile_h5 -> Star catalog index file in h5 format, which records the center pointing of each sky area, the pixel coordinates, the triangle invariants and the asterism indices of the stars.
    Outputs:
        fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in deg
        pixel_width_estimate -> [float] Pixel width of camera in [deg]
    """
    import multiprocessing as mp
    from functools import partial

    # set the number of feasible cores
    number_of_cores = mp.cpu_count() - 1
    
    # Read the star catalog index file
    fp_radecs,stars_xy,stars_invariants,stars_asterisms = StarCatalog.read_h5_indices(catalogfile_h5)

    n = len(fp_radecs)
    params = (camera_xy,camera_asterisms,camera_invariant_tree,pixel_width)
    # Divide tasks across multiple cores
    args = [(fp_radecs[i],stars_xy[i],stars_asterisms[i],stars_invariants[i]) for i in range(n)]
    func = partial(find_transform_mp, params)
    pool = mp.Pool(number_of_cores)
    reslist = pool.imap_unordered(func,args)
    pool.close()

    for res in reslist:
        if res:  # You can set other condition here
            pool.terminate()
            break
    pool.join()    

    return res 

def find_transform_mp(params,arg):
    """
    Accompanying function for get_orientation_mp.

    Inputs:
        params -> Parameters that do not need to be allocated to multiple cores
        arg -> Variables that need to be allocated to multiple cores
    Outputs:  
        fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in [deg]
        pixel_width_estimate -> [float] Pixel width of camera in [deg]
    """
    camera_xy,camera_asterisms,camera_invariant_tree,pixel_width = params
    camera = (camera_xy,camera_asterisms,camera_invariant_tree)
    fp_radecs_i,stars_xy_i,stars_asterisms_i,stars_invariants_i = arg

    ra_c,dec_c = fp_radecs_i
    catalog_xy = stars_xy_i
    catalog_asterisms = stars_asterisms_i
    catalog_invariant_tree = KDTree(stars_invariants_i)

    stars = (catalog_xy,catalog_asterisms,catalog_invariant_tree)
    # Align the sources in camera and the stars in catalog, and establish the mapping relationship.
    try:
        transf, (pixels_camera_match, pixels_catalog_match),_s,_d = find_transform_tree(camera,stars)
        # Roughly calibrate the center pointing of the camera
        pixels_cc_affine = matrix_transform([0,0],transf.params)  
        pixels_cc_affine_x,pixels_cc_affine_y = pixels_cc_affine[:,0],pixels_cc_affine[:,1]
        wcs = wcs_trans(pixel_width,fp_radecs_i)
        cc_radec_estimate = wcs.pixel_to_world(pixels_cc_affine_x,pixels_cc_affine_y)
        fp_radec = fp_ra,fp_dec = cc_radec_estimate.ra.deg[0],cc_radec_estimate.dec.deg[0]
        pixel_width_estimate = pixel_width * transf.scale
        res = (fp_radec,pixel_width_estimate)
    except:
        res = None
    return res