from starcatalogquery.wcs import wcs_trans
from starcatalogquery.catalog_index import find_healpix_level
from scipy.spatial import KDTree
import multiprocessing as mp
from functools import partial

from .astroalign import find_transform_tree,matrix_transform

PIXEL_WIDTH_GUESS = 0.01 # Initial guess value of pixel width in deg

def generate_args_list(catalog_data):
    """
    Generates a list of arguments for parallel processing of star catalog blind matching.

    Inputs:
        catalog_data -> [dict] A dictionary with levels 'K1' to 'K11'.
                           Each level has a list of tuples (xy, invariants, asterisms).
    Outputs:
        args_list -> [list] A list of tuples. Each tuple prepared for parallel processing.
    """
    args_list = []
    fp_radecs = catalog_data.get('fp_radecs', [])

    # Iterate through levels K1 to K11
    for i in range(1, 12):
        level_key = f"K{i}"
        level_data = catalog_data.get(level_key, [])

        # Create argument tuples for parallel processing
        for fp_radec, partition_data in zip(fp_radecs, level_data):
            stars_xy_i, stars_invariants_i, stars_asterisms_i = partition_data
            arg = (fp_radec, stars_xy_i, stars_asterisms_i, stars_invariants_i)
            args_list.append(arg)

    return args_list

def get_orientation_mp(camera_xy,camera_asterisms,camera_invariant_tree,fov,hashed_data):
    """
    Obtain the center pointing and pixel width of the camera through blind matching of star maps with the multi-core parallel computing.

    Inputs:
        camera_xy -> [2d array, n*2] Pixel coordinates of sources in camera.
        camera_asterisms -> [2d array, m*3 or m*4] Indices of sources that correspond to each invariant.
        camera_invariant_tree -> [scipy.spatial.KDTree] kd-tree constructed from the geometric invariants.
        fov -> [2-ele tuple] Field of view in deg, such as (2,2). If not None, try to find the HEALPix level using fov.
        hashed_data -> Hashed data in h5-formatted file, which records
            - center pointing of each sky area
            - pixel coordinates of stars
            - geometric invariants
            - asterism indices of the stars
    Outputs:
        fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in deg
        pixel_width_estimate -> [float] Pixel width of camera in [deg]
    """
    # set the number of feasible cores
    number_of_cores = mp.cpu_count() - 1

    params = (camera_xy,camera_asterisms,camera_invariant_tree)

    # Divide tasks across multiple cores
    args = generate_args_list(hashed_data)

    # try to find the HEALPix level
    if fov is not None:
        level, nside, npix, pixel_size = find_healpix_level(min(fov))
        start = (int(level[1:]) - 1) * 48
        end = start + 48
        args = args[start:end]

    # Function to be executed in parallel
    func = partial(find_transform_mp, params)

    with mp.Pool(number_of_cores) as pool:
        # imap_unordered allows processing as soon as any worker is free
        result_objects = pool.imap_unordered(func, args)
        try:
            for result in result_objects:
                if result:  # Assuming result is not None if successful
                    return result  # Return the successful result immediately  
        finally:
            pool.terminate()  # Terminate all processes immediately
            pool.join()  # Wait for the worker processes to terminate

def find_transform_mp(params,arg):
    """
    Accompanying function for get_orientation_mp.

    Inputs:
        params -> Variables that do not need to be allocated to multiple cores
        arg -> Variables that need to be allocated to multiple cores
    Outputs:  
        fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in [deg]
        pixel_width_estimate -> [float] Pixel width of camera in [deg]
    """
    camera_xy,camera_asterisms,camera_invariant_tree = params
    camera = (camera_xy,camera_asterisms,camera_invariant_tree)
    fp_radecs_i,stars_xy_i,stars_asterisms_i,stars_invariants_i = arg

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

        # Construct the WCS transform between the pixel coordinates and celestial coordinates
        wcs = wcs_trans(PIXEL_WIDTH_GUESS,fp_radecs_i)

        # Estimate the center pointing
        cc_radec_estimate = wcs.pixel_to_world(pixels_cc_affine_x,pixels_cc_affine_y)
        fp_radec = cc_radec_estimate.ra.deg[0],cc_radec_estimate.dec.deg[0]

        # Estimate the pixel width
        pixel_width_estimate = PIXEL_WIDTH_GUESS * transf.scale
        res = (fp_radec,pixel_width_estimate)
    except:
        res = None
    return res