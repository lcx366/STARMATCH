import numpy as np
from numpy.linalg import norm,det
from scipy.spatial import KDTree
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from pathlib import Path
from skimage.transform import PiecewiseAffineTransform,PolynomialTransform
from starcatalogquery.invariantfeatures import calculate_invariantfeatures
from GPy.kern import RBF
from GPy.models import GPRegression
import statsmodels.api as sm

from .orientation import get_orientation_mp
from .astroalign import find_transform_tree,matrix_transform
from .preprocessing import lowess_smooth,iqr_outliers
from .distortion import distortion_model
from .plot import show_image

# Maximum number of stars to extract from each tile.
MAX_NUM_PER_TILE = 5

# Maximum number of sources used to execute the star map matching.
MAX_CONTROL_POINTS = 30

# Pixel distance tolerance to assume two points are the same for the primary and secondary affine transformation.
PIXEL_TOLS = (20,3)

def photometric_model(F, C):
    """
    Photometric model to calculate the apparent magnitude M of a celestial object
    based on its flux F and the magnitude constant C.

    The relationship follows the classical formula for magnitude:
    M = C - 2.5 * log10(F)

    In this model:
    - M represents the **apparent magnitude** of the object (how bright it appears from Earth).
    - C is the **magnitude constant**, which depends on the system of measurements used.
    - F is the **radiative flux**, typically measured in units like watts per square meter, representing the energy received from the object per unit area.

    The formula is derived from the logarithmic scale used in astronomy for brightness, where a difference of 5 magnitudes corresponds to a factor of 100 in brightness.
    Specifically, each decrease in magnitude by 1 corresponds to an increase in brightness by approximately 2.512 times.

    Usage:
        >>> M = photometric_model(F, C)
    Inputs:
        F -> [float, array-like] The radiative flux of the celestial objects.
        C -> [float] The magnitude constant, typically related to the reference flux.
    Outputs:
        M -> [float] The apparent magnitude of the objects.
    """
    return C - 2.5 * np.log10(F)

def photometric_robust_linear_fit(F, M):
    """
    Perform a robust linear fit to estimate the magnitude constant C and its uncertainty
    based on observed fluxes F and apparent magnitudes M.

    Inputs:
        F -> [array-like,float] The radiative flux of the celestial objects.
        M -> [array-like,float] The apparent magnitudes of the celestial objects.
    Outputs:
        params -> [float] The estimated magnitude constant C.
        params_err -> [float] The standard error (uncertainty) of the estimated C.
    """

    # Calculate the magnitude constant C for each observation.
    # According to the photometric model: M = C - 2.5 * log10(F)
    C = M + 2.5 * np.log10(F)

    # Prepare the design matrix X for regression.
    # Since we are estimating a constant C (intercept only), X is an array of ones.
    X = np.ones_like(F)[:, None]

    # Create a Robust Linear Model (RLM) using statsmodels.
    # The RLM is used instead of ordinary least squares to reduce the influence of outliers.
    rlm_model = sm.RLM(C, X)

    # Fit the model to estimate the parameter(s).
    rlm_results = rlm_model.fit()

    # Extract the estimated parameter (C) and its standard error (uncertainty).
    params, = rlm_results.params
    params_err, = rlm_results.bse

    return params, params_err

def radec_res_rms(wcs,xy,catalog_df):
    """
    1. Given the WCS transformation, convert the pixel coordinates of sources to celestial coordinates.
    2. Compare cthe calculated celestial coordinates with those in star catalog, and compute the residuals and RMS of Ra and Dec components respectively.

    Usage:
        >>> radec_res,radec_rms = radec_res_rms(wcs,xy,catalog_df)
    Inputs:
        wcs -> [Object of class WCS] WCS transformation system
        xy -> [2d array (n*2)] Pixel coordinates of sources
        catalog_df -> Star catalog in form of pandas.DataFrame
    Outputs:
        radec_res -> [2d array(n*2)] Residual of Ra and Dec in arcseconds
        radec_rms -> [array(2 elements)] RMS of Ra and Dec in arcseconds
    """
    # Convert the pixel coordinates of stars to celestial coordinates
    ra_estimate,dec_estimate = wcs.pixel_to_world_values(xy[:,0],xy[:,1])
    # Calculate the residual of Ra and Dec components
    radec_res = catalog_df[['ra','dec']].values - np.stack([ra_estimate,dec_estimate]).T
    # Metric correction for residual of Ra components
    radec_res[:,0] *=  np.cos(np.deg2rad(dec_estimate))
    radec_res *= 3600 # Convert degrees to arcseconds 
    # Calculate the RMS of Ra and Dec components
    radec_rms = np.sqrt(np.mean(radec_res**2,axis=0))

    return radec_res,radec_rms

class ResultContainer(object):
    """
    Class ResultContainer.
    Group and package the calculation results.

    Attributes:
        _description -> [str] Description of the results.
        catalog_df -> [pandas.DataFrame] DataFrame containing the matched stars data.
        Additional attributes can be dynamically added from the info dictionary.
    """
    def __init__(self, info):
        """
        Initialize the ResultContainer instance with the provided information dictionary.

        Inputs:
            info -> [dict] Dictionary containing calculation results and their descriptions.
        """
        self.__dict__.update(info)

    def __repr__(self):
        """
        Return a string representation of the ResultContainer instance, showing the description.

        Returns:
            str : Formatted string with key attributes of the ResultContainer instance.
        """
        mag_rms_str = "mag_rms = {:.2f}".format(self.mag_rms) if hasattr(self, 'mag_rms') else ""
        return "<ResultContainer object: {:s} xy_rms = [{:.2f}, {:.2f}] radec_rms = [{:.4e}, {:.4e}] {:s}>".format(
            self._description,*self.xy_rms,*self.radec_rms,mag_rms_str)

    def to_csv(self,path_res='csv/starmatch.csv'):
        """
        Save data frame of the calculation results to a csv-formatted file.

        Usage:
            >>> path_res = sources.affined_results.to_csv(path_res)
            >>> path_res = sources.matched_results.to_csv(path_res)
            >>> path_res = sources.calibrated_results.to_csv(path_res)
        Inputs:
            path_res -> [str,optional,default=None] Path to save the csv-formatted file
        Outputs:
            path_res -> [str] Path of the csv-formatted file
        """
        # Ensure the directory exists
        Path(path_res).parent.mkdir(parents=True, exist_ok=True)

        # Save the dataframe to a csv-formatted file
        self.catalog_df.to_csv(path_res)

class StarMatch(object):
    """
    Class StarMatch.
    Generate an instance of class Sources as an entrance to star map matching and astronomical calibration.
    """
    def from_sources(xy_raw,camera_params,flux_raw=None,mode_invariants='triangles',distortion=None):
        """
        Generate an instance of the class Sources as an entry point to star map matching and astronomical calibration.

        Usage:
            >>> # Example 1: No distortion is considered
            >>> from starmatch import StarMatch
            >>> # Configure the FOV[deg], pixel width[deg], resolution of the camera
            >>> camera_params = {'fov':(2,2),'pixel_width':0.002,'res':(1024,1024)}
            >>> # We use the top 30 brightest sources to compute the triangle or quad geometric invariants.
            >>> # For the case where distortion pre-correction is not considered
            >>> sources1 = StarMatch.from_sources(xy_raw,camera_params,flux_raw=flux,mode_invariants='triangles')
            >>> # Example 2: 'Brown–Conrady' distortion is considered
            >>> from starmatch.classes import Distortion
            >>> model = 'Brown–Conrady'
            >>> coeffs = [[-1e-4,1e-4],[1e-3,1e-3,1e-4,1e-5]]
            >>> dc = [0.1,0.1]
            >>> distortion_scale = 128
            >>> distortion = Distortion(model,coeffs,dc,distortion_scale)
            >>> sources2 = StarMatch.from_sources(xy_raw,camera_params,flux_raw=flux,mode_invariants='triangles',distortion=distortion)
        Inputs:
            xy_raw -> [2d array] Pixel coordinates of sources
            camera_params -> [dict] The necessary parameters of the camera, such as {'fov':(2,2),'pixel_width':0.02,'res':(1024,1024)}
            where 'fov' and 'pixel_width' are in [deg], and 'res' represents the resolution of the camera.
            flux_raw -> [array,optional,default=None] Flux(Grayscale value) of sources. If None, skip the calculation of point source apparent magnitude.
            mode_invariants -> [str] Mode of geometric invariants to use. Available options are 'triangles' or 'quads'.
            distortion -> [Object of class Distortion, optional, default=None] Distortion model to use. If None, no distortion is applied.
        Outputs:
            sources -> Instance of class Sources, which includes the following attributes:
                xy_raw -> [2d array, n*2] Pixel coordinates of sources
                flux_raw -> [array] Flux(Grayscale value) of sources
                xy -> [2d array, n*2] Pixel coordinates of sources truncated by max_control_points
                flux -> [array] Flux(Grayscale value) of sources truncated by max_control_points
                invariants -> [2d array, n*2] Array of (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of sources.
                asterisms -> [2d array n*3] Array of the indices of sources that correspond to each invariant triangle.
                kdtree -> [Instance of class scipy.spatial.KDTree] 2d-tree for quick nearest-neighbor lookup
                max_control_points -> [int] Maximum number of sources used to execute the star map matching
                _fov -> [2-ele tuple] Camera field of view in [deg], such as (2,2)
                _pixel_width -> [float] Camera pixel size in [deg]
                _res -> [tuple of int] Camera resolution
        Note:
            The sources need to be sorted in descending order according to their fulx in advance.
        """
        # Extract camera parameters
        fov = camera_params.get('fov', None)
        pixel_width = camera_params.get('pixel_width', None)
        res = camera_params['res']

        # Apply distortion correction if a distortion model is provided
        if distortion is not None: xy_raw = distortion.apply(xy_raw,1)

        # Truncate the number of control points if necessary
        n = len(xy_raw)
        max_control_points = min(MAX_CONTROL_POINTS, n)
        xy = xy_raw[:max_control_points]

        if flux_raw is not None:
            flux = flux_raw[:max_control_points]
        else:
            flux = flux_raw

        # Calculate geometric invariants and construct KDTree
        invariants,asterisms,kdtree = calculate_invariantfeatures(xy,mode_invariants)

        if mode_invariants == 'triangles':
            min_matches = 6
        elif mode_invariants == 'quads':
            min_matches = 4

        # Create a dictionary of source information
        info = {
            'xy_raw': xy_raw,
            'flux_raw': flux_raw,
            'xy': xy,
            'flux': flux,
            'invariants': invariants,
            'asterisms': asterisms,
            'kdtree': kdtree,
            'max_control_points': max_control_points,
            '_fov': fov,
            '_pixel_width': pixel_width,
            '_res': np.array(res),
            '_mode_invariants': mode_invariants,
            '_min_matches': min_matches,
            '_pixel_tols': PIXEL_TOLS
        }

        return Sources(info)      

class Sources(object):
    """
    Class Sources represents a collection of star sources with their attributes and methods for star map matching and astronomical calibration.
    It includes pixel coordinates, flux values, geometric invariants, and methods for computing these invariants and estimating the camera's center pointing.

    Attributes:
        xy_raw -> [2d array (n, 2)] Pixel coordinates of sources.
        flux_raw -> [array] Flux (grayscale value) of sources.
        xy -> [2d array (n, 2)] Pixel coordinates of sources truncated by max_control_points.
        flux -> [array] Flux (grayscale value) of sources truncated by max_control_points.
        invariants -> [2d array (m, 2) or (m, 4)] Geometric invariants of sources.
        Every three sources form a 2-component invariant, and every four sources form a 4-component invariant.
        _mode_invariants -> [str] Mode of invariant of sources. 'triangles' and 'quads' are available.
        asterisms -> [2d array (m, 3) or (m, 4)] Indices of sources that correspond to each invariant triangle.
        kdtree -> [Object KDTree] 2D-tree of triangle invariants or 4D-tree of quad invariants for quick nearest-neighbor lookup.
        max_control_points -> [int] Maximum number of sources used to execute the star map matching.
        _fov -> [2-ele tuple] Field of view of the camera in degrees, such as (2,2).
        _pixel_width -> [float] Pixel width of the camera in degrees.
        _res -> [tuple of int] Resolution of the camera, such as [1024,1024].
        _wcs -> [Object WCS] World Coordinate System that defines the projection relationship between pixel coordinates and celestial coordinates.
        fp_radec_affine -> [2-ele array] Center pointing for 4-parameters affine transformations(similarity transformations).
        affine_matrix -> [2d array (3, 3)] Affine matrix in the form of [[λcosα, λsinα, λtx], [-λsinα, λcosα, λty], [0, 0, 1]].
        _affine_translation -> [2-ele array] Affine translation params in the form of [tx, ty].
        _affine_rotation -> [float] Rotation angle of affine transformations in radians.
        _affine_scale -> [float] Scale coefficient of affine transformations.
        _L -> [int] Normalized length scale, used for normalizing pixel coordinates.
        affined_results -> [Object ResultContainer] Affined results, where a small number of sources are used for initiating star map matching.
        The results include
            - xy -> [2d array (n, 2)] Affined pixel coordinates of sources.
            - xy_res -> [2d array (n, 2)] Residuals of xy.
            - xy_rms -> [2-ele array] RMS of xy.
            - mag_res -> [float] Residuals of magnitudes of sources.
            - mag_rms -> [float] RMS of magnitudes.
            - C -> [float] Magnitudes constant.
            - C_sigma -> [float] Uncertainty of magnitudes constant.
            - catalog_df -> [pandas.DataFrame] DataFrame of matched stars.
            - _description -> [str] Results description.
            - pixels_camera_match -> [2d array (n, 2)] Pixel coordinates of matched sources.
            - radec_res -> [2d array (n, 2)] Residuals of celestial coordinates.
            - radec_rms -> [2-ele array] RMS of celestial coordinates.
        matched_results -> [Object ResultContainer] Matched results, where a large number of sources are used for enhancing star map matching, similar to affined_results.
        calibrated_results -> [Object ResultContainer] Calibrated results, where distortion is corrected based on the matched_results, similar to affined_results..

    Methods:
        invariantfeatures:
            Calculate geometric invariants of sources and construct 2D-Tree from these invariants.
        center_pointing:
            Obtain the center pointing of the camera through blind matching of star maps with multi-core parallel computing.
        align:
            Find the mapping model between the sources and the stars.
        apply:
            Apply the mapping model to unknown sources.
        fp_calibrate:
            Calibrate the orientation of the camera center.
        show_distortion:
            Show the distortion of the camera.
        show_starmatch:
            Mark the matching stars on the original image.
    """

    def __init__(self, info):
        """
        Initialize the Sources instance with the provided information dictionary.

        Inputs:
            info -> [dict] Dictionary containing source attributes and their values.
        """
        self.__dict__.update(info)

    def __repr__(self):
        """
        Return a string representation of the Sources instance, including key attributes.

        Returns:
            str : Formatted string with key attributes of the Sources instance.
        """
        return "<Sources object: max_control_points = {:d}, res = {:}>".format(self.max_control_points,self._res)

    def invariantfeatures(self,max_control_points=None):
        """
        Computes geometric invariant features for a set of source points by generating unique triangles or quads.
        These features are then used to build a KDTree for efficient matching.

        Usage:
            >>> new_sources = sources.invariantfeatures(max_control_points)
        Inputs:
            max_control_points -> [int,optional,default=None] Maximum number of sources used to calculate invariant features.
            If None, use all sources.
        Outputs:
            Updated Object Sources
        """
        info = self.__dict__.copy()
        n = len(self.xy_raw)
        if max_control_points is None or max_control_points > n: max_control_points = n
        if self.max_control_points != max_control_points:  
            xy = self.xy_raw[:max_control_points]
            invariants,asterisms,kdtree = calculate_invariantfeatures(xy,self._mode_invariants)
            info.update({'xy':xy,'invariants':invariants,'asterisms':asterisms,'kdtree':kdtree,'max_control_points':max_control_points})
        return Sources(info)

    def center_pointing(self,sc_simplified_hashed):
        """
        Estimate the center pointing of the camera through blind matching over star maps with the multi-core parallel computing.

        Usage:
            >>> k_min = 1 # Minimum HEALPix hierarchy level.
            >>> k_max = 6 # Maximum HEALPix hierarchy level.
            >>> mode_invariants = 'triangles' # alternative is 'quads'
            >>> # Generate a h5-formatted star catalog geometric invariants hashed file.
            >>> simplified_catalog.h5_hashes(k_min,k_max,mode_invariants) # Traverse from level 'K1' to level 'K11'
            >>> # Read the hashed file
            >>> simplified_catalog.read_h5_hashes()
            >>> fp_radec,pixel_width_estimate,fov_estimate = sources.center_pointing(simplified_catalog)
        Inputs:
            sc_simplified_hashed -> [H5HashesData] An instance of H5HashesData containing the geometric invariants data.
        Outputs:
            fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in [deg]  
            pixel_width_estimate -> [float] Pixel width of camera in [deg]
            fov_estimate -> [2-ele tuple] FOV of camera in [deg]
        """
        # Check the mode of invariant features for both sources and star catalogs.
        mode_invariants_sources = self._mode_invariants
        mode_invariants_catalogs = sc_simplified_hashed.mode_invariants
        if mode_invariants_sources != mode_invariants_catalogs:
            raise Exception("The mode of the invariant feature of ss is '{mode_invariants_sources}', while that of star catalog is '{mode_invariants_catalogs}'. The two are inconsistent.")

        hashed_data = sc_simplified_hashed.hashed_data
        simplified_catalog = sc_simplified_hashed.sc_simplified

        fp_radec,pixel_width_estimate = get_orientation_mp(self.xy,self.asterisms,self.kdtree,self._fov,self._res,self._mode_invariants,self._pixel_tols,self._min_matches,simplified_catalog,hashed_data)

        fov_estimate = pixel_width_estimate * self._res
        self._pixel_width = pixel_width_estimate
        self._fov = fov_estimate

        return fp_radec,pixel_width_estimate,fov_estimate

    def align(self,fp_radec,simplified_catalog,L=150,distortion_calibrate=None,astrometry_corrections={},outlier_remove='lowess'):
        """
        Given the approximate center pointing, find the mapping model between the sources in image and the stars in catalogs.

        Usage:
            >>> fp_radec = [141.8,-2] # The approximate center pointing [Ra,Dec] in [deg]
            >>> astrometry_corrections = {'t':'2019-02-26T20:11:14.347','proper-motion':None,'aberration':(0.55952273, -1.17780654,  7.50324956),'parallax':None}
            >>> sources.align(fp_radec,simplified_catalog,distortion_calibrate='gpr',astrometry_corrections=astrometry_corrections)
        Inputs:
            fp_radec -> [tuple of float] Approximate center pointing of in form of [Ra,Dec] in deg
            simplified_catalog -> [Object of class StarCatalogSimplified] A basic catalog created from the reduced catalog by applying magnitude truncation and proper motion correction, suitable for quick lookups.
            L -> [int,optional,default=150] The number of pixels in a unit length. It controls the tolerance of the 3D-Tree composed of (x,y,mag) for sources of camera and catalog
            For example, if we set the tolerance to 0.2 and L to 150(default), it means the difference within 0.2*150=30 for pixel coordinates and 0.2 for magnitude is the correct matching,
            distortion_calibrate -> [str,optional,default=None] If not None, the distortion correction will be employed. Available options are
                - 'gpr': nonparametric Gaussian Process Regression(GPR).
                - 'piecewise-affine': The transform is based on a Delaunay triangulation of the points to form a mesh. Each triangle is used to find a local affine transform.
                - 'polynomial': 2D polynomial transformation with the following form

                                X = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i ))
                                Y = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i ))

            astrometry_corrections -> [dict] Dictionary specifying the types of corrections to apply.
                - 't' -> [str] Observation time in UTC, such as '2019-02-26T20:11:14.347'.
                   It specifies the time at which corrections are to be applied.
                - 'proper-motion' -> [None] If present, apply proper motion correction.
                   This term corrects for the motion of stars across the sky due to their velocities.
                - 'aberration' -> [tuple] Aberration correction parameters. Observer's velocity relative to Earth's center (vx, vy, vz) in km/s.
                   This term corrects for the apparent shift in star positions due to the motion of the observer.
                - 'parallax' -> [None] If present, apply parallax correction.
                   This term corrects for the apparent shift in star positions due to the change in observer's viewpoint as the Earth orbits the Sun.
                - 'deflection' -> [None] If present, apply light deflection correction.
                   This term corrects for the bending of light from stars due to the gravitational field of the Sun, based on general relativity.
            outlier_remove -> [str] Method of outlier removal. Available options are:
                - 'lowess' -> Identifies outliers with the method of LOWESS (Locally Weighted Scatterplot Smoothing). Here, LOWESS uses a weighted **linear regression** by default.
                - 'iqr' -> Identifies outliers with the method of Interquartile Range (IQR).
        Outputs:
            self : Updated instance with alignment and calibration results.
        """
        fov,pixel_width,res = self._fov,self._pixel_width,self._res
        pixel_tols, min_matches = self._pixel_tols, self._min_matches
        fov_min,fov_max = min(fov),max(fov)
        search_radius = 1.06*fov_max
        pixels_camera,flux_camera = self.xy,self.flux

        # Query Star Catalog around the fiducial point.
        stars = simplified_catalog.search_cone(fp_radec,search_radius,fov_min,max_num_per_tile=MAX_NUM_PER_TILE,astrometry_corrections=astrometry_corrections)
        stars.pixel_xy(pixel_width) # Calculate the pixel coordinates of stars
        stars.invariantfeatures(self._mode_invariants) # Calculate the triangle invariants and constructs a 2D Tree of stars; and records the asterism indices for each triangle.
        wcs = stars.wcs # Object of WCS transformation

        # Align sources from the camera and from the star catalog
        camera_tuple = (self.xy,self.asterisms,self.kdtree)
        catalog_tuple = (stars.xy,stars.asterisms,stars.kdtree)
        transf, (pixels_camera_match, pixels_catalog_match),_s,_d = find_transform_tree(camera_tuple,catalog_tuple,pixel_tols[0],min_matches)

        # Roughly calibrate the center pointing of the camera
        pixels_cc_affine = matrix_transform([0,0],transf.params)
        fp_radec_affine = wcs.pixel_to_world_values(pixels_cc_affine[:,0],pixels_cc_affine[:,1])
        fp_radec_affine = np.hstack(fp_radec_affine)

        # Re-calculate the affine transform by the updated center pointing of the camera
        if norm(pixels_cc_affine) > min(res)/10:
            stars = simplified_catalog.search_cone(fp_radec_affine,search_radius,fov_min,max_num_per_tile=MAX_NUM_PER_TILE,astrometry_corrections=astrometry_corrections)
        else:
            stars.center = fp_radec_affine

        stars.pixel_xy(pixel_width) 
        stars.invariantfeatures(self._mode_invariants)
        catalog_df = stars.df
        wcs = stars.wcs

        catalog_tuple = (stars.xy,stars.asterisms,stars.kdtree)
        transf, (pixels_camera_match, pixels_catalog_match),_s,_d = find_transform_tree(camera_tuple,catalog_tuple,pixel_tols[1],min_matches*2)
        affine_matrix = transf.params
        affine_translation = transf.translation
        affine_rotation = transf.rotation # in radians
        affine_scale = np.sqrt(det(affine_matrix[:2, :2])) # For similarity transform, it is equal to transf.scale

        _, ind_catalog_match, ind_catalog = np.intersect1d(pixels_catalog_match[:,0],stars.xy[:,0],return_indices=True)
        pixels_camera_match = pixels_camera_match[ind_catalog_match]
        pixels_catalog_match = pixels_catalog_match[ind_catalog_match]
        catalog_df_affine = catalog_df.loc[ind_catalog]

        xy_affine = pixels_camera_affine = matrix_transform(pixels_camera_match,affine_matrix)
        xy_res_affine = pixels_catalog_match - pixels_camera_affine
        xy_rms_affine = np.sqrt(np.mean(xy_res_affine**2,axis=0))
        catalog_df_affine[['x_camera', 'y_camera']] = pixels_camera_match
        catalog_df_affine[['dx','dy']] = xy_res_affine
        # Calculate the celestial coordinates of sources, and compare with the star catalog
        # Calculate the residual and RMS of the Ra and Dec components respectively.
        radec_res_affine,radec_rms_affine = radec_res_rms(wcs,xy_affine,catalog_df_affine)
        catalog_df_affine[['dRa', 'dDec']] = radec_res_affine

        info_affined_results_photometric = {}
        if flux_camera is not None:
            flux_camera_match = flux_camera[_s][ind_catalog_match]

            # Parameters of the photometric model (magnitude constant) are fitted using the robust linear least squares method.
            C_affine, C_sigma_affine = photometric_robust_linear_fit(flux_camera_match, catalog_df_affine['mag'])
            
            # Calculate the magnitude residual
            mag_res_affine = catalog_df_affine['mag'].values - photometric_model(flux_camera_match, C_affine)
            # Calculate RMS of the magnitude residual
            mag_rms_affine = np.sqrt(np.mean(mag_res_affine ** 2))

            catalog_df_affine['dmag'] = mag_res_affine
            dict_values = mag_res_affine, mag_rms_affine, C_affine, C_sigma_affine
            dict_keys = 'mag_res', 'mag_rms', 'C', 'C_sigma'
            info_affined_results_photometric = dict(zip(dict_keys, dict_values))

        # Group and package the calculation results for affined sources(pixel coordinates and flux).
        description = 'affined results'
        dict_values = xy_affine,xy_res_affine,xy_rms_affine,catalog_df_affine,description,pixels_camera_match,radec_res_affine,radec_rms_affine
        dict_keys = 'xy','xy_res','xy_rms','catalog_df','_description','pixels_camera_match','radec_res','radec_rms'
        info_affined_results_astrometric = dict(zip(dict_keys, dict_values))
        info_affined_results = info_affined_results_astrometric | info_affined_results_photometric
        affined_results = ResultContainer(info_affined_results)

        # Basic results for affine transformation
        dict_values = affined_results,wcs,fp_radec_affine,affine_matrix,affine_translation,affine_rotation,affine_scale,L,pixel_width,fov
        dict_keys = 'affined_results','_wcs','fp_radec_affine','affine_matrix','_affine_translation','_affine_rotation','_affine_scale','_L','_pixel_width','_fov'
        info_update = dict(zip(dict_keys, dict_values))
        self.__dict__.update(info_update)

        # This part replaces the sources of the affine matching with the sources of the 3D-Tree matching and performs calculations similar to the previous part.
        # Apply the affine matrix and the magnitude constant to all sources in camera image, then build a dimensionless 3D-Tree for camera and starcatalog 
        pixels_camera_affine = matrix_transform(self.xy_raw,affine_matrix)
        stars = simplified_catalog.search_cone(fp_radec_affine, search_radius / 1.5, fov_min/2, max_num_per_tile=MAX_NUM_PER_TILE,astrometry_corrections=astrometry_corrections)
        stars.pixel_xy(pixel_width) 
        catalog_df = stars.df

        if flux_camera is not None:
            source_xymag = np.hstack([pixels_camera_affine / L, photometric_model(self.flux_raw,C_affine)[:, None]])
            catalog_xymag = np.hstack([stars.xy/L,stars.mag[:,None]])
        else:
            source_xymag = pixels_camera_affine / L
            catalog_xymag = stars.xy/L

        source_xymag_tree = KDTree(source_xymag)
        catalog_xymag_tree = KDTree(catalog_xymag)
        matches_list = source_xymag_tree.query_ball_tree(catalog_xymag_tree, r=0.2)

        # Initialize the matching lists
        matches, matches_index_source, matches_index_catalog = [], [], []

        # Find the matching pairs using enumerate and list comprehension
        for i, (t1, t2_list) in enumerate(zip(source_xymag, matches_list)):
            if len(t2_list) == 1:
                matches_index_source.append(i)
                matches_index_catalog.append(t2_list[0])
                t2 = catalog_xymag_tree.data[t2_list[0]]
                matches.append(list(zip(t1, t2)))

        pixels_camera_match = self.xy_raw[matches_index_source]
        catalog_df = catalog_df.loc[matches_index_catalog]
        matches = np.array(matches)
        matches_source,matches_catalog = matches[:,:,0],matches[:, :, 1]
        matches_res = matches_catalog - matches_source

        xy_L_source = matches_source[:,:2]
        xy_L_catalog = matches_catalog[:,:2]
        xy_res_L = matches_res[:,:2]

        xy_source = xy_L_source * L
        xy_catalog = xy_L_catalog * L
        xy_res = xy_res_L * L

        # Detecting anomalies caused by mismatches
        if outlier_remove == 'lowess':
            flag_outliers = lowess_smooth(xy_L_source, xy_res)
        elif outlier_remove == 'iqr':
            flag_outliers = iqr_outliers(xy_res)
        else:
            raise ValueError('outlier_remove must be either "lowess" or "iqr".')
        flag_inliers = ~flag_outliers

        # Remove outliers
        xy_res = xy_res[flag_inliers]
        catalog_df = catalog_df[flag_inliers]
        pixels_camera_match = pixels_camera_match[flag_inliers]
        xy_source = xy_source[flag_inliers]
        xy_catalog = xy_catalog[flag_inliers]

        xy_rms = np.sqrt(np.mean(xy_res**2,axis=0))
        catalog_df[['x_camera', 'y_camera']] = pixels_camera_match
        catalog_df[['dx','dy']] = xy_res
        radec_res,radec_rms = radec_res_rms(wcs,xy_source,catalog_df)
        catalog_df[['dRa', 'dDec']] = radec_res

        info_matched_results_photometric = {}
        if flux_camera is not None:
            # The parameters of the photometric model (magnitude constant) were fitted using the least squares method.
            C, C_var = curve_fit(photometric_model,self.flux_raw[matches_index_source][flag_inliers],catalog_df['mag'],p0=[C_affine])
            # Calculate the magnitude residual
            mag_res = catalog_df['mag'].values - photometric_model(self.flux_raw[matches_index_source][flag_inliers], C)
            # Calculate RMS of the magnitude residual
            mag_rms = np.sqrt(np.mean(mag_res ** 2))
            # Calculate the uncertainty in magnitude constant
            C_sigma = np.sqrt(C_var.item())

            catalog_df['dmag'] = mag_res
            dict_values = mag_res, mag_rms, C, C_sigma
            dict_keys = 'mag_res', 'mag_rms', 'C', 'C_sigma'
            info_matched_results_photometric = dict(zip(dict_keys, dict_values))

        # Group and package the calculation results for matched sources(pixel coordinates and flux).
        description = 'matched results'
        dict_values = xy_source,xy_res,xy_rms,catalog_df,description,pixels_camera_match,radec_res,radec_rms
        dict_keys = 'xy','xy_res','xy_rms','catalog_df','_description','pixels_camera_match','radec_res','radec_rms'
        info_matched_results_astrometric = dict(zip(dict_keys, dict_values))
        info_matched_results = info_matched_results_astrometric | info_matched_results_photometric
        matched_results = ResultContainer(info_matched_results)

        pixel_width_estimate = pixel_width*affine_scale
        fov_estimate = pixel_width_estimate * self._res

        dict_values = matched_results,pixel_width_estimate,fov_estimate
        dict_keys = 'matched_results','pixel_width_estimate','fov_estimate'
        info_update = dict(zip(dict_keys, dict_values))
        self.__dict__.update(info_update)

        xy_calibrated = xy_source
        xy_res_calibrated = xy_res
        xy_rms_calibrated = xy_rms
        catalog_df_calibrated = catalog_df
        radec_res_calibrated = radec_res
        radec_rms_calibrated = radec_rms
        info_calibrated_results_photometric = info_matched_results_photometric

        # Distortion correction
        if distortion_calibrate is not None:
            if distortion_calibrate in ['piecewise-affine','polynomial']:
                if distortion_calibrate == 'piecewise-affine':
                    tform = PiecewiseAffineTransform()
                elif distortion_calibrate == 'polynomial':
                    tform = PolynomialTransform()
                tform.estimate(xy_source, xy_catalog)
                xy_calibrated = tform(xy_source)
                xy_res_calibrated = xy_catalog - xy_calibrated
                xy_rms_calibrated = np.sqrt(np.mean(xy_res_calibrated**2,axis=0))

                catalog_df_calibrated = catalog_df.copy()
                catalog_df_calibrated[['dx', 'dy']] = xy_res_calibrated
                radec_res_calibrated,radec_rms_calibrated = radec_res_rms(wcs,xy_calibrated,catalog_df_calibrated)
                catalog_df_calibrated[['dRa', 'dDec']] = radec_res_calibrated

                self.__dict__.update({'_tform':tform})

            elif distortion_calibrate == 'gpr':
                L_GPR = 512  # Used to scale coordinates to avoid large squared distances between points
                xy_L = xy_source/L_GPR # Normalize the coordinates
                U,V = xy_res[None,:].T

                kerx = RBF(2) # Define 2D kernel for distortion along x axis
                mx_L = GPRegression(xy_L,U,kerx) # Create simple GP model
                mx_L.optimize() # Calculate the maximum likelihood solution of hyperparameters

                kery = RBF(2) # Define 2D kernel for distortion along y axis
                my_L = GPRegression(xy_L,V,kery)
                my_L.optimize()

                meanx_L,varx_L = mx_L.predict(xy_L) # fitted distortion along x axis
                meany_L,vary_L = my_L.predict(xy_L) # fitted distortion along y axis

                distortion_fitted = np.hstack([meanx_L,meany_L])

                xy_calibrated = xy_source + distortion_fitted
                xy_res_calibrated = xy_res - distortion_fitted
                xy_rms_calibrated = np.sqrt(np.mean(xy_res_calibrated**2,axis=0))

                catalog_df_calibrated = catalog_df.copy()
                catalog_df_calibrated[['dx','dy']] = xy_res_calibrated
                radec_res_calibrated,radec_rms_calibrated = radec_res_rms(wcs,xy_calibrated,catalog_df_calibrated)

                self.__dict__.update({'_mx_L':mx_L,'_my_L':my_L,'_L_GPR':L_GPR})
            else:
                raise Exception("Unrecognized calibration method. Only 'gpr', 'piecewise-affine', and 'polynomial' are supported.")

        info_calibrated_results_astrometric = {
            'xy': xy_calibrated,
            'xy_res': xy_res_calibrated,
            'xy_rms': xy_rms_calibrated,
            'catalog_df': catalog_df_calibrated,
            'radec_res': radec_res_calibrated,
            'radec_rms': radec_rms_calibrated,
            '_description': f"calibrated results method = '{distortion_calibrate}'",
            '_method': distortion_calibrate
        }
        info_calibrated_results = info_calibrated_results_astrometric | info_calibrated_results_photometric
        self.calibrated_results = ResultContainer(info_calibrated_results)   

    def apply(self,xy_target,flux_target=None):
        """
        Apply the mapping model to unknown sources.

        Usage:
            >>> radec1,M1 = sources1.apply([-400,-400],5e3) 
            >>> radec2,M2 = sources2.apply([[-400,-400],[500,500]],[1e4,5e3])
        Inputs:
            xy_target -> [array] Pixel coordinates of unknown sources
            flux_target -> [array] Flux(Grayscale value) of unknown sources
        Outputs:
            radec -> [array] Estimated celestial coordinates (Ra, Dec) of the target sources.
            M_affine -> [array, optional] Apparent magnitudes estimated from the affine model.
            M_matched -> [array, optional] Apparent magnitudes estimated from the matched results.
        """
        # Apply the similarity transform (affine matrix)
        xy_target_affine = matrix_transform(xy_target,self.affine_matrix)

        # Apply distortion calibration if available
        calibrated_xy = xy_target_affine  # Default to affine transformation

        if hasattr(self,'calibrated_results'):
            # If a distortion calibration model exists
            if self.calibrated_results._method == 'gpr':
                meanx_L,varx_L = self._mx_L.predict(xy_target_affine/self._L_GPR)
                meany_L,vary_L = self._my_L.predict(xy_target_affine/self._L_GPR)
                mean_xy = np.hstack([meanx_L,meany_L])
                calibrated_xy = xy_target_affine + mean_xy
            elif self.calibrated_results._method in ['piecewise-affine','polynomial']:
                calibrated_xy = self._tform(xy_target_affine)
                
        # Convert pixel coordinates to celestial coordinates
        radec = self._wcs.pixel_to_world_values(calibrated_xy[:,0],calibrated_xy[:,1])
        radec = np.stack(radec).T

        # Photometric
        if flux_target is not None:
            M_affine = photometric_model(flux_target,self.affined_results.C)
            M_matchd = photometric_model(flux_target,self.matched_results.C)
            return radec,M_affine,M_matchd
        else:
            return radec

    def fp_calibrate(self):
        """
        Calibrate the orientation of the camera center

        Usage:
            >>> sources.fp_calibrate()
            >>> print(sources.fp_radec_calibrated)
        Outputs:
            self -> Updated instance with calibrated center pointing.
        """
        fp_radec_calibrated = self.apply([0,0])
        self.fp_radec_calibrated = fp_radec_calibrated[0]

        return self
            
    def show_distortion(self,mode='vector',fig_file=None):
        """
        Show the distortion of the camera.

        Usage:
            >>> sources.show_distortion('vector')
            >>> sources.show_distortion('contourf')
        Inputs:
            mode -> [str,optional,default='vector'] The way of the distortion is displayed. Avaliable options include 'vector' plot and 'contourf' plot.  
            fig_file -> [str] Path to save the distortion map
        outputs:
            Distortion map displayed or saved to file.
        """

        x,y = self.matched_results.xy.T
        u,v = self.matched_results.xy_res.T

        plt.clf()
        ylim,xlim = self._res/2

        # 'vector' plot mode
        if mode == 'vector':
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(x, y,facecolors="None", color='m',s=10,alpha=0.5) 
            plt.quiver(x, y, u, v, color='b', units='xy')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim(-xlim, xlim)
            ax.set_ylim(-ylim, ylim)
            ax.tick_params(axis='both', which='major', labelsize=7)

        # 'contourf' plot mode
        elif mode == 'contourf':
            xx = np.linspace(-xlim, xlim, 50)
            yy = np.linspace(-ylim, ylim, 50)
            X, Y = np.meshgrid(xx, yy)

            U = griddata(self.matched_results.xy, u, (X, Y), method='nearest')
            V = griddata(self.matched_results.xy, v, (X, Y), method='nearest')

            fig = plt.figure(dpi=300)
            if xlim > 1.2*ylim:
                axes = fig.subplots(2,1)
            else:
                axes = fig.subplots(1,2)

            mean = (U,V)
            for i in range(2):
                if mean[i].max() * mean[i].min() < 0:
                    abs_Z_max = np.abs(mean[i]).max()
                    Z_levels = np.linspace(-abs_Z_max, abs_Z_max, 51)
                    cs = axes[i].contourf(X, Y, mean[i], levels=Z_levels, extend='both', cmap=plt.cm.RdBu_r)
                else:
                    cs = axes[i].contourf(X, Y, mean[i], extend='both')

                divider = make_axes_locatable(axes[i])  
                cax = divider.append_axes("right", size="5%", pad=0.05) 
                cb = fig.colorbar(cs, ax=axes[i],cax=cax,format='%4.1f')
                cb.ax.tick_params(labelsize=7)

                axes[i].set_aspect('equal', adjustable='box')
                axes[i].scatter(x, y,facecolors="None", color='m',s=10,alpha=0.5) 
                axes[i].set_xlim(-xlim, xlim)
                axes[i].set_ylim(-ylim, ylim)
                axes[i].tick_params(axis='both', which='major', labelsize=7)

            plt.subplots_adjust(wspace=0.3)
        
        if fig_file is None:
            plt.show() 
        else:
            Path(fig_file).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(fig_file,bbox_inches='tight')      

    def show_starmatch(self,image_raw,offset,fig_out=None):
        """
        Mark the matching stars on the original image

        Usage:
            >>> sources.show_starmatch(image_raw,offset)
        Inputs:
            image_raw -> [2d array of float] Raw grayscale image with origin at bottom left corner point
            offset -> [array of float] Pixel coordinates of the center of image
            fig_out -> [str,optional,default=None] The path to save the output image.
        Outputs:
            starmatch map  
        """
        xy = self.matched_results.pixels_camera_match + offset
        text = self.matched_results.catalog_df.index.astype(str)
        if fig_out is None:
            plot_kwargs = {'mark':(xy,'s','red',text)}
        else:
            plot_kwargs = {'mark':(xy,'s','red',text),'fig_path':fig_out}
        show_image(image_raw,origin='lower',**plot_kwargs)    
        
class Distortion(object):
    """
    Class Distortion

    Attributes:
        - model -> [str] Type of distortion model. Available distortion models include 'RadialStandard', 'RadialDivision', 'Tangential', 'Brown–Conrady'.
        - coeffs -> [list] Coefficients associated to distortion models
        - dc -> [list of two elements] Pixel coordinates of the distortion center
        - distortion_scale -> [int] The length scale of the distortion model, that is, the number of pixels per unit length

    Methods:
        - apply : Compute the distortion-corrected pixel coordinates
        - sketchmap : Sketch the vector plot of distirtion
    """    

    def __init__(self,model,coeffs,dc,distortion_scale=1):  
        """
        Initialize the distortion model

        Usage:
            >>> from starmatch.classes import Distortion
            >>> model = 'Brown–Conrady' # Type of distortion model
            >>> coeffs = [[-1e-4,1e-4],[1e-3,1e-3,1e-4,1e-5]] # Coefficients associated to distortion models
            >>> dc = [0.1,0.1] # Pixel coordinates of the distortion center
            >>> distortion_scale = 128 # The length scale of the distortion model, that is, the number of pixels per unit length
            >>> distortion = Distortion(model,coeffs,dc,distortion_scale) # Establish a distortion model
        Inputs:
            model -> [str] Type of distortion model. Available distortion models include 
                1. 'RadialStandard': Standard Radial Distortion Model(SRDM)
                2. 'RadialDivision': Division-mode Radial Distortion Model(DRDM)
                3. 'Tangential': Tangential Distortion Model(also known as the de-centering distortion)
                4. 'Brown–Conrady': Brown-Conrady Distortion Model(BCDM)
            coeffs -> [list] Coefficients associated to distortion models
            dc -> [list of two elements] Pixel coordinates of the distortion center
            distortion_scale -> [int,optional,default=1] The length scale of the distortion model, that is, the number of pixels per unit length
        Outputs:
            distortion -> Instance of class Distortion

        Note: 
        Considering that the distortion model involves the power of the distance from the distortion center to pixels, the length scale is introduced here. 
        Associaated to the length scale, the normalized model coefficients is also necessary.

        For more details on distortion models, please refer to 
            1. https://en.wikipedia.org/wiki/Distortion_(optics)
            2. https://www.imatest.com/docs/distortion-methods-and-modules/
            3. https://www.imatest.com/support/docs/pre-5-2/geometric-calibration-deprecated/distortion-models/ 
        """
        model_list = ['RadialStandard','RadialDivision','Tangential','Brown–Conrady']
        if model in model_list:
            info = {'model':model,'coeffs':coeffs,'dc':dc,'distortion_scale':distortion_scale}
        else:    
            raise Exception('Distortion model should be one of {:s}'.format(str(model_list)))

        self.__dict__.update(info)

    def __repr__(self):
        """
        Return a string representation of the Distortion instance.

        Returns:
            str : Formatted string with key attributes of the Distortion instance.
        """

        return f"<Distortion object: model='{self.model}', distortion_scale={self.distortion_scale}>"

    def apply(self,pixels_xy,pixel_scale=1):
        """
        Compute the distortion-corrected pixel coordinates.

        Usage:
            >>> # Calculate the distortion-corrected pixel coordinates
            >>> pixels_xy = [367,125]
            >>> pixels_XY = distortion.apply(pixels_xy)
            >>>
            >>> # Calculate the distortion-corrected pixel coordinates at normalized scale
            >>> pixels_xy = [[5.734375, 1.953125], [1.109375, -0.875]]
            >>> pixel_scale = 128
            >>> pixels_XY = distortion.apply(pixels_xy,pixel_scale)
        Inputs:
            pixels_xy -> [list or array of two float elements, or 2d array] Pixel coordinates to calculate distortion
            pixel_scale -> [int,optional,default=1] The length scale of the normalized pixel coordinates, that is, the number of pixels per unit length
        Outputs:
            pixels_XY -> [array of two float elements, or 2d array(n*2)] Distortion-corrected pixel coordinates    
        """
        ratio = pixel_scale/self.distortion_scale

        # Convert the normalized pixel coordinates to the distortion-scaled pixel coordinates
        pixels_xy = np.array(pixels_xy)
        if pixels_xy.ndim == 1: 
            pixels_xy_ = pixels_xy[None,:] * ratio
        else:
            pixels_xy_ = pixels_xy * ratio
                
        if self.model == 'RadialStandard':
            pixels_XY_ = distortion_model.RadialStandard(self.coeffs,self.dc,pixels_xy_)
        elif self.model == 'RadialDivision':   
            pixels_XY_ = distortion_model.RadialDivision(self.coeffs,self.dc,pixels_xy_) 
        elif self.model == 'Tangential':   
            pixels_XY_ = distortion_model.Tangential(self.coeffs,self.dc,pixels_xy_)     
        elif self.model == 'Brown–Conrady':   
            pixels_XY_ = distortion_model.Brown_Conrady(self.coeffs,self.dc,pixels_xy_) 

        # Restore to the normalized pixel coordinates
        pixels_XY = pixels_XY_ / ratio  
        
        return pixels_XY

    def sketchmap(self,xlim,ylim,pixel_scale=1,mode='vector'):   
        """
        Sketch the vector plot or contour plot of distortion.

        Usage:
            >>> xlim,ylim = 512,512
            >>> distortion.sketchmap(xlim,ylim)
            >>> # For normalized pixel coordinates
            >>> xlim,ylim = 2,2
            >>> pixel_scale = 256
            >>> distortion.sketchmap(xlim,ylim,pixel_scale)
        Inputs:
            xlim -> [float] The x-direction boundary of the sketch map
            ylim -> [float] The y-direction boundary of the sketch map
            pixel_scale -> [int,optional,default=1] The length scale of the normalized pixel coordinates, that is, the number of pixels per unit length   
            mode -> [str,optional,default='vector'] The way of the distortion is displayed. Avaliable options include 'vector' plot and 'contourf' plot.   
        """ 
        ratio = pixel_scale/self.distortion_scale

        x = np.linspace(-xlim, xlim, 20)
        y = np.linspace(-ylim, ylim, 20)
            
        # Create meshgrid
        X, Y = np.meshgrid(x, y)
        xy = np.stack([X,Y]).transpose(1,2,0).reshape(-1,2)

        # Convert the normalized pixel coordinates to the distortion-scaled pixel coordinates
        xy_ = xy * ratio
        vx_,vy_ = (xy_ - self.apply(xy_,self.distortion_scale)).T # distortion 

        # Restore to the normalized scale
        vxy = vx,vy = (vx_/ratio).reshape(X.shape),(vy_/ratio).reshape(Y.shape)

        plt.clf()

        if mode == 'vector':
            fig, ax = plt.subplots(dpi=300)
            ax.quiver(X,Y,vx,vy,color="C0") # scale=20
            ax.set_aspect('equal', adjustable='box')
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.tick_params(axis='both', which='major', labelsize=7)

        elif mode == 'contourf':
            
            fig = plt.figure(dpi=300)

            if xlim > 1.2*ylim:
                axes = fig.subplots(2,1)
            else:
                axes = fig.subplots(1,2)
            
            for i in range(2):
                if vxy[i].max()*vxy[i].min() < 0:
                    abs_Z_max = np.abs(vxy[i]).max()     
                    Z_levels = np.linspace(-abs_Z_max,abs_Z_max, 51)
                    cs = axes[i].contourf(X, Y, vxy[i],levels = Z_levels,extend='both',cmap=plt.cm.RdBu_r) 
                else:    
                    cs = axes[i].contourf(X, Y, vxy[i],extend='both') 

                divider = make_axes_locatable(axes[i])  
                cax = divider.append_axes("right", size="5%", pad=0.05)    
                cb = fig.colorbar(cs, ax=axes[i],cax=cax,format='%4.1f')
                cb.ax.tick_params(labelsize=7)
                axes[i].set_aspect('equal', adjustable='box')
                axes[i].set_xlim(-xlim, xlim)
                axes[i].set_ylim(-ylim, ylim)
                axes[i].tick_params(axis='both', which='major', labelsize=7)
            plt.subplots_adjust(wspace=0.4)
        
        plt.show() 

