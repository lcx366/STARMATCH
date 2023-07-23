import numpy as np
from numpy.linalg import norm
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import GPy,copy
from pathlib import Path

from .orientation import get_orientation,get_orientation_mp
from .astroalign import invariantfeatures,find_transform_tree,matrix_transform
from .distortion import distortion_model
from .plot import show_image

def solidangle_ratio(fov,r):
    """
    Calculate the ratio of the solid angles spanned by the square and the cone.

    Usage:
        >>> ratio = Solid_angles_ratio(8,10)
    Inputs:
        fov -> [float] FOV of a camera in [deg] that determines a square
        r -> [float] Angular radius in [deg] that determines a cone
    Outputs:
        ratio -> [float] The ratio of the solid angles spanned by the square and the cone    
    """
    fov_rad = np.deg2rad(fov)
    r_rad = np.deg2rad(r)
    Omega_square = 4*np.arcsin(np.tan(fov_rad/2)**2) # solid angles spanned by the square
    Omega_cone = 4*np.pi*np.sin(r_rad/2)**2 # solid angles spanned by the cone
    return Omega_square/Omega_cone

def radec_res_rms(wcs,xy,catalog_df):
    """
    1. Given the WCS transformation, convert the pixel coordinates of stars to celestial coordinates.
    2. Compare cthe calculated celestial coordinates with those in star catalog, then compute the residual and RMS of Ra and Dec components respectively.

    Usage:
        >>> radec_res,radec_rms = radec_res_rms(wcs,xy,catalog_df)
    Inputs:
        wcs -> Object of class WCS, which defines the WCS transformation
        xy -> [list(2 elements) or 2d array(n*2)] Pixel coordinates of stars
        catalog_df -> Star catalog in form of pandas dataframe
    Outputs:
        radec_res -> [2d array(n*2)] Residual of Ra and Dec components
        radec_rms -> [array(2 elements)] RMS of Ra and Dec components
    """
    # Convert the pixel coordinates of stars to celestial coordinates
    radec_estimate = wcs.pixel_to_world(xy[:,0],xy[:,1])
    ra_estimate,dec_estimate = radec_estimate.ra.deg,radec_estimate.dec.deg
    # Calculate the residual of Ra and Dec components
    radec_res = np.array(catalog_df[['ra','dec']] - np.stack([ra_estimate,dec_estimate]).T)
    # Metric correction for residual of Ra components
    radec_res[:,0] *=  np.cos(np.deg2rad(dec_estimate))
    # Calculate the RMS of Ra and Dec components
    radec_rms = np.sqrt(np.mean(radec_res**2,axis=0))

    return radec_res,radec_rms

class Constructor(object):
    """
    Class Constructor, mainly used to group and package the calculation results.
    """
    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return self._description

    def to_csv(self,path_catalog=None):
        """
        Save catalog_df to a csv file.

        Usage:
            >>> path_catalog = sources.affine_results.to_csv()
            >>> path_catalog = sources.match_results.to_csv(path_catalog)
            >>> path_catalog = sources.calibrate_results.to_csv(path_catalog)

        Inputs:
            path_catalog -> [str,optional,default=None] Path to save the csv file

        Outputs:
            path_catalog -> [str] Path of the csv file    
        """
        df = self.catalog_df

        if path_catalog is None: path_catalog = 'csv/starmatch.csv' 
        Path(path_catalog).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path_catalog) # Save the pandas dataframe to a csv-formatted file    

        return path_catalog

class StarMatch(object):
    """
    Class StarMatch, mainly used to generate an instance of class Sources as an entrance to star map matching and astronomical calibration.
    """

    def from_sources(xy_raw,flux_raw,camera_params,max_control_points=20,distortion=None):
        """
        Generate an instance of class Sources using pixel coordinates and flux(grayscale value) of sources, as is the entrance to star map matching and astronomical calibration.

        Usage:
            >>> # Example 1: No distortion is considered
            >>> from starmatch import StarMatch
            >>> # set the FOV[deg], pixel width[deg], resolution of the camera
            >>> camera_params = {'fov':8,'pixel_width':0.01,'res':(1024,1024)}
            >>> # We use the first 20 points to compute the triangle invariants involved in affine transformation. 
            >>> # Usually these points are those sources that are brightest.
            >>> sources1 = StarMatch.from_sources(xy,flux,camera_params,20) # No distortion corrected 
            >>>
            >>> # Example 2: 'Brown–Conrady' distortion is considered
            >>> from starmatch.classes import Distortion
            >>> from starmatch import StarMatch
            >>> model = 'Brown–Conrady'
            >>> coeffs = [[-1e-4,1e-4],[1e-3,1e-3,1e-4,1e-5]]
            >>> dc = [0.1,0.1]
            >>> distortion_scale = 128
            >>> distortion = Distortion(model,coeffs,dc,distortion_scale)
            >>> sources2 = StarMatch.from_sources(xy,flux,camera_params,20,distortion)
        Inputs:
            xy_raw -> [2d array] Pixel coordinates of sources
            flux_raw -> [array] Flux(Grayscale value) of sources
            camera_params -> [dict] The necessary parameters of the camera, such as {'fov':8,'pixel_width':0.01,'res':(1024,1024)}
            where 'fov' and 'pixel_width' are in [deg], and 'res' represents the resolution of the camera.
            max_control_points -> [int,optional,default=20] Maximum number of sources used to execute the star map matching
            distortion -> Instance of class Distortion, which defines the distortion model
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
                _fov -> [float] Camera field of view in [deg]
                _pixel_width -> [float] Camera pixel size in [deg]
                _res -> [tuple of int] Camera resolution
        Note:
            When truncated by max_control_points, the sources need to be sorted from largest to smallest according to their fulx in advance.        
        """
        fov,pixel_width,res = camera_params['fov'],camera_params['pixel_width'],camera_params['res']

        # Distortion correction
        if distortion is not None: xy_raw = distortion.apply(xy_raw,1)
               
        n = len(xy_raw)
        if max_control_points > n: max_control_points = n
        xy = xy_raw[:max_control_points]
        flux = flux_raw[:max_control_points]
        invariants,asterisms,kdtree = invariantfeatures(xy) # equivalent to invariants,asterisms,kdtree = invariantfeatures(xy_L)

        dict_values = xy_raw,flux_raw,xy,flux,invariants,asterisms,kdtree,max_control_points,fov,pixel_width,np.array(res)
        dict_keys = 'xy_raw','flux_raw','xy','flux','invariants','asterisms','kdtree','max_control_points','_fov','_pixel_width','_res'
        info = dict(zip(dict_keys, dict_values))

        return Sources(info)      

class Sources(object):
    """
    Class Sources

    Attributes:
        - xy_raw -> [2d array(n*2)] Pixel coordinates of sources
        - flux_raw -> [array] Flux(Grayscale value) of sources
        - xy -> [2d array(n*2)] Pixel coordinates of sources truncated by max_control_points
        - flux -> [array] Flux(Grayscale value) of sources truncated by max_control_points
        - invariants -> [2d array(m*2)] Array of (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of sources.
        - asterisms -> [2d array(m*3)] Array of the indices of sources that correspond to invariant triangles.
        - kdtree -> [Instance of class scipy.spatial.KDTree] 2d-tree of triangel invariants for quick nearest-neighbor lookup
        - max_control_points -> [int] Maximum number of sources used to execute the star map matching
        - _fov -> [float] Field of View of camera in [deg]
        - _pixel_width -> [float] Pixel width of camera in [deg]
        - _res -> [tuple of int] Resolution of camera
        - _wcs -> object of class WCS, which defines the projection relationship between pixel coordinates and celestial coordinates
        - fp_radec_affine -> [array] Center pointing for affine transformations
        - affine_matrix -> [2d array] Affine matrix in form of [[λcosα,λsinα,λtx],[−λsinα,λcosα,λty],[0,0,1]]
        - _affine_translation -> [array] Affine translation in form of [tx,ty]
        - _affine_rotation -> [float] Rotation angle of affine transformations α in [rad]
        - _affine_scale -> [float] Scale coefficient of affine transformations λ
        - _L -> [int] The normalized length scale. 
        - affine_results -> Instance of class Constructor, whose attributes include 
            xy -> Affined pixel coordinates of affine-matched sources
            xy_res -> Residuals of xy
            xy_rms -> RMS of xy
            mag_res -> Residuals of magnitudes for affine-matched sources
            mag_rms -> RMS of magnitudes for affine-matched sources
            C -> Magtitudes Constant for affine-matched sources
            C_sigma -> Uncertainty of magtitudes Constant
            catalog_df -> Pandas dataframe of matched stars for affine-matched sources
            _description -> results description
            pixels_camera_match -> Pixel coordinates of affine-matched sources
            radec_res -> Residuals of celestial coordinates for affine-matched sources
            radec_rms -> RMS of celestial coordinates for affine-matched sources
        - match_results -> Instance of class Constructor, whose attributes include 
            xy ->  Affined pixel coordinates of 3D-Tree-matched sources
            xy_res -> Residuals of xy
            xy_rms -> RMS of xy
            mag_res -> Residuals of magnitudes for 3D-Tree-matched sources
            mag_rms -> RMS of magnitudes for 3D-Tree-matched sources
            C -> Magtitudes Constant for 3D-Tree-matched sources
            C_sigma -> Uncertainty of magtitudes Constant
            catalog_df -> Pandas dataframe of matched stars for 3D-Tree-matched sources
            _description -> results description
            pixels_camera_match -> Pixel coordinates of 3D-Tree-matched sources
            radec_res -> Residuals of celestial coordinates for 3D-Tree-matched sources
            radec_rms -> RMS of celestial coordinates for 3D-Tree-matched sources
        - _mx_L -> Normalized GPR model in distortion of the x component
        - _my_L -> Normalized GPR model in distortion of the y component  
        - calibrated_results -> Instance of class Constructor, whose attributes include 
            xy -> Affined pixel coordinates for distortion-corrected sources
            xy_res -> Residuals of xy
            xy_rms -> RMS of xy
            mag_res -> Same to that of match_results
            mag_rms -> Same to that of match_results
            C -> Same to that of match_results
            C_sigma -> Same to that of match_results
            catalog_df -> Pandas dataframe of matched stars for distortion-corrected sources
            _description -> results description
            pixels_camera_match -> Pixel coordinates of distortion-corrected sources
            radec_res -> Residuals of celestial coordinates for distortion-corrected sources
            radec_rms -> RMS of celestial coordinates for distortion-corrected sources
 
    Methods:
        - invariantfeatures : Given the max_control_points, refresh the invariants, and 2D-Tree.   
        - center_pointing : Obtain the center pointing of the camera through blind matching of star maps.
        - center_pointing_mp : Obtain the center pointing of the camera through blind matching of star maps with the multi-core parallel computing.
        - align : Given the approximate center pointing of the camera, find the mapping model between the source image and the star map. 
        - apply : Apply the mapping model(affine model or calibration model) to unknown sources
        - fp_calibrate : Calibrate the orientation of the camera center
        - show_distortion : Show the distortion of the camera in modes of 'vector' or 'contourf'.
        - show_starmatch : Mark the matching stars on the original image
    """   
    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class Sources'

    def invariantfeatures(self,max_control_points):
        """
        Given the max_control_points, refresh the following calculations and update the self:
        1. Calculate the unique invariants (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of centroids.
        2. Construct the 2D Tree from the the unique invariants.
        3. Record an array of the indices of centroids that correspond to each invariant.

        Usage:
            >>> # update the sources
            >>> sources.invariantfeatures(max_control_points)
        Inputs:
            max_control_points -> [int] Maximum number of sources used to re-execute invariantfeatures
        Outputs:
            Updated self
        """
        n = len(self.xy_raw)
        if max_control_points > n: max_control_points = n
        if self.max_control_points != max_control_points:  

            xy = self.xy_raw[:max_control_points]
            invariants,asterisms,kdtree = invariantfeatures(xy)

            self.info.update({'xy':xy,'invariants':invariants,'asterisms':asterisms,'kdtree':kdtree,'max_control_points':max_control_points})
            self.xy,self.invariants,self.asterisms,self.kdtree,self.max_control_points = xy,invariants,asterisms,kdtree,max_control_points

        return self   

    def center_pointing(self,simplified_catalog,max_control_points=30):
        """
        Estimate the center pointing and pixel width of the camera through blind matching over star maps.

        Usage:
            >>> fp_radec,pixel_width = sources1.center_pointing(simplified_catalog)
        Inputs:
            simplified_catalog -> Instance of class StarCatalogSimplified
            max_control_points -> [int,optional,default=30] Maxinum number of the stars sorted by brightness used for a sky area.
        Outputs:
            fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in [deg]  
            pixel_width -> [float] Pixel width of camera in [deg]
        """

        indices_h5,mcp_ratio = simplified_catalog.h5_incices(self._fov,self._pixel_width,max_control_points)
        max_control_points = round(max_control_points*mcp_ratio)
        self.invariantfeatures(max_control_points) 
        fp_radec,pixel_width_estimate = get_orientation(self.xy,self.asterisms,self.kdtree,self._pixel_width,indices_h5)

        return fp_radec,pixel_width_estimate

    def center_pointing_mp(self,simplified_catalog,max_control_points=30):
        """
        Estimate the center pointing of the camera through blind matching over star maps with the multi-core parallel computing.

        Usage:
            >>> fp_radec,pixel_width = sources1.center_pointing_mp(simplified_catalog)
        Inputs:
            simplified_catalog -> Instance of class StarCatalogSimplified
            max_control_points -> [int,optional,default=30] Maxinum number of the stars sorted by brightness used for a sky area.
        Outputs:
            fp_radec -> [tuple of float] Center pointing of the camera in form of [Ra,Dec] in [deg]  
            pixel_width -> [float] Pixel width of camera in [deg]
        """
        indices_h5,mcp_ratio = simplified_catalog.h5_incices(self._fov,self._pixel_width,max_control_points)
        max_control_points = round(max_control_points*mcp_ratio)  
        self.invariantfeatures(max_control_points)
        fp_radec,pixel_width_estimate = get_orientation_mp(self.xy,self.asterisms,self.kdtree,self._pixel_width,indices_h5)

        return fp_radec,pixel_width_estimate

    def align(self,fp_radec,simplified_catalog,L=32,calibrate=False):
        """
        Given the approximate center pointing of the camera, find the mapping model between the source image and the star map. 
        The mapping model includes two types: the affine model without considering the geometric distortion and the affine model + distortion model.
        Usage:
            >>> # Example 1: with no distortion on pixel coordinates of sources
            >>> from starmatch import StarMatch
            >>> camera_params = {'fov':8,'pixel_width':0.01,'res':(1024,1024)}
            >>> sources1 = StarMatch.from_sources(xy,flux,camera_params)
            >>> # Load the simplified star catalog
            >>> from starcatalogquery import StarCatalog
            >>> dir_from_simplified = 'starcatalogs/simplified/hygv3/res5/mag9.0/epoch2022.0/'
            >>> hygv3_simplified = StarCatalog.load('simplified','hygv3',5,dir_from_simplified)
            >>> # star map matching
            >>> fp_radec = [202,31] # [Ra,Dec] in [deg]
            >>> sources1.align(fp_radec,hygv3_simplified)
            >>> # Example 2: with 'Brown–Conrady' distortion on pixel coordinates of sources
            >>> from starmatch.classes import Distortion
            >>> from starmatch import StarMatch
            >>> # Make a distortion correction
            >>> sources2 = StarMatch.from_sources(xy,flux,camera_params,20,distortion)
            >>> sources2.align(fp_radec,hygv3_simplified,L=32,calibrate=True)
        Inputs:
            fp_radec -> [tuple of float] Approximate center pointing of the camera in form of [Ra,Dec] in deg
            simplified_catalog -> [Instance of class StarCatalogSimplified] The minimalist star catalog, which only includes the position and apparent magnitude of stars at a specific epoch
            L -> [int,optional,default=32] The number of pixels in a unit length. It controls the tolerance of the 3D-Tree composed of (x,y,mag) for sources of camera and catalog
            For example, if we set the tolerance to 0.2(default), it means the difference within 0.2*32=6.4 for pixel coordinates and 0.2 for magnitude is the correct matching, 
            Another use is associated with the GPR-baed distortion correction. Experience has shown that normalizing the pixel coordinates by L is beneficial to the operation of GPR.
            calibrate -> [bool,optional,default=False] Whether to perform a distortion correction. If False, no distortion correction is applied. If True, the distortion correction will be employed using the nonparametric Gaussian Process Regression(GPR) method.
        Outputs:
            updated self      
        """
        fov,pixel_width,res = self._fov,self._pixel_width,self._res
        search_radius = 0.75*fov
        ratio = solidangle_ratio(fov,search_radius)

        max_control_points = int(self.max_control_points/ratio)
        pixels_camera,flux_camera = self.xy,self.flux

        # Query Star Catalog around the fiducial point.
        stars = simplified_catalog.search_cone(fp_radec,search_radius,max_control_points)
        stars.pixel_xy(pixel_width) # Calculate the pixel coordinates of stars
        stars.invariantfeatures() # Calculate the triangle invariants and constructs a 2D Tree of stars; and records the asterism indices for each triangle.
        wcs = stars.wcs # Object of WCS transformation

        # Align sources from the camera and from the star catalog
        camera_tuple = (self.xy,self.asterisms,self.kdtree)
        catalog_tuple = (stars.xy,stars.asterisms,stars.kdtree)
        transf, (pixels_camera_match, pixels_catalog_match),_s,_d = find_transform_tree(camera_tuple,catalog_tuple)

        # Roughly calibrate the center pointing of the camera
        pixels_cc_affine = matrix_transform([0,0],transf.params)  
        pixels_cc_affine_x,pixels_cc_affine_y = pixels_cc_affine[:,0],pixels_cc_affine[:,1]
        cc_radec_estimate = wcs.pixel_to_world(pixels_cc_affine_x,pixels_cc_affine_y)

        # Re-calculate the affine transform by the updated center pointing of the camera
        fp_radec_affine = cc_radec_estimate.ra.deg[0],cc_radec_estimate.dec.deg[0]
        if norm([pixels_cc_affine_x,pixels_cc_affine_y]) > min(res)/10:
            stars = simplified_catalog.search_cone(fp_radec_affine,search_radius,max_control_points)
        else:
            stars.center = fp_radec_affine

        stars.pixel_xy(pixel_width) 
        stars.invariantfeatures() 
        pixels_catalog = stars.xy
        catalog_df = stars.df
        wcs = stars.wcs

        catalog_tuple = (stars.xy,stars.asterisms,stars.kdtree)
        transf, (pixels_camera_match, pixels_catalog_match),_s,_d = find_transform_tree(camera_tuple,catalog_tuple)
        affine_matrix = transf.params
        affine_translation = transf.translation
        affine_rotation = transf.rotation # in radians
        affine_scale = transf.scale

        _, ind_catalog_match, ind_catalog = np.intersect1d(pixels_catalog_match[:,0],stars.xy[:,0],return_indices=True)
        pixels_camera_match = pixels_camera_match[ind_catalog_match]
        pixels_catalog_match = pixels_catalog_match[ind_catalog_match]
        catalog_df_affine = catalog_df.loc[ind_catalog]

        flux_camera_match = flux_camera[_s][ind_catalog_match]
        xy_affine = pixels_camera_affine = matrix_transform(pixels_camera_match,affine_matrix)
        xy_res_affine = pixels_catalog_match - pixels_camera_affine
        xy_rms_affine = np.sqrt(np.mean(xy_res_affine**2,axis=0))
        catalog_df_affine[['dx','dy']] = xy_res_affine 

        # M = C - 2.5log10(F), where C is an undetermined magnitude constant
        C_affine = (2.5*np.log10(flux_camera_match) + catalog_df_affine['mag']).mean()
        mag_res_affine = np.array(catalog_df_affine['mag']) -(C_affine - 2.5*np.log10(flux_camera_match))
        mag_rms_affine = np.sqrt(np.mean(mag_res_affine**2))
        n = len(mag_res_affine)
        C_sigma_affine = np.sqrt(np.dot(mag_res_affine,mag_res_affine)/((n-1)*n))

        catalog_df_affine['dmag'] = mag_res_affine 
        catalog_df_affine[['x_camera','y_camera']] = pixels_camera_match

        # Calculate the celestial coordinates of sources, and compare with the star catalog
        # Calculate the residual and RMS of the Ra and Dec components respectively.
        radec_res_affine,radec_rms_affine = radec_res_rms(wcs,xy_affine,catalog_df_affine)

        # Group and package the calculation results for affined sources(pixel coordinates and flux).
        description = 'Constructor of affine results'
        dict_values = xy_affine,xy_res_affine,xy_rms_affine,mag_res_affine,mag_rms_affine,C_affine,C_sigma_affine,catalog_df_affine,description,pixels_camera_match,radec_res_affine,radec_rms_affine
        dict_keys = 'xy','xy_res','xy_rms','mag_res','mag_rms','C','C_sigma','catalog_df','_description','pixels_camera_match','radec_res','radec_rms'
        info_affine_results = dict(zip(dict_keys, dict_values))
        affine_results = Constructor(info_affine_results)
        self.info['affine_results'] = affine_results
        self.affine_results = affine_results

        # Basic results for affine transformation
        dict_values = wcs,fp_radec_affine,affine_matrix,affine_translation,affine_rotation,affine_scale,L,pixel_width,fov
        dict_keys = '_wcs','fp_radec_affine','affine_matrix','_affine_translation','_affine_rotation','_affine_scale','_L','_pixel_width','_fov'
        info_update = dict(zip(dict_keys, dict_values))
        self.info.update(info_update)
        self._wcs = wcs
        self.fp_radec_affine = fp_radec_affine
        self.affine_matrix = affine_matrix
        self._affine_translation = affine_translation
        self._affine_rotation = affine_rotation
        self._affine_scale = affine_scale
        self._L = L

        # This part replaces the sources of the affine matching with the sources of the 3D-Tree matching and performs calculations similar to the previous part.
        # Apply the affine matrix and the magnitude constant to all sources in camera image, then build a dimensionless 3D-Tree for camera and starcatalog 
        pixels_camera_affine = matrix_transform(self.xy_raw,affine_matrix)
        source_xymag = np.hstack([pixels_camera_affine/L,(C_affine - 2.5*np.log10(self.flux_raw))[:,None]])
        # In order to make the stars in catalog cover sources as much as possible, the number of stars in search area is expanded to twice that of sources
        stars = simplified_catalog.search_cone(fp_radec_affine,search_radius,2*len(self.xy_raw))
        stars.pixel_xy(pixel_width) 
        catalog_df = stars.df
        catalog_xymag = np.hstack([stars.xy/L,stars.mag[:,None]])
        source_xymag_tree = KDTree(source_xymag)
        catalog_xymag_tree = KDTree(catalog_xymag)
        matches_list = source_xymag_tree.query_ball_tree(catalog_xymag_tree, r=0.2)

        # Find the matching pairs
        matches,matches_index_source,matches_index_catalog = [],[],[]
        i = 0
        for t1, t2_list in zip(source_xymag, matches_list):
            if len(t2_list) == 1: 
                matches_index_source.append(i)
                matches_index_catalog.append(t2_list[0])
                t2 = catalog_xymag_tree.data[t2_list[0]]
                matches.append(list(zip(t1, t2)))
            i += 1    

        pixels_camera_match = self.xy_raw[matches_index_source]
        catalog_df = catalog_df.loc[matches_index_catalog]
        matches = np.array(matches)
        matches_res = (matches[:,:,1] - matches[:,:,0])

        xy_L = matches[:,:,0][:,:2]
        xy_res_L = matches_res[:,:2]
        mag_res = matches_res[:,2]
        xy = xy_L * L
        xy_res = xy_res_L * L
        xy_rms = np.sqrt(np.mean(xy_res**2,axis=0))
        mag_rms = np.sqrt(np.mean(mag_res**2))

        C = (2.5*np.log10(self.flux_raw[matches_index_source]) + catalog_df['mag']).mean()
        n = len(mag_res)
        C_sigma = np.sqrt(np.dot(mag_res,mag_res)/((n-1)*n))

        catalog_df[['dx','dy']] = xy_res
        catalog_df['dmag'] = mag_res
        catalog_df[['x_camera','y_camera']] = pixels_camera_match
        radec_res,radec_rms = radec_res_rms(wcs,xy,catalog_df)

        # Group and package the calculation results for matched sources(pixel coordinates and flux).
        description = 'Constructor of matchd results'
        dict_values = xy,xy_res,xy_rms,mag_res,mag_rms,C,C_sigma,catalog_df,description,pixels_camera_match,radec_res,radec_rms
        dict_keys = 'xy','xy_res','xy_rms','mag_res','mag_rms','C','C_sigma','catalog_df','_description','pixels_camera_match','radec_res','radec_rms'
        info_match_results = dict(zip(dict_keys, dict_values))
        match_results = Constructor(info_match_results)
        self.info['match_results'] = match_results
        self.match_results = match_results

        pixel_width_estimate = pixel_width*affine_scale
        fov_estimate = pixel_width_estimate * self._res
        self.pixel_width_estimate = pixel_width_estimate
        self.fov_estimate = fov_estimate
        self.info.update(dict(zip(['pixel_width_estimate','fov_estmate'],[pixel_width_estimate,fov_estimate])))

        # Distortion correction
        if calibrate:
            # Normalize the coordinates
            xy_L = xy/L
            U_L = xy_res[:,0][:,None]/L
            V_L = xy_res[:,1][:,None]/L

            # Define kernel
            kerx = GPy.kern.RBF(2)
            # Create simple GP model
            mx_L = GPy.models.GPRegression(xy_L,U_L,kerx)
            mx_L.optimize() 
            mx_L.optimize_restarts(num_restarts = 100,verbose=False)

            kery = GPy.kern.RBF(2)
            my_L = GPy.models.GPRegression(xy_L,V_L,kery)
            my_L.optimize() 
            my_L.optimize_restarts(num_restarts = 100,verbose=False)

            self.info.update({'_mx_L':mx_L,'_my_L':my_L})
            self._mx_L,self._my_L = mx_L,my_L

            meanx_L,varx_L = self._mx_L.predict(xy_L)
            meany_L,vary_L = self._my_L.predict(xy_L)

            distortion_fitted = np.hstack([meanx_L,meany_L]) * L

            xy_calibrated = xy + distortion_fitted
            xy_res_calibrated = xy_res - distortion_fitted
            xy_rms_calibrated = np.sqrt(np.mean(xy_res_calibrated**2,axis=0))

            catalog_df_calibrated = catalog_df.copy()
            catalog_df_calibrated[['dx','dy']] = xy_res_calibrated
            radec_res_calibrated,radec_rms_calibrated = radec_res_rms(wcs,xy_calibrated,catalog_df_calibrated)

            # Group and package the calculation results for calibrated sources(pixel coordinates).
            calibrate_results = copy.deepcopy(match_results) 
            calibrate_results.xy = xy_calibrated
            calibrate_results.xy_res = xy_res_calibrated
            calibrate_results.xy_rms = xy_rms_calibrated
            calibrate_results.catalog_df = catalog_df_calibrated
            calibrate_results.radec_res = radec_res_calibrated
            calibrate_results.radec_rms = radec_rms_calibrated
            calibrate_results._description = 'Constructor of calibrated results'
            info_update = {'radec_res':radec_res_calibrated,'radec_rms':radec_rms_calibrated}
            calibrate_results.info.update(info_update)

            self.info['calibrate_results'] = calibrate_results
            self.calibrate_results = calibrate_results

        return self   

    def apply(self,xy_target,flux_target=None):
        """
        Apply the mapping model(affine model or calibration model) to unknown sources

        Usage:
            >>> radec1,M1 = sources1.apply([-400,-400],5e3) 
            >>> radec2,M2 = sources2.apply([[-400,-400],[500,500]],[1e4,5e3])

        Inputs:
            xy_target -> [array] Pixel coordinates of unknown sources
            flux_target -> [array] Flux(Grayscale value) of unknown sources
        """

        # Affine
        xy_target_affine = matrix_transform(xy_target,self.affine_matrix)

        # If a distortion calibration model exists
        if hasattr(self,'_mx_L'):
            meanx_L,varx_L = self._mx_L.predict(xy_target_affine/self._L)
            meany_L,vary_L = self._my_L.predict(xy_target_affine/self._L)
            mean_xy = np.hstack([meanx_L,meany_L]) * self._L
            calibrated_xy = xy_target_affine + mean_xy 
        else:
            calibrated_xy = xy_target_affine
                
        # Convert pixel coordinates to celestial coordinates
        radec_estimate = self._wcs.pixel_to_world(calibrated_xy[:,0],calibrated_xy[:,1])
        ra,dec = radec_estimate.ra.deg,radec_estimate.dec.deg
        radec = np.stack([ra,dec]).T

        # Photometric
        if flux_target is not None:
            M_affine = self.affine_results.C - 2.5*np.log10(flux_target)
            M_matchd = self.match_results.C - 2.5*np.log10(flux_target)
            return radec,M_affine,M_matchd
        else:
            return radec    

    def fp_calibrate(self):
        """
        Calibrate the orientation of the camera center

        Usage:
            >>> sources1.fp_calibrate()
            >>> print(sources1.fp_radec_calibrated)
        """
        fp_radec_calibrated = self.apply([0,0])
        info_update = {'fp_radec_calibrated':fp_radec_calibrated}
        self.info.update(info_update)
        self.fp_radec_calibrated = fp_radec_calibrated

        return self
            
    def show_distortion(self,mode,fig_file=None):
        """
        Show the distortion of the camera in modes of vector or contourf.

        Usage:
            >>> sources2.show_distortion('vector')
            >>> sources2.show_distortion('contourf')
        Inputs:
            mode -> [str] How distortions are displayed, including 'vector' plots or 'contourf' plots    
            fig_file -> [str] Path to save the distortion map
        outputs:
            distortion map
        """

        x,y = self.match_results.xy.T
        u,v = self.match_results.xy_res.T

        plt.clf()

        ylim,xlim = self._res/2
        # 'vector' plots
        if mode == 'vector':
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(x, y,facecolors="None", color='m',s=10,alpha=0.5) 
            plt.quiver(x, y, u, v, color='b', units='xy')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim(-xlim, xlim)
            ax.set_ylim(-ylim, ylim)
            ax.tick_params(axis='both', which='major', labelsize=7)

        # 'contourf' plots
        elif mode == 'contourf':   
            if not hasattr(self,'_mx_L'): raise Exception('No distortion correction is found')
            xx = np.linspace(-xlim, xlim, 50)
            yy = np.linspace(-ylim, ylim, 50)
            X, Y = np.meshgrid(xx, yy)

            xx_L = xx/self._L
            yy_L = yy/self._L
            X_L, Y_L = np.meshgrid(xx_L, yy_L)

            mesh_L = np.stack([X_L,Y_L]).transpose(1,2,0).reshape(-1,2)
            meanx_L,varx_L = self._mx_L.predict(mesh_L)
            meany_L,vary_L = self._my_L.predict(mesh_L)
            meanx_L = meanx_L.reshape(50,50)
            meany_L = meany_L.reshape(50,50)
            mean_L = np.array([meanx_L,meany_L])
            mean = mean_L * self._L

            fig = plt.figure(dpi=300)
            axes = fig.subplots(1,2)

            for i in range(2):

                axes[i].set_aspect('equal', adjustable='box')

                if mean[i].max()*mean[i].min() < 0:
                    abs_Z_max = np.abs(mean[i]).max()     
                    Z_levels = np.linspace(-abs_Z_max,abs_Z_max, 51)
                    cs = axes[i].contourf(X, Y, mean[i],levels = Z_levels,extend='both',cmap=plt.cm.RdBu_r) 
                else:    
                    cs = axes[i].contourf(X, Y, mean[i],extend='both') 

                cb = fig.colorbar(cs, ax=axes[i],format='%4.1f',pad=0.07,shrink=0.5)
                cb.ax.tick_params(labelsize=7)

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
            sources1.show_starmatch(self,image_raw,offset)

        Inputs:
            image_raw -> [2d array of float] Raw grayscale image with origin at bottom left corner point
            offset -> [array of float] Pixel coordinates of the center of image
            fig_out -> [str,optional,default=None] The path to save the output image.
        Outputs:
            starmatch map  
        """
        xy = self.match_results.pixels_camera_match + offset
        text = self.match_results.catalog_df.index.astype(str)
        if fig_out is None:
            plot_kwargs = {'mark':(xy,'s','red',text)}
        else:
            plot_kwargs = {'mark':(xy,'s','red',text),'figname':fig_out}
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

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class Distortion'

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

    def sketchmap(self,xlim,ylim,pixel_scale=1):   
        """
        Sketch the vector plot of distortion

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
        """ 
        plt.clf()
        fig, ax = plt.subplots(dpi=300)

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
        vx,vy = vx_/ratio,vy_/ratio

        ax.quiver(X,Y,vx.reshape(X.shape),vy.reshape(Y.shape),color="C0") # scale=20
        ax.set_aspect('equal', adjustable='box')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.tick_params(axis='both', which='major', labelsize=7)
        plt.show()
