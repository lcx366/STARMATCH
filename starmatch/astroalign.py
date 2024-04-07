
import numpy as np
from numpy.linalg import norm
from scipy.spatial import KDTree
from itertools import combinations
from skimage.transform import estimate_transform,matrix_transform

PIXEL_TOL = 2 # Pixel distance tolerance to assume two points are the same.
MIN_MATCHES = 5 # Minimum number of triangle matches to accept a transformation.

class _MatchTransform:
    """
    A class to find the best 2D similarity transform that maps a set of source points
    to a set of target points.
    """
    def __init__(self, source, target):
        """
        Initializes the _MatchTransform with source and target points.

        Inputs:
            source -> [array-like] Array of source points with shape of (n,2).
            target -> [array-like] Array of target points with shape of (n,2).
        """
        self.source = source
        self.target = target

    def fit(self, data):
        """
        Return the best 2D similarity transform from the points given in data.

        data: N sets of similar corresponding triangles.
            3 indices for a triangle in ref
            and the 3 indices for the corresponding triangle in target;
            arranged in a (N, 3, 2) array.
        """
        """
        Estimates the best 2D similarity transform from the source points to the target points
        given a set of corresponding triangles.

        Inputs:
            data -> [array-like] An array of shape (N, 3, 2) containing N sets of triangles, where each triangle
                is represented by 3 indices. The first dimension corresponds to the source triangles, and the second to the target triangles.
        Outputs:
            A transform object that contains the estimated similarity transformation.
        """
        # Reshape data to separate source and destination indices.
        d1, d2, d3 = data.shape
        s, d = data.reshape(d1 * d2, d3).T
        # Estimate the transformation using the "similarity" type.
        approx_t = estimate_transform("similarity", self.source[s], self.target[d])
        return approx_t

    def get_error(self, data, approx_t):
        """
        Calculates the error for each set of corresponding triangles after applying the estimated transformation.

        Inputs:
            data -> [array-like] An array of shape (N, 2, 3) as in the fit method.
            approx_t -> The estimated transformation object returned by the fit method.
        Outputs:
            error -> [array-like] An array of maximum residual errors for each set of triangles.
        """
        # Reshape data to separate source and destination indices for error calculation
        d1, d2, d3 = data.shape
        s, d = data.reshape(d1 * d2, d3).T
        # Calculate residuals after transformation and reshape to original triangle groups
        resid = approx_t.residuals(self.source[s], self.target[d]).reshape(d1, d2)
        # Determine the maximum error for each group of triangles
        error = resid.max(axis=1)
        return error

def unique_triangles(points,NUM_NEAREST_NEIGHBORS=5):
    """
    Generates unique triangles from a set of points by finding the nearest neighbors for each point,
    then calculating the ratios of the sides of these triangles. This helps in identifying invariant features
    of triangles that can be used for matching purposes.

    Inputs:
        points -> [array-like] A numpy array of shape (n, 2), representing n points in 2D space.
        NUM_NEAREST_NEIGHBORS -> [int,optional,default=5] Number of nearest neighbors to consider for triangle formation.
    Outputs:
        inv_uniq -> [array-like] A numpy array containing the ratios [L3/L2, L2/L1] for each unique triangle, where L3 is the
      longest side, L2 is the middle, and L1 is the shortest side of the triangle.
        triang_vrtx_uniq -> [array-like] A numpy array of vertex indices [A, B, C] for each triangle, sorted such that A is
      opposite L3, B is opposite L1, and C is opposite L2.
    """
    # Construct a KDTree for efficient nearest neighbor search.
    tree = KDTree(points)

    # Use a set to track unique triangle combinations based on their vertex indices.
    triangles = set()

    # For each point, query the KDTree to find its k nearest neighbors, including itself.
    for i in range(len(points)):
        # Generate all possible combinations of 3 points (triangles) from these neighbors.
        _, indices = tree.query(points[i], k=NUM_NEAREST_NEIGHBORS+1)
        for combo in combinations(indices, 3):
            # Add each triangle to the set, ensuring uniqueness.
            triangles.add(tuple(sorted(combo)))

    # Lists to hold the output ratios and vertex indices.
    inv_uniq = []
    triang_vrtx_uniq = []

    # Process each unique triangle.
    for tri in triangles:
        # Extract the coordinates of the triangle's vertices.
        tri_points = points[list(tri)]
        # Calculate the lengths of each side of the triangle and associate them with the
        # index of the vertex opposite each side. This creates a mapping of side lengths to vertices.
        sides = [
            (norm(tri_points[0] - tri_points[1]), tri[2]),  # Side opposite to vertex C
            (norm(tri_points[0] - tri_points[2]), tri[1]),  # Side opposite to vertex B
            (norm(tri_points[1] - tri_points[2]), tri[0])   # Side opposite to vertex A
        ]
        # Sort the sides by their lengths to identify L1, L2, L3, along with their corresponding opposite vertices.
        sides.sort(key=lambda x: x[0])

        # Calculate the ratios of the side lengths: L3/L2 and L2/L1.
        ratios = [sides[2][0] / sides[1][0], sides[1][0] / sides[0][0]]  # L3/L2, L2/L1
        inv_uniq.append(ratios)
        
        # Sort the vertices according to the problem specification: A opposite L3, B opposite L1, and C opposite L2.
        # This is done by retrieving the indices of the vertices in the order of the side lengths they are opposite to.
        indices_sorted_by_length = [sides[2][1], sides[0][1], sides[1][1]]
        triang_vrtx_uniq.append(indices_sorted_by_length)

    return np.array(inv_uniq), np.array(triang_vrtx_uniq)  

def invariantfeatures(sources):
    """

    Inputs:
        sources -> [2d array] The pixel coordinates of sources
    Outputs:
        invariants -> [2d array, n*2] Array of (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of sources.
        asterisms -> [2d array n*3] Array of the indices of sources that correspond to each invariant.
        kdtree -> [Instance of class scipy.spatial.KDTree] kd-tree for quick nearest-neighbor lookup
    """
    """
    Computes invariant features for a set of source points by generating unique triangles and calculating
    their side length ratios. These features are then used to build a KDTree for efficient matching.

    Usage:
        invariants,asterisms,kdtree = invariantfeatures(sources,max_control_points)
    Inputs:
        sources -> [array-like] A 2D numpy array representing the pixel coordinates of source points.
    Outputs:
        invariants -> [array-like] A 2D numpy array containing the invariant ratios (L2/L1, L1/L0) for each triangle, where L2 >= L1 >= L0 are the sides of the triangle.
        asterisms -> [array-like] A 2D numpy array containing the indices of source points that form each triangle.
        kdtree -> An instance of scipy.spatial.KDTree built from the invariants for quick nearest-neighbor lookup.
    """
    invariants,asterisms = unique_triangles(sources)
    kdtree = KDTree(invariants) # Build KDTree for the invariants.
    return invariants,asterisms,kdtree

def find_transform(source, target, max_control_points=20):
    """
    This function aims to find a Affine Transform that best maps the source points to the target points. 

    Inputs:
        source -> [array-like] A 2D array representing an iterable of (x, y) coordinates of source control points.
        target [array-like] A 2D array representing an iterable of (x, y) coordinates of target control points.
        max_control_points -> [int,optional,default=20] The maximum number of control points to use for finding the transformation.
    Outputs:
        best_t -> The transformation object with transformation parameters - rotation, translation, and scale.
        (source_controlp[s], target_controlp[d]) -> [tuple] Arrays of corresponding coordinates in the source and target.
    Raises:
        ValueError: If fewer than 3 stars are found in either the source or target.
        MaxIterError: If no transformation is found.
    """
    # Attempt to convert source and target to numpy arrays, limited by max_control_points
    try:
        source_controlp = np.array(source)[:max_control_points]
    except Exception:
        raise TypeError("Input type for source not supported.")

    try:
        target_controlp = np.array(target)[:max_control_points]
    except Exception:
        raise TypeError("Input type for target not supported.")

    # Ensure there are enough reference points in both source and target
    if len(source_controlp) < 3 or len(target_controlp) < 3:
        raise ValueError("Not enough reference stars in source or target image; minimum required is 3.")

    # Generate invariants (unique features) and asterisms (star patterns) for both source and target
    source_invariants, source_asterisms, source_invariant_tree = invariantfeatures(source_controlp)
    target_invariants, target_asterisms, target_invariant_tree = invariantfeatures(target_controlp)

    # Find matches between source and target invariants within a certain radius
    matches_list = source_invariant_tree.query_ball_tree(target_invariant_tree, r=0.1)

    # Process matches to find corresponding triangles (asterisms) in source and target
    matches = []
    for t1, t2_list in zip(source_asterisms, matches_list):
        for t2 in target_asterisms[t2_list]:
            matches.append(list(zip(t1, t2)))
    matches = np.array(matches)
    n_invariants = len(matches)

    # Attempt to fit a transformation model using the matches
    inv_model = _MatchTransform(source_controlp, target_controlp)
    
    # Decide on the fitting approach based on the number of control points and matches
    if (len(source_controlp) == 3 or len(target_controlp) == 3) and n_invariants == 1:
        # Directly fit the model if only one match is found
        best_t = inv_model.fit(matches)
        # Assume all indices are inliers since there's only one match
        inlier_ind = np.arange(n_invariants)
    else:
        # Use RANSAC to find the best model while excluding outliers
        best_t, inlier_ind = _ransac(matches, inv_model, PIXEL_TOL, MIN_MATCHES) 

    # Flatten the inlier matches to a 2D array for processing
    inlier_matches_flat = matches[inlier_ind].reshape(-1, 2)

    # Create a set of unique pairs to ensure each combination is evaluated once
    unique_pairs = set(tuple(pair) for pair in inlier_matches_flat)

    # Initialize a dictionary to track the best (lowest error) target index for each source index
    best_matches_dict = {}
    for source_idx, target_idx in unique_pairs:
        # Transform the source point using the estimated transformation parameters
        predicted_target = matrix_transform(source_controlp[source_idx], best_t.params)
        # Calculate the reprojection error for this source-target pair
        reprojection_error = norm(predicted_target - target_controlp[target_idx])

        # Update the dictionary only if this source index is not already present,
        # or if the current pair has a lower error than the previously stored one.
        if source_idx not in best_matches_dict or reprojection_error < best_matches_dict[source_idx][1]:
            best_matches_dict[source_idx] = (target_idx, reprojection_error) 

    # Extract the best source-target pairs from the dictionary
    best_source_indices, best_target_indices_errors = zip(*best_matches_dict.items())
    best_target_indices = [target_idx for target_idx, _ in best_target_indices_errors]
    best_source_indices = list(best_source_indices)

    # Return the transformation object and the arrays of best matching source and target points
    return best_t, (source_controlp[best_source_indices], target_controlp[best_target_indices])

def find_transform_tree(source, target):
    """
    Estimates an Affine Transform that best maps source points to target points. 
    This function uses control points from both the source and target to find the transformation parameters, 
    including rotation, translation, and scaling using RANSAC or direct fitting based on the number of matches.

    Inputs:
        source -> [tuple] A tuple containing a 2D array of (x, y) coordinates of source control points, source asterisms, and a source invariant KDTree.
        target -> [tuple] A tuple containing a 2D array of (x, y) coordinates of target control points, target asterisms, and a target invariant KDTree.
    Outputs:
        best_t -> The transformation object with transformation parameters - rotation, translation, and scale.
        matched_pairs -> [tuple of arrays] Arrays of corresponding coordinates in the source and target that match based on the estimated transformation.
        best_source_indices -> [list] List of indices in the source control points that are part of the best matches.
        best_target_indices -> [list] List of indices in the target control points that correspond to the best_source_indices based on the estimated transformation.
    Raises:
        ValueError: If fewer than 3 stars are found in either the source or target.
        MaxIterError: If no transformation is found.
    """
    source_controlp,source_asterisms,source_invariant_tree = source
    target_controlp,target_asterisms,target_invariant_tree = target 

    # Ensure there are enough reference points in both source and target
    if len(source_controlp) < 3 or len(target_controlp) < 3:
        raise ValueError("Not enough reference stars in source or target image; minimum required is 3.")

    # Find matches between source and target invariants within a certain radius
    matches_list = source_invariant_tree.query_ball_tree(target_invariant_tree, r=0.1)

    # Process matches to find corresponding triangles (asterisms) in source and target
    matches = []
    for t1, t2_list in zip(source_asterisms, matches_list):
        for t2 in target_asterisms[t2_list]:
            matches.append(list(zip(t1, t2)))
    matches = np.array(matches)
    n_invariants = len(matches)

    # Attempt to fit a transformation model using the matches
    inv_model = _MatchTransform(source_controlp, target_controlp)
    
    # Decide on the fitting approach based on the number of control points and matches
    if (len(source_controlp) == 3 or len(target_controlp) == 3) and n_invariants == 1:
        # Directly fit the model if only one match is found
        best_t = inv_model.fit(matches)
        # Assume all indices are inliers since there's only one match
        inlier_ind = np.arange(n_invariants)
    else:
        # Use RANSAC to find the best model while excluding outliers
        best_t, inlier_ind = _ransac(matches, inv_model, PIXEL_TOL, MIN_MATCHES) 

    # Flatten the inlier matches to a 2D array for processing
    inlier_matches_flat = matches[inlier_ind].reshape(-1, 2)

    # Create a set of unique pairs to ensure each combination is evaluated once
    unique_pairs = set(tuple(pair) for pair in inlier_matches_flat)

    # Initialize a dictionary to track the best (lowest error) target index for each source index
    best_matches_dict = {}
    for source_idx, target_idx in unique_pairs:
        # Transform the source point using the estimated transformation parameters
        predicted_target = matrix_transform(source_controlp[source_idx], best_t.params)
        # Calculate the reprojection error for this source-target pair
        reprojection_error = norm(predicted_target - target_controlp[target_idx])

        # Update the dictionary only if this source index is not already present,
        # or if the current pair has a lower error than the previously stored one.
        if source_idx not in best_matches_dict or reprojection_error < best_matches_dict[source_idx][1]:
            best_matches_dict[source_idx] = (target_idx, reprojection_error) 

    # Extract the best source-target pairs from the dictionary
    best_source_indices, best_target_indices_errors = zip(*best_matches_dict.items())
    best_target_indices = [target_idx for target_idx, _ in best_target_indices_errors]
    best_source_indices = list(best_source_indices)

    return best_t, (source_controlp[best_source_indices], target_controlp[best_target_indices]),best_source_indices,best_target_indices

def _ransac(data, model, thresh, min_matches):
    """
    Fit model parameters to the given dataset using the RANSAC algorithm. 

    Inputs:
        data -> [array-like] An array of matched points, where each match consists of corresponding indices from two datasets.
        model -> [_MatchTransform object] A model that supports a fit method to estimate its parameters from a dataset and a get_error 
        method to compute the error of fitting the data to the model.
        thresh -> [float] A threshold value for determining when a data point fits the model. Data points with an error below 
        this threshold are considered inliers.
        min_matches -> [int] The minimum number of data points required to fit the model. If a set of data points results in a 
        model fit that has fewer inliers than this threshold, the model is discarded.
    Outputs:
        bestfit -> [SimilarityTransform object] The model parameters that best fit the data, or None if no suitable model is found.
        best_inlier_idxs -> [array-like] Indices of the data points that are considered inliers to the best fitting model.
    Raises:
        Raises an exception if no acceptable model is found after iterating through the dataset.
    """
    good_fit = None
    n_data = data.shape[0]
    all_idxs = np.arange(n_data)
    np.random.shuffle(all_idxs)

    for iter_i in range(n_data):
        # Randomly sample one point (minimal sample set to fit the model) in each iteration
        maybe_idxs = all_idxs[iter_i:iter_i+1]
        maybeinliers = data[maybe_idxs]
        maybemodel = model.fit(maybeinliers)
        
        # Compute errors for all data points using the tentative model
        err = model.get_error(data, maybemodel)

        # Select inliers based on the error threshold
        inliers = data[err < thresh]

        # Update the best model if the current iteration yields more inliers
        if len(inliers) >= min_matches:
            # Fit the model to all inliers
            good_fit = model.fit(inliers)
            break       

    if good_fit is None:
        raise Exception("List of matching triangles exhausted before an acceptable transformation was found")

    # Fit the model to the final set of inliers
    err = model.get_error(data, good_fit)

    better_data = data[err < thresh]
    better_fit = model.fit(better_data)
    better_inlier_idxs = np.arange(n_data)[err < thresh] # Indices of inliers for the best model

    return better_fit, better_inlier_idxs
