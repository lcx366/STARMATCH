
# The basic implementation of this script comes from Astroalign: A Python module for astrological image registration.
# We acknowledge the author's work and contributions.

import numpy as np
from numpy.linalg import norm
from skimage.transform import estimate_transform,matrix_transform
class _MatchTransform:
    """
    A class to find the best 2D similarity transform that maps a set of source points
    to a set of target points.
    """
    def __init__(self, source, target, ttpte):
        """
        Initializes the _MatchTransform with source and target points.

        Inputs:
            source -> [array-like] Array of source points with shape of (n,2).
            target -> [array-like] Array of target points with shape of (n,2).
            ttpte -> [str] Type of transform. Available options are 'similarity' and 'affine'.

        The similarity transformation has the following form in 2D:
        X = a0 * x - b0 * y + a1 = s * x * cos(rotation) - s * y * sin(rotation) + a1
        Y = b0 * x + a0 * y + b1 = s * x * sin(rotation) + s * y * cos(rotation) + b1

        M = [[a0  -b0  a1]
             [b0  a0  b1]
             [0   0    1]]

        The 6-parameters affine transformation has the following form in 2D:
        X = a0 * x + a1 * y + a2 = s * x * [cos(rotation) + tan(shear_y) * sin(rotation)] - s * y * [tan(shear_x) * cos(rotation) + sin(rotation)] + translation_x
        Y = b0 * x + b1 * y + b2 = s * x * [sin(rotation) - tan(shear_y) * cos(rotation)] - s * y * [tan(shear_x) * sin(rotation) - cos(rotation)] + translation_y

        M = [[a0  a1  a2]
             [b0  b1  b2]
             [0   0    1]]

        s = sqrt(det(M[:2,:2]))
        """
        self.source = source
        self.target = target
        self.ttpte = ttpte

    def fit(self, data):
        """
        Estimates the best 2D transform from the source points to the target points given a set of corresponding triangles/quads.

        Inputs:
            data -> [array-like] An array of shape (N, M, 2) containing N sets of triangles/quads,
            where each triangle is represented by 3 indices, and each quad is represented by 4 indices.
            The first dimension corresponds to the source triangles/quads, and the second to the target triangles/quads.
        Outputs:
            A transform object that contains the estimated transformation.
        """
        flattened_arr = data.reshape(-1, 2)
        s,d = np.unique(flattened_arr, axis=0).T
        # If duplicate points are not removed, the number of duplicate points is equivalent to the weight of the least squares.
        # s,d = flattened_arr.T # Uncomment to keep the duplicate points.
        # Estimate the transformation
        approx_t = estimate_transform(self.ttpte, self.source[s], self.target[d])
        return approx_t

    def get_error(self, data, approx_t):
        """
        Calculates the error for each set of corresponding triangles after applying the estimated transformation.

        Inputs:
            data -> [array-like] An array of shape (N, M, 2) as in the fit method.
            approx_t -> The estimated transformation object returned by the fit method.
        Outputs:
            error -> [array-like] An array of maximum residual errors for each set of triangles/quads.
        """
        # Reshape data to separate source and destination indices for error calculation
        d1, d2, d3 = data.shape
        s, d = data.reshape(d1 * d2, d3).T
        # Calculate residuals after transformation and reshape to original triangle/quad groups
        resid = approx_t.residuals(self.source[s], self.target[d]).reshape(d1, d2)
        # Determine the maximum error for each group of triangles/quads
        error = resid.max(axis=1)
        return error

def find_transform_tree(source, target, pixel_tol, min_matches, ttpte='similarity'):
    """
    Estimates a similarity transform that best maps source points to target points,
    including rotation, translation, and scaling, using RANSAC or direct fitting based on the number of matches.

    Inputs:
        source -> [tuple] A tuple containing
            - a 2D array of (x, y) coordinates of source points
            - source asterisms
            - a source invariant features KDTree.
        target -> [tuple] A tuple containing
            - a 2D array of (x, y) coordinates of target points
            - target asterisms
            - target invariant features KDTree.
        min_matches -> [int] Minimum number of triangle or quad matches to accept a transformation.    
        ttpte -> [str,optional,default='similarity'] Type of transform. Available options are 'similarity' and 'affine'.
    Outputs:
        best_t -> The transformation object with transformation parameters - rotation, translation, and scale.
        matched_pairs -> [tuple of arrays] Arrays of corresponding coordinates in the source and target that match based on the estimated transformation.
        best_source_indices -> [list] List of indices in the source points that are part of the best matches.
        best_target_indices -> [list] List of indices in the target points that correspond to the best_source_indices.
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
    matches_list = source_invariant_tree.query_ball_tree(target_invariant_tree, r=0.01)

    # Process matches to find corresponding asterisms in source and target
    matches = []
    for t1, t2_list in zip(source_asterisms, matches_list):
        if t2_list:  # Ensure t2_list is not empty
            for t2 in target_asterisms[t2_list]:
                matches.append(list(zip(t1, t2)))

    if not matches:
        raise MaxIterError("No transformation found; insufficient matches.")

    matches = np.array(matches)
    n_invariants = len(matches)

    # Attempt to fit a transformation model using the matches
    inv_model = _MatchTransform(source_controlp, target_controlp, ttpte)
    
    # Decide on the fitting approach based on the number of control points and matches
    if (len(source_controlp) == 3 or len(target_controlp) == 3) and n_invariants == 1:
        # Directly fit the model if only one match is found
        best_t = inv_model.fit(matches)
        # Assume all indices are inliers since there's only one match
        inlier_ind = np.arange(n_invariants)
    else:
        # Use RANSAC to find the best model while excluding outliers
        best_t, inlier_ind = _ransac(matches, inv_model, pixel_tol, min_matches)

    # Flatten the inlier matches to a 2D array for processing
    inlier_matches_flat = matches[inlier_ind].reshape(-1, 2)

    # Create a set of unique pairs to ensure each combination is evaluated once
    unique_pairs = set(map(tuple, inlier_matches_flat))

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
        raise Exception("List of matching patterns exhausted before an acceptable transformation was found.")
    # Fit the model to the final set of inliers
    err = model.get_error(data, good_fit)

    best_data = data[err < thresh]
    best_fit = model.fit(best_data)
    best_inlier_idxs = np.arange(n_data)[err < thresh] # Indices of inliers for the best model

    return best_fit, best_inlier_idxs