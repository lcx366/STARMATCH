# Maximum number of stars to extract from each HEALPix tile.
MAX_NUM_PER_TILE = 5

# Number of nearest stars to select from the catalog (not detected sources) for neighbor searches.
NUM_NEAREST_NEIGHBORS = 10

# Pixel distance tolerances for considering two points identical during blind, primary, and secondary affine transformations.
PIXEL_TOLS = (60, 20, 3)

# Minimum number of inliers for RANSAC in initial astrometric matching: (triangles, quads)
MIN_MATCHES = (8, 2)

