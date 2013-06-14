from .utils import read_vtk
from .geometric_moments_orig_m import geometric_moments_orig
from .zernike_m import zernike
from .feature_extraction_m import feature_extraction
from .reformat_zernike_m import reformat_zernike

def extract_features(volume,facets,order,num_facets,num_vertices) :
    V,F,N = volume, facets, order
    G = zernike3d.geometric_moments_orig(V,F,N,num_facets,num_vertices)
    Z = zernike3d.zernike(G,N)
    D = zernike3d.feature_extraction(Z,N)
    return D
