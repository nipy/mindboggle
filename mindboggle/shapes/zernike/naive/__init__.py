import numpy
import numpy.linalg
import scipy
import scipy.misc

from .geometric_moments_orig_m import geometric_moments_orig
from .zernike_m import zernike
from .feature_extraction_m import feature_extraction
from .reformat_zernike_m import reformat_zernike
from .trinomial_matrix_m import trinomial_matrix
from .D_CV_orig_m import D_CV_orig
from .Dabc_orig_m import Dabc_orig
from .factorial_precalc_m import factorial_precalc
from .D_SG_orig_m import D_SG_orig
from .D_SG_orig_part_m import D_SG_orig_part
from .mon_comb_m import mon_comb
from .Qklnu_m import Qklnu
from .Yljm_m import Yljm
from .trinomial_m import trinomial

from scipy.io import loadmat as load
from numpy import zeros, sum

from math import sqrt
from scipy.misc import comb as nchoosek
from numpy import power, sqrt
from .trinomial_m import trinomial
from numpy import zeros, power

class NaivePipeline(object) :
    geometric_moments = geometric_moments_orig
    zernike_moments = zernike
    feature_extraction = feature_extraction
    
    zeros = numpy.zeros
    mod = numpy.mod
    conj = numpy.conj
    power = numpy.power
    flipud = numpy.flipud
    nan = numpy.NaN
    argwhere = numpy.argwhere
    transpose = numpy.transpose
    where = numpy.where
    cat = numpy.concatenate
    norm = numpy.linalg.norm
    array = numpy.array
    sum = numpy.sum
    det = numpy.linalg.det
    hstack = numpy.hstack
    imag = numpy.imag
    real = numpy.real
    conj = numpy.conj
    pi = numpy.pi

    sqrt = scipy.sqrt
    factorial = scipy.misc.factorial
    nchoosek = scipy.misc.comb

def extract_features(volume,facets,order,num_facets,num_vertices) :
    V,F,N = volume, facets, order
    pl = NaivePipeline()
    G = pl.geometric_moments(V,F,N,num_facets,num_vertices)
    Z = pl.zernike(G,N)
    D = pl.feature_extraction(Z,N)
    return D
