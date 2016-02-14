import numpy
import numpy.linalg
import scipy
import scipy.misc
import time
#import profilehooks
import copy_reg
import types
import itertools

from .Dabc_orig_m import iter_rev_sum as Dabc_orig
from .D_CV_orig_m import D_CV_orig
from .D_SG_orig_m import D_SG_orig
from .D_SG_orig_part_m import D_SG_orig_part
from .factorial_precalc_m import factorial_precalc
from .feature_extraction_m import feature_extraction
from .geometric_moments_orig_m import mproc as geometric_moments_orig
from .mon_comb_m import mon_comb
from .Qklnu_m import Qklnu
from .reformat_zernike_m import reformat_zernike
from .trinomial_matrix_m import trinomial_matrix
from .trinomial_m import trinomial
from .Yljm_m import Yljm
from .zernike_m import zernike
from .zk_demo_m import zk_demo

class MultiprocPipeline(object) :
    geometric_moments_orig = geometric_moments_orig
    trinomial_matrix = trinomial_matrix
    trinomial = trinomial
    D_CV_orig = D_CV_orig
    mon_comb = mon_comb
    #Dabc_orig = Dabc_orig
    D_SG_orig = D_SG_orig
    D_SG_orig_part = D_SG_orig_part
    factorial_precalc = factorial_precalc
    Yljm = Yljm
    Qklnu = Qklnu
    zernike = zernike
    feature_extraction = feature_extraction

    pi = numpy.pi
    NaN = numpy.NaN

    def __init__(self) :
        self.tic_time = None

    def Dabc_orig(*args,**dargs) : return Dabc_orig(*args,**dargs) # you have to do this for "import ... as ..." functions that are run asynch
    def display(self,astring,*args) :
        if len(args) > 0 : astring = astring.format(*args)
        print(astring)
    def zeros(self,*args,**dargs) : return numpy.zeros(args,**dargs)
    def rng(self,*args) :
        if len(args) == 1 : return xrange(args[0]+1)
        elif len(args) == 2 : return xrange(args[0],args[1]+1)
        elif len(args) == 3 : return xrange(args[0],args[1]+1,args[2])
        else : raise Exception()
    def rng_prod(self,*args,**dargs) :
        args = tuple([ self.rng(*a) for a in args ])
        return itertools.product(*args,**dargs)
    def factorial(self,*args,**dargs) : return scipy.misc.factorial(*args,**dargs)
    def tic(self) :
        self.tic_time = time.time()
    def toc(self) :
        if self.tic_time is None : raise Exception()
        t = self.tic_time
        self.tic_time = None
        return time.time() - t
    def mod(self,*args,**dargs) : return numpy.mod(*args,**dargs)
    def power(self,*args,**dargs) : return numpy.power(*args,**dargs)
    def det(self,*args,**dargs) : return numpy.linalg.det(*args,**dargs)
    def cat(self,axis,*args) :
        assert len(args) >= 2
        ndim = args[0].ndim
        assert all([ a.ndim == ndim for a in args])
        if axis >= ndim : args = tuple([ numpy.expand_dims(a,axis) for a in args ])
        return numpy.concatenate(args,axis=axis)
    def sqrt(self,*args,**dargs) : return scipy.sqrt(*args,**dargs)
    def nchoosek(self,*args,**dargs) : return scipy.misc.comb(*args,**dargs)
    def floor(self,*args,**dargs) : return numpy.floor(*args,**dargs).astype(int)
    def conj(self,*args,**dargs) : return numpy.conj(*args,**dargs)
    def sum(self,*args,**dargs) : return numpy.sum(*args,**dargs)
    def real(self,*args,**dargs) : return numpy.real(*args,**dargs)
    def imag(self,*args,**dargs) : return numpy.imag(*args,**dargs)
    def norm(self,arr,ord=None) :
        if ord is not None : return numpy.linalg.norm(arr,ord=ord)
        return numpy.linalg.norm(arr)
    def rot90(self,arr,k) :
        if len(arr.shape) == 1 and k == 2 : return numpy.flipud(arr)
        else : raise NotImplementedError()
    def find(self,bool_arr) : return bool_arr
    def size(self,arr,dim=None) :
        if dim is None : return arr.shape
        else : return arr.shape[dim]
    #@profilehooks.profile(filename='demo.prfl')
    def demo(self,V,F,ZMvtk) : return zk_demo(self,V,F,ZMvtk)


def unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try : func = cls.__dict__[func_name]
        except KeyError : pass
        else : break
    return func.__get__(obj, cls)


def pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return unpickle_method, (func_name, obj, cls)


copy_reg.pickle(types.MethodType, pickle_method, unpickle_method)
