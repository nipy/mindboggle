import mindboggle.guts
import numpy
import os.path
import scipy.io
import glob
import functools

from . import MultiprocPipeline as Pipeline
 
def read_vtk(filename) :
    faces, lines, indices, points, npoints, depths, name, input_vtk = mindboggle.guts.read_vtk(input_vtk)
    return numpy.array(points), numpy.array(faces)

head,tail = os.path.split(__file__)
DATA_DIR = os.path.join(head,'test_data')
ALLOWED_ERROR = 1e-8

FNS_AND_ARGS = [ #('D_CV_orig',              ('facet','num_vertices','N','X','K','tri_matrix',), ('C','Vol',) ),
                 #('trinomial',              ('i','j','k',),                                     ('t',) ),
                 #('trinomial_matrix',       ('N',),                                             ('tri_matrix',) ),
                 #('mon_comb',               ('tri_matrix','vertex','N',),                       ('c',) ),
                 #('D_SG_orig_part',         ('num_facets','i','j','k','C','D','S','Vol','F',),  ('S',) ),
                 #('D_SG_orig',              ('num_facets','i','N','C','D','Vol','F',),          ('G',) ),
                 #('Dabc_orig',              ('C','N',),                                         ('D',) ), #slow
                 #('factorial_precalc',      ('N',),                                             ('F',) ),
                 #('Yljm',                   ('l','j','m',),                                     ('y',) ),
                 #('Qklnu',                  ('k','l','nu',),                                    ('q',) ),
                 #('zernike',                ('G','N',),                                         ('Z',) ),
                 #('feature_extraction',     ('Z','N',),                                         ('Descriptors',) ),
                 #('geometric_moments_orig', ('X','K','N','num_facets','num_vertices',),         ('G',) ), #slow
                 ('demo',                   ('V','F','ZMvtk'),                                   ('Descriptors',) ),
                 #('reformat_zernike',      ('Z','N',),                                         ('ZM,') ),
                 ]

def test() :
    pl = Pipeline()
    for fnname, inargs, outargs in FNS_AND_ARGS :
        globname = '{}-*.mat'.format(fnname)
        fn = getattr(pl,fnname)
        files = sorted( glob.glob( os.path.join(DATA_DIR,globname) ) )
        if len(files) == 0 : raise Exception(globname)
        for filename in files :
            p = functools.partial(run, fn, filename, inargs, outargs)
            p.description = filename
            yield (p, )

def run(fn,filename,inarg_names,outarg_names) :
    matfile = scipy.io.loadmat( filename,squeeze_me=True )
    inargs = tuple([ matfile[a] for a in inarg_names ])
    outargs = tuple([ matfile[a] for a in outarg_names ])
    new_outargs = fn(*inargs)
    if not isinstance(new_outargs,tuple) : new_outargs = tuple([new_outargs])
    for n,a,aa in zip(outarg_names, outargs, new_outargs) :
        if isinstance(a,numpy.ndarray) and isinstance(aa,numpy.ndarray) : assert a.shape == aa.shape, 'Output {}: {} != {}'.format(n,a.shape,aa.shape)
        err_max = numpy.max( numpy.abs( a-aa ) )
        assert err_max < ALLOWED_ERROR, 'Output {}: Error ({}) > ALLOWED_ERROR ({})'.format(n,err_max, ALLOWED_ERROR)
        #assert a.dtype == aa.dtype, '{} != {}'.format(a.dtype,aa.dtype)
