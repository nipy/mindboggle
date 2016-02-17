from . import zernike
#from .test.multiproc import MultiprocPipeline
from mindboggle.mio.vtks import read_vtk
import numpy as np

import argparse
import logging
#import profilehooks

def example1():
    #    >>> # Example 1: simple cube (decimation results in a Segmentation Fault):
    #    >>> from mindboggle.shapes.zernike import zernike_moments
    #    >>> points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    #    >>> faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    #    >>> order = 3
    #    >>> scale_input = True
    #    >>> zernike_moments(points, faces, order, scale_input)
    #    [0.0918881492369654,
    #     0.09357431096617608,
    #     0.04309029164656885,
    #     0.06466432586854755,
    #     0.03820155248327533,
    #     0.04138011726544602]

    print "Example 1: Simple cube"
    points = np.array([[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]])
    faces = np.array([[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]])
    result = zernike(points, faces, order=3, scale_input=True)
    assert np.allclose(result, np.array([0.0918881492369654, 0.09357431096617608, 0.04309029164656885,
                                         0.06466432586854755, 0.03820155248327533, 0.04138011726544602]))
    
def main():
    zernike_fn = zernike

    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', nargs='?', default=None, const='debug', choices=['debug', 'info', 'warning', 'error', 'critical']) 
    parser.add_argument('vtk_file', nargs='?', default=None)
    parser.add_argument('-o', '--order', type=int, default=3)
    parser.add_argument('-p', '--profile', nargs='?', default=None, const='stdout')
    parser.add_argument('-t', '--timecall', default=False, action='store_true')
    parser.add_argument('-v', '--validate', default=False, action='store_true')
    ns = parser.parse_args()

    if ns.debug is not None:
        logging.basicConfig(level=getattr(logging, ns.debug.upper()))

#    if ns.profile is not None:
#        filename = ns.profile
#        if ns.profile == 'stdout':
#            filename = None
#        zernike_fn = profilehooks.profile(zernike_fn, immediate=False, filename=filename)

#    if ns.timecall:
#        zernike_fn = profilehooks.timecall(zernike_fn)

    if ns.vtk_file is not None:
        points, indices, lines, faces, depths, scalar_names, npoints, \
            input_vtk = read_vtk(ns.vtk_file)
        print(len(faces), len(points))
        X = zernike_fn(points, faces, order=ns.order, scale_input=True)
        if ns.validate:
            Y = zernike_fn(points, faces, order=ns.order, scale_input=True, pl_cls=MultiprocPipeline)
            assert np.allclose(X, Y)
    else:
        example1()


if __name__ == '__main__':
    main()
