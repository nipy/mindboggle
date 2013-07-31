from mindboggle.shapes.zernike.multiproc import MultiprocPipeline as Pipeline
from mindboggle.utils.io_vtk import read_vtk
from argparse import ArgumentParser
from numpy import array
#from profilehooks import timecall

def main() :

    parser = ArgumentParser()
    parser.add_argument('vtk_file')
    ns = parser.parse_args()
    faces, lines, indices, points, npoints, scalars, scalar_names, input_vtk = read_vtk(ns.vtk_file)

    points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]

    zernike_moments( array(points), array(faces) )

#@timecall
def zernike_moments(V,F) :

    pl = Pipeline()
    num_vertices = 3
    num_facets = pl.size(F,0)
    N = 20 # original MATLAB code claims it will only work up to order 20
    G = pl.geometric_moments_orig(V,F,N,num_facets,num_vertices) 
    Z = pl.zernike(G,N)                                          
    Descriptors = pl.feature_extraction(Z,N)                     
    return Descriptors

if __name__ == '__main__' : main()
