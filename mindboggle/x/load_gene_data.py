import os
import numpy as np
from scipy.io import loadmat
from mindboggle.utils.io_vtk import read_vtk
from mindboggle.utils.compute import point_distance

G = loadmat('/Users/arno/Dropbox/MB/data/allen/H0351.2002.mat')
values = G['expr']
gene_mni = G['mni']
path = os.environ['MINDBOGGLE_DATA']
genes = []
gene_values = []
for i in range(25):
    print(i)
    input_vtk = os.path.join(path, 'allen', 'labels_traveldepth' + str(i) + '.vtk')
    if os.path.exists(input_vtk):
        faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(input_vtk)

        I = [i for i,x in enumerate(depths) if x>-1]
        gene = 0
        value = 0
        print(len(I))
        points2 = np.array(points)
        points2 = points2[I]
        for point in points2:
            mind, minI = point_distance(point, gene_mni)
            if np.max(values[minI]) > value:
                gene = minI
                value = np.max(values[minI])

        genes.append(gene)
        gene_values.append(value)

    else:
        genes.append(0)
        gene_values.append(0)


