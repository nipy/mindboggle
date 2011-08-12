from nibabel import load
import numpy as np

def hausdorff_distance(a, b):
    dist = []
    for pt1 in a:
        dist2b = np.sqrt(np.sum((b-pt1)*(b-pt1), axis=1))
        dist.append(min(dist2b))
    return np.mean(np.array(dist))

m1 = load('m1.nii.gz')
m2 = load('m2.nii.gz')

voxels1 = np.nonzero(m1.get_data())
voxels1 = np.concatenate([voxels1[0][:,None], voxels1[1][:,None], voxels1[2][:,None], np.ones((voxels1[0].shape[0],1))], axis=1)
m1xyz = np.dot(m1.get_affine(), voxels1.T).T[:,0:3]
    
voxels2 = np.nonzero(m2.get_data())
voxels2 = np.concatenate([voxels2[0][:,None], voxels2[1][:,None], voxels2[2][:,None], np.ones((voxels2[0].shape[0],1))], axis=1)
m2xyz = np.dot(m2.get_affine(), voxels2.T).T[:,0:3]


dist12 = hausdorff_distance(m1xyz, m2xyz)
dist21 = hausdorff_distance(m2xyz, m1xyz)

print dist12, dist21, (dist12+dist21)/2
