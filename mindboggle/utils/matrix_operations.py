#!/usr/bin/python
"""
Operations on volume matrices.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#------------------------------------------------------------------------------
# Find all neighbors from faces
#------------------------------------------------------------------------------
def flip_axes(input_file, flipx=True, flipy=True, flipz=False):
    """
    Flip image volume data about x, y, and/or z axes
    without regard to the file header. The default is to flip
    x and y to rotate the image volume by 180 degrees.

    Parameters
    ----------
    input_file : image volume
    flipx : Boolean
        flip about x-axis?
    flipy : Boolean
        flip about y-axis?
    flipz : Boolean
        flip about z-axis?

    Returns
    -------
    out_file : string
        output image volume

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load image volume
    data = nb.load(input_file).get_data()
    lenx, leny, lenz = np.shape(data)
    data_new = np.zeros((lenx, leny, lenz))

    # Flip x
    if flipx:
        for x in range(lenx):
            data_new[lenx-1-x,:,:] = data[x,:,:]

    # Flip y
    if flipy:
        for y in range(leny):
            data_new[:,leny-1-y,:] = data[:,y,:]

    # Flip z
    if flipz:
        for z in range(lenz):
            data_new[:,:,lenz-1-z] = data[z,:,:]

    # Save output
    out_file = 'reorient_' + os.path.basename(input_file)
    img = nb.Nifti1Image(data_new, np.eye(4))
    img.to_filename(out_file)

    return out_file

#------------------------------------------------------------------------------
# Example
#------------------------------------------------------------------------------
if __name__ == "__main__":

    from mindboggle.utils.matrix_operations import flip_axes

    input_file = 'PhStr_MG_20120801.nii.gz'
    flip_axes(input_file)
