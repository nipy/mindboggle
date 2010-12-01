#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Ray Razlighi on 2010-09-08.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os

import nibabel as nb
import numpy as np

img = nb.load(os.path.join('/Users/ray/Research/Images/ColumbiaScans/Head2_2_Registered_To_2_1.nii.gz'))

print(img.get_shape())
print(img.get_data_dtype())
print(img.get_affine().shape)


data = img.get_data()
print(data.shape)
print(type(data))


hdr = img.get_header()
print(hdr)

def main():
	pass


if __name__ == '__main__':
	main()



