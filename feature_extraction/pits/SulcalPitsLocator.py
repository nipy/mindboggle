#!/usr/bin/env python
# encoding: utf-8
"""
SulcalPitsLocator.py
This scripts (1) reads the Sulcal Pits given my Dr. Kiho Im's code, (2) finds the coordinates of them from the 
FreeSurfer vertices (lh.sulc.asc), (3) creates a blank nifti image with a bright dots at the locations of the Sulcal Pits.
The output nifti image will have the same size as the input image, which makes it possible to be overlaied with the 
original image for visualization purposes.

Created by Ray Razlighi on 2010-09-08.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import getopt
import os
import re

import nibabel as nb
import numpy as np


def PitsReader(VerticesFileName, PitsLocationsFileName):  #Extract the SulcalPits coordinates from cortical mesh or vertices
	
	SulcalPitsOnVertices = np.loadtxt(PitsLocationsFileName)
	Vertices = np.loadtxt(VerticesFileName)
	VerticesCounter = 0
	SulcalPitsCoordinatesList = []
	
	for x in SulcalPitsOnVertices:
		if x != 0:
			SulcalPitsCoordinatesList.append(Vertices[VerticesCounter, 1:4])
		VerticesCounter = VerticesCounter + 1
	
	SulcalPitsCoordinates = np.vstack(SulcalPitsCoordinatesList)
	return SulcalPitsCoordinates


def NiftiPitsFileGenerator(Points, NiftifileHandle, OutputNiftiFileName):
	
	NiftiFileSize = NiftifileHandle.get_shape()
	NiftiFileAffine = NiftifileHandle.get_affine()
	NiftiFileHeader = NiftifileHandle.get_header()
	Data = np.zeros(NiftiFileSize, dtype=np.int16)
		
	for x in Points:
		Data[(x[0]+round(NiftiFileSize[0]/2)), (x[1]+round(NiftiFileSize[1]/2)), (x[2]+round(NiftiFileSize[2]/2))] = 6500
		
	OutputNiftiFileHandle = nb.Nifti1Image(Data, NiftiFileAffine, header=NiftiFileHeader)
	OutputNiftiFileHandle.to_filename(OutputNiftiFileName)




def main(argv=None):
	
	if sys.argv is None:
		print 'Usage: SulcalePitsLocator <Image Name> <Sulcal Pits Directory>'
	   	sys.exit(0)
	if sys.argv[1]=='help':
		print 'Usage: SulcalePitsLocator <Image Name> <Sulcal Pits Directory>'
		sys.exit(0)
	
	ImageName = sys.argv[1]
	SulcalPitsDirectory = sys.argv[2]
	print('Extracting coresponding Image from:'+ImageName)
	print('Extracting coresponding Sulcal Pits from:'+SulcalPitsDirectory)
	
	SulcalPitsCoordinates1 = PitsReader( SulcalPitsDirectory+'lh.sulc.asc', SulcalPitsDirectory+'lh.sulc.f10ma20md5.wss.dset' )
	print('number of Sulacel Pits located:', len(SulcalPitsCoordinates1))
	
	Image1_Handle = nb.load(ImageName)
	ImageName = re.sub(r'.nii.gz$', '', ImageName)
	ImageName = re.sub(r'.nii$', '', ImageName)
	
	OutputNiftiFileName1 = ImageName+'_SulcalPits.nii.gz'
	NiftiPitsFileGenerator(SulcalPitsCoordinates1, Image1_Handle, OutputNiftiFileName1)
	print('Writing the output Nifti file in'+OutputNiftiFileName1)
	


if __name__ == '__main__':
	main()





