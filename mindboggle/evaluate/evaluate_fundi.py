#!/usr/bin/python
	
def findNeighboringLabels(ind,Faces,labels):

	nlabels = [1000,1000,1000,1000]

	for i in range(len(Faces)):
		for j in range(0,3):
			if Faces[i,j] == ind:
				
				for k in range(0,3):
					found = 0
					nextFree = -1
					for s in range(0,4):
						if nlabels[s] == 1000:
							if nextFree < 0:
								nextFree = s
						
						if nlabels[s] == labels[Faces[i,k]]:
							found = 1
					if found == 0:
						nlabels[nextFree] = labels[Faces[i,k]]
	
	nlabels = sorted(nlabels)
	
	return nlabels
	
	
def extractLabelFundi(B):


	import numpy as np
	
	n_vertices = len(B)
	
	LF = np.zeros(n_vertices)
	
	lab1 = 0
	lab2 = 0
	fnum = 0

	for i in range(n_vertices):
		if B[i,0] < 1000:
		
		
			for j in range(1,46):
		
				if j == 1:		
					lab1 = 12
					lab2 = 28
					fnum = 1
										
				elif j == 2:
					lab1 = 3
					lab2 = 28
					fnum = 2

				elif j == 3:
					lab1 = 3
					lab2 = 18
					fnum = 3
			
				elif j == 4:
					lab1 = 3
					lab2 = 24
					fnum = 4

				elif j == 5:
					lab1 = 24
					lab2 = 28
					fnum = 4

				elif j == 6:
					lab1 = 18
					lab2 = 24
					fnum = 4

				elif j == 7:
					lab1 = 22
					lab2 = 24
					fnum = 5
	
				elif j == 8:
					lab1 = 22
					lab2 = 29
					fnum = 6

				elif j == 9:
					lab1 = 22
					lab2 = 31
					fnum = 6
			
				elif j == 10:
					lab1 = 29
					lab2 = 31
					fnum = 7

				elif j == 11:
					lab1 = 8
					lab2 = 29
					fnum = 7
					
				elif j == 12:
					lab1 = 8
					lab2 = 31
					fnum = 8
					
				elif j == 13:
					lab1 = 30
					lab2 = 31
					fnum = 9
					
				elif j == 14:
					lab1 = 8
					lab2 = 11
					fnum = 10
					
				elif j == 15:
					lab1 = 11
					lab2 = 29
					fnum = 10
					
				elif j == 16:
					lab1 = 11
					lab2 = 15
					fnum = 11
					
				elif j == 17:
					lab1 = 9
					lab2 = 11
					fnum = 11
					
				elif j == 18:
					lab1 = 15
					lab2 = 30
					fnum = 12
					
				elif j == 19:
					lab1 = 9
					lab2 = 15
					fnum = 13
					
				elif j == 20:
					lab1 = 30
					lab2 = 35
					fnum = 14
					
				elif j == 21:
					lab1 = 34
					lab2 = 35
					fnum = 14
					
				elif j == 22:
					lab1 = 12
					lab2 = 35
					fnum = 14
					
				elif j == 23:
					lab1 = 2
					lab2 = 35
					fnum = 14

				elif j == 24:
					lab1 = 24
					lab2 = 35
					fnum = 14
					
				elif j == 25:
					lab1 = 22
					lab2 = 35
					fnum = 14
					
				elif j == 26:
					lab1 = 31
					lab2 = 35
					fnum = 14
					
				elif j == 27:
					lab1 = 30
					lab2 = 34
					fnum = 15
					
				elif j == 28:
					lab1 = 2
					lab2 = 14
					fnum = 16
					
				elif j == 29:
					lab1 = 2
					lab2 = 28
					fnum = 16
					
				elif j == 30:
					lab1 = 2
					lab2 = 17
					fnum = 16
					
				elif j == 31:
					lab1 = 17
					lab2 = 25
					fnum = 16
					
				elif j == 32:
					lab1 = 17
					lab2 = 28
					fnum = 17
					
				elif j == 33:
					lab1 = 5
					lab2 = 25
					fnum = 18
					
				elif j == 34:
					lab1 = 13
					lab2 = 25
					fnum = 19
					
				elif j == 35:
					lab1 = 2
					lab2 = 13
					fnum = 19
					
				elif j == 36:
					lab1 = 14
					lab2 = 28
					fnum = 20
					
				elif j == 37:
					lab1 = 2
					lab2 = 4
					fnum = 21
					
				elif j == 38:
					lab1 = 12
					lab2 = 18
					fnum = 22
					
				elif j == 39:
					lab1 = 3
					lab2 = 12
					fnum = 22
					
				elif j == 40:
					lab1 = 12
					lab2 = 14
					fnum = 23
					
				elif j == 41:
					lab1 = 7
					lab2 = 9
					fnum = 24

				elif j == 42:
					lab1 = 7
					lab2 = 11
					fnum = 24
					
				elif j == 43:
					lab1 = 6
					lab2 = 7
					fnum = 25

				elif j == 44:
					lab1 = 7
					lab2 = 16
					fnum = 25
					
				elif j == 45:
					lab1 = 7
					lab2 = 13
					fnum = 25

				
				if B[i,0] == lab1 or B[i,1] == lab1 or B[i,2] == lab1:
					if B[i,1] == lab2 or B[i,2] == lab2 or B[i,3] == lab2:						
						LF[i] = fnum						
						#print('LOC3 {0}, {1}, {2}'.format(i,B[i],LF[i]))
				
							
	return LF
		
									
																
def extractLabelBoundaries(n_vertices,Faces,labels,folds):
	
	import numpy as np
	
	B = np.zeros((n_vertices,4))
	B = B + 1000
	
	counter = 0	
	for i in range(len(Faces)):
		
		for j in range(0,3):
			currVert = Faces[i,j]
	
			if folds[currVert] > 0.0:
				i1 = Faces[i,0]
				i2 = Faces[i,1]
				i3 = Faces[i,2]
				
				nlabels = [labels[i1],labels[i2],labels[i3]]
				

				for k in range(0,3):
					found = 0
					nextFree = -1
					for s in range(0,4):
					
						if B[currVert,s] == 1000:
							if nextFree < 0:
								nextFree = s
						
						if nlabels[k] == B[currVert,s]:
							found = 1
					if found == 0:
						B[currVert,nextFree] = nlabels[k]

				B[currVert] = sorted(B[currVert])
				
				if B[currVert,1] == 1000:
					B[currVert] = [1000,1000,1000,1000]
																							
	return B																										
			
def compEucDist(a,b):

	import math

	dist = pow(a[0] - b[0],2) + pow(a[1] - b[1],2) + pow(a[2] - b[2], 2)
	dist = math.sqrt(dist)		
	
	return dist

												
																		
def findDistToClosestFundusPoint(currVertex,fundi,vertices):
	minDist = 10000
	
	for i in range(len(fundi)):
		if fundi[i] > 1.0:
			currDist = compEucDist(currVertex,vertices[i])
			if currDist < minDist:
				minDist = currDist

	return minDist																																																																																																																								
																																																																																																																																																																																																																																											
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																					

def compFundusPointDists(LF,fundi,vertices):

	import numpy as np
	
	n_vertices = len(vertices)
	
	D = np.zeros((n_vertices,max(LF)+1))
	D = D - 1

	allDist = 0
	numDist = 0

	for i in range(n_vertices):
	#for i in range(0,6000):
		if LF[i] > 0:
			D[i,0] = findDistToClosestFundusPoint(vertices[i],fundi,vertices)
			D[i,LF[i]] = D[i,0]
			
			allDist += D[i,0]
			numDist += 1 			
			
		if i % 1000 == 0 and numDist > 0:
			pctDone = 100.0 * float(i) / float(n_vertices)
			meanDist = allDist / numDist
			print('Done: {0} pct, mean dist {1}, {2}'.format(pctDone,meanDist,i))

	
	meanDist = allDist / numDist
	print('Done: 100 pct, mean dist {0}'.format(meanDist))
	return D																																					
																																																										

def main():

	import re, string, sys, os
	
	mb_path = os.path.abspath('/Users/yrjo/yrjo_work/mindboggle2/fundus_evaluation/mindboggle')
	sys.path.append(mb_path)
	vtk_python_path = os.path.abspath('/Users/yrjo/yrjo_work/VTK_bin/Wrapping/Python')
	sys.path.append(vtk_python_path)
	vtk_bin_path = os.path.abspath('/Users/yrjo/yrjo_work/VTK_bin/bin')
	sys.path.append(vtk_bin_path)

	
	import vtk
	import numpy as np
	import mindboggle
	from mindboggle.utils.io_vtk import load_scalar, write_scalars
	from mindboggle.utils.io_vtk import write_scalars
	
	print(sys.argv[0])
	print('Version 20120908a')
	
	print('***')
	print('Input fundi:')
	print(sys.argv[1])

	print('***')	
	print('Input folds:')
	print(sys.argv[2])
	
	print('***')
	print('Input labels:')
	print(sys.argv[3])
	print('***')
	
	
	fundi_file = sys.argv[1]
	folds_file = sys.argv[2]
	labels_file = sys.argv[3]	
	
	#fundi_file = '/Users/yrjo/yrjo_work/mindboggle2/fundus_evaluation/results/features/_hemi_lh_subject_HLN-12-3/fundi.vtk'
	#folds_file = '/Users/yrjo/yrjo_work/mindboggle2/fundus_evaluation/results/features/_hemi_lh_subject_HLN-12-3/folds.vtk'
	#labels_file = '/Users/yrjo/yrjo_work/mindboggle2/fundus_evaluation/results/labels/_hemi_lh_subject_HLN-12-3/lh.pial.labels.max.vtk'
	
	vertices, Faces, fundi = load_scalar(fundi_file, return_arrays=1)
	vertices, Faces, folds = load_scalar(folds_file, return_arrays=1)	
	vertices, Faces, labels = load_scalar(labels_file, return_arrays=1)	
	
	
	n_vertices = len(vertices)
			
	B = extractLabelBoundaries(n_vertices,Faces,labels,folds)
															
	LF = extractLabelFundi(B)
	
	D = compFundusPointDists(LF,fundi,vertices)	
		
	#write_scalars(vtk_file, vertices, len(), Faces, D, LUT_names=[])
	TEST = write_scalars('test.vtk', vertices, len(vertices), Faces, D)
	
														
	return 0

if __name__ == '__main__':
	main()
