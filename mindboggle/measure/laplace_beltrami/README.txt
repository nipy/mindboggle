For Laplace Beltrami operator:
	Module of interest - laplace_beltrami.py

For Computing dissimilarity between shapes:
	Module of interest - shape_analysis.py
	
Under the hood:
	vtk_operations.py - extracts nodes and meshes from vtk file
	compute_weights.py - constructs affinity matrix based on graph extracted above
	kernels.py - contains functions necessary for computing the weight between vertices
	graph_operations.py - computes graph Laplacian based on affinity matrix computed above
	
Modules are pretty well documented and should be self-explanatory.

Of note, though:
(1)	LBO computation will take a long time for graphs comprising >6500 nodes.
	This problem may be addressed by fixing the call to scipy.sparse.linalg.eigen.arpack.eigsh()
	At present, the algorithm is asked to find the eignevalues of smallest magnitude. 
	This is known to be computationally intensive.
	The various solutions which exist include changing the tolerance threshold (which will produce inaccurate results)
	or using the shift-invert mode, which in theory should work, though I've been having difficulties getting consistent/correct results. 
	Perhaps there are other pieces of code which we can use for this problem. 

(2) There is currently (in the code base) no method of computing the volume of a 3D surface (as found, for example, in a vtk file).
	There is probably a simple way to calculate it.
	In any event, it is not necessary for running the code, as it only appears in computing the normalized weighted spectral distance metric.
	Therefore the code can still compute Reuter's Shape-DNA and Konukoglu's unnormalized weighted spectral distance.

(3) When doing some preliminary tests, I've been using cotangent_kernel for 3D surfaces, and rbf_kernel for 1D curves.
