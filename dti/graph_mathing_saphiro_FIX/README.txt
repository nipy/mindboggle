WEIGHTED GRAPH MATCHING METHOD.

The code implement the weighted graph matching algorithm proposed in:

Larry S Shapiro and J Michael Brady "Feature-based correspondence: an eigenvector approach" Image and vision computing, 10 (1992) 283-288.

However the algorithm does not perform well because:
- the method is very sensitive to noise.
- one of the key step of the method is the reorientation of the eigenvectors to obtain the same orientation between the two graphs. Unfortunately the reorientation method is described in a PhD thesis I could not find.