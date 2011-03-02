
#include "fmarch_v2.h"
/* Allow non-constant speed term, but slower */

double compute_dist(Fvector3d *v0, Fvector3d *v1)
{
  double x, y, z;

  x = v1->x - v0->x;
  y = v1->y - v0->y;
  z = v1->z - v0->z;

  return(sqrt(x*x+y*y+z*z));
}

void GetTwoVertex(Surface *surface, int faceID, int vertex0ID, int *vertex1ID, int *vertex2ID)
{
	if(	surface->faces[faceID].v[0] == vertex0ID)
	{
		*vertex1ID = surface->faces[faceID].v[1];
		*vertex2ID = surface->faces[faceID].v[2];
	}
	else if(	surface->faces[faceID].v[1] == vertex0ID)
	{
		*vertex1ID = surface->faces[faceID].v[0];
		*vertex2ID = surface->faces[faceID].v[2];
	}
	else if(	surface->faces[faceID].v[2] == vertex0ID)
	{
		*vertex1ID = surface->faces[faceID].v[1];
		*vertex2ID = surface->faces[faceID].v[0];
	}
	else
	{
		printf("Error: can not find the vertex in the face!\n");
	}

	 return;
}

double *FastMarchMesh(Surface *mesh, int *contour, int numinitvert, double *F,int terminator)
{
  double dummtT;
  double *T;  /* final distances of each vertex to the contour */
  double tempvalue, newvalue;
  unsigned char *label;      /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */  
  int *BackPointer;  /* backpointer to the narrowband heap */
  Xheap H;           /* Narrrow band heap */
  XheapElement he;
  int heapindex;
  double currentF;

  int i, j, k, NN, NN1, NN2, VN, VId, neighVId, n0, n1;
  int m;

  //VN = XUGetVertexTableVN(mesh);
  VN = mesh->vertexNum;
  T = (double *) malloc(sizeof(double)*VN);
  label = (unsigned char *) malloc(sizeof(unsigned char)*VN);
  BackPointer = (int *) malloc(sizeof(int)*VN);

  /* Initialize heap structure */
  H = xhInitEmpty();

  /* Initialization for marching inwards and outwards simultaneously */
  for (i=0; i<VN; i++) {
    label[i] = (unsigned char) FAWAY;  /*All points are labelled as FAR AWAY */
    T[i] = INFINITY;
  }

  for (i=0; i<numinitvert; i++) {
    label[contour[i]] = (unsigned char) ALIVE; /* all vertices on contour are alive */
    T[contour[i]] = 0; /* all vertices on contour have distance zero */
  }

  /* Initialize NarrowBand Heap */

  for (i=0; i<VN; i++) {

    if (label[i] == (unsigned char)ALIVE) {

      /* put its neighbors into the NarrowBand */ 
      //NN1 = XUGetNeighborEN(mesh, i);
		NN1 = mesh->neighbors[i].size();
      for (j=0; j<NN1; j++) {
        //neighVId = XUGetVertexNeighborVId(mesh, i, j);
		  neighVId = mesh->neighbors[i][j];

        if (label[neighVId] == (unsigned char)FAWAY) {
          /* compute the distance to this point and add it to the heap */ 
          label[neighVId] = (unsigned char)NBAND;
          /* Check all possible triangles surrounding vertex neighVId1 */
          /* Note: Only ALIVE points contribute to the distance computation */
          //NN2 = XUGetNeighborEN(mesh, neighVId);
/*
		  NN2 = mesh->neighbors[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
			n0 = mesh->neighbors[neighVId][k];
            n1 = mesh->neighbors[neighVId][(k+1)%NN2];

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);
            if (tempvalue < newvalue)  newvalue = tempvalue;
          }
*/ 
		  NN2 = mesh->adjacentfaces[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
		   
//			n0 = mesh->neighbors[neighVId][k];
//          n1 = mesh->neighbors[neighVId][(k+1)%NN2];

			GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][k], neighVId, &n0, &n1);

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);

            if (tempvalue < newvalue)  newvalue = tempvalue;
          }

          T[neighVId] = newvalue;
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);

        } /* end if */
      } /* end for */    
    } /* end if */
  } /* end for */

  /* End of Initialization */
 
  /* Begin Fast Marching to get the unsigned distance function inwords
   * and outwards simultaneously
   * since points inside and outside the contour won't interfere with
   * each other
   */
  //printf("   "); 
  m = 0;
  while(!xhIsEmpty(H)){ /* There are still points not yet accepted */
    he = xhRemove(H); /* Label the point with smallest value among all
                         NarrowBand points as ALIVE */
    m++;
    VId = he.id; 

    T[VId] = he.value;

    label[VId] = (unsigned char)ALIVE;
    if ( T[terminator]*2 < he.value )
	  break; 
    
    /* Update its neighbor */
    /* Put FARAWAY neighbors into NarrowBand, Recompute values at
     * NarrowBand neighbors,
     * Keep ALIVE (Accepted) neighbor unchanged
     */
//    NN1 = XUGetNeighborEN(mesh, VId);
	NN1 = mesh->neighbors[VId].size();
    for (i=0; i<NN1; i++) {
      //neighVId = XUGetVertexNeighborVId(mesh, VId, i);
		neighVId = mesh->neighbors[VId][i];
      
	  /* Don't change ALIVE neighbors */
      if (label[neighVId] != (unsigned char)ALIVE) { 
/*
//       NN2 = XUGetNeighborEN(mesh, neighVId);
		NN2 = mesh->neighbors[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {
//          n0 = XUGetVertexNeighborVId(mesh, neighVId, j);
//          n1 = XUGetVertexNeighborVId(mesh, neighVId, (j+1)%NN2);
		  n0 = mesh->neighbors[neighVId][j];
		  n1 = mesh->neighbors[neighVId][(j+1)%NN2];
*/
		NN2 = mesh->adjacentfaces[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {

		GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][j], neighVId, &n0, &n1);

		  tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				BackPointer, &H);
          if (tempvalue < newvalue)  newvalue = tempvalue;
        } /* end for */

        
        /* If it was a FARAWAY point, add it to the NarrowBand Heap;
         * otherwise, just update its value
         * using the backpointer
         */
        if(label[neighVId] == (PGbyte)NBAND)
          xhChangeValue(BackPointer[neighVId],newvalue, H);
        else{
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);
          label[neighVId] = (PGbyte)NBAND;
        }
          
      } /* end if ALIVE */
    } /* end updating neighbors */

  } /* end of marching loop */   
             
  free(label);
  free(BackPointer);
  xhDestroy(H);

return T;
}


double *SearchConnectingVertexUsingFastMarching(Surface *mesh, int startVertexID, double *F, int *faceSulci, int geodesicThresh, int sulcusLabel, int *vertexID)
{
  double dummtT;
  double *T;  /* final distances of each vertex to the contour */
  double tempvalue, newvalue;
  unsigned char *label;      /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */  
  int *BackPointer;  /* backpointer to the narrowband heap */
  Xheap H;           /* Narrrow band heap */
  XheapElement he;
  int heapindex;
  double currentF;

  int i, j, k, NN, NN1, NN2, VN, VId, neighVId, n0, n1;
  int m;

  //VN = XUGetVertexTableVN(mesh);
  VN = mesh->vertexNum;
  T = (double *) malloc(sizeof(double)*VN);
  label = (unsigned char *) malloc(sizeof(unsigned char)*VN);
  BackPointer = (int *) malloc(sizeof(int)*VN);

  /* Initialize heap structure */
  H = xhInitEmpty();

  /* Initialization for marching inwards and outwards simultaneously */
  for (i=0; i<VN; i++) {
    label[i] = (unsigned char) FAWAY;  /*All points are labelled as FAR AWAY */
    T[i] = INFINITY;
  }

  /*
  for (i=0; i<numinitvert; i++) {
    label[contour[i]] = (unsigned char) ALIVE; // all vertices on contour are alive 
    T[contour[i]] = 0; // all vertices on contour have distance zero 
  }
*/
	 label[startVertexID] = (unsigned char)ALIVE;
	 T[startVertexID] = 0.0;


  /* Initialize NarrowBand Heap */

  for (i=0; i<VN; i++) {

    if (label[i] == (unsigned char)ALIVE) {

      /* put its neighbors into the NarrowBand */ 
      //NN1 = XUGetNeighborEN(mesh, i);
		NN1 = mesh->neighbors[i].size();
      for (j=0; j<NN1; j++) {
        //neighVId = XUGetVertexNeighborVId(mesh, i, j);
		  neighVId = mesh->neighbors[i][j];

        if (label[neighVId] == (unsigned char)FAWAY) {
          /* compute the distance to this point and add it to the heap */ 
          label[neighVId] = (unsigned char)NBAND;
          /* Check all possible triangles surrounding vertex neighVId1 */
          /* Note: Only ALIVE points contribute to the distance computation */
          //NN2 = XUGetNeighborEN(mesh, neighVId);
/*
		  NN2 = mesh->neighbors[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
			n0 = mesh->neighbors[neighVId][k];
            n1 = mesh->neighbors[neighVId][(k+1)%NN2];

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);
            if (tempvalue < newvalue)  newvalue = tempvalue;
          }
*/ 
		  NN2 = mesh->adjacentfaces[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
		   
//			n0 = mesh->neighbors[neighVId][k];
//          n1 = mesh->neighbors[neighVId][(k+1)%NN2];

			GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][k], neighVId, &n0, &n1);

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);

            if (tempvalue < newvalue)  newvalue = tempvalue;
          }

          T[neighVId] = newvalue;
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);

        } /* end if */
      } /* end for */    
    } /* end if */
  } /* end for */

  /* End of Initialization */
 
  /* Begin Fast Marching to get the unsigned distance function inwords
   * and outwards simultaneously
   * since points inside and outside the contour won't interfere with
   * each other
   */
	bool findVertex = false;

  //printf("   "); 
  m = 0;

  while(!xhIsEmpty(H)){ /* There are still points not yet accepted */
    he = xhRemove(H); /* Label the point with smallest value among all
                         NarrowBand points as ALIVE */
    m++;
    VId = he.id; 

    T[VId] = he.value;

    label[VId] = (unsigned char)ALIVE;
//    if ( T[terminator]*2 < he.value )
//	  break; 
	
	if(he.value >= 1.0*geodesicThresh)
		break;

	for(i=0; i<mesh->adjacentfaces[VId].size(); i++)
	{
		if(faceSulci[mesh->adjacentfaces[VId][i]] >= 0 && faceSulci[mesh->adjacentfaces[VId][i]] != sulcusLabel)
		{
			*vertexID = VId;
			findVertex = true;
			break;
		}
	}
    
	if(findVertex)	break;

    /* Update its neighbor */
    /* Put FARAWAY neighbors into NarrowBand, Recompute values at
     * NarrowBand neighbors,
     * Keep ALIVE (Accepted) neighbor unchanged
     */
//    NN1 = XUGetNeighborEN(mesh, VId);
	NN1 = mesh->neighbors[VId].size();
    for (i=0; i<NN1; i++) {
      //neighVId = XUGetVertexNeighborVId(mesh, VId, i);
		neighVId = mesh->neighbors[VId][i];
      
	  /* Don't change ALIVE neighbors */
      if (label[neighVId] != (unsigned char)ALIVE) { 
/*
//       NN2 = XUGetNeighborEN(mesh, neighVId);
		NN2 = mesh->neighbors[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {
//          n0 = XUGetVertexNeighborVId(mesh, neighVId, j);
//          n1 = XUGetVertexNeighborVId(mesh, neighVId, (j+1)%NN2);
		  n0 = mesh->neighbors[neighVId][j];
		  n1 = mesh->neighbors[neighVId][(j+1)%NN2];
*/
		NN2 = mesh->adjacentfaces[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {

		GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][j], neighVId, &n0, &n1);

		  tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				BackPointer, &H);
          if (tempvalue < newvalue)  newvalue = tempvalue;
        } /* end for */

        
        /* If it was a FARAWAY point, add it to the NarrowBand Heap;
         * otherwise, just update its value
         * using the backpointer
         */
        if(label[neighVId] == (PGbyte)NBAND)
          xhChangeValue(BackPointer[neighVId],newvalue, H);
        else{
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);
          label[neighVId] = (PGbyte)NBAND;
        }
          
      } /* end if ALIVE */
    } /* end updating neighbors */

  } /* end of marching loop */   
             
  free(label);
  free(BackPointer);
  xhDestroy(H);

return T;
}



double *FastMarchMeshHalfVertex(Surface *mesh, VtxOnFace start, VtxOnFace end, double *F)
{
  double dummtT;
  double *T;  /* final distances of each vertex to the contour */
  double tempvalue, newvalue;
  unsigned char *label;      /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */  
  int *BackPointer;  /* backpointer to the narrowband heap */
  Xheap H;           /* Narrrow band heap */
  XheapElement he;
  int heapindex;
  double currentF;

  int i, j, k, NN, NN1, NN2, VN, VId, neighVId, n0, n1;
  int m;

  //VN = XUGetVertexTableVN(mesh);
  VN = mesh->vertexNum;
  T = (double *) malloc(sizeof(double)*VN);
  label = (unsigned char *) malloc(sizeof(unsigned char)*VN);
  BackPointer = (int *) malloc(sizeof(int)*VN);

  /* Initialize heap structure */
  H = xhInitEmpty();

  /* Initialization for marching inwards and outwards simultaneously */
  for (i=0; i<VN; i++) {
    label[i] = (unsigned char) FAWAY;  /*All points are labelled as FAR AWAY */
    T[i] = INFINITY;
  }

  if(start.lambda == 0.0)
  {
	  label[start.ID] = (unsigned char)ALIVE;
	  T[start.ID] = 0.0;
  }
  else
  {
	  Fvector3d halfVertex;

	  halfVertex.x = (1.0-start.lambda)*mesh->vertex[mesh->edges[start.ID].vertexID[0]].x + start.lambda*mesh->vertex[mesh->edges[start.ID].vertexID[1]].x;
	  halfVertex.y = (1.0-start.lambda)*mesh->vertex[mesh->edges[start.ID].vertexID[0]].y + start.lambda*mesh->vertex[mesh->edges[start.ID].vertexID[1]].y;
	  halfVertex.z = (1.0-start.lambda)*mesh->vertex[mesh->edges[start.ID].vertexID[0]].z + start.lambda*mesh->vertex[mesh->edges[start.ID].vertexID[1]].z;

	  int tempFaceID;

	  for(i=0; i<2; i++)
	  {
		  tempFaceID = mesh->edges[start.ID].faceID[i];
		  
		  if(tempFaceID >= 0)
		  {
			  for(j=0; j<3; j++)
			  {
					label[mesh->faces[tempFaceID].v[j]] = (unsigned char)ALIVE;
					T[mesh->faces[tempFaceID].v[j]] = compute_dist( &halfVertex, &(mesh->vertex[mesh->faces[tempFaceID].v[j]]) );
				}
		  }
	  }
  }



  /* Initialize NarrowBand Heap */

  for (i=0; i<VN; i++) {

    if (label[i] == (unsigned char)ALIVE) {

      /* put its neighbors into the NarrowBand */ 
      //NN1 = XUGetNeighborEN(mesh, i);
		NN1 = mesh->neighbors[i].size();
      for (j=0; j<NN1; j++) {
        //neighVId = XUGetVertexNeighborVId(mesh, i, j);
		  neighVId = mesh->neighbors[i][j];

        if (label[neighVId] == (unsigned char)FAWAY) {
          /* compute the distance to this point and add it to the heap */ 
          label[neighVId] = (unsigned char)NBAND;
          /* Check all possible triangles surrounding vertex neighVId1 */
          /* Note: Only ALIVE points contribute to the distance computation */
          //NN2 = XUGetNeighborEN(mesh, neighVId);
/*
		  NN2 = mesh->neighbors[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
			n0 = mesh->neighbors[neighVId][k];
            n1 = mesh->neighbors[neighVId][(k+1)%NN2];

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);
            if (tempvalue < newvalue)  newvalue = tempvalue;
          }
*/ 
		  NN2 = mesh->adjacentfaces[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
		   
//			n0 = mesh->neighbors[neighVId][k];
//          n1 = mesh->neighbors[neighVId][(k+1)%NN2];

			GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][k], neighVId, &n0, &n1);

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);

            if (tempvalue < newvalue)  newvalue = tempvalue;
          }

          T[neighVId] = newvalue;
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);

        } /* end if */
      } /* end for */    
    } /* end if */
  } /* end for */

  /* End of Initialization */
 
  /* Begin Fast Marching to get the unsigned distance function inwords
   * and outwards simultaneously
   * since points inside and outside the contour won't interfere with
   * each other
   */
	bool findVertex = false;

  //printf("   "); 
  m = 0;

  int vertexID;
  while(!xhIsEmpty(H))
  { /* There are still points not yet accepted */
    he = xhRemove(H); /* Label the point with smallest value among all
                         NarrowBand points as ALIVE */
    m++;
    VId = he.id; 

    T[VId] = he.value;

    label[VId] = (unsigned char)ALIVE;
	
	if(end.lambda == 0.0)
	{
		if ( T[end.ID]*2 < he.value )
			break;
	}
	else if(T[mesh->edges[end.ID].vertexID[0]] < INFINITY && T[mesh->edges[end.ID].vertexID[1]] < INFINITY)
		break;

    /* Update its neighbor */
    /* Put FARAWAY neighbors into NarrowBand, Recompute values at
     * NarrowBand neighbors,
     * Keep ALIVE (Accepted) neighbor unchanged
     */
//    NN1 = XUGetNeighborEN(mesh, VId);
	NN1 = mesh->neighbors[VId].size();
    for (i=0; i<NN1; i++) {
      //neighVId = XUGetVertexNeighborVId(mesh, VId, i);
		neighVId = mesh->neighbors[VId][i];
      
	  /* Don't change ALIVE neighbors */
      if (label[neighVId] != (unsigned char)ALIVE) { 
/*
//       NN2 = XUGetNeighborEN(mesh, neighVId);
		NN2 = mesh->neighbors[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {
//          n0 = XUGetVertexNeighborVId(mesh, neighVId, j);
//          n1 = XUGetVertexNeighborVId(mesh, neighVId, (j+1)%NN2);
		  n0 = mesh->neighbors[neighVId][j];
		  n1 = mesh->neighbors[neighVId][(j+1)%NN2];
*/
		NN2 = mesh->adjacentfaces[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {

		GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][j], neighVId, &n0, &n1);

		  tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				BackPointer, &H);
          if (tempvalue < newvalue)  newvalue = tempvalue;
        } /* end for */

        
        /* If it was a FARAWAY point, add it to the NarrowBand Heap;
         * otherwise, just update its value
         * using the backpointer
         */
        if(label[neighVId] == (PGbyte)NBAND)
          xhChangeValue(BackPointer[neighVId],newvalue, H);
        else{
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);
          label[neighVId] = (PGbyte)NBAND;
        }
          
      } /* end if ALIVE */
    } /* end updating neighbors */

  } /* end of marching loop */   
             
  free(label);
  free(BackPointer);
  xhDestroy(H);

return T;
}



double *SearchConnectingHalfVertexUsingFastMarching(Surface *mesh, Segment *segment, VtxOnFace start, double *F, 
													int *edgeSulci,Fvector3d *edgeHalfVertex, int geodesicThresh, int sulcusLabel, VtxOnFace *end)
{
  double dummtT;
  double *T;  /* final distances of each vertex to the contour */
  double tempvalue, newvalue;
  unsigned char *label;      /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */  
  int *BackPointer;  /* backpointer to the narrowband heap */
  Xheap H;           /* Narrrow band heap */
  XheapElement he;
  int heapindex;
  double currentF;

  int i, j, k, NN, NN1, NN2, VN, VId, neighVId, n0, n1;
  int m;

  //VN = XUGetVertexTableVN(mesh);
  VN = mesh->vertexNum;
  T = (double *) malloc(sizeof(double)*VN);
  label = (unsigned char *) malloc(sizeof(unsigned char)*VN);
  BackPointer = (int *) malloc(sizeof(int)*VN);

  /* Initialize heap structure */
  H = xhInitEmpty();

  /* Initialization for marching inwards and outwards simultaneously */
  for (i=0; i<VN; i++) {
    label[i] = (unsigned char) FAWAY;  /*All points are labelled as FAR AWAY */
    T[i] = INFINITY;
  }

  /*
  for (i=0; i<numinitvert; i++) {
    label[contour[i]] = (unsigned char) ALIVE; // all vertices on contour are alive 
    T[contour[i]] = 0; // all vertices on contour have distance zero 
  }
*/
  if(start.lambda == 0.0)
  {
	  label[start.ID] = (unsigned char)ALIVE;
	  T[start.ID] = 0.0;
  }
  else
  {
	  Fvector3d halfVertex;

	  halfVertex.x = (1.0-start.lambda)*mesh->vertex[mesh->edges[start.ID].vertexID[0]].x + start.lambda*mesh->vertex[mesh->edges[start.ID].vertexID[1]].x;
	  halfVertex.y = (1.0-start.lambda)*mesh->vertex[mesh->edges[start.ID].vertexID[0]].y + start.lambda*mesh->vertex[mesh->edges[start.ID].vertexID[1]].y;
	  halfVertex.z = (1.0-start.lambda)*mesh->vertex[mesh->edges[start.ID].vertexID[0]].z + start.lambda*mesh->vertex[mesh->edges[start.ID].vertexID[1]].z;


	  int tempFaceID;

	  for(i=0; i<2; i++)
	  {
		  tempFaceID = mesh->edges[start.ID].faceID[i];
		  
		  if(tempFaceID >= 0)
		  {
			  for(j=0; j<3; j++)
			  {
					label[mesh->faces[tempFaceID].v[j]] = (unsigned char)ALIVE;
					T[mesh->faces[tempFaceID].v[j]] = compute_dist( &halfVertex, &(mesh->vertex[mesh->faces[tempFaceID].v[j]]) );
				}
		  }
	  }

  }



  /* Initialize NarrowBand Heap */

  for (i=0; i<VN; i++) {

    if (label[i] == (unsigned char)ALIVE) {

      /* put its neighbors into the NarrowBand */ 
      //NN1 = XUGetNeighborEN(mesh, i);
		NN1 = mesh->neighbors[i].size();
      for (j=0; j<NN1; j++) {
        //neighVId = XUGetVertexNeighborVId(mesh, i, j);
		  neighVId = mesh->neighbors[i][j];

        if (label[neighVId] == (unsigned char)FAWAY) {
          /* compute the distance to this point and add it to the heap */ 
          label[neighVId] = (unsigned char)NBAND;
          /* Check all possible triangles surrounding vertex neighVId1 */
          /* Note: Only ALIVE points contribute to the distance computation */
          //NN2 = XUGetNeighborEN(mesh, neighVId);
/*
		  NN2 = mesh->neighbors[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
			n0 = mesh->neighbors[neighVId][k];
            n1 = mesh->neighbors[neighVId][(k+1)%NN2];

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);
            if (tempvalue < newvalue)  newvalue = tempvalue;
          }
*/ 
		  NN2 = mesh->adjacentfaces[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
		   
//			n0 = mesh->neighbors[neighVId][k];
//          n1 = mesh->neighbors[neighVId][(k+1)%NN2];

			GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][k], neighVId, &n0, &n1);

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);

            if (tempvalue < newvalue)  newvalue = tempvalue;
          }

          T[neighVId] = newvalue;
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);

        } /* end if */
      } /* end for */    
    } /* end if */
  } /* end for */

  /* End of Initialization */
 
  /* Begin Fast Marching to get the unsigned distance function inwords
   * and outwards simultaneously
   * since points inside and outside the contour won't interfere with
   * each other
   */
	bool findVertex = false;

  //printf("   "); 
  m = 0;

  int vertexID;
  while(!xhIsEmpty(H))
  { /* There are still points not yet accepted */
    he = xhRemove(H); /* Label the point with smallest value among all
                         NarrowBand points as ALIVE */
    m++;
    VId = he.id; 

    T[VId] = he.value;

    label[VId] = (unsigned char)ALIVE;
//    if ( T[terminator]*2 < he.value )
//	  break; 
	
	if(he.value > 1.0*geodesicThresh)
		break;

	if(!findVertex)
	{
		for(i=0; i<mesh->connectedges[VId].size(); i++)
		{
			if(edgeSulci[mesh->connectedges[VId][i]] >= 0 && edgeSulci[mesh->connectedges[VId][i]] != sulcusLabel)
			{
				vertexID = VId;
				findVertex = true;
				break;
			}
		}
	}
    
	if(findVertex)	
	{
		bool visitAllVertex = true;

		for(i=0; i<mesh->neighbors[vertexID].size(); i++)
		{
			if(T[mesh->neighbors[vertexID][i]] == INFINITY)
			{
				visitAllVertex = false;
				break;
			}
		}

		if(visitAllVertex)
		{
			float minGeoDist = 100000.0;

			for(i=0; i<mesh->connectedges[vertexID].size(); i++)
			{
				if(edgeSulci[mesh->connectedges[vertexID][i]] >= 0 && edgeSulci[mesh->connectedges[vertexID][i]] != sulcusLabel)
				{
					float dist1 = compute_dist( &(mesh->vertex[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[0]]), &(mesh->vertex[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[1]]) );
					float dist2 = compute_dist( &(edgeHalfVertex[mesh->connectedges[vertexID][i]]), &(mesh->vertex[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[0]]) );
					float lamda =  dist2/dist1;

					float geoDist = (1.0-lamda)*T[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[0]] + (lamda)*T[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[1]];

					if(geoDist < minGeoDist)
					{
						minGeoDist = geoDist;
						end->lambda = lamda;
						end->ID = mesh->connectedges[vertexID][i];
					}
				}
			}
			break;
		}
	}
    /* Update its neighbor */
    /* Put FARAWAY neighbors into NarrowBand, Recompute values at
     * NarrowBand neighbors,
     * Keep ALIVE (Accepted) neighbor unchanged
     */
//    NN1 = XUGetNeighborEN(mesh, VId);
	NN1 = mesh->neighbors[VId].size();
    for (i=0; i<NN1; i++) {
      //neighVId = XUGetVertexNeighborVId(mesh, VId, i);
		neighVId = mesh->neighbors[VId][i];
      
	  /* Don't change ALIVE neighbors */
      if (label[neighVId] != (unsigned char)ALIVE) { 
/*
//       NN2 = XUGetNeighborEN(mesh, neighVId);
		NN2 = mesh->neighbors[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {
//          n0 = XUGetVertexNeighborVId(mesh, neighVId, j);
//          n1 = XUGetVertexNeighborVId(mesh, neighVId, (j+1)%NN2);
		  n0 = mesh->neighbors[neighVId][j];
		  n1 = mesh->neighbors[neighVId][(j+1)%NN2];
*/
		NN2 = mesh->adjacentfaces[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {

		GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][j], neighVId, &n0, &n1);

		  tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				BackPointer, &H);
          if (tempvalue < newvalue)  newvalue = tempvalue;
        } /* end for */

        
        /* If it was a FARAWAY point, add it to the NarrowBand Heap;
         * otherwise, just update its value
         * using the backpointer
         */
        if(label[neighVId] == (PGbyte)NBAND)
          xhChangeValue(BackPointer[neighVId],newvalue, H);
        else{
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);
          label[neighVId] = (PGbyte)NBAND;
        }
          
      } /* end if ALIVE */
    } /* end updating neighbors */

  } /* end of marching loop */   
             
  free(label);
  free(BackPointer);
  xhDestroy(H);

return T;
}




double *SearchConnectingHalfVertexAndVertex(Surface *mesh, Segment *segment, SulciTrackingOut *out, VtxOnFace start, double *F, 
													int *edgeSulci,int *vertexSulci, Fvector3d *edgeHalfVertex, int geodesicThresh, int sulcusLabel, VtxOnFace *end)
{
  double dummtT;
  double *T;  /* final distances of each vertex to the contour */
  double tempvalue, newvalue;
  unsigned char *label;      /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */  
  int *BackPointer;  /* backpointer to the narrowband heap */
  Xheap H;           /* Narrrow band heap */
  XheapElement he;
  int heapindex;
  double currentF;

  int i, j, k, NN, NN1, NN2, VN, VId, neighVId, n0, n1;
  int m;
  Fvector3d halfVertex;
  float dist;

  VN = mesh->vertexNum;
  T = (double *) malloc(sizeof(double)*VN);
  label = (unsigned char *) malloc(sizeof(unsigned char)*VN);
  BackPointer = (int *) malloc(sizeof(int)*VN);

  /* Initialize heap structure */
  H = xhInitEmpty();

  /* Initialization for marching inwards and outwards simultaneously */
  for (i=0; i<VN; i++) {
    label[i] = (unsigned char) FAWAY;  /*All points are labelled as FAR AWAY */
    T[i] = INFINITY;
  }

  /*
  for (i=0; i<numinitvert; i++) {
    label[contour[i]] = (unsigned char) ALIVE; // all vertices on contour are alive 
    T[contour[i]] = 0; // all vertices on contour have distance zero 
  }
*/
  if(start.lambda == 0.0)
  {
	  label[start.ID] = (unsigned char)ALIVE;
	  T[start.ID] = 0.0;
  }
  else
  {
	  halfVertex.x = edgeHalfVertex[start.ID].x;
	  halfVertex.y = edgeHalfVertex[start.ID].y;
	  halfVertex.z = edgeHalfVertex[start.ID].z;

	  int tempFaceID;

	  for(i=0; i<2; i++)
	  {
		  tempFaceID = mesh->edges[start.ID].faceID[i];
		  
		  if(tempFaceID >= 0)
		  {
			  for(j=0; j<3; j++)
			  {
					label[mesh->faces[tempFaceID].v[j]] = (unsigned char)ALIVE;
					T[mesh->faces[tempFaceID].v[j]] = compute_dist( &halfVertex, &(mesh->vertex[mesh->faces[tempFaceID].v[j]]) );
				}
		  }
	  }
  }



  /* Initialize NarrowBand Heap */

  for (i=0; i<VN; i++) {

    if (label[i] == (unsigned char)ALIVE) {

      /* put its neighbors into the NarrowBand */ 
      //NN1 = XUGetNeighborEN(mesh, i);
		NN1 = mesh->neighbors[i].size();
      for (j=0; j<NN1; j++) {
        //neighVId = XUGetVertexNeighborVId(mesh, i, j);
		  neighVId = mesh->neighbors[i][j];

        if (label[neighVId] == (unsigned char)FAWAY) {
          /* compute the distance to this point and add it to the heap */ 
          label[neighVId] = (unsigned char)NBAND;
          /* Check all possible triangles surrounding vertex neighVId1 */
          /* Note: Only ALIVE points contribute to the distance computation */
          //NN2 = XUGetNeighborEN(mesh, neighVId);
/*
		  NN2 = mesh->neighbors[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
			n0 = mesh->neighbors[neighVId][k];
            n1 = mesh->neighbors[neighVId][(k+1)%NN2];

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);
            if (tempvalue < newvalue)  newvalue = tempvalue;
          }
*/ 
		  NN2 = mesh->adjacentfaces[neighVId].size();
          newvalue = INFINITY;
		  currentF = F[neighVId];

          for (k=0; k<NN2; k++) {
//			n0 = XUGetVertexNeighborVId(mesh, neighVId, k);
//			n1 = XUGetVertexNeighborVId(mesh, neighVId, (k+1)%NN2);
		   
//			n0 = mesh->neighbors[neighVId][k];
//          n1 = mesh->neighbors[neighVId][(k+1)%NN2];

			GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][k], neighVId, &n0, &n1);

            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				   BackPointer, &H);

            if (tempvalue < newvalue)  newvalue = tempvalue;
          }

          T[neighVId] = newvalue;
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);

        } /* end if */
      } /* end for */    
    } /* end if */
  } /* end for */

  /* End of Initialization */
 
  /* Begin Fast Marching to get the unsigned distance function inwords
   * and outwards simultaneously
   * since points inside and outside the contour won't interfere with
   * each other
   */


  //printf("   "); 
  m = 0;
  bool findVertex = false;
  int vertexID;
  while(!xhIsEmpty(H))
  { /* There are still points not yet accepted */
    he = xhRemove(H); /* Label the point with smallest value among all
                         NarrowBand points as ALIVE */
    m++;
    VId = he.id; 

    T[VId] = he.value;

    label[VId] = (unsigned char)ALIVE;
	
//	if(he.value > 1.0*geodesicThresh)
//		break;

	 if(start.lambda == 0.0)
	{
		dist = MyDistance(mesh->vertex[start.ID], mesh->vertex[he.id]);

		if(dist > geodesicThresh)
			break;
	}
	else
	{
		dist = MyDistance(halfVertex, mesh->vertex[he.id]);

		if(dist > geodesicThresh)
			break;
	}

	if(!findVertex)
	{
		if(vertexSulci[VId] >=0 && out->sulci[vertexSulci[VId]].color != sulcusLabel)
		{
			vertexID = VId;
			findVertex = true;
		}
		else
		{
			// edge
			for(i=0; i<mesh->adjacentedges[VId].size(); i++)
			{
				if(edgeSulci[mesh->adjacentedges[VId][i]] >= 0 && out->sulci[edgeSulci[mesh->adjacentedges[VId][i]]].color != sulcusLabel)
				{
					vertexID = VId;
					findVertex = true;
					break;
				}
			}
		}
	}
    
	if(findVertex)	
	{
		bool visitAllVertex = true;

		for(i=0; i<mesh->neighbors[vertexID].size(); i++)
		{
			if(T[mesh->neighbors[vertexID][i]] == INFINITY)
			{
				visitAllVertex = false;
				break;
			}
		}

		if(visitAllVertex)
		{
			float minGeoDist = 100000.0;
			bool stop = false;
			
			// vertex
			if(vertexSulci[vertexID] >= 0)
			{
				end->lambda = 0.0;
				end->ID = vertexID;
				stop = true;
				break;
			}

			// edge
			for(i=0; i<mesh->connectedges[vertexID].size(); i++)
			{
				if(edgeSulci[mesh->connectedges[vertexID][i]] >= 0 && out->sulci[edgeSulci[mesh->connectedges[vertexID][i]]].color != sulcusLabel)
				{
					float dist1 = compute_dist( &(mesh->vertex[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[0]]), &(mesh->vertex[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[1]]) );
					float dist2 = compute_dist( &(edgeHalfVertex[mesh->connectedges[vertexID][i]]), &(mesh->vertex[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[0]]) );
					float lamda =  dist2/dist1;

					float geoDist = (1.0-lamda)*T[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[0]] + (lamda)*T[mesh->edges[mesh->connectedges[vertexID][i]].vertexID[1]];

					if(geoDist < minGeoDist)
					{
						minGeoDist = geoDist;
						end->lambda = lamda;
						end->ID = mesh->connectedges[vertexID][i];
						stop = true;
					}
				}
			}
			if(stop)
				break;

			break;
		}
	}
    /* Update its neighbor */
    /* Put FARAWAY neighbors into NarrowBand, Recompute values at
     * NarrowBand neighbors,
     * Keep ALIVE (Accepted) neighbor unchanged
     */
//    NN1 = XUGetNeighborEN(mesh, VId);
	NN1 = mesh->neighbors[VId].size();
    for (i=0; i<NN1; i++) {
      //neighVId = XUGetVertexNeighborVId(mesh, VId, i);
		neighVId = mesh->neighbors[VId][i];
      
	  /* Don't change ALIVE neighbors */
      if (label[neighVId] != (unsigned char)ALIVE) { 
/*
//       NN2 = XUGetNeighborEN(mesh, neighVId);
		NN2 = mesh->neighbors[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {
//          n0 = XUGetVertexNeighborVId(mesh, neighVId, j);
//          n1 = XUGetVertexNeighborVId(mesh, neighVId, (j+1)%NN2);
		  n0 = mesh->neighbors[neighVId][j];
		  n1 = mesh->neighbors[neighVId][(j+1)%NN2];
*/
		NN2 = mesh->adjacentfaces[neighVId].size();
        newvalue = INFINITY;
		currentF = F[neighVId];
        for (j=0; j<NN2; j++) {

		GetTwoVertex(mesh, mesh->adjacentfaces[neighVId][j], neighVId, &n0, &n1);

		  tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label,currentF,
				BackPointer, &H);
          if (tempvalue < newvalue)  newvalue = tempvalue;
        } /* end for */

        
        /* If it was a FARAWAY point, add it to the NarrowBand Heap;
         * otherwise, just update its value
         * using the backpointer
         */
        if(label[neighVId] == (PGbyte)NBAND)
          xhChangeValue(BackPointer[neighVId],newvalue, H);
        else{
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);
          label[neighVId] = (PGbyte)NBAND;
        }
          
      } /* end if ALIVE */
    } /* end updating neighbors */

  } /* end of marching loop */   
             
  free(label);
  free(BackPointer);
  xhDestroy(H);

return T;
}





/* vIDa, vIDb & vIDc are the three vertices making up the triangle */
/* mesh is the mesh to compute the distances on */
/* T is the current transit time values */
/* label is the label map */
/* compute the new T for vIDc based on the current T values of vIDa & vIDb */
/* Unfold the mesh whenever C is an obtuse angle */
/* see Kimmel and Sethian for derivation of these formulae */
double ReCompute(int vIDc, int vIDa, int vIDb, Surface *mesh, double *T,
		 unsigned char *label, double Fc, int *Backpointer, Xheap *pH) 
{
  double u, a, b, c;
  Fvector3d *Va, *Vb, *Vc;
  unsigned char la, lb;
  double Ta, Tb;
  double t1, t2, t;
  int tmpint;
  double tmpdouble;

  int tmpvIDa;
  double tmpa;
  double tmpTa;
  unsigned char tmpla;

  /* Define some auxillary variables */
  int NN, neighVId, i, UNotID, UID, P1ID, P2ID;

  /* SQUARES of lengthes (Captital Letters) */
  double P1A, P1B, P1C, P2A, P2B, P2C, UA, UB, UC, P1P2, P1U, P2U; 
  double cos12A, sin12A, cos12B, sin12B, cos12C, sin12C, cos12U, sin12U;
  double AC, BC, AB;
	
  Fvector3d *P1, *P2, *U;
  Fvector3d vP1U, vP2U;
 
  int iters; 
	 
  la = label[vIDa];
  lb = label[vIDb];

  if ( (la != (unsigned char) ALIVE) && (lb != (unsigned char) ALIVE) ) 
    return INFINITY;
  
  if(Fc == 0) return INFINITY;

  Ta = T[vIDa];
  Tb = T[vIDb];
  
  if(la == (unsigned char) ALIVE){
    if(Ta == INFINITY) return INFINITY;
  }
  else Ta = INFINITY;
  
  if(lb == (unsigned char) ALIVE){
    if(Tb == INFINITY) return INFINITY;
  }
  else Tb = INFINITY;
/*
  Va = XUGetVertexPtr(mesh, vIDa);
  Vb = XUGetVertexPtr(mesh, vIDb);
  Vc = XUGetVertexPtr(mesh, vIDc);
*/
  Va = &(mesh->vertex[vIDa]);
  Vb = &(mesh->vertex[vIDb]);
  Vc = &(mesh->vertex[vIDc]);

  a = compute_dist(Vb, Vc);
  b = compute_dist(Va, Vc);
  
  if(Ta > Tb){
    tmpTa = Ta;
    Ta = Tb;
    Tb = tmpTa;
    
    tmpla = la;
    la = lb;
    lb = tmpla;
    
    tmpa = a;
    a = b;
    b = tmpa;
    
    tmpvIDa = vIDa;
    vIDa = vIDb;
    vIDb = tmpvIDa;
/*    
    Va = XUGetVertexPtr(mesh, vIDa);
    Vb = XUGetVertexPtr(mesh, vIDb);
*/	
	Va = &(mesh->vertex[vIDa]);
	Vb = &(mesh->vertex[vIDb]);
  }
  
  c =  compute_dist(Va, Vb);
  AC = b*b; 
  BC = a*a;
  AB = c*c;
  if(AB < AC+BC)
  {/* Acute Triangle */
    if(la == ALIVE && lb == ALIVE){
      t = ComputeTInAcute(Ta, Tb, a,  b,  c, Fc);
      return t;
    }
    else return (MIN(Ta+b/Fc, Tb+a/Fc));  /* Infinity + sth = Infinity */
  }

  /* Otherwise, perform unfolding */
  
  /*Initialization */
  P1ID = vIDa;
  P2ID = vIDb;
  UNotID = vIDc;
  
  P1A = 0; P1B = AB; P1C = AC;
  P2B = 0; P2A = P1B; P2C = BC;
  P1P2 = P1B;
  
  cos12A = 1; sin12A = 0; cos12B = 1; sin12B = 0; 
  cos12C = (P1P2 + P2C - P1C)/(2*sqrt(P1P2*P2C));
  sin12C = (1-cos12C*cos12C); /* Notice: Square of sine */
  
  int leftVertexID, rightVertexID;

  /* Now iteratively unfolding */
  iters = 0;
  while(iters <10){
    /* Find the newly unfolded vertex ID */
    //NN = XUGetNeighborEN(mesh, P1ID);
	  //NN =  mesh->adjacentedges[P1ID].size();
	  NN =  mesh->neighbors[P1ID].size();
    for (i=0; i<NN; i++) {
//      neighVId = XUGetVertexNeighborVId(mesh, P1ID, i);
		neighVId = mesh->neighbors[P1ID][i];
      if(neighVId == P2ID)
	  {
/*	
	//if(XUGetVertexNeighborVId(mesh, P1ID, (i+NN-1)%NN) == UNotID){
	if(mesh->neighbors[P1ID][(i+NN-1)%NN] == UNotID)
	{
	  //UID = XUGetVertexNeighborVId(mesh, P1ID, (i+1)%NN);
		UID = mesh->neighbors[P1ID][(i+1)%NN];
	  break;
	}
	//if(XUGetVertexNeighborVId(mesh, P1ID, (i+1)%NN) == UNotID){
	if(mesh->neighbors[P1ID][(i+1)%NN] == UNotID)
	{
	  //UID = XUGetVertexNeighborVId(mesh, P1ID, (i+NN-1)%NN);
		UID = mesh->neighbors[P1ID][(i+NN-1)%NN];
	  break;
	}
*/		


		  leftVertexID = -1;
		  rightVertexID = -1;
		GetLeftRightVertex(mesh, P1ID, P2ID, &leftVertexID, &rightVertexID);

		if(leftVertexID == UNotID)
		{
			UID = rightVertexID;
			break;
		}
		if(rightVertexID == UNotID)
		{	
			UID = leftVertexID;
			break;
		}

      }
    }

/*    
    P1 = XUGetVertexPtr(mesh, P1ID); 
    P2 = XUGetVertexPtr(mesh, P2ID); 
    U = XUGetVertexPtr(mesh, UID);
*/
	if(UID < 0)
	{
		printf("UID < 0 in ReCompute\n");
		break;
	}

	P1 = &(mesh->vertex[P1ID]); 
    P2 = &(mesh->vertex[P2ID]); 
    U = &(mesh->vertex[UID]);

//    vP1U = vectsubptr(P1, U);
//    vP2U = vectsubptr(P2, U);
	vP1U.x = P1->x - U->x;
	vP1U.y = P1->y - U->y;
	vP1U.z = P1->z - U->z;
	vP2U.x = P2->x - U->x;
	vP2U.y = P2->y - U->y;
	vP2U.z = P2->z - U->z;
    
    P1U = (vP1U.x * vP1U.x + vP1U.y * vP1U.y + vP1U.z * vP1U.z);
    P2U = (vP2U.x * vP2U.x + vP2U.y * vP2U.y + vP2U.z * vP2U.z);
    
    cos12U = (P1P2 + P2U - P1U)/(2*sqrt(P1P2*P2U));
    sin12U = (1-cos12U*cos12U); /* Notice: Square of sine */
    
    /* Now compute three lengthes (squared) */
    UA = P2U + P2A - 2*sqrt(P2U*P2A)*(cos12A*cos12U - sqrt(sin12A*sin12U));
    UB = P2U + P2B - 2*sqrt(P2U*P2B)*(cos12B*cos12U - sqrt(sin12B*sin12U));
    UC = P2U + P2C - 2*sqrt(P2U*P2C)*(cos12C*cos12U - sqrt(sin12C*sin12U));
    
    /* Now Judge Which Side to continue unfolding */
    if(UA > (UC + AC))
	{/* Unfold along P1U */
      UNotID = P2ID;
      P2ID = UID;
      P1P2 = P1U;
      P2A = UA; P2B = UB; P2C = UC;
    }
	else if(UB > (UC + BC))
	{ /* Unfold along P2U */
      UNotID = P1ID;
      P1ID = UID;
      P1P2 = P2U;
      P1A = UA; P1B = UB; P1C = UC;
    }
	else
	{ /* Stop Unfolding and compute T*/
      /* Compute the actual lengthes */
      UC = sqrt(UC);
      UA = sqrt(UA);
      UB = sqrt(UB);
      if(label[UID] == (PGbyte) ALIVE)
	  {
	if(la == ALIVE){
	  t1 = ComputeTInAcute(Ta, T[UID], UC, b, UA, Fc);
	  if(t1 < Ta){ /* Reset A to NBAND and put into Heap */
	    xhInsert(Ta, vIDa, &(Backpointer[vIDa]), *pH);
	    label[vIDa] = (PGbyte)NBAND;
	  }
	}else t1 = INFINITY;
        if(lb == ALIVE){
	  t2 = ComputeTInAcute(Tb, T[UID], UC, a, UB, Fc);
	  if(t2 < Tb){ /* Reset B to NBAND and put into Heap */
	    xhInsert(Tb, vIDb, &(Backpointer[vIDb]), *pH);
	    label[vIDb] = (PGbyte)NBAND;
	  }
	}else t2 = INFINITY;
	return MIN(t1,t2);
      }	
      else return (MIN(Ta+b/Fc,Tb+a/Fc)); 
    } 
    
    /* Update angles */
    cos12A = (P1P2 + P2A - P1A)/(2*sqrt(P1P2*P2A));
    if(P2B != 0.0)
      cos12B = (P1P2 + P2B - P1B)/(2*sqrt(P1P2*P2B));
    cos12C = (P1P2 + P2C - P1C)/(2*sqrt(P1P2*P2C));
    
    sin12A = 1 - cos12A*cos12A;
    sin12B = 1 - cos12B*cos12B;
    sin12C = 1 - cos12C*cos12C;
    
    iters++;
  }/* End of while loop */

  return (MIN(Ta+b/Fc, Tb+a/Fc));
}

double ComputeTInAcute(double Ta, double Tb, double a, double b, double c, double Fc)
{
  double t1, t2, t, CD, costheta;
  double aa,bb,cc,u,tmp;
  
  costheta = (a*a+b*b-c*c)/(2*a*b); 
  
  u = Tb - Ta; 
  
  Fc = 1/Fc; /* Inverted here ! */

  aa = a*a + b*b -2*a*b*costheta; 
  bb = 2*b*u*(a*costheta-b);
  cc = b*b*(u*u-a*a*Fc*Fc*(1-costheta*costheta));
  
  tmp = bb*bb -4*aa*cc;
  if(tmp < 0) return (MIN(b*Fc+Ta,a*Fc+Tb));

  tmp = sqrt(tmp);
 
  t1 = (-bb + tmp)/(2*aa);  
  t2 = (-bb - tmp)/(2*aa);  
  t = MAX(t1,t2); 
  CD = (b*(t-u))/t;
  
  if ( (u<t) && (a*costheta<CD) && (CD<(a/costheta)) )
    return (t+Ta);  
  else
    return (MIN(b*Fc+Ta,a*Fc+Tb));
  
}
