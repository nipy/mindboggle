#include "pathfinder.h"

void WriteSurfaceAttributeMarchingSpeed(char *fname, Surface *surface, double *marchingSpeed)
{
	int i, j;

	FILE *fp;
	fp = fopen(fname,"wt");
		
	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"vtk output\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET POLYDATA\n");
	fprintf(fp,"POINTS %d float\n",surface->vertexNum);

	// vertex
	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%f %f %f\n",surface->vertex[i].x, surface->vertex[i].y, surface->vertex[i].z);
	} 

	// face
	fprintf(fp,"POLYGONS %d %d\n",surface->faceNum,surface->faceNum*4);

	for(i=0;i<surface->faceNum;i++)
	{
		fprintf(fp,"3 %d %d %d\n",surface->faces[i].v[0],surface->faces[i].v[1],surface->faces[i].v[2]);
	}

	// principal curvature 1
	fprintf(fp,"POINT_DATA %d\n",surface->vertexNum);
	fprintf(fp,"SCALARS pcurv1 float\n");
	fprintf(fp,"LOOKUP_TABLE pcurv1Table\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%f\n",surface->curv1[i]);
	}
		
	// principal curvature 2
	fprintf(fp,"SCALARS pcurv2 float\n");
	fprintf(fp,"LOOKUP_TABLE pcurv2Table\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%f\n",surface->curv2[i]);
	}

	// marching speed
	fprintf(fp,"SCALARS marchingSpeed float\n");
	fprintf(fp,"LOOKUP_TABLE marchingSpeedTable\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%f\n",marchingSpeed[i]);
	}

	// curvature derivative 0 
	fprintf(fp,"SCALARS dcurv0 float\n");
	fprintf(fp,"LOOKUP_TABLE dcurv0Table\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%f\n",surface->dcurv[0][i]);
	}

	// curvature derivative 1 
	fprintf(fp,"SCALARS dcurv1 float\n");
	fprintf(fp,"LOOKUP_TABLE dcurv1Table\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%f\n",surface->dcurv[1][i]);
	}

	// curvature derivative 2 
	fprintf(fp,"SCALARS dcurv2 float\n");
	fprintf(fp,"LOOKUP_TABLE dcurv2Table\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%f\n",surface->dcurv[2][i]);
	}

	// curvature derivative 3 
	fprintf(fp,"SCALARS dcurv3 float\n");
	fprintf(fp,"LOOKUP_TABLE dcurv3Table\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%f\n",surface->dcurv[3][i]);
	}

	// principal direction 1
	fprintf(fp,"VECTORS pdir1 float\n");

	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%f %f %f\n",surface->pdir1[i].x, surface->pdir1[i].y, surface->pdir1[i].z);
	} 

	// principal direction 2
	fprintf(fp,"VECTORS pdir2 float\n");

	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%f %f %f\n",surface->pdir2[i].x, surface->pdir2[i].y, surface->pdir2[i].z);
	} 

	// normal
	fprintf(fp,"NORMALS Normals float\n");

	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%f %f %f\n",surface->normal[i].x, surface->normal[i].y, surface->normal[i].z);
	} 

	fclose(fp);

	return;
}



void PathFinder(Surface *pmesh, double* marchingSpeed, int* VtxID, int N, VtxOnFace* PathVtxList, int* pathLen)
{
  int i,j,k,l,jm;
  int Vnum;
  int nBlkEdge = 0;
  int BlkEdge[32];
  int nBlkVtx = 0;
  int BlkVtx[32];
  VtxOnFace CurrentPosition,NextPosition,dummyPosition;
  int DONE;
  int OnEdge,Obtuse;
  int Found,specialVtx;
  int BackUpVtx;
  double* GeoDist;
  double* UnitSpeed;
  double distdummy;
  float* dummyF;
  double maxGrad,dummyGrad,dummyDouble;
  int nbrVtxNum,nbrPolyNum,nbrEdgeNum,dummyEdge;
  int nbrVtxID[NBR_MAX_NUM],nbrPolyID[NBR_MAX_NUM],nbrEdgeID[NBR_MAX_NUM];
  int PathLenSave;
  double elapse_time;
  Fvector3d V0, V1, V2;
  double lambda, b0, b1, d0;
  int CurrentTri = -1;

//  Vnum = XUGetVertexTableVN(pmesh);
  Vnum = pmesh->vertexNum;
  GeoDist = (double*)malloc(Vnum*sizeof(double));
  UnitSpeed = (double*)malloc(Vnum*sizeof(double));
  Found = 0;

  for (k = 0; k < Vnum; k ++)
    UnitSpeed[k] = 1.0;

  *pathLen = 1;
  PathVtxList[0].ID = *(VtxID+N-1);
  PathVtxList[0].lambda = 0.0;
  
  for (j=N-2;j>=0;j--)
    {
      PathLenSave = *pathLen;

      GeoDist = FastMarchMesh(pmesh, VtxID+j, 1,marchingSpeed,*(VtxID+j+1));

	  WriteSurfaceAttributeMarchingSpeed("GeoDist.vtk", pmesh, GeoDist);

      CurrentPosition.ID = *(VtxID+j+1);
      CurrentPosition.lambda = 0.0;
      OnEdge = 0;

	  int count = 0;
      while (CurrentPosition.ID != *(VtxID+j) && *pathLen < MAXPATHLEN)
		{

			count++;
		  if(!OnEdge)
			GetNbrV2(CurrentPosition.ID,pmesh,nbrVtxID,nbrPolyID,nbrEdgeID,&nbrVtxNum,&nbrPolyNum,&nbrEdgeNum,1);
		  else
			{
			  nbrEdgeNum = 4;
			  nbrVtxNum = 5;
			  nbrVtxID[0] = 0;
//			  nbrPolyID[0] = pmesh->edge_table.edge[CurrentPosition.ID].P[0];
//			  nbrPolyID[1] = pmesh->edge_table.edge[CurrentPosition.ID].P[1];
			  nbrPolyID[0] = pmesh->edges[CurrentPosition.ID].faceID[0];
			  nbrPolyID[1] = pmesh->edges[CurrentPosition.ID].faceID[1];
			  i = 0;jm = 1;
			  for (l = 0;l<3;l++){
//				dummyEdge = pmesh->poly_table.poly[nbrPolyID[0]].elist[l];
				dummyEdge = pmesh->faceEdges[nbrPolyID[0]][l];
				  if(dummyEdge!=CurrentPosition.ID)
				  {
//					nbrEdgeID[i] = pmesh->poly_table.poly[nbrPolyID[0]].elist[l];
					nbrEdgeID[i] = pmesh->faceEdges[nbrPolyID[0]][l];
					i++;
				  }
//				if(!IsInSet(nbrVtxID,jm,pmesh->edge_table.edge[dummyEdge].V[0]))
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[0]))
				  {
//					nbrVtxID[jm] = pmesh->edge_table.edge[dummyEdge].V[0];
					nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[0];
					jm++;
				  }
//				if(!IsInSet(nbrVtxID,jm,pmesh->edge_table.edge[dummyEdge].V[1]))
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[1]))
				  {
//					nbrVtxID[jm] = pmesh->edge_table.edge[dummyEdge].V[1];
					  nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[1];
					jm++;
				  }
			  }
			  for (l = 0;l<3;l++){
//				dummyEdge = pmesh->poly_table.poly[nbrPolyID[1]].elist[l];
				dummyEdge = pmesh->faceEdges[nbrPolyID[1]][l];
				if(dummyEdge!=CurrentPosition.ID)
				  {
//					nbrEdgeID[i] = pmesh->poly_table.poly[nbrPolyID[1]].elist[l];
					nbrEdgeID[i] = pmesh->faceEdges[nbrPolyID[1]][l];
					i++;
				  }
//				if(!IsInSet(nbrVtxID,jm,pmesh->edge_table.edge[dummyEdge].V[0]))
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[0]))				  {
//					nbrVtxID[jm] = pmesh->edge_table.edge[dummyEdge].V[0];
					  nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[0];
					jm++;
				  }
//				if(!IsInSet(nbrVtxID,jm,pmesh->edge_table.edge[dummyEdge].V[1]))
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[1]))
				  {
					  nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[1];
					jm++;
				  }
			  }
			}
		  // Found the neighboring vertices and edges of current position. If current position
		//	 is a vertex, the neighboring information is obtained using GetNbrV2. The information
		//	 includes neighboring vertices, edges, and faces. If current position is on an edge
		//	 the neighboring information is obtained by first checking the TWO neighboring faces 
		//	 of the edge, then put the edges of the triangles (faces) other than THIS one into 
		//	 the neighboring edge list, all vertices forming the two triangles are in the 
		//	 neighboring vertex list. //

		  // From the current position and its neighboring information to find the next position
		//	 of the curve. Steps:
		//	 1. Find the point in the neighboring vertices whose distance value is the smallest.
		//	    This is a backup incase minimal gradient can not be found or is illegal;
		//	 2. For each neighboring edge, compute the smallest gradient;
		//	 3. Decide whether to take the result in 1 or the result in 2.
		  //
		  
		  distdummy = InitialGuess (pmesh, CurrentPosition, GeoDist, nbrVtxID, nbrVtxNum, &dummyPosition, nBlkVtx, BlkVtx);
		  for (k = 0; k < nbrEdgeNum; k ++){
			if (IsInSet(BlkEdge,nBlkEdge,nbrEdgeID[k]))
			  continue;
			dummyGrad = Gradient(pmesh, nbrEdgeID[k], CurrentPosition, GeoDist, &NextPosition);
			//			if (NextPosition.lambda == 0 && IsInSet(BlkVtx,nBlkVtx,NextPosition.ID))
			//			continue; //
			if ( dummyGrad > distdummy){
			  dummyPosition.ID = NextPosition.ID;
			  dummyPosition.lambda = NextPosition.lambda;
			  distdummy = dummyGrad;
			}
		  }

		  // check to avoid going from an edge to its own vtx, and touring all three vertices of a triangle //
		  if ((CurrentPosition.lambda != 0) && (dummyPosition.lambda == 0)) {
			// going from an edge to a vtx, possible of going from an edge to its own vtx //
			if (IsOnEdge(pmesh, dummyPosition.ID, CurrentPosition.ID))
			  DONE = 0;  // bad case //
			else 
			  DONE = 1;  // good case //

			while (!DONE){
			  CurrentPosition.lambda = PathVtxList[(*pathLen)-1].lambda;
			  CurrentPosition.ID = PathVtxList[(*pathLen)-1].ID;
			  if ((CurrentPosition.lambda == 0) || !IsOnEdge(pmesh, dummyPosition.ID, CurrentPosition.ID))
				DONE = 1;
			  else
				(*pathLen) --;
			}
		  }  
		  else if ((CurrentPosition.lambda == 0) && (dummyPosition.lambda == 0)) {
			// going from a vtx to a vtx, possible of touring all three vertices of a triangle //
			GetNbrV2(dummyPosition.ID,pmesh,nbrVtxID,nbrPolyID,nbrEdgeID,&nbrVtxNum,&nbrPolyNum,&nbrEdgeNum,1);
			if ((PathVtxList[(*pathLen)-2].lambda== 0) && IsInSet(nbrVtxID,nbrVtxNum,PathVtxList[(*pathLen)-2].ID))
			  (*pathLen)--;
		  }		 
		  // determine the block edges //
		  nBlkEdge = DetermineblkEdge(pmesh, PathVtxList[(*pathLen)-1], dummyPosition, BlkEdge);
		  nBlkVtx = DetermineblkVtx(pmesh, PathVtxList[(*pathLen)-1], dummyPosition, BlkVtx);

		  PathVtxList[*pathLen].ID = dummyPosition.ID;
		  PathVtxList[*pathLen].lambda = dummyPosition.lambda;
		  (*pathLen) ++;
		  
		  CurrentPosition.ID = dummyPosition.ID;
		  CurrentPosition.lambda =  dummyPosition.lambda;
		  if (CurrentPosition.lambda == 0.0)
			OnEdge = 0;
		  else
			OnEdge = 1;

		  printf("count = %d: ID = %d, lambda = %f\n",count,CurrentPosition.ID, CurrentPosition.lambda);
		}

	}
  free(GeoDist);
  return;
}


bool CheckToStart(Surface *surface, VtxOnFace start, VtxOnFace cp)
{
	int i;

	if(start.lambda == 0.0)
	{
		if(cp.lambda == 0.0)
		{
			if(cp.ID == start.ID)
				return true;
			else
				return false;
		}
		else
		{
			for(i=0; i<surface->adjacentedges[start.ID].size(); i++)
			{
				if(cp.ID == surface->adjacentedges[start.ID][i])
				{
					return true;
				}
			}
			return false;
		}
	}
	else
	{
		if(cp.lambda == 0.0)
		{
			if(IsVertexInFace(surface, surface->edges[start.ID].faceID[0], cp.ID) || IsVertexInFace(surface, surface->edges[start.ID].faceID[1], cp.ID))
				return true;
		}
		else
		{
			if(IsVertexInFace(surface, surface->edges[start.ID].faceID[0], surface->edges[cp.ID].vertexID[0]) && IsVertexInFace(surface, surface->edges[start.ID].faceID[0], surface->edges[cp.ID].vertexID[1]) )
				return true;
			if(IsVertexInFace(surface, surface->edges[start.ID].faceID[1], surface->edges[cp.ID].vertexID[0]) && IsVertexInFace(surface, surface->edges[start.ID].faceID[1], surface->edges[cp.ID].vertexID[1]) )
				return true;
		}
	}

	return false;
}

bool PathFinderHalfVertex(Surface *pmesh, double* marchingSpeed, VtxOnFace start, VtxOnFace end, VtxOnFace* PathVtxList, int* pathLen)
{
  int i,j,k,l,jm;
  int Vnum;
  int nBlkEdge = 0;
  int BlkEdge[32];
  int nBlkVtx = 0;
  int BlkVtx[32];
  VtxOnFace CurrentPosition,NextPosition,dummyPosition;
  int DONE;
  int OnEdge,Obtuse;
  int Found,specialVtx;
  int BackUpVtx;
  double* GeoDist;
  double distdummy;
  float* dummyF;
  double maxGrad,dummyGrad,dummyDouble;
  int nbrVtxNum,nbrPolyNum,nbrEdgeNum,dummyEdge;
  int nbrVtxID[NBR_MAX_NUM],nbrPolyID[NBR_MAX_NUM],nbrEdgeID[NBR_MAX_NUM];
  double elapse_time;
  Fvector3d V0, V1, V2;
  double lambda, b0, b1, d0;
  int CurrentTri = -1;

  Vnum = pmesh->vertexNum;
  GeoDist = (double*)malloc(Vnum*sizeof(double));
  Found = 0;

  *pathLen = 1;
  PathVtxList[0].ID = end.ID;
  PathVtxList[0].lambda = end.lambda;
  

  GeoDist = FastMarchMeshHalfVertex(pmesh, start, end, marchingSpeed);

// WriteSurfaceAttributeMarchingSpeed("GeoDist.vtk", pmesh, GeoDist);

	  CurrentPosition.ID = end.ID;
	  CurrentPosition.lambda = end.lambda;
	  if(CurrentPosition.lambda == 0.0)  OnEdge = 0;
	  else OnEdge = 1;
		
	  int count = 0;
		bool endlessLoop = false;

	  while (!CheckToStart(pmesh, start,  CurrentPosition) && *pathLen < MAXPATHLEN)
		{
			if(count == 1000)
			{
				printf("Error: endless loop!\n");
				endlessLoop = true;
				break;
			}
			
			count++;
		  if(!OnEdge)
			GetNbrV2(CurrentPosition.ID,pmesh,nbrVtxID,nbrPolyID,nbrEdgeID,&nbrVtxNum,&nbrPolyNum,&nbrEdgeNum,1);
		  else
			{
			  nbrEdgeNum = 4;
			  nbrVtxNum = 5;
			  nbrVtxID[0] = 0;

			  nbrPolyID[0] = pmesh->edges[CurrentPosition.ID].faceID[0];
			  nbrPolyID[1] = pmesh->edges[CurrentPosition.ID].faceID[1];
			  i = 0;jm = 1;
			  for (l = 0;l<3;l++){
				dummyEdge = pmesh->faceEdges[nbrPolyID[0]][l];
				  if(dummyEdge!=CurrentPosition.ID)
				  {
					nbrEdgeID[i] = pmesh->faceEdges[nbrPolyID[0]][l];
					i++;
				  }
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[0]))
				  {
					nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[0];
					jm++;
				  }
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[1]))
				  {
					  nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[1];
					jm++;
				  }
			  }
			  for (l = 0;l<3;l++){
				dummyEdge = pmesh->faceEdges[nbrPolyID[1]][l];
				if(dummyEdge!=CurrentPosition.ID)
				  {
					nbrEdgeID[i] = pmesh->faceEdges[nbrPolyID[1]][l];
					i++;
				  }
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[0]))				 
				  {
					  nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[0];
					jm++;
				  }
				  if(!IsInSet(nbrVtxID,jm,pmesh->edges[dummyEdge].vertexID[1]))
				  {
					  nbrVtxID[jm] = pmesh->edges[dummyEdge].vertexID[1];
					jm++;
				  }
			  }
			}
		  // Found the neighboring vertices and edges of current position. If current position
		//	 is a vertex, the neighboring information is obtained using GetNbrV2. The information
		//	 includes neighboring vertices, edges, and faces. If current position is on an edge
		//	 the neighboring information is obtained by first checking the TWO neighboring faces 
		//	 of the edge, then put the edges of the triangles (faces) other than THIS one into 
		//	 the neighboring edge list, all vertices forming the two triangles are in the 
		//	 neighboring vertex list. //

		  // From the current position and its neighboring information to find the next position
		//	 of the curve. Steps:
		//	 1. Find the point in the neighboring vertices whose distance value is the smallest.
		//	    This is a backup incase minimal gradient can not be found or is illegal;
		//	 2. For each neighboring edge, compute the smallest gradient;
		//	 3. Decide whether to take the result in 1 or the result in 2.
		  //
		  
		  distdummy = InitialGuess (pmesh, CurrentPosition, GeoDist, nbrVtxID, nbrVtxNum, &dummyPosition, nBlkVtx, BlkVtx);
		  for (k = 0; k < nbrEdgeNum; k ++){
			if (IsInSet(BlkEdge,nBlkEdge,nbrEdgeID[k]))
			  continue;
			dummyGrad = Gradient(pmesh, nbrEdgeID[k], CurrentPosition, GeoDist, &NextPosition);

			if ( dummyGrad > distdummy){
			  dummyPosition.ID = NextPosition.ID;
			  dummyPosition.lambda = NextPosition.lambda;
			  distdummy = dummyGrad;
			}
		  }

		  // check to avoid going from an edge to its own vtx, and touring all three vertices of a triangle //
		  if ((CurrentPosition.lambda != 0) && (dummyPosition.lambda == 0)) {
			// going from an edge to a vtx, possible of going from an edge to its own vtx //
			if (IsOnEdge(pmesh, dummyPosition.ID, CurrentPosition.ID))
			  DONE = 0;  // bad case //
			else 
			  DONE = 1;  // good case //

			while (!DONE){
			  CurrentPosition.lambda = PathVtxList[(*pathLen)-1].lambda;
			  CurrentPosition.ID = PathVtxList[(*pathLen)-1].ID;
			  if ((CurrentPosition.lambda == 0) || !IsOnEdge(pmesh, dummyPosition.ID, CurrentPosition.ID))
				DONE = 1;
			  else
				(*pathLen) --;
			}
		  }  
		  else if ((CurrentPosition.lambda == 0) && (dummyPosition.lambda == 0)) {
			// going from a vtx to a vtx, possible of touring all three vertices of a triangle //
			GetNbrV2(dummyPosition.ID,pmesh,nbrVtxID,nbrPolyID,nbrEdgeID,&nbrVtxNum,&nbrPolyNum,&nbrEdgeNum,1);
			if ((PathVtxList[(*pathLen)-2].lambda== 0) && IsInSet(nbrVtxID,nbrVtxNum,PathVtxList[(*pathLen)-2].ID))
			  (*pathLen)--;
		  }		 
		  // determine the block edges //
		  nBlkEdge = DetermineblkEdge(pmesh, PathVtxList[(*pathLen)-1], dummyPosition, BlkEdge);
		  nBlkVtx = DetermineblkVtx(pmesh, PathVtxList[(*pathLen)-1], dummyPosition, BlkVtx);

		  PathVtxList[*pathLen].ID = dummyPosition.ID;
		  PathVtxList[*pathLen].lambda = dummyPosition.lambda;
		  (*pathLen) ++;
		  
		  CurrentPosition.ID = dummyPosition.ID;
		  CurrentPosition.lambda =  dummyPosition.lambda;
		  if (CurrentPosition.lambda == 0.0)
			OnEdge = 0;
		  else
			OnEdge = 1;

		  // check endless loop
			if((*pathLen) > 3)
			{
				if(PathVtxList[(*pathLen)-1].ID == PathVtxList[(*pathLen)-3].ID)
				{
					printf("Error: endless loop!\n");
					endlessLoop = true;
					break;
				}
			}

		  printf("count = %d: ID = %d, lambda = %f\n",count,CurrentPosition.ID, CurrentPosition.lambda);
		}

	free(GeoDist);
	
	return endlessLoop;
}

double CalculateLambdaGradient(double b0, double b1, double b2,double A, double B, double C,double *Grad)
{
  double lambda;
  double g1,g2;
  double b12,b02,b01,alpha;
  double x0,y0;

  b12 = b1-b2;
  b02 = b0-b2;
  b01 = b0-b1;
  g1 = b01/B;
  g2 = b02/A;
  x0 = (B*B+C*C-A*A)/(2*C);
  y0 = sqrt(B*B-x0*x0);

  lambda = (B*B*b12 + C*x0*b01)/(C*x0*b12 + C*C*b01);
  *Grad = (b01 + lambda * b12)/sqrt((lambda * C - x0)*(lambda * C - x0)+y0*y0);
  
  if ((g1 > *Grad) && (g1 >= g2)){
	*Grad = g1;
	return(0.0);
  }
  if ((g2 > *Grad) && (g2 >= g1)){
	*Grad = g2;
	return(1.0);
  }

  if (lambda < 0.0 || lambda > 1.0){
	if (g1 > g2){
	  *Grad = g1;
	  lambda = 0.0;
	}
	else{
	  *Grad = g2;
	  lambda = 1.0;
	}
  }
  return(lambda);
}

double Gradient(Surface* pmesh,int EdgeID,VtxOnFace CurrentPosition,double *GeoDist,VtxOnFace *Position)
{
  Fvector3d V0,V1,V2,Vd;
  double b0,b1,b2,b,Gradient,d,lambda;
  double Gradient2,lambda2;
  int ID0,ID1,ID2;
  double A, B, C, D;
  
  VtxOnFace EndPosition;

  // variables for unfolding 
  int PolyID0,PolyID1;
  int dummyID;
  int Obtuse;
  int P1ID,P2ID,UID,UNotID;
  int i,j,iters,NN,neighVId;

  Fvector3d *P1, *P2, *U;
  Fvector3d vP1U, vP2U;

  double P1A, P1B, P1C, P2A, P2B, P2C, UA, UB, UC, P1P2, P1U, P2U; 
  double cos12A, sin12A, cos12B, sin12B, cos12C, sin12C, cos12U, sin12U;
  double AC, BC, AB;

  double cos_alpha1,cos_alpha2,cos_beta,sin_beta,sin_gamma;
  double CosGamma, SinGamma, CosAlpha, SinAlpha, CosDelta, SinDelta, CosCAU;
  double CP;

  Obtuse = IsObtuse(pmesh,CurrentPosition, EdgeID);  

  b0 = FindValueAndPosition(pmesh, CurrentPosition, GeoDist, &V0);
  
//  ID1 = pmesh->edge_table.edge[EdgeID].V[0];
//  ID2 = pmesh->edge_table.edge[EdgeID].V[1];
  ID1 = pmesh->edges[EdgeID].vertexID[0];
  ID2 = pmesh->edges[EdgeID].vertexID[1];

//  V1.x = pmesh->vertex_table.vertex[ID1].x;
 // V1.y = pmesh->vertex_table.vertex[ID1].y;
 // V1.z = pmesh->vertex_table.vertex[ID1].z;
  //V2.x = pmesh->vertex_table.vertex[ID2].x;
  //V2.y = pmesh->vertex_table.vertex[ID2].y;
  //V2.z = pmesh->vertex_table.vertex[ID2].z;

  V1.x = pmesh->vertex[ID1].x;
  V1.y = pmesh->vertex[ID1].y;
  V1.z = pmesh->vertex[ID1].z;
  V2.x = pmesh->vertex[ID2].x;
  V2.y = pmesh->vertex[ID2].y;
  V2.z = pmesh->vertex[ID2].z;

  b1 = GeoDist[ID1];
  b2 = GeoDist[ID2];
  
  A = MyDistance(V0,V2);
  B = MyDistance(V0,V1);
  C = MyDistance(V1,V2);

  if (Obtuse == 0){
	lambda = CalculateLambdaGradient(b0,b1,b2,A,B,C,&Gradient);
	if(lambda == 0.0)
	{
//	  Position->ID = pmesh->edge_table.edge[EdgeID].V[0];
		Position->ID = pmesh->edges[EdgeID].vertexID[0];
	  Position->lambda = 0.0;
	}
	else if(lambda == 1.0)
	  {
//		Position->ID = pmesh->edge_table.edge[EdgeID].V[1];
		  Position->ID = pmesh->edges[EdgeID].vertexID[1];
		Position->lambda = 0.0;
	  }
	else 
	  {
		Position->ID = EdgeID;
		Position->lambda = lambda;
	  }
	
	return(Gradient);
  }

  else { // If the triangle is an obtuse triangle, perform unfolding here 
	// The edge AB is always in the right order, meaning A is always the start point of the edge 
	if (CurrentPosition.lambda == 0)
	  UNotID = CurrentPosition.ID;
//	else if ((pmesh->edge_table.edge[CurrentPosition.ID].V[0] != ID1) && (pmesh->edge_table.edge[CurrentPosition.ID].V[0] != ID2)){
	else if ((pmesh->edges[CurrentPosition.ID].vertexID[0] != ID1) && (pmesh->edges[CurrentPosition.ID].vertexID[0] != ID2)){
//	  UNotID = pmesh->edge_table.edge[CurrentPosition.ID].V[0];
//	  dummyID = pmesh->edge_table.edge[CurrentPosition.ID].V[1];
		UNotID = pmesh->edges[CurrentPosition.ID].vertexID[0];
		dummyID = pmesh->edges[CurrentPosition.ID].vertexID[1];
	}
	else{ 
//	  UNotID = pmesh->edge_table.edge[CurrentPosition.ID].V[1];
//	  dummyID = pmesh->edge_table.edge[CurrentPosition.ID].V[0];
	  UNotID = pmesh->edges[CurrentPosition.ID].vertexID[1];
	  dummyID = pmesh->edges[CurrentPosition.ID].vertexID[0];
	}

	P1ID = ID1;
	P2ID = ID2;
	
//	Vd.x = pmesh->vertex_table.vertex[UNotID].x;
//	Vd.y = pmesh->vertex_table.vertex[UNotID].y;
//	Vd.z = pmesh->vertex_table.vertex[UNotID].z;
	Vd.x = pmesh->vertex[UNotID].x;
	Vd.y = pmesh->vertex[UNotID].y;
	Vd.z = pmesh->vertex[UNotID].z;
	D = MyDistance(Vd,V1);
	
	AC = MyDistanceSquare(Vd,V1);
	BC = MyDistanceSquare(Vd,V2);
	P1A = 0; P1B = C*C; P1C = D*D;
	P2B = 0; P2A = P1B; P2C = MyDistanceSquare(Vd,V2);
	P1P2 = P1B;
    
	cos12A = 1; sin12A = 0; cos12B = 1; sin12B = 0; 
	cos12C = (P1P2 + P2C - P1C)/(2*sqrt(P1P2*P2C));
	sin12C = (1-cos12C*cos12C); // Notice: Square of sine 
    
	// Now iteratively unfolding 
	iters = 0;
	while(iters <10){
	  
	  // Find the newly unfolded vertex ID 
	  
//	  NN = XUGetNeighborEN(pmesh, P1ID);
	  NN = pmesh->neighbors[P1ID].size();
	  for (i=0; i<NN; i++) {
//		neighVId = XUGetVertexNeighborVId(pmesh, P1ID, i);
		neighVId = pmesh->neighbors[P1ID][i];
		if(neighVId == P2ID)
		{
			/*
		  //if(XUGetVertexNeighborVId(pmesh, P1ID, (i+NN-1)%NN) == UNotID){
			if(pmesh->neighbors[P1ID][(i+NN-1)%NN] == UNotID)
			{
//			UID = XUGetVertexNeighborVId(pmesh, P1ID, (i+1)%NN);
			UID = pmesh->neighbors[P1ID][(i+1)%NN];
			b = GeoDist[UID];
			break;
			}
//		  if(XUGetVertexNeighborVId(pmesh, P1ID, (i+1)%NN) == UNotID){
		  if(pmesh->neighbors[P1ID][(i+1)%NN] == UNotID)
		  {
//			UID = XUGetVertexNeighborVId(pmesh, P1ID, (i+NN-1)%NN);
			UID = pmesh->neighbors[P1ID][(i+NN-1)%NN];
			b = GeoDist[UID];
			break;
		  }
		  */
			
			int leftVertexID, rightVertexID;

			leftVertexID = rightVertexID = -1;

			GetLeftRightVertex(pmesh, P1ID, P2ID, &leftVertexID, &rightVertexID);

			if(leftVertexID == UNotID)
			{
				UID = rightVertexID;
				if(UID < 0)
				{
					printf("UID < 0 in Gradient\n");
					UID = leftVertexID;
				}
		
				b = GeoDist[UID];
				break;
			}
			if(rightVertexID == UNotID)
			{	
				UID = leftVertexID;
				if(UID < 0)
				{
					printf("UID < 0 in Gradient\n");
					UID = rightVertexID;
				}

				b = GeoDist[UID];
				break;
			}
		}
	  }
//
//	  P1 = XUGetVertexPtr(pmesh,P1ID); 
//	  P2 = XUGetVertexPtr(pmesh,P2ID); 
//	  U = XUGetVertexPtr(pmesh,UID);
//

	  P1 = &(pmesh->vertex[P1ID]); 
	  P2 = &(pmesh->vertex[P2ID]); 
	  U = &(pmesh->vertex[UID]);
    
//	  vP1U = vectsubptr(P1, U);
//	  vP2U = vectsubptr(P2, U);
	  vP1U.x = P1->x - U->x;
	  vP1U.y = P1->y - U->y;
	  vP1U.z = P1->z - U->z;
	  vP2U.x = P2->x - U->x;
	  vP2U.y = P2->y - U->y;
	  vP2U.z = P2->z - U->z;

 
	  P1U = (vP1U.x * vP1U.x + vP1U.y * vP1U.y + vP1U.z * vP1U.z);
	  P2U = (vP2U.x * vP2U.x + vP2U.y * vP2U.y + vP2U.z * vP2U.z);
	  
	  cos12U = (P1P2 + P2U - P1U)/(2*sqrt(P1P2*P2U));
	  sin12U = (1-cos12U*cos12U); // Notice: Square of sine //
      
	  // Now compute three lengthes (squared) //
	  UA = P2U + P2A - 2*sqrt(P2U*P2A)*(cos12A*cos12U - sqrt(sin12A*sin12U));
	  UB = P2U + P2B - 2*sqrt(P2U*P2B)*(cos12B*cos12U - sqrt(sin12B*sin12U));
	  UC = P2U + P2C - 2*sqrt(P2U*P2C)*(cos12C*cos12U - sqrt(sin12C*sin12U));
	  
	 // Now Judge Which Side to continue unfolding //
	  if(UA > (UC + AC)){// Unfold along P1U //
		UNotID = P2ID;
		P2ID = UID;
		P1P2 = P1U;
		P2A = UA; P2B = UB; P2C = UC;
	  }else if(UB > (UC + BC)){ // Unfold along P2U //
		UNotID = P1ID;
		P1ID = UID;
		P1P2 = P2U;
		P1A = UA; P1B = UB; P1C = UC;
	  }else{ // Stop Unfolding and compute T//
		// Compute the actual lengthes ///
		UC = sqrt(UC);
		UA = sqrt(UA);
		UB = sqrt(UB);
		break;
	  } 
	  
	  // Update angles //
	  cos12A = (P1P2 + P2A - P1A)/(2*sqrt(P1P2*P2A));
	  if(P2B != 0.0)
		cos12B = (P1P2 + P2B - P1B)/(2*sqrt(P1P2*P2B));
	  cos12C = (P1P2 + P2C - P1C)/(2*sqrt(P1P2*P2C));
	  
	  sin12A = 1 - cos12A*cos12A;
	  sin12B = 1 - cos12B*cos12B;
	  sin12C = 1 - cos12C*cos12C;
      
	  iters++;
	}// End of while loop //

	// After unfolding, there are three cases where we need to consider        //
	// For each case, there is a different procedure to compute the next point //
	if (CurrentPosition.lambda == 0){ // case one: the easiest one //
	  	lambda = CalculateLambdaGradient(b0,b1,b,UC,B,UA,&Gradient);
		lambda2 = CalculateLambdaGradient(b0,b,b2,A,UC,UB,&Gradient2);
		if (Gradient >= Gradient2){
		  CosGamma = (UA*UA+C*C-UB*UB)/(2*UA*C);
		  SinGamma = sqrt(1-CosGamma*CosGamma);
		  CosAlpha = (B*B+C*C-A*A)/(2*B*C);
		  SinAlpha = sqrt(1-CosAlpha*CosAlpha);
		  CosCAU = CosGamma*CosAlpha -SinGamma*SinAlpha;
		  CP = sqrt(B*B + lambda*UA*lambda*UA - 2*B*lambda*UA*CosCAU);
		  CosDelta = (B*B+CP*CP-lambda*UA*lambda*UA)/(2*B*CP);
		  SinDelta = sqrt(1-CosDelta*CosDelta);
		  lambda = B*SinDelta/((SinAlpha*CosDelta + CosAlpha*SinDelta)*C);
		  if (lambda == 0){
//			Position->ID = pmesh->edge_table.edge[EdgeID].V[0];
			  Position->ID = pmesh->edges[EdgeID].vertexID[0];
			Position->lambda = 0;
		  }
		  else {
			Position->ID = EdgeID;
			Position->lambda = lambda;
		  }
		  return(Gradient);
		}

		else {
		  lambda = 1-lambda2;
		  CosGamma = (UB*UB+C*C-UA*UA)/(2*UB*C);
		  SinGamma = sqrt(1-CosGamma*CosGamma);
		  CosAlpha = (A*A+C*C-B*B)/(2*A*C);
		  SinAlpha = sqrt(1-CosAlpha*CosAlpha);
		  CosCAU = CosGamma*CosAlpha -SinGamma*SinAlpha;
		  CP = sqrt(A*A + lambda*UB*lambda*UB - 2*A*lambda*UB*CosCAU);
		  CosDelta = (A*A+CP*CP-lambda*UB*lambda*UB)/(2*A*CP);
		  SinDelta = sqrt(1-CosDelta*CosDelta);
		  lambda = 1 - A*SinDelta/((SinAlpha*CosDelta + CosAlpha*SinDelta)*C);
		  if (lambda == 1){
//			Position->ID = pmesh->edge_table.edge[EdgeID].V[1];
			  Position->ID = pmesh->edges[EdgeID].vertexID[1];
			Position->lambda = 0;
		  }
		  else {
			Position->ID = EdgeID;
			Position->lambda = lambda;
		  }
		  return(Gradient);
		}
	} // The easist case //

	else { // If current point is on an edge //
	  // Need to compute UC2 //
	  if (dummyID == ID1) { // C' is on AC //
		CosGamma = (UB*UB+C*C-UA*UA)/(2*UB*C);
		SinGamma = sqrt(1-CosGamma*CosGamma);
		CosAlpha = (A*A+C*C-B*B)/(2*A*C);
		SinAlpha = sqrt(1-CosAlpha*CosAlpha);
		CosCAU = CosGamma*CosAlpha -SinGamma*SinAlpha;
		UC = sqrt(A*A+UB*UB-2*A*UB*CosAlpha);
	  }
	  else if (dummyID == ID2) { // C' is on BC //
		CosGamma = (UA*UA+C*C-UB*UB)/(2*UA*C);
		SinGamma = sqrt(1-CosGamma*CosGamma);
		CosAlpha = (B*B+C*C-A*A)/(2*B*C);
		SinAlpha = sqrt(1-CosAlpha*CosAlpha);
		CosCAU = CosGamma*CosAlpha -SinGamma*SinAlpha;
		UC = sqrt(B*B+UA*UA-2*B*UA*CosAlpha);
	  }
	  lambda = CalculateLambdaGradient(b0,b1,b,UC,B,UA,&Gradient);
	  lambda2 = CalculateLambdaGradient(b0,b,b2,A,UC,UB,&Gradient2);
	  if (Gradient >= Gradient2){
		CosGamma = (UA*UA+C*C-UB*UB)/(2*UA*C);
		SinGamma = sqrt(1-CosGamma*CosGamma);
		CosAlpha = (B*B+C*C-A*A)/(2*B*C);
		SinAlpha = sqrt(1-CosAlpha*CosAlpha);
		CosCAU = CosGamma*CosAlpha -SinGamma*SinAlpha;
		CP = sqrt(B*B + lambda*UA*lambda*UA - 2*B*lambda*UA*CosCAU);
		CosDelta = (B*B+CP*CP-lambda*UA*lambda*UA)/(2*B*CP);
		SinDelta = sqrt(1-CosDelta*CosDelta);
		lambda = B*SinDelta/((SinAlpha*CosDelta + CosAlpha*SinDelta)*C);
		if (lambda == 0){
//		  Position->ID = pmesh->edge_table.edge[EdgeID].V[0];
		  Position->ID = pmesh->edges[EdgeID].vertexID[0];
		  Position->lambda = 0;
		}
		else {
		  Position->ID = EdgeID;
		  Position->lambda = lambda;
		}
		return(Gradient);
		}
	  
	  else {
		lambda = 1-lambda2;
		CosGamma = (UB*UB+C*C-UA*UA)/(2*UB*C);
		SinGamma = sqrt(1-CosGamma*CosGamma);
		CosAlpha = (A*A+C*C-B*B)/(2*A*C);
		SinAlpha = sqrt(1-CosAlpha*CosAlpha);
		CosCAU = CosGamma*CosAlpha -SinGamma*SinAlpha;
		CP = sqrt(A*A + lambda*UB*lambda*UB - 2*A*lambda*UB*CosCAU);
		CosDelta = (A*A+CP*CP-lambda*UB*lambda*UB)/(2*A*CP);
		SinDelta = sqrt(1-CosDelta*CosDelta);
		lambda = 1 - A*SinDelta/((SinAlpha*CosDelta + CosAlpha*SinDelta)*C);
		if (lambda == 1){
//			Position->ID = pmesh->edge_table.edge[EdgeID].V[1];
			Position->ID = pmesh->edges[EdgeID].vertexID[1];
			Position->lambda = 0;
		}
		else {
		  Position->ID = EdgeID;
		  Position->lambda = lambda;
		}
		  return(Gradient);
	  }
	} 
  }
return(0);
}
//
double InitialGuess (Surface *pmesh, VtxOnFace CurrentPosition, double* GeoDist, 
					 int* nbrVtxID, int nbrVtxNum, VtxOnFace *dummyPosition, int nBV, int* BV)
{

  double b0, d0;
  Fvector3d V0;
  double distdummy = 0.0;
  int k, l;

  b0 = FindValueAndPosition(pmesh, CurrentPosition, GeoDist, &V0);
  
  for (k = 1; k < nbrVtxNum; k ++){
	if (!IsInSet(BV, nBV, nbrVtxID[k])){
//	  d0 = MyDistance(V0, (Vector)pmesh->vertex_table.vertex[nbrVtxID[k]]);
	  d0 = MyDistance(V0, pmesh->vertex[nbrVtxID[k]]);
	  distdummy = (b0 - GeoDist[nbrVtxID[k]])/d0;
	  dummyPosition->ID = nbrVtxID[k];
	  break;
	}
  }
  l = k + 1;
  for (k = l; k < nbrVtxNum; k ++){
	if (IsInSet(BV, nBV, nbrVtxID[k]))
	  continue;
//	d0 = MyDistance(V0, (Vector)pmesh->vertex_table.vertex[nbrVtxID[k]]);
	d0 = MyDistance(V0, pmesh->vertex[nbrVtxID[k]]);
	d0 = (b0 - GeoDist[nbrVtxID[k]])/d0;
	if (d0 >  distdummy){
	  dummyPosition->ID = nbrVtxID[k];
	  distdummy = d0;
	}
  }
  dummyPosition->lambda = 0.0;
  return(distdummy);
}
//
int DetermineblkVtx(Surface *pmesh, VtxOnFace v0, VtxOnFace v1, int* BV)
{
  int nBV;
  int id0, id1;
  int j,k,l;
  int nbrVtxNum0,nbrVtxNum1,nbrPolyNum,nbrEdgeNum,dummyEdge;
  int nbrVtxID0[NBR_MAX_NUM],nbrVtxID1[NBR_MAX_NUM],nbrPolyID[NBR_MAX_NUM],nbrEdgeID[NBR_MAX_NUM];
  
  if (v0.lambda == 0 && v1.lambda == 0){ // vtx to vtx 
	nBV = 1;
	BV[0] = v0.ID;
  }
  else if (v0.lambda == 0 && v1.lambda !=0 ){ // vtx to edge
	nBV = 1;
	BV[0] = v0.ID;
  }
  else if (v0.lambda != 0 && v1.lambda == 0){ // edge to vtx 
	nBV = 2;
//	BV[0] = pmesh->edge_table.edge[v0.ID].V[0];
//	BV[1] = pmesh->edge_table.edge[v0.ID].V[1];
	BV[0] = pmesh->edges[v0.ID].vertexID[0];
	BV[1] = pmesh->edges[v0.ID].vertexID[1];
  }
  else { // edge to edge 
	nBV = 2;
	//BV[0] = pmesh->edge_table.edge[v0.ID].V[0];
	//BV[1] = pmesh->edge_table.edge[v0.ID].V[1];
	BV[0] = pmesh->edges[v0.ID].vertexID[0];
	BV[1] = pmesh->edges[v0.ID].vertexID[1];
  }
  return(nBV);
}


int DetermineblkEdge(Surface *pmesh, VtxOnFace v0, VtxOnFace v1, int* BE)
{
  int nBE;
  if (v0.lambda == 0 && v1.lambda == 0){ // vtx to vtx
	nBE = 2;
	FindTwoEdges(pmesh, v0, v1, BE);
  }
  else if (v0.lambda == 0 && v1.lambda !=0 ){ // vtx to edge 
	nBE = 2;
	FindTwoEdges(pmesh, v0, v1, BE);
  }
  else if (v0.lambda != 0 && v1.lambda == 0){ // edge to vtx 
	nBE = 1;
	BE[0] = v0.ID;
  }
  else { // edge to edge 
	nBE = 2;
	BE[0] = v0.ID;
	BE[1] = FindAnotherEdge(pmesh, v0.ID, v1.ID);
  }
  return(nBE);//
}

void FindTwoEdges(Surface *pmesh, VtxOnFace v0, VtxOnFace v1, int *BE)
{
  int p;
  int j,k,l;
  int nbrVtxNum,nbrPolyNum,nbrEdgeNum,dummyEdge;
  int nbrVtxID[NBR_MAX_NUM],nbrPolyID[NBR_MAX_NUM],nbrEdgeID[NBR_MAX_NUM];

  if (v1.lambda != 0){  //vtx to edge 
	p = FindTriangle(pmesh, v0, v1);
	//if (pmesh->poly_table.poly[p].elist[0] == v1.ID){
	if (pmesh->faceEdges[p][0] == v1.ID)
	{
//	  BE[0] =  pmesh->poly_table.poly[p].elist[1];
//	  BE[1] =  pmesh->poly_table.poly[p].elist[2];
	  BE[0] =  pmesh->faceEdges[p][1];
	  BE[1] =  pmesh->faceEdges[p][2];

	  return;
	}
//	if (pmesh->poly_table.poly[p].elist[1] == v1.ID){
	if (pmesh->faceEdges[p][1] == v1.ID){
//	  BE[0] =  pmesh->poly_table.poly[p].elist[0];
//	  BE[1] =  pmesh->poly_table.poly[p].elist[2];
  	  BE[0] =  pmesh->faceEdges[p][0];
	  BE[1] =  pmesh->faceEdges[p][2];

	  return;
	}
//	if (pmesh->poly_table.poly[p].elist[2] == v1.ID){
	if (pmesh->faceEdges[p][2] == v1.ID){
	  BE[0] =  pmesh->faceEdges[p][1];
	  BE[1] =  pmesh->faceEdges[p][0];
	  return;
	}
  }
  else { // vtx to vtx 
	GetNbrV2(v0.ID, pmesh,nbrVtxID,nbrPolyID,nbrEdgeID,&nbrVtxNum,&nbrPolyNum,&nbrEdgeNum,1);	
	l = 0;
	for (k = 0; k < nbrPolyNum; k++){
	  if (VertexOnPolygon(pmesh, v1.ID, nbrPolyID[k])){
		for (j = 0; j < 3; j++){
		  //p = pmesh->poly_table.poly[nbrPolyID[k]].elist[j];
		  p = pmesh->faceEdges[nbrPolyID[k]][j];
		  if (IsOnEdge(pmesh, v0.ID, p) && !IsOnEdge(pmesh, v1.ID, p)){
			BE[l] = p;
			l ++;
			if (l == 2)
			  return;
		  }
		}
	  }
	}
  }
}

int FindAnotherEdge(Surface *pmesh, int e0, int e1)
{
  int e2,k;
  int nbrVtxNum,nbrPolyNum,nbrEdgeNum,dummyEdge;
  int nbrVtxID[NBR_MAX_NUM],nbrPolyID[NBR_MAX_NUM],nbrEdgeID[NBR_MAX_NUM];
  int vId, vId0, vId1, vId2;

  // First find the vtx ID of the intersection of e0 and e1 
//  vId0 = pmesh->edge_table.edge[e0].V[0];
 // vId1 = pmesh->edge_table.edge[e0].V[1];
 // vId  = pmesh->edge_table.edge[e1].V[0];
 // vId2 = pmesh->edge_table.edge[e1].V[1]; 
  vId0 = pmesh->edges[e0].vertexID[0];
  vId1 = pmesh->edges[e0].vertexID[1];
  vId  =  pmesh->edges[e1].vertexID[0];
  vId2 =  pmesh->edges[e1].vertexID[1];

  if (vId == vId0) {
	vId0 = vId1;
	vId1 = vId2;
  }
  else if(vId == vId1) 
	vId1 = vId2;
  else if (vId2 == vId0) {
	vId0 = vId;
	vId = vId2;
  }
  else {
	vId1 = vId;
	vId = vId2;
  }
  
  GetNbrV2(vId, pmesh,nbrVtxID,nbrPolyID,nbrEdgeID,&nbrVtxNum,&nbrPolyNum,&nbrEdgeNum,1);
  for (k = 0; k < nbrEdgeNum; k++) {
	//if (pmesh->edge_table.edge[nbrEdgeID[k]].V[0] == vId0 && pmesh->edge_table.edge[nbrEdgeID[k]].V[1] == vId1)
	  if (pmesh->edges[nbrEdgeID[k]].vertexID[0] == vId0 && pmesh->edges[nbrEdgeID[k]].vertexID[1] == vId1)
	  return(nbrEdgeID[k]);
	//if (pmesh->edge_table.edge[nbrEdgeID[k]].V[0] == vId1 && pmesh->edge_table.edge[nbrEdgeID[k]].V[1] == vId0)
	if (pmesh->edges[nbrEdgeID[k]].vertexID[0] == vId1 && pmesh->edges[nbrEdgeID[k]].vertexID[1] == vId0)
	  return(nbrEdgeID[k]);
  }
  return(-1);
}
