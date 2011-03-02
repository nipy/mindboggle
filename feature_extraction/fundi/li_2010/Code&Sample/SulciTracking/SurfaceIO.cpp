#include "SurfaceIO.h"


void ReadSurfaceAttribute(char *fname, Surface *surface)
{
	FILE *fp;
	char s[20];
	int pi[3];
	float pf[3];
	int m, n, junk;
	int vertexNum, faceNum;
	int neighborsNum, adjacentfacesNum;

	fp = fopen(fname,"rt");

	if (!fp)
	{
		printf("error in reading file %s\n",fname);
			
		exit(1);
	}
	
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"POINTS"));

	// read vertex
	fscanf(fp,"%d ",&vertexNum);
	fscanf(fp,"%s ",s);

	if(vertexNum <= 0)
	{
		printf("error: vertexNum = %d\n",vertexNum);
				
		exit(1);
	}

	surface->vertexNum = vertexNum;	
	surface->vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*vertexNum);

	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f %f %f ",pf,pf+1,pf+2);
		surface->vertex[n].x = pf[0];
		surface->vertex[n].y = pf[1];
		surface->vertex[n].z = pf[2];
	}

	printf("surface->vertexNum  = %d\n",surface->vertexNum);

	// read face
	fscanf(fp,"%s ",s);
	fscanf(fp,"%d ",&faceNum);
	fscanf(fp,"%d ",&junk);
	
	if(faceNum <= 0)
	{
		printf("error: faceNum = %d\n",faceNum);
				
		exit(1);
	}

	surface->faceNum = faceNum;
	surface->faces = (Face *)malloc(sizeof(Face)*faceNum);
	
	for(n=0; n<faceNum; n++)
	{
		fscanf(fp,"%d ",pi);
		fscanf(fp,"%d %d %d ",pi,pi+1,pi+2);
		surface->faces[n].v[0] = pi[0];
		surface->faces[n].v[1] = pi[1];
		surface->faces[n].v[2] = pi[2];
	}
	printf("surface->faceNum = %d\n",surface->faceNum);

	surface->curv1 = (float *)malloc(sizeof(float)*vertexNum);
	surface->curv2 = (float *)malloc(sizeof(float)*vertexNum);

	// read principal curvature 1
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"pcurv1Table"));
	
	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f",pf);
		surface->curv1[n] = pf[0];
	}

	// read principal curvature 2
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"pcurv2Table"));
	
	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f",pf);
		surface->curv2[n] = pf[0];
	}


	for(int i=0; i<4; i++)
		surface->dcurv[i] = (float *)malloc(sizeof(float)*vertexNum);

	// read curvature derivative 0
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"dcurv0Table"));
	
	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f",pf);
		surface->dcurv[0][n] = pf[0];
	}

	// read curvature derivative 1
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"dcurv1Table"));
	
	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f",pf);
		surface->dcurv[1][n] = pf[0];
	}

	// read curvature derivative 2
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"dcurv2Table"));
	
	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f",pf);
		surface->dcurv[2][n] = pf[0];
	}

	// read curvature derivative 3
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"dcurv3Table"));
	
	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f",pf);
		surface->dcurv[3][n] = pf[0];
	}


	surface->pdir1 = (Fvector3d *)malloc(sizeof(Fvector3d)*vertexNum);
	surface->pdir2 = (Fvector3d *)malloc(sizeof(Fvector3d)*vertexNum);

	// read principal direction 1
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"pdir1"));
	
	fscanf(fp,"%s ",s);

	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f %f %f ",pf,pf+1,pf+2);
		surface->pdir1[n].x = pf[0];
		surface->pdir1[n].y = pf[1];
		surface->pdir1[n].z = pf[2];
	}
	
	// read principal direction 2
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"pdir2"));
	
	fscanf(fp,"%s ",s);

	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f %f %f ",pf,pf+1,pf+2);
		surface->pdir2[n].x = pf[0];
		surface->pdir2[n].y = pf[1];
		surface->pdir2[n].z = pf[2];
	}

	surface->normal = (Fvector3d *)malloc(sizeof(Fvector3d)*vertexNum);


	// read normal
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"Normals"));
	
	fscanf(fp,"%s ",s);

	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%f %f %f ",pf,pf+1,pf+2);
		surface->normal[n].x = pf[0];
		surface->normal[n].y = pf[1];
		surface->normal[n].z = pf[2];
	}

	
	surface->neighbors.resize(vertexNum);
	// read neighbors
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"neighbors"));
	
	surface->numMaxNeighbors = 0;

	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%d ",&neighborsNum);
		
		if(neighborsNum > surface->numMaxNeighbors)
			surface->numMaxNeighbors  = neighborsNum;
		
		for(m=0; m<neighborsNum; m++)
		{
			fscanf(fp,"%d ",pi);
			surface->neighbors[n].push_back(pi[0]);
		}
	}

	surface->adjacentfaces.resize(vertexNum);
	// read adjacentfaces
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"adjacentfaces"));
	
	surface->numMaxAdjacentfaces = 0;

	for(n=0; n<vertexNum; n++)
	{
		fscanf(fp,"%d ",&adjacentfacesNum);
		
		if(adjacentfacesNum > surface->numMaxAdjacentfaces)
			surface->numMaxAdjacentfaces  = adjacentfacesNum;

		for(m=0; m<adjacentfacesNum; m++)
		{
			fscanf(fp,"%d ",pi);
			surface->adjacentfaces[n].push_back(pi[0]);
		}
	}

	// read across_edge
	do
	{	
		fscanf(fp,"%s ",s);
	}while(strcmp(s,"across_edge"));
	
	surface->across_edge = (Face *)malloc(sizeof(Face)*faceNum);
	
	for(n=0; n<faceNum; n++)
	{
		fscanf(fp,"%d ",pi);
		fscanf(fp,"%d %d %d ",pi,pi+1,pi+2);
		surface->across_edge[n].v[0] = pi[0];
		surface->across_edge[n].v[1] = pi[1];
		surface->across_edge[n].v[2] = pi[2];
	}

	fclose(fp);

	printf("surface->numMaxNeighbors = %d\n", surface->numMaxNeighbors);
	printf("surface->numMaxAdjacentfaces= %d\n", surface->numMaxAdjacentfaces);

	return;
}

void WriteSurfaceAttribute(char *fname, Surface *surface, bool saveEdge)
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

	// neighbors
	fprintf(fp,"neighbors\n");
		
	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%d ",surface->neighbors[i].size());
			
		for(j=0; j<surface->neighbors[i].size(); j++)
			fprintf(fp,"%d ",surface->neighbors[i][j]);
			
		fprintf(fp,"\n");
	} 

	// adjacentfaces
	fprintf(fp,"adjacentfaces\n");

	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%d ",surface->adjacentfaces[i].size());
			
		for(j=0; j<surface->adjacentfaces[i].size(); j++)
			fprintf(fp,"%d ",surface->adjacentfaces[i][j]);
			
		fprintf(fp,"\n");
	}
		
	// across_edge
	fprintf(fp,"across_edge\n");

	for(i=0;i<surface->faceNum;i++)
	{	  
		fprintf(fp,"3 %d %d %d\n",surface->across_edge[i].v[0], surface->across_edge[i].v[1], surface->across_edge[i].v[2]);
	}
		
	// edge
	if(saveEdge)
	{
		fprintf(fp,"edge\n");
		for(i=0; i<surface->edges.size(); i++)
			fprintf(fp,"%d %d %d %d\n",surface->edges[i].vertexID[0], surface->edges[i].vertexID[1], surface->edges[i].faceID[0], surface->edges[i].faceID[1]);

		fprintf(fp,"adjacentedges\n");
		for(i=0;i<surface->vertexNum;i++)
		{	  
			fprintf(fp,"%d ",surface->adjacentedges[i].size());
			
			for(j=0; j<surface->adjacentedges[i].size(); j++)
				fprintf(fp,"%d ",surface->adjacentedges[i][j]);
			
			fprintf(fp,"\n");
		} 
	}

	fclose(fp);

	return;
}


void WriteSurfaceAttributeLookupTable(char *fname, Surface *surface, int *label, int *lookupTable)
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

	// lookupTable
	fprintf(fp,"SCALARS lookupTable int\n");
	fprintf(fp,"LOOKUP_TABLE lookupTable\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%d\n",lookupTable[label[i]]);
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

	// neighbors
	fprintf(fp,"neighbors\n");
		
	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%d ",surface->neighbors[i].size());
			
		for(j=0; j<surface->neighbors[i].size(); j++)
			fprintf(fp,"%d ",surface->neighbors[i][j]);
			
		fprintf(fp,"\n");
	} 

	// adjacentfaces
	fprintf(fp,"adjacentfaces\n");

	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%d ",surface->adjacentfaces[i].size());
			
		for(j=0; j<surface->adjacentfaces[i].size(); j++)
			fprintf(fp,"%d ",surface->adjacentfaces[i][j]);
			
		fprintf(fp,"\n");
	}
		
	// across_edge
	fprintf(fp,"across_edge\n");

	for(i=0;i<surface->faceNum;i++)
	{	  
		fprintf(fp,"3 %d %d %d\n",surface->across_edge[i].v[0], surface->across_edge[i].v[1], surface->across_edge[i].v[2]);
	}
		
	fclose(fp);

	return;
}

void WriteSurfaceAttributeColor(char *fname, Surface *surface, int *label)
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

	// lookupTable
	fprintf(fp,"SCALARS lookupTable int\n");
	fprintf(fp,"LOOKUP_TABLE lookupTable\n");

	for(i=0;i<surface->vertexNum;i++)
	{
		fprintf(fp,"%d\n",label[i]);
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

	// neighbors
	fprintf(fp,"neighbors\n");
		
	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%d ",surface->neighbors[i].size());
			
		for(j=0; j<surface->neighbors[i].size(); j++)
			fprintf(fp,"%d ",surface->neighbors[i][j]);
			
		fprintf(fp,"\n");
	} 

	// adjacentfaces
	fprintf(fp,"adjacentfaces\n");

	for(i=0;i<surface->vertexNum;i++)
	{	  
		fprintf(fp,"%d ",surface->adjacentfaces[i].size());
			
		for(j=0; j<surface->adjacentfaces[i].size(); j++)
			fprintf(fp,"%d ",surface->adjacentfaces[i][j]);
			
		fprintf(fp,"\n");
	}
		
	// across_edge
	fprintf(fp,"across_edge\n");

	for(i=0;i<surface->faceNum;i++)
	{	  
		fprintf(fp,"3 %d %d %d\n",surface->across_edge[i].v[0], surface->across_edge[i].v[1], surface->across_edge[i].v[2]);
	}
		
	fclose(fp);

	return;
}



bool IsTwoVertexInFace(Surface *surface, int vertex1, int vertex2)
{
	for(int i=0; i<surface->adjacentfaces[vertex1].size(); i++)
	{
		if( surface->faces[surface->adjacentfaces[vertex1][i]].v[0] == vertex2 
			|| surface->faces[surface->adjacentfaces[vertex1][i]].v[1] == vertex2
			|| surface->faces[surface->adjacentfaces[vertex1][i]].v[2] == vertex2 )
		{
			return true;
		}
	}

    return false;
}



void ChangeInt(int *v0, int *v1)
{
	int temp;

	temp = *v0;	*v0 = *v1;	*v1 = temp;

	return;
}


void FindAnotherVertex(Surface *surface, int vertexID, int passedVertexID, int passedFaceID, int *nextVertexID)
{
	for(int j=0; j<surface->adjacentfaces[vertexID].size(); j++)
	{
		if(surface->faces[surface->adjacentfaces[vertexID][j]].v[0] == passedVertexID && surface->adjacentfaces[vertexID][j] != passedFaceID)
		{
			if(surface->faces[surface->adjacentfaces[vertexID][j]].v[1] == vertexID)
				*nextVertexID = surface->faces[surface->adjacentfaces[vertexID][j]].v[2];
			else
				*nextVertexID = surface->faces[surface->adjacentfaces[vertexID][j]].v[1];
			return;
		}
		else 	if(surface->faces[surface->adjacentfaces[vertexID][j]].v[1] == passedVertexID && surface->adjacentfaces[vertexID][j] != passedFaceID)
		{
			if(surface->faces[surface->adjacentfaces[vertexID][j]].v[0] == vertexID)
				*nextVertexID = surface->faces[surface->adjacentfaces[vertexID][j]].v[2];
			else
				*nextVertexID = surface->faces[surface->adjacentfaces[vertexID][j]].v[0];
			return;
		}
		else 	if(surface->faces[surface->adjacentfaces[vertexID][j]].v[2] == passedVertexID && surface->adjacentfaces[vertexID][j] != passedFaceID)
		{
			if(surface->faces[surface->adjacentfaces[vertexID][j]].v[1] == vertexID)
				*nextVertexID = surface->faces[surface->adjacentfaces[vertexID][j]].v[0];
			else
				*nextVertexID = surface->faces[surface->adjacentfaces[vertexID][j]].v[1];
			return;
		}
	}

	printf("Error: cannot find the next vertex!\n");

	return;
}

void GetSurfaceEdgesAttribute(Surface *surface)
{
	int i, j, k, n;
	int temp;
	bool existEdge;
	bool findNext;
	Edge tempEdge;
	int tempFace;
	int tempVertex;
	int iter;

	int neighborEdgeNum = 0;
	int currentEdgeNum = 0;
	int nVertexNum;
	int currentVertex;
	int existingNeighborNum;

	int *tempNeighbors;

	tempNeighbors = (int *)malloc(sizeof(int)*surface->numMaxNeighbors);

	// reorder neighbors for each vertex

	/*
	for(i=0; i<surface->vertexNum; i++)
	{
		 nVertexNum = surface->neighbors[i].size(); 
		for(j=0; j<nVertexNum-1; j++)
		{
			if( !IsTwoVertexInFace(surface,surface->neighbors[i][j], surface->neighbors[i][j+1]) )
			{
				for(k=j+2; k<nVertexNum; k++)
				{
					if( IsTwoVertexInFace(surface,surface->neighbors[i][j], surface->neighbors[i][k]) )
					{
						temp = surface->neighbors[i][k];
						surface->neighbors[i][k] = surface->neighbors[i][j+1];
						surface->neighbors[i][j+1] = temp;
						break;
					}
				}
			}
		}
	}
*/
	// check close
	for(i=0; i<surface->vertexNum; i++)
	{
		if(surface->neighbors[i].size() != surface->adjacentfaces[i].size())
			printf("i = %d, neighbors = %d, adjacentfaces = %d\n",i,surface->neighbors[i].size(),surface->adjacentfaces[i].size());
	}

	/*
	// check whether neighbor vertex is neighbor
	int count = 0;
	for(i=0; i<surface->vertexNum; i++)
	{
		if(surface->neighbors[i].size() == surface->adjacentfaces[i].size())
		for(j=0; j<surface->neighbors[i].size(); j++)
		{
			if(!IsTwoVertexInFace(surface,surface->neighbors[i][j], surface->neighbors[i][(j+1)%surface->neighbors[i].size()]))
			{
				printf("i = %d, j = %d\n",i,j);
				count++;
			}
		}
	}
*/

	surface->connectedges.resize(surface->vertexNum);

	// add vertex to edge
	for(i=0; i<surface->vertexNum; i++)
	{
		for(j=0; j<surface->neighbors[i].size(); j++)
		{
			tempEdge.vertexID[0] = i;
			tempEdge.vertexID[1] = surface->neighbors[i][j];
			
			tempEdge.faceID[0] = -1;
			tempEdge.faceID[1] = -1;

			// check whether the edge exists
			existEdge = false;
			for(n=0; n<surface->connectedges[tempEdge.vertexID[1]].size(); n++)
			{
				if( surface->edges[surface->connectedges[tempEdge.vertexID[1]][n]].vertexID[0] == tempEdge.vertexID[0]
				|| surface->edges[surface->connectedges[tempEdge.vertexID[1]][n]].vertexID[1] == tempEdge.vertexID[0] )
				{
					existEdge = true;
					break;
				}
			}
			
			// add temp edge to edgelist
			if(!existEdge)
			{
				if( surface->edges.size() >0)
					currentEdgeNum++;
				
				surface->edges.push_back(tempEdge);
				surface->connectedges[i].push_back(currentEdgeNum);
				surface->connectedges[tempEdge.vertexID[1]].push_back(currentEdgeNum);
			}
		}
	}
	
	for(i=0; i<surface->vertexNum; i++)
	{
		if(surface->connectedges[i].size() != surface->neighbors[i].size())
			printf("Error: surface->connectedges[%d].size() = %d, surface->neighbors[%d].size() = %d",i,surface->connectedges[i].size(),i,surface->neighbors[i].size());
	}

	surface->edgeNum = surface->edges.size();
	printf("surface->edgeNum  = %d\n",surface->edges.size());
	
	
	surface->adjacentedges.resize(surface->vertexNum);

	// find adjacent edges for each vertex
	for(i=0; i<surface->vertexNum; i++)
	{
		for(j=0; j<surface->adjacentfaces[i].size(); j++)
		{
			int faceID = surface->adjacentfaces[i][j];

			int v0, v1, v2;

			v0 = surface->faces[faceID].v[0];
			v1 = surface->faces[faceID].v[1];
			v2 = surface->faces[faceID].v[2];

			if(v1 == i )
				ChangeInt(&v0, &v1);
			else if(v2 == i)
				ChangeInt(&v0, &v2);
			
			for(k=0; k<surface->connectedges[v1].size(); k++)
			{
				if( v2 == surface->edges[surface->connectedges[v1][k]].vertexID[0] 
				|| v2 == surface->edges[surface->connectedges[v1][k]].vertexID[1] )
				{
					surface->adjacentedges[i].push_back(surface->connectedges[v1][k]);
					break;
				}
			}
		}
	}

	for(i=0; i<surface->vertexNum; i++)
	{
		if(surface->adjacentedges[i].size() != surface->adjacentfaces[i].size())
			printf("Error: surface->adjacentedges[%d].size() = %d, adjacentfaces[%d].size() = %d",i,surface->adjacentedges[i].size(),i,surface->adjacentfaces[i].size());
	}

	// add face to edge
	for(i=0; i<surface->edges.size(); i++)
	{
		int tempVertexID = surface->edges[i].vertexID[0];
		
		int count = 0;
		for(j=0; j<surface->adjacentfaces[tempVertexID].size(); j++)
		{
			int tempFaceID =  surface->adjacentfaces[tempVertexID][j];
			
			if( surface->faces[tempFaceID].v[0] == surface->edges[i].vertexID[1] 
			|| surface->faces[tempFaceID].v[1] == surface->edges[i].vertexID[1] 
			|| surface->faces[tempFaceID].v[2] == surface->edges[i].vertexID[1] )
			{
				surface->edges[i].faceID[count] = tempFaceID;
				count++;
			}
			
			if(count == 2)
				break;
		}
	}

	// 3 edges in each face
	surface->faceEdges.resize(surface->faceNum);

	for(i=0; i<surface->edgeNum; i++)
	{
		for(j=0; j<2; j++)
		{
			if(surface->edges[i].faceID[j] >= 0)
				surface->faceEdges[surface->edges[i].faceID[j]].push_back(i);
		}
	}

	for(i=0; i<surface->faceNum; i++)
		if(surface->faceEdges[i].size() != 3)
			printf("Error: surface->faceEdges[%d].size() = %d",i,surface->faceEdges[i].size());

	return;
}



bool IsVertexInFace(Surface *surface, int faceID, int vertexID)
{
	if(surface->faces[faceID].v[0] == vertexID || surface->faces[faceID].v[1] == vertexID || surface->faces[faceID].v[2] == vertexID )
		return true;
	else 
		return false;
}

void GetOneVertex(Surface *surface, int faceID, int vertexID0, int vertexID1, int *vertex)
{
	for(int i=0; i<3; i++)
	{
		if(surface->faces[faceID].v[i] != vertexID0 && surface->faces[faceID].v[i] != vertexID1)
		{
			*vertex = surface->faces[faceID].v[i];
			break;
		}
	}
	
	return;
}

void GetLeftRightVertex(Surface *surface, int vertexID, int nVertexID, int *leftVertexID, int *rightVertexID)
{
	bool findLeft = false;

	for(int i=0; i<surface->adjacentfaces[vertexID].size(); i++)
	{
		if( IsVertexInFace(surface, surface->adjacentfaces[vertexID][i], nVertexID) )
		{
			int vertex;
			
			GetOneVertex(surface, surface->adjacentfaces[vertexID][i], vertexID, nVertexID, &vertex);

			if( !findLeft)
			{
				*leftVertexID = vertex;
				findLeft = true;
			}
			else
			{
				*rightVertexID = vertex;
				
				return;
			}
		}
	}

	return;
}

void FindEdgeID(Surface *surface, int vertex1ID, int vertex2ID, int *edgeID)
{
	*edgeID = -1;

	for(int i=0; i<surface->connectedges[vertex1ID].size(); i++)
	{
		if( surface->edges[surface->connectedges[vertex1ID][i]].vertexID[0] == vertex2ID 
			|| surface->edges[surface->connectedges[vertex1ID][i]].vertexID[1] == vertex2ID)
		{
			*edgeID = surface->connectedges[vertex1ID][i];
			return;
		}
	}

	for(int i=0; i<surface->connectedges[vertex2ID].size(); i++)
	{
		if( surface->edges[surface->connectedges[vertex2ID][i]].vertexID[0] == vertex1ID 
			|| surface->edges[surface->connectedges[vertex2ID][i]].vertexID[1] == vertex1ID)
		{
			*edgeID = surface->connectedges[vertex2ID][i];
			return;
		}
	}

	return;
}

void NormalizeFvector3d(Fvector3d *vector)
{
	float norm;

	norm = sqrt(vector->x*vector->x + vector->y*vector->y + vector->z*vector->z);

	if(norm > 0)
	{
		vector->x /= norm;
		vector->y /= norm;
		vector->z /= norm;
	}

	return;
}

Fvector3d ScalarMultiplyFvector3d(float scalar, Fvector3d vectorIn)
{
	Fvector3d vectorOut;

	vectorOut.x = scalar * vectorIn.x;
	vectorOut.y = scalar * vectorIn.y;
	vectorOut.z = scalar * vectorIn.z;

	return vectorOut;
}

float Fvector3dDOTFvector3d(Fvector3d vector1, Fvector3d vector2)
{
	return(vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z);
}

Fvector3d Fvector3dADDFvector3d(Fvector3d vector1, Fvector3d vector2)
{
	Fvector3d vector;

	vector.x = vector1.x + vector2.x;
	vector.y = vector1.y + vector2.y;
	vector.z = vector1.z + vector2.z;
	
	return vector;
}

Fvector3d Fvector3dMINUSFvector3d(Fvector3d vector1, Fvector3d vector2)
{
	Fvector3d vector;

	vector.x = vector1.x - vector2.x;
	vector.y = vector1.y - vector2.y;
	vector.z = vector1.z - vector2.z;
	
	return vector;
}

float DistanceBetweenTwoVertices(Fvector3d vertex1, Fvector3d vertex2)
{
	float distance;

	distance = sqrt( (vertex1.x - vertex2.x)*(vertex1.x - vertex2.x) 
		+ (vertex1.y - vertex2.y)*(vertex1.y - vertex2.y) + (vertex1.z - vertex2.z)*(vertex1.z - vertex2.z) );

	return distance;
}

int FindShareEdge(Surface *surface, int faceID1, int faceID2)
{
	for(int i=0; i<surface->faceEdges[faceID1].size(); i++)
	{
		for(int j=0; j<surface->faceEdges[faceID2].size(); j++)
		{
			if(surface->faceEdges[faceID1][i] == surface->faceEdges[faceID2][j])
				return (surface->faceEdges[faceID1][i]);
		}
	}

	return -1;
}

