#include "TriMesh.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
		int i, j;

        char *fnameIn = argv[1];
		char *fnameOut = argv[2];

		if(argc < 3)
		{
			printf("Usage: %s input.ply output.vtk\n", argv[0]);

			exit(1);
		}

		TriMesh *m = TriMesh::read(fnameIn);
        
		if (!m)
		{
			printf("error in reading file %s\n",fnameIn);
			
			exit(1);
		}

		m->need_faces();
		
		m->need_normals();
		
		m->need_neighbors();
		
		m->need_adjacentfaces();
		
		m->need_curvatures();
		
		m->need_dcurv();

		m->need_across_edge();

		// write vtk file
		FILE *fp;
		fp = fopen(fnameOut,"wt");
		
		fprintf(fp,"# vtk DataFile Version 3.0\n");
		fprintf(fp,"vtk output\n");
		fprintf(fp,"ASCII\n");
		fprintf(fp,"DATASET POLYDATA\n");
		fprintf(fp,"POINTS %d float\n",m->vertices.size());

		// vertex
		for(i=0;i<m->vertices.size();i++)
		{	  
			fprintf(fp,"%f %f %f\n",m->vertices[i][0],m->vertices[i][1],m->vertices[i][2]);
		} 

		// face
		fprintf(fp,"POLYGONS %d %d\n",m->faces.size(),m->faces.size()*4);

		for(i=0;i<m->faces.size();i++)
		{
			fprintf(fp,"3 %d %d %d\n",m->faces[i][0],m->faces[i][1],m->faces[i][2]);
		}

		// principal curvature 1
		fprintf(fp,"POINT_DATA %d\n",m->vertices.size());
		fprintf(fp,"SCALARS pcurv1 float\n");
		fprintf(fp,"LOOKUP_TABLE pcurv1Table\n");

		for(i=0;i<m->vertices.size();i++)
		{
			fprintf(fp,"%f\n",m->curv1[i]);
		}
		
		// principal curvature 2
		fprintf(fp,"SCALARS pcurv2 float\n");
		fprintf(fp,"LOOKUP_TABLE pcurv2Table\n");

		for(i=0;i<m->vertices.size();i++)
		{
			fprintf(fp,"%f\n",m->curv2[i]);
		}

		// curvature derivative 0 
		fprintf(fp,"SCALARS dcurv0 float\n");
		fprintf(fp,"LOOKUP_TABLE dcurv0Table\n");

		for(i=0;i<m->vertices.size();i++)
		{
			fprintf(fp,"%f\n",m->dcurv[i][0]);
		}

		// curvature derivative 1 
		fprintf(fp,"SCALARS dcurv1 float\n");
		fprintf(fp,"LOOKUP_TABLE dcurv1Table\n");

		for(i=0;i<m->vertices.size();i++)
		{
			fprintf(fp,"%f\n",m->dcurv[i][1]);
		}

		// curvature derivative 2 
		fprintf(fp,"SCALARS dcurv2 float\n");
		fprintf(fp,"LOOKUP_TABLE dcurv2Table\n");

		for(i=0;i<m->vertices.size();i++)
		{
			fprintf(fp,"%f\n",m->dcurv[i][2]);
		}

		// curvature derivative 3 
		fprintf(fp,"SCALARS dcurv3 float\n");
		fprintf(fp,"LOOKUP_TABLE dcurv3Table\n");

		for(i=0;i<m->vertices.size();i++)
		{
			fprintf(fp,"%f\n",m->dcurv[i][3]);
		}

		// principal direction 1
		fprintf(fp,"VECTORS pdir1 float\n");

		for(i=0;i<m->vertices.size();i++)
		{	  
			fprintf(fp,"%f %f %f\n",m->pdir1[i][0],m->pdir1[i][1],m->pdir1[i][2]);
		} 

		// principal direction 2
		fprintf(fp,"VECTORS pdir2 float\n");

		for(i=0;i<m->vertices.size();i++)
		{	  
			fprintf(fp,"%f %f %f\n",m->pdir2[i][0],m->pdir2[i][1],m->pdir2[i][2]);
		} 

		// normal
		fprintf(fp,"NORMALS Normals float\n");

		for(i=0;i<m->vertices.size();i++)
		{	  
			fprintf(fp,"%f %f %f\n",m->normals[i][0],m->normals[i][1],m->normals[i][2]);
		} 

		// neighbors
		fprintf(fp,"neighbors\n");
		
		for(i=0;i<m->vertices.size();i++)
		{	  
			fprintf(fp,"%d ",m->neighbors[i].size());
			
			for(j=0; j<m->neighbors[i].size(); j++)
				fprintf(fp,"%d ",m->neighbors[i][j]);
			
			fprintf(fp,"\n");
		} 

		// adjacentfaces
		fprintf(fp,"adjacentfaces\n");

		for(i=0;i<m->vertices.size();i++)
		{	  
			fprintf(fp,"%d ",m->adjacentfaces[i].size());
			
			for(j=0; j<m->adjacentfaces[i].size(); j++)
				fprintf(fp,"%d ",m->adjacentfaces[i][j]);
			
			fprintf(fp,"\n");
		}
		
		// across_edge
		fprintf(fp,"across_edge\n");

		for(i=0;i<m->faces.size();i++)
		{	  
			fprintf(fp,"3 %d %d %d\n",m->across_edge[i][0],m->across_edge[i][1],m->across_edge[i][2]);
		}

		fclose(fp);

		return 0;
}
