#include "MyUtil.h"

/*
void GetNbrV2(int CurrentVtxID, Surface *surface,int *VtxID,int *PolyID,int *EdgeID,int *VNum,int *PNum,int *Enum,int Layer)
// This routine is used to find the nearest Layer layers of neighbor of vertex CurrentVtxID.
//   The neighboring vertex IDs are stored in VtxID and the neighboring polygon IDs are in PolyID.
 //  The VtxID is not ordered except the immidiate neighbors. Neither is PolyID.
{

  int i,j,k,l;
  int Vtx0,Poly0,Edge0;
  int VtxNum;
  int VtxAdded;
  int v0,v1,p0,p1;
  Neighbor nb;

  Vtx0 = 0;Poly0 = 0;Edge0 = 0;
  VtxAdded = 0;
  VtxID[Vtx0] = CurrentVtxID; 
  VtxNum = 1;
 
  for (i = 0;i<Layer;i++)
    {
      VtxAdded = 0;
      for (j = Vtx0;j<VtxNum;j++)
		{
		  nb = XUGetNeighbor(pmesh,VtxID[j]);
		  for (k = 0;k<nb.EN;k++)
			{
			  v0 = pmesh->edge_table.edge[nb.elist[k]].V[0];
			  v1 = pmesh->edge_table.edge[nb.elist[k]].V[1];
			  p0 = pmesh->edge_table.edge[nb.elist[k]].P[0];
			  p1 = pmesh->edge_table.edge[nb.elist[k]].P[1];
			  
			  if(!IsInSet(VtxID,VtxNum+VtxAdded,v0))
				{
				  VtxID[VtxNum+VtxAdded] = v0;
				  VtxAdded++;
				}
			  if(!IsInSet(VtxID,VtxNum+VtxAdded,v1))
				{
				  VtxID[VtxNum+VtxAdded] = v1;
				  VtxAdded++;
				}
			  
			  if(!IsInSet(PolyID,Poly0,p0))
				{
				  PolyID[Poly0++] = p0;
				  for (l = 0;l<3;l++)
					{
					  if(pmesh->edge_table.edge[pmesh->poly_table.poly[p0].elist[l]].V[0]!=CurrentVtxID && 
						 pmesh->edge_table.edge[pmesh->poly_table.poly[p0].elist[l]].V[1]!=CurrentVtxID &&
						 !IsInSet(EdgeID,Edge0,pmesh->poly_table.poly[p0].elist[l]))
						EdgeID[Edge0++] = pmesh->poly_table.poly[p0].elist[l];
					}	
				  
				}
			  if(!IsInSet(PolyID,Poly0,p1))
				{
				  PolyID[Poly0++] = p1;
				  for (l = 0;l<3;l++)
					{
					  if(pmesh->edge_table.edge[pmesh->poly_table.poly[p1].elist[l]].V[0]!=CurrentVtxID && 
						 pmesh->edge_table.edge[pmesh->poly_table.poly[p1].elist[l]].V[1]!=CurrentVtxID &&
						 !IsInSet(EdgeID,Edge0,pmesh->poly_table.poly[p1].elist[l]))
						EdgeID[Edge0++] = pmesh->poly_table.poly[p1].elist[l];
					}	
				  
				}
			  
			}
		}
      Vtx0 = VtxNum;
      VtxNum += VtxAdded;
    }
	
  *VNum = VtxNum;
  *PNum = Poly0;
  *Enum = Edge0;  
  
  return;
}
*/
void GetNbrV2(int CurrentVtxID, Surface *surface,int *VtxID,int *PolyID,int *EdgeID,int *VNum,int *PNum,int *Enum,int Layer)
{
	int i;

	*VNum = surface->neighbors[CurrentVtxID].size();

	for(i=0; i<surface->neighbors[CurrentVtxID].size(); i++)
	{
		VtxID[i] = surface->neighbors[CurrentVtxID][i];
	}


	*PNum = surface->adjacentfaces[CurrentVtxID].size();

	for(i = 0; i<surface->adjacentfaces[CurrentVtxID].size(); i++)
	{
		PolyID[i] = surface->adjacentfaces[CurrentVtxID][i];
	}


	*Enum = surface->adjacentedges[CurrentVtxID].size();

	for(i = 0; i<surface->adjacentedges[CurrentVtxID].size(); i++)
	{
		EdgeID[i] = surface->adjacentedges[CurrentVtxID][i];
	}

	return;

}


int IsInSet(int *set,int N,int ID)
{
	int k;
  
	if (N==0)
		return(0);
  
	for (k = 1;k<N;k++)
    {
		if(ID==set[k])
			return(1);
    }
	
	return(0);
}


int IsObtuse(Surface *surface, VtxOnFace Position, int EdgeID)
{
  
	Fvector3d V0, V1, V2;
	int ID1, ID2;
	double a, b, c;
 
	V1.x = surface->vertex[surface->edges[EdgeID].vertexID[0]].x;
	V1.y = surface->vertex[surface->edges[EdgeID].vertexID[0]].y;
	V1.z = surface->vertex[surface->edges[EdgeID].vertexID[0]].z;
	V2.x = surface->vertex[surface->edges[EdgeID].vertexID[1]].x;
	V2.y = surface->vertex[surface->edges[EdgeID].vertexID[1]].y;
	V2.z = surface->vertex[surface->edges[EdgeID].vertexID[1]].z;
	
	if(Position.lambda == 0.0)
    {
		V0.x = surface->vertex[Position.ID].x;
		V0.y = surface->vertex[Position.ID].y;
		V0.z = surface->vertex[Position.ID].z;
    }
	else
    {
		ID1 = surface->edges[Position.ID].vertexID[0];
		ID2 = surface->edges[Position.ID].vertexID[1];
		
		if(ID1 == surface->edges[EdgeID].vertexID[0] || ID1 == surface->edges[EdgeID].vertexID[1])
		{
			V0.x = surface->vertex[ID2].x;
			V0.y = surface->vertex[ID2].y;
			V0.z = surface->vertex[ID2].z;
		}
		else
		{
			V0.x = surface->vertex[ID1].x;
			V0.y = surface->vertex[ID1].y;
			V0.z = surface->vertex[ID1].z;
		}
	}
	
	a = MyDistanceSquare(V0,V1);
	b = MyDistanceSquare(V0,V2);
	c = MyDistanceSquare(V1,V2);
	
	if(c>a+b)
		return 1;
	else
		return 0;
}

int VertexOnPolygon(Surface *surface,int VtxID,int PolyID)
{
	int i, j, k;
	
	for(i=0; i<3; i++)
	{
		if(surface->faces[PolyID].v[i] == VtxID)
			return 1;
	}

	return 0;

/*
	//for(i = 0;i<pmesh->poly_table.poly[PolyID].EN;i++)
	for(i = 0;i<surface->faceEdges[PolyID].size();i++)
    {
		 j = surface->faceEdges[PolyID][i];
		 if(VtxID == surface->edges[j].vertexID[0] || VtxID == surface->edges[j].vertexID[1])
			return(1);
    }
	
	return (0);
*/
}

int EdgeOnPolygon(Surface *surface,int EdgeID,int PolyID)
{
	int i;
	
	for (i = 0;i<3;i++)
		if(EdgeID == surface->faceEdges[PolyID][i])
			return(1);
	
	return(0);
}

double FindValueAndPosition(Surface *surface, VtxOnFace CurrentPosition, double *GeoDist, Fvector3d *V0)
{
	int ID0, ID1, ID2;
	Fvector3d V1, V2;
	double b0;

	if(CurrentPosition.lambda == 0.0)
    {
		b0 = GeoDist[CurrentPosition.ID];
		ID0 = CurrentPosition.ID;

		V0->x = surface->vertex[ID0].x;
		V0->y = surface->vertex[ID0].y;
		V0->z = surface->vertex[ID0].z;
	}
	else
    {

		ID1 = surface->edges[CurrentPosition.ID].vertexID[0];
		ID2 = surface->edges[CurrentPosition.ID].vertexID[1];

		V1.x = surface->vertex[ID1].x;
		V1.y = surface->vertex[ID1].y;
		V1.z = surface->vertex[ID1].z;
		V2.x = surface->vertex[ID2].x;
		V2.y = surface->vertex[ID2].y;
		V2.z = surface->vertex[ID2].z;

		b0 = GeoDist[ID1]+CurrentPosition.lambda*(GeoDist[ID2]-GeoDist[ID1]);
		V0->x = V1.x+CurrentPosition.lambda*(V2.x-V1.x);
		V0->y = V1.y+CurrentPosition.lambda*(V2.y-V1.y);
		V0->z = V1.z+CurrentPosition.lambda*(V2.z-V1.z);
	}
	
	return(b0);
}


int FindTriangle(Surface* surface, VtxOnFace V1, VtxOnFace V2)
{
	int P0, P1;
	
	if (V1.lambda == 0)
	{
		//P0 = pmesh->edge_table.edge[V2.ID].P[0];
		P0 = surface->edges[V2.ID].faceID[0];
        
		if (VertexOnPolygon(surface, V1.ID, P0))
			return(P0);
		else
			//return(pmesh->edge_table.edge[V2.ID].P[1]);
			return(surface->edges[V2.ID].faceID[1]);
	}
	else 
	{
		//P0 = pmesh->edge_table.edge[V1.ID].P[0];
		P0 = surface->edges[V1.ID].faceID[0];
		if (V2.lambda == 0)
		{
			if (VertexOnPolygon(surface, V2.ID, P0))
				return(P0);
			else
				//return(pmesh->edge_table.edge[V1.ID].P[1]);
				return(surface->edges[V1.ID].faceID[1]);
		}
		else
		{
			if (EdgeOnPolygon(surface, V2.ID, P0))
				return(P0);
			else
				//return(pmesh->edge_table.edge[V1.ID].P[1]);
				return(surface->edges[V1.ID].faceID[1]);
		}
	}
}


double MyDistanceSquare(Fvector3d V1,Fvector3d V2)
{
	double d;
  
	d = (V1.x-V2.x)*(V1.x-V2.x)+(V1.y-V2.y)*(V1.y-V2.y)+(V1.z-V2.z)*(V1.z-V2.z);
	
	return(d);
}

double MyDistance(Fvector3d V1,Fvector3d V2)
{
	double d;
  
	d = (V1.x-V2.x)*(V1.x-V2.x)+(V1.y-V2.y)*(V1.y-V2.y)+(V1.z-V2.z)*(V1.z-V2.z);
  
	d = sqrt(d);
  
	return(d);
}


int IsOnEdge(Surface *surface, int vId, int eId)
{
	//if (vId == pmesh->edge_table.edge[eId].V[0] || vId == pmesh->edge_table.edge[eId].V[1])
	if (vId == surface->edges[eId].vertexID[0] || vId == surface->edges[eId].vertexID[1])
		return(1);
	else 
		return(0);
}
