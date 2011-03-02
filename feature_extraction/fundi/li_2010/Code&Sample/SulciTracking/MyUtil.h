#ifndef UTIL_H
#define UTIL_H

#include "SurfaceIO.h"
//#include "fmarch_v2.h"

#define NBR_MAX_NUM 5000
#define SIGN(a) (a>=0?1:-1)
#define MAXSULCINUM 50
#define MAXPATHLEN 10000
#define MAXBUFLEN 100
#define TINY 1e-6
#define PI 3.14159265358979
 
typedef struct 
{
  int ID;
  double lambda;
} VtxOnFace;

void GetNbrV2(int CurrentVtxID, Surface *surface,int *VtxID,int *PolyID,int *EdgeID,int *VtxNum,int *PolyNum,int *EdgeNum,int Layer);

int IsInSet(int *set,int N,int ID);

double MyDistance(Fvector3d V1,Fvector3d V2);

double MyDistanceSquare(Fvector3d V1,Fvector3d V2);

int IsObtuse(Surface *surface,VtxOnFace Position, int EdgeID);

int IsOnEdge(Surface *surface, int vId, int eId);

int VertexOnPolygon(Surface *surface, int VtxID, int PolyID);

int EdgeOnPolygon(Surface *surface, int VtxID, int PolyID);

double FindValueAndPosition(Surface *surface, VtxOnFace CurrentPosition, double * GeoDist, Fvector3d * V0);

int FindTriangle(Surface *surface, VtxOnFace V1, VtxOnFace V2);

#endif