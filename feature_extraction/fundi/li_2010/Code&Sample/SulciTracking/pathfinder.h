
#include "MyUtil.h"
#include "fmarch_v2.h"

void PathFinder(Surface *pmesh, double* marchingSpeed, int* VtxID, int N, VtxOnFace* PathVtxList, int* pathLen);

bool PathFinderHalfVertex(Surface *pmesh, double* marchingSpeed, VtxOnFace start, VtxOnFace end, VtxOnFace* PathVtxList, int* pathLen);

double Gradient(Surface* pmesh,int EdgeID,VtxOnFace CurrentPosition,double *GeoDist,VtxOnFace *Position);

double InitialGuess (Surface*, VtxOnFace , double* , int* , int , VtxOnFace *, int, int *);

int DetermineblkEdge(Surface *pmesh, VtxOnFace v0, VtxOnFace v1, int* BE);

int DetermineblkVtx(Surface *pmesh, VtxOnFace v0, VtxOnFace v1, int* BV);

int FindAnotherEdge(Surface *pmesh, int e0, int e1);

void FindTwoEdges(Surface *pmesh, VtxOnFace v0, VtxOnFace v1, int *BE);