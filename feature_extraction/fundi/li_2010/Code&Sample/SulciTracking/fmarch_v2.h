#ifndef FMARCHING_V2
#define FMARCHING_V2

#include <math.h>
/*
#include "pgstd.h"
#include "pgpm.h"
#include "pgdx.h"
#include "dxio.h"
#include "pgutil.h"
*/
#include "heap.h"
#include "MyUtil.h"
//#include "vector.h"

#define ALIVE 1
#define NBAND 2
#define FAWAY 3 
#define MIN(a,b) (a <b? a: b)
#define MAX(a,b) (a >b? a: b)
#define ABS(a) (a>0? a: -(a))
#define INFINITY 1000000000 

/* function declarations */
static void usage(void);
double *FastMarchMesh(Surface *surface, int *contour, int numinitvert, double *F,int terminator);
double *FastMarchMeshHalfVertex(Surface *mesh, VtxOnFace start, VtxOnFace end, double *F);

double ReCompute(int vIDc, int vIDa, int vIDb, Surface *surface, double *T, unsigned char *label, double Fc, int *BackPointer, Xheap *pH);
double ComputeTInAcute(double,double,double,double,double, double);
double *SearchConnectingVertexUsingFastMarching(Surface *mesh, int startVertexID, double *F, int *faceSulci, int geodesicThresh, int sulcusLabel, int *vertexID);
double *SearchConnectingHalfVertexUsingFastMarching(Surface *mesh, Segment *segment, VtxOnFace start, double *F, 
													int *edgeSulci,Fvector3d *edgeHalfVertex, int geodesicThresh, int sulcusLabel, VtxOnFace *end);
double *SearchConnectingHalfVertexUsingFastMarchingAfterFastMarching(Surface *mesh, Segment *segment, VtxOnFace start, double *F, 
													int *edgeSulci,int *vertexSulci, Fvector3d *edgeHalfVertex, int geodesicThresh, int sulcusLabel, VtxOnFace *end);

double *SearchConnectingHalfVertexAndVertex(Surface *mesh, Segment *segment, SulciTrackingOut *out, VtxOnFace start, double *F, 
													int *edgeSulci,int *vertexSulci, Fvector3d *edgeHalfVertex, int geodesicThresh, int sulcusLabel, VtxOnFace *end);

#endif
