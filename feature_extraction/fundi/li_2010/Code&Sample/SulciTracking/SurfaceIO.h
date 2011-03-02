#ifndef SURFACEIO_H
#define SURFACEIO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include "mvcd.h"
using namespace std;
/*
typedef struct 
{
	float x;
	float y;
	float z;
} Fvector3d;
*/
typedef struct 
{
	int v[3];
} Face;

typedef struct
{
	int edgeID;
	int vertexID[2];
	int faceID[2];
} Edge;

// Structure for surface storage space for the algorithm 
typedef struct{
	int vertexNum;
	int faceNum;
	int edgeNum;
	unsigned char *isSulci; // Used to flag potential sulci vertices
	Fvector3d *vertex; 
	Face *faces; // the vertex ID in each face
	vector<Edge> edges; 
    Fvector3d *normal;
	Fvector3d *pdir1, *pdir2;
	float *curv1, *curv2;
	float *depth;
	float *dcurv[4];
	vector< vector<int> > neighbors; 		//   For each vertex, all neighboring vertices
	vector< vector<int> > adjacentfaces; 	//	 For each vertex, all neighboring faces
	vector< vector<int> > adjacentedges;  //  For each vertex, all neighboring edges
	vector< vector<int> > connectedges;  //  For each vertex, all connected edges
	Face *across_edge; // For each face, 3 neighboring faces that share the same edge
	vector< vector<int> > faceEdges; 		//   // the edge ID in each face
	int numMaxNeighbors;
	int numMaxAdjacentfaces;
	int numSulciStartVertex;	 // a variable which contains the number of possible starting vertices for sulci 
	int numSulciVertex;		 // A variable which contains the number of possible sulci vertices 
} Surface;

typedef struct 
{  
	int  vertexNum;
	int  faceNum;
	float highCurv;	// the upper threshold on the principal curvature 1
	float lowCurv;		// the lower threshold on the principal curvature 1
	float highDepth; // the upper threshold on the principal curvature 1
	float lowDepth;   // the lower threshold on the principal curvature 1
} Params;

typedef struct 
{
	int vertexID;
	int pValue;  // value of the sulci vertex: 2 if a sulci start vertex else 1 if an ordinary sulci vertex
                 // A sulci start vertex is not necessarily a vertex where a sulci begins but just the vertex at 
                 // which we begin to search for a sulcus
	float curv1;
	float curv2;
	float depth;
	Fvector3d vertex;
	Fvector3d pdir1;
	Fvector3d pdir2;
	Fvector3d normal;
	float dcurv[4];
	vector<int> neighbors;
	vector<int> adjacentfaces;
	char flag;   // to flag the sulci start vertex if they are already processed 
} SulciVertex;

typedef struct
{
	int n;			// number of vertices in the sulcus
	int len;
	int color;     
	char junction;  // =1 if the sulcus has a junction  
	char trash;		// Linked list postprocessing K-keep D-drop 
	Fvector3d *vertex;
	float *curv1; // the curvature 1
	float *depth; // the geodeisc depth
//	float *dcurv[4];
	Fvector3d *pdir2; // the pricipal direction for sulci following
	int *segmentID;
	int *vertexID;
	int *faceID;
	int *edgeID;
	bool *isJunction; // 1: vertexJunction;
	bool *isCenter;
} Sulci;

typedef struct
{
	int n;
	vector<Fvector3d> start;
	vector<Fvector3d> end;
	vector<Fvector3d> pdirS;
	vector<Fvector3d> pdirE;
	vector<unsigned char> isStartSegment; // the start and end segment 
	vector<bool> isStrictSegment; 
	vector<bool> endIsCenter;
	vector<float> curvS;
	vector<float> curvE;
	vector<float> edgeS;
	vector<float> edgeE;
	vector<int> face;
	vector<int> segment2Sulci;
	vector< vector<int> > Face2Segment;
	vector< vector<int> > adjacentSegment;
} Segment;


typedef struct
{
	int numberOfSulci;    // the number of sulci found
	vector< vector<int> > color2Sulci;
	Sulci *sulci;          // a pointer to a array of numberOfSulci sulci
	Surface surface; // a pointer to the surface 
} SulciTrackingOut;

void ReadSurfaceAttribute(char *fname, Surface *surface);
void WriteSurfaceAttribute(char *fname, Surface *surface, bool saveEdge);
void WriteSurfaceAttributeLookupTable(char *fname, Surface *surface, int *label, int *lookupTable);
void WriteSurfaceAttributeColor(char *fname, Surface *surface, int *label);

void GetSurfaceEdgesAttribute(Surface *surface);
void FindEdgeID(Surface *surface, int vertex1, int vertex2, int *edgeID);
bool IsVertexInFace(Surface *surface, int faceID, int vertexID);
void GetOneVertex(Surface *surface, int faceID, int vertexID0, int vertexID1, int *vertex);
void GetLeftRightVertex(Surface *surface, int vertexID, int nVertexID, int *leftVertexID, int *rightVertexID);

void NormalizeFvector3d(Fvector3d *vector);
Fvector3d ScalarMultiplyFvector3d(float scalar, Fvector3d vectorIn);
float Fvector3dDOTFvector3d(Fvector3d vector1, Fvector3d vector2);
Fvector3d Fvector3dADDFvector3d(Fvector3d vector1, Fvector3d vector2);
Fvector3d Fvector3dMINUSFvector3d(Fvector3d vector1, Fvector3d vector2);
float DistanceBetweenTwoVertices(Fvector3d vertex1, Fvector3d vertex2);
int FindShareEdge(Surface *surface, int faceID1, int faceID2);

#endif