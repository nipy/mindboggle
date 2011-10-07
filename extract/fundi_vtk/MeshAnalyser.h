/* ********************************************************************
 * MeshAnalyser
 * Author: Joachim Giard
 * Date (last update): 2 Nov 2009
 * 
 * Class containing various methods to apply on 3D meshes. The 
 * meshanalyser object takes a vtkPolyData object as input.
 * 
 * *******************************************************************/



#ifndef MESHANALYSER_H_
#define MESHANALYSER_H_

#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include "vtkCellArray.h"
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkBox.h>
#include <vtkBoundingBox.h>
#include <vtkObjectFactory.h>
#include <vector>

#define pi 3.14159265358979323846
#define sign(a)(a<0)?-1:1
#define norm(a)sqrt(pow(a[0],2)+pow(a[1],2))
using namespace std;

class MeshAnalyser
{
public:

	//constructor and destructor
	MeshAnalyser(vtkPolyData* mesh);
	MeshAnalyser(char* fileName);
	~MeshAnalyser();
	
	//SetMesh(char* arg1)
	//SetMesh(vtkPolyData* arg1)
	//
	//Read a polydata and make it the mesh to treat by this instance 
	//of meshAnalyser
	//
	//arg1 = file name *.vtk OR a polydata
	//Fill this->mesh with arg1
	void SetMesh(char* fileName);
	void SetMesh(vtkPolyData* mesh);
	
	//WriteIntoFile(char* arg1)
	//WriteIntoFile(char* arg1, char* arg2)
	//WriteIntoFile(char* arg1, vtkDataArray* arg2)
	//
	//Write the current mesh contained in this->mesh with or without 
	//a scalars field which could be defined in a vector from the 
	//exterior or a text refering to an instance vector.
	//
	//arg1 = file name *.vtk
	//arg2 = property as a dataArray OR as one of these property:
	//	"geoDist" = this->geoDistRing: geodesic distance
	//	"depth"	= this->depth: travel depth
	//	"curv" = this->curv: mean curvature
	//	"gCurv" = this->gCurv: gaussian curvature
	//	"test" = this->test: a test vector for new methods
	//	"surf" = this->pointSurf: the area of the point neighborhood
	//	"voronoi" = this->voronoiBin: the index of the voronoi zone
	void WriteIntoFile(char* fileName);
	void WriteIntoFile(char* fileName, char* prop);
	void WriteIntoFile(char* fileName, vtkDataArray* propExt);
	
	//ComputeTravelDepth(bool arg1)
	//
	//Compute the travel depth, i.e., the shortest distance between the 
	//mesh surface and its convex hull without going through the mesh 
	//interior.
	//
	//arg1 = bool (true -> depth normalized between 0 and 1)
	//Compute Travel Depth and fill this->depth
	void ComputeTravelDepth(bool norm);

	void ComputeTravelDepth(bool norm, bool sphere);

	
	//GeoDistRing(vtkIdType arg1, double arg2)
	//GeoDistRing(vtkIdType arg1, double arg2, double arg3)
	//
	//Compute a geodiesic distance with or without mesh simplification
	//using a front propagation method
	//
	//arg1 = starting point index
	//arg2 = maximal distance computed (approximation of the infinity)
	//arg3 = decimation factor (< 1) OR number of points for an approximation using a decimated mesh
	//compute geodesic distance from arg1 to all points and fill this->geoDistRing with distances,
	//this->geoCenter with arg1 and this->inRing with indexes of points closer than arg2
	void GeoDistRing(vtkIdType stPoint, double maxDist);
	void GeoDistRing(vtkIdType stPoint, double maxDist, double approx);
	
	//GetPointNeighbors(vtkIdType arg1, vtkIdList* arg2)
	//
	//return the 1-ring neighbors of a vertex
	//
	//arg1 = a point id
	//arg2 = an empty list of ids
	//fill arg2 with all the neighbors of arg1
	void GetPointNeighbors(vtkIdType id, vtkIdList* neighbors);

	//the same with the mesh as a parameter
	void GetPointNeighbors2(vtkIdType id, vtkIdList* neighbors, vtkPolyData* Mesh);

	//Compute the roughnes. Falta explicar!!
	void ComputeRoughness();

	//ComputeCurvature(double arg1)
	//
	//Compute the mean curvature alone using laplacian filtering and the
	//dot product between normal and displacement vectors. 
	//Faster than ComputeBothCurvatures.
	//
	//arg1 = resolution (0 <= arg1 <= 1) the more it's small, the more it's local
	//Compute the local surface curvature and fill this->curv
	void ComputeCurvature(double res);
	
	//ComputeBothCurvatures(double arg1)
	//
	//Compute the mean and the gaussian curvature using ratio between surfaces
	//of the original one and the one displaced in the normal direction (mean)
	//and between the original one and the laplacian filtered one (gaussian).
	//
	//arg1 = ray of computatation. Warning: the smoothing of the curvature
	//is done in a euclidean neighborhood, can make appear some errors.
	//The more the ray is small, the more it's local
	//Compute the local surface mean and gaussian curvature and fill this->curv and this->gCurv
	void ComputeBothCurvatures(double ray);
	
	//Compute the v=curv*rough
	void ComputeMask();

	void ComputeParamCorsini();


	void ComputeCorsini(MeshAnalyser* ma2);

	//compute thetaxy.  Find closest neighbor
	void ComputeSigmaCorsini(vtkPolyData* wMesh,vtkPolyData* Mesh, MeshAnalyser* ma2, MeshAnalyser* ma3);

	//compute thetaxy. Points equal index
	void ComputeSigmaCorsini2(MeshAnalyser* ma2);
	//Simplify(double arg1)
	//
	//Reduce the number of vertices and save the new mesh in a new polydata.
	//
	//arg1 = decimation factor (< 1) OR number of resulting points
	//Decimate this->mesh. this->simpl is filled with the result  
	void Simplify(double factor);

	//ComputeVoronoi
	//
	//Compute the geodesic voronoi regions on the surface given e list 
	//of centroids (slow and naive implementation)
	//
	//arg1 = a list of the indices of centroids
	//arg2 = maximal distance. points too far from centroids are labeled with -1
	//Compute the voronoi partition and fill this->voronoiBin with the indices
	// of the closest centroid for each point of the mesh
	vtkDoubleArray* ComputeVoronoi(vtkIdList* centroidsList);



	//detect prong points in a mesh. Calculates the sum of geodesic distances of each
	//point to the others points of the mesh.
	void ProngsDetection();
	
	void CenterData(vtkPolyData *p);

	//vector containg travel depth values (ComputeTravelDepth)
	vtkDoubleArray* derDepth;

	//vector containg travel depth values (ComputeTravelDepth)
	vtkDoubleArray* depth;

	vtkDoubleArray* VoronoiPatchesColors2;


	//Methods to get instance variables
	vtkDoubleArray* GetVoronoiPatchesColors2(){return this->VoronoiPatchesColors2;}
	vtkDoubleArray* GetGeoDistArray(){return this->geoDistArray;}
	vtkDoubleArray* GetGeoDistRing(){return this->geoDistRing;}
	vtkIdList* GetPointsInRing(){return this->inRing;}
	int GetNumberOfPoints(){return this->nbPoints;}
	vtkPolyData* GetMesh(){return this->mesh;}
	vtkDoubleArray* GetPointSurface(){return this->pointSurf;}
	vtkDoubleArray* GetTravelDepth(){return this->depth;}
	vtkDoubleArray* GetCurvature(){return this->curv;}
	vtkDoubleArray* GetRoughness(){return this->rough;}
	vtkDoubleArray* GetMask(){return this->mask;}
	vtkDoubleArray* GetMu(){return this->mu;}
	vtkDoubleArray* GetTheta(){return this->theta;}
	vtkDoubleArray* GetThetaXY(){return this->thetaxy;}
	vtkDoubleArray* GetL(){return this->L;}
	vtkDoubleArray* GetC(){return this->C;}
	vtkDoubleArray* GetS(){return this->S;}
	vtkDoubleArray* GetLmsdm(){return this->lmsdm;}
	vtkDoubleArray* GetMsdm(){return this->msdm;}
	vtkDoubleArray* GetDerDepth(){return this->derDepth;}
	vtkBoundingBox* GetBBox(){return this->BBox;}
	vtkDoubleArray* GetGaussianCurvature(){return this->gCurv;}
	vtkPolyData* GetSimpleMesh(){return this->simpl;}
	vtkIdList* GetVoronoiBin(){return this->voronoiBin;}
	int GetDiag(){return this->diag;}

protected:
	//same as GeoDistRing but on this->simpl. Used in GeoDistRing(vtkIdType arg1, double arg2, double arg3)
	void GeoDistRingSimple(vtkIdType stPoint, double maxDist);
	//same as GetPointNeighbors but on this->simpl. Used in GeoDistRing(vtkIdType arg1, double arg2, double arg3)
	void GetPointNeighborsSimple(vtkIdType id, vtkIdList* neighbors);
		
	//ComputeNormals
	//Compute the normals of the mesh and fill this->normals
	void ComputeNormals();
		
	//ComputePointSurface()
	//Compute surface corresponding to each point and fill this->pointSurf
	//Warning: works only for triangular meshes
	void ComputePointSurface();
	
	//ComputePointSurfaceSimple()
	//Compute surface corresponding to each point of this->simpl and fill this->pointSurfSimple
	void ComputePointSurfaceSimple();
	


private:

	//Analysed mesh (constructor or SetMesh)
	vtkPolyData* mesh;
	
	//Simplified mesh (Simplify)
	vtkPolyData* simpl;
	
	//vector containg point surface areas (ComputePointSurface)
	vtkDoubleArray* pointSurf;
	
	//vector containg this->simpl point surface areas (ComputePointSurfaceSimple)
	vtkDoubleArray* pointSurfSimple;
	
	//vector containg geodesic distances (GeoDistRing)
	vtkDoubleArray* geoDistRing;
	
	//vector containg geodesic distances for the simplified mesh (GeoDistRingSimple)
	vtkDoubleArray* geoDistRingSimple;
	
	//vector containg curvature values (ComputeCurvature or computeBothCurvatures)
	vtkDoubleArray* curv;
	
	//vector containg gaussian curvature values (ComputeBothCurvatures)
	vtkDoubleArray* gCurv;
	
	//vector containg gaussian roughness values (ComputeRoughness)
	vtkDoubleArray* rough;

	//vector containg mask v (ComputeMask)
	vtkDoubleArray* mask;

	//vector of the averages of mesh x
	vtkDoubleArray* mu;

	//vector of the standard deviation of mesh x
	vtkDoubleArray* theta;

	//vector of the joint standard deviation of mesh x and y
	vtkDoubleArray* thetaxy;

	//vector containing curvature comparison
	vtkDoubleArray* L;

	//vector containing contrast comparison
	vtkDoubleArray* C;

	//vector containing structure comparison
	vtkDoubleArray* S;

	//vector containing local distance measure lmsdm
	vtkDoubleArray* lmsdm;

	//vector containing local distance measure msdm (metrique de Corsini)
	vtkDoubleArray* msdm;

	//vector with the sum of geodesic distances from each point to all the others points of the mesh for every point.
	vtkDoubleArray* geoDistArray;

	//vector of the normals (ComputeNormals)
	vtkDataArray* normals;
		
	//square matrix containing geodesic distances between all pairs of points (GeoDistRingSimple)
	std::vector <float> allDist;

	//a test vector
	vtkDoubleArray* test;
	
	//list of points at a distance smaller than a parameter (GeoDistRing)
	vtkIdList* inRing;
	
	//list of points at a distance smaller than a parameter for the decimated mesh (GeoDistRingSimple)
	vtkIdList* inRingSimple;
	
	//starting point for current geodesic distances (GeoDistRing)
	vtkIdType geoCenter;
	
	//list of the predecessor of the current point for the geodesic distance with a front propagation (GeoDistRing);
	vtkIdList* predGeo;
	
	//Bounding Box
	vtkBoundingBox* BBox;

	//number of points in this->mesh
	int nbPoints;
	
	//diagonal of this->mesh
	int diag;

	//number of points in this->simpl
	int nbPointsSimple;
	
	//total surface of this->mesh (ComputePointSurface)
	double totalSurface;
	
	//total surface of this->simpl (ComputePointSurfaceSimple)
	double totalSurfaceSimple;
	
	//current approximation factor for this->simpl
	double approxFactor;
	
	//list of point of this->simpl which are close from this->mesh. length N x nbPoints
	// (geoDistRing(arg1,arg2,arg3))
	vtkIdList* close;

	//Indices of the voronoi bin in which the vertex is contained (ComputeVoronoi)
	vtkIdList* voronoiBin;

};

#endif
