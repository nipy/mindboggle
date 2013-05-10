/* ********************************************************************
 * MeshAnalyser
 *
 * Copyright 2009 Joachim Giard
 * Université catholique de Louvain, Belgium
 * Apache License, Version 2.0
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
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vector>
#include <vtkCellLocator.h>


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
    //	"1color" = a single color
    // else = color by the index of the point
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

    //ComputeTravelDepth(bool arg1, vtkPolyData* arg2)
    //
    //Compute the travel depth, i.e., the shortest distance between the
    //mesh surface and a reference mesh without going through the mesh
    //interior.
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //arg2 = vtkPolyData: reference mesh
    //Compute Travel Depth and fill this->depth
    void ComputeTravelDepth(bool norm, vtkPolyData* refMesh);

    //ComputeTravelDepth(bool arg1)
    //
    //Compute the travel depth, i.e., the shortest distance between the
    //mesh surface and a reference mesh (the closed mesh) without going
    //through the mesh interior.
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //Compute Travel Depth and fill this->depth
    void ComputeTravelDepthFromClosed(bool norm);


    //ComputeGeodesicDepth(bool arg1)
    //
    //Compute the geodesic depth, i.e., the shortest distance between the
    //mesh surface and a reference mesh going along the surface mesh
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //Compute Geodesic Depth and fill this->geoDepth
    void ComputeGeodesicDepth(bool norm);


    //ComputeGeodesicDepth(bool arg1, vtkPolyData* arg2)
    //
    //Compute the geodesic depth, i.e., the shortest distance between the
    //mesh surface and its convex hull going along the surface mesh
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //arg2 = vtkPolyData: reference mesh
    //Compute Geodesic Depth and fill this->geoDepth
    void ComputeGeodesicDepth(bool norm, vtkPolyData *pq);


    //ComputeGeodesicDepth(bool arg1)
    //
    //Compute the geodesic depth, i.e., the shortest distance between the
    //mesh surface and iand a reference mesh (the closed mesh)
    //going along the surface mesh
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //Compute Geodesic Depth and fill this->geoDepth
    void ComputeGeodesicDepthFromClosed(bool norm);



    //ComputeEuclideanDepth(bool arg1)
    //
    //Compute the Euclidean depth, i.e., the shortest Euclidean distance between the
    //mesh surface and its convex hull
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //Compute Travel Depth and fill this->euclideanDepth
    void ComputeEuclideanDepth(bool norm);

    //ComputeEuclideanDepth(bool arg1, vtkPolyData* arg2)
    //
    //Compute the Euclidean depth, i.e., the shortest Euclidean distance between the
    //mesh surface and a reference mesh
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //arg2 = vtkPolyData: reference mesh
    //Compute Travel Depth and fill this->euclideanDepth
    void ComputeEuclideanDepth(bool norm, vtkPolyData* refMesh);

    //ComputeEuclideanDepth(bool arg1)
    //
    //Compute the Euclidean depth, i.e., the shortest Euclidean distance between the
    //mesh surface and a reference mesh (the closed mesh)
    //
    //arg1 = bool (true -> depth normalized between 0 and 1)
    //Compute Travel Depth and fill this->euclideanDepth
    void ComputeEuclideanDepthFromClosed(bool norm);


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

    //ComputeCurvature(double arg1)
    //
    //Compute the mean curvature alone using laplacian filtering and the
    //dot product between normal and displacement vectors.
    //Faster than ComputeBothCurvatures.
    //
    //arg1 = resolution (0 <= arg1 <= 1) the more it's small, the more it's local
    //Compute the local surface curvature and fill this->curv
    void ComputeCurvature(double res, int nbIt);

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

    //ComputePrincipalCurvatures()
    //
    //Compute the principal curvature, the mean and the gaussian curvature the deviation of the normal vectors.
    //
    //The procipal curvatures are computed by comparing the normal vectors (dot product) of the 1-ring neighbors
    //of the each point. The more the normal are divergent, the more the curvature value is low,
    // The more the normal are convergent, the more the curvature value is high.
    // Return the vector field of the minimal curvature direction
    //Compute the local surface mean and gaussian curvature and fill this->curv1, this->curv2,
    //this->curv and this->gCurv
    vtkDoubleArray *ComputePrincipalCurvatures(double nebSize);

    //Simplify(double arg1)
    //
    //Reduce the number of vertices and save the new mesh in a new polydata.
    //
    //arg1 = decimation factor (< 1) OR number of resulting points
    //Decimate this->mesh. this->simpl is filled with the result
    void Simplify(double factor);

    //ComputeClosedMesh
    //
    //Compute a mesh that correspond to a morphological closing
    //of the interior volume of the mesh. Stores it to this->closedMesh
    void ComputeClosedMesh(double kernelSize);

    //ComputeClosedMesh
    //
    //Compute a mesh that estimates a morphological closing
    //of the interior volume of the mesh. Stores it to this->closedMesh
    void ComputeClosedMeshFast();

    //Methods to get instance variables
    vtkDoubleArray* GetGeoDistRing(){return this->geoDistRing;}
    vtkIdList* GetPointsInRing(){return this->inRing;}
    int GetNumberOfPoints(){return this->nbPoints;}
    vtkPolyData* GetMesh(){return this->mesh;}
    vtkDoubleArray* GetPointSurface(){return this->pointSurf;}
    vtkDoubleArray* GetTravelDepth(){return this->depth;}
    vtkDoubleArray* GetEuclideanDepth(){return this->euclideanDepth;}
    vtkDoubleArray* GetCurvature(){return this->curv;}
    vtkDoubleArray* GetCurvature1(){return this->curv1;}
    vtkDoubleArray* GetCurvature2(){return this->curv2;}
    vtkDoubleArray* GetGaussianCurvature(){return this->gCurv;}
    vtkPolyData* GetSimpleMesh(){return this->simpl;}
    vtkDataArray* GetNormals(){return this->normals;}

    //ComputeHistogram
    //
    //Compute an histogram for a particular surface variable
    //
    //arg1 = variable name (the same name as for WriteIntoFile)
    //arg2 = the number of bins in the histogram
    //
    //Compute the number of occurance of the variable comprised in the arg2 bins
    //covering [min,max] the min and the max of the variable. Print out the result.
    void ComputeHistogram(char *prop, const int nbBins);
    void ComputeHistogram(vtkDataArray *data, int nbBins);

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

    //Initialize
    //
    //The common part of the constructors
    void Initialize();

    double IsIntersecting(double point1[3], double point2[2]);


    //Analysed mesh (constructor or SetMesh)
    vtkPolyData* mesh;

    //Simplified mesh (Simplify)
    vtkPolyData* simpl;

    //vector containing travel depth values (ComputeTravelDepth)
    vtkDoubleArray* depth;

    //vector containing geodesic depth values (ComputeGeodesicDepth)
    vtkDoubleArray* geoDepth;

    //vector containing point surface areas (ComputePointSurface)
    vtkDoubleArray* pointSurf;

    //vector containing this->simpl point surface areas (ComputePointSurfaceSimple)
    vtkDoubleArray* pointSurfSimple;

    //vector containing geodesic distances (GeoDistRing)
    vtkDoubleArray* geoDistRing;

    //vector containing geodesic distances for the simplified mesh (GeoDistRingSimple)
    vtkDoubleArray* geoDistRingSimple;

    //vector containing curvature values (ComputeCurvature or computeBothCurvatures or computePrincipalCurvatures)
    vtkDoubleArray* curv;

    //vector containing the first principal curvature values (ComputeCurvature or computeBothCurvatures or computePrincipalCurvatures)
    vtkDoubleArray* curv1;

    //vector containing the second principal curvature values (ComputeCurvature or computeBothCurvature or computePrincipalCurvaturess)
    vtkDoubleArray* curv2;

    //vector containing gaussian curvature values (ComputeBothCurvatures or computePrincipalCurvatures)
    vtkDoubleArray* gCurv;

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

    //number of points in this->mesh
    int nbPoints;

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

    //A mesh corresponding to the morphological closing
    //of the interior volume of the original mesh.
    vtkPolyData* closedMesh;

    //vector containing Euclidean depth values (ComputeEuclideanDepth)
    vtkDoubleArray* euclideanDepth;

    vtkCellLocator* meshLocator;
    vtkPolyData* medialSurface;

};

#endif
