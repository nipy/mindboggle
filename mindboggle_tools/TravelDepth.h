// *********************************************************
//      Copyright (c) Universite catholique de Louvain.
//      All rights reserved
//
//      Place du Levant, 2
//      1348 Louvain-la-Neuve
//      Belgium
//      Tel : +32 10 47 23 00
//
// *********************************************************

#ifndef TRAVELDEPTH_H
#define TRAVELDEPTH_H

#include "MeshAnalyser.h"

class TravelDepth
{
public:
    TravelDepth(vtkPolyData* mesh);
    TravelDepth(char* fileName);
    ~TravelDepth();

    void ComputeDepth();

    vtkDoubleArray* GetDepth();

    void WriteIntoFile(char* fileName);

private:

    void Initialize();

    void ComputeConvexHull();
    void AffectDepthToVisiblePoints();
    void DoublePropagationLoop();

    void GeodesicPropagation();
    void EuclideanPushPropagation();
    void EuclideanPullPropagation();
    void CompleteGeodesicPropagation();

    void DecreaseReferenceSet(vtkIdList *changed);
    void RefillReferenceSet();

    void ReachMaxConfidence(int i);

    void GetPointNeighbors(vtkIdType id, vtkIdList *neighbors);

    bool UpdateDepth(vtkIdType idParent, vtkIdType idChild, double valueDiff);
    void UpdateConfidence(vtkIdList* idsToCheck, vtkIdList* idsNotToUpdate);

    bool CheckIntersection(double x1, double y1, double z1, double x2, double y2, double z2);


    vtkDoubleArray* m_depth;
    vtkDoubleArray* m_normalizedDepth;

    vtkPolyData* m_mesh;
    vtkPolyData* m_hull;

    MeshAnalyser* m_meshAnalyser;

    vtkCellLocator *m_hullLocator;
    vtkCellLocator *m_meshLocator;
    vtkPointLocator *m_meshPointLocator;

    vtkDoubleArray* m_confidence;
    vtkIdList* m_lowConfidenceIds;
    vtkIdList* m_referenceIds;

    int m_maxIt;
    int m_maxConfidence;
    int m_maxBound;



};


#endif // TRAVELDEPTH_H
