/* ********************************************************************
 * MeshAnalyser
 *
 * Copyright 2009 Joachim Giard
 * Université catholique de Louvain, Belgium
 * Apache License, Version 2.0
 *
 * Class containing methods to compute medial surfaces of brains.
 *
 * *******************************************************************/

#include <vtkPolyData.h>
#include <vtkPointLocator.h>
#include <vtkOBBTree.h>

#include <vector>

using namespace std;

#ifndef MEDIALSURFACECOMPUTER_H
#define MEDIALSURFACECOMPUTER_H

class MedialSurfaceComputer
{
public:
    MedialSurfaceComputer(vtkPolyData *originalMesh, vtkPolyData *fundi);
    ~MedialSurfaceComputer();

    void WriteIntoFile(char* fileName);

private:
    vtkPolyData* m_medialSurface;

    vtkPolyData* m_originalFundi;

    vtkPolyData* m_originalMesh;
    vtkPointLocator* m_originalPointLocator;
    vtkOBBTree* m_originalCellLocator;


    vtkPolyData* m_processedFundi;

    vector<vtkPolyData*> m_layers;

    vtkDataArray* m_normals;
    vtkDoubleArray* m_pitsNormals;

    vtkPoints* m_candidatePoints;
    vtkPolyData* m_candidatePolyData;

    vtkDataArray* m_metric;
    vtkPointLocator* m_candidateLocator;

    int m_nbLayers;

    vtkIdList* m_bins;

    void ProcessFundi();

    void BuildMedialSurface();
    void AffectNormalsToFundi();
    void BuildMeshFromLayers();
    void FindCandidatePoints();
    void FilterPointPosition(vtkPolyData* mesh);
    void AnchorPointToCandidate(double point[3]);
    void AnchorAllPointsToCandidates();

    void SelectNextPoint(double point[3], vtkPolyData *prevLayer, int i);

    void GetPointNeighbors(vtkPolyData* pd, vtkIdType id, vtkIdList *neighbors);

    bool ProjectInNormalDirection(vtkIdType pointId, double interPoint[]);

    void FrontVoronoi();

};




#endif // MEDIALSURFACECOMPUTER_H
