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


#ifndef MEDIALSURFACECOMPUTER_H
#define MEDIALSURFACECOMPUTER_H

#include "MeshAnalyser.h"

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
    MeshAnalyser* m_originalMeshAnalyser;

    vtkPolyData* m_processedFundi;

    vector<vtkPolyData*> m_layers;

    vtkDoubleArray* m_pitsNorm;

    int m_nbLayers;

    void ProcessFundi();

    void BuildMedialSurface();
    void AffectNormalsToFundi();
    void BuildMeshFromLayers();

    void SelectNextPoint(double point[3], MeshAnalyser *mal0, int i);


};




#endif // MEDIALSURFACECOMPUTER_H
