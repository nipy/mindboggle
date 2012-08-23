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

#ifndef POINTAREACOMPUTER_H
#define POINTAREACOMPUTER_H

#include "vtkPolyData.h"

class PointAreaComputer
{
public:
    // Constructor.
    PointAreaComputer(vtkPolyData* mesh);
    // destructor
    ~PointAreaComputer();

    //Compute the voronoi area of each point
    void ComputeArea();

    vtkDoubleArray* GetArea();

    void WriteIntoFile(char* fileName);

private:
    vtkDoubleArray *ComputeThridArea();
    vtkDoubleArray *ComputeVoronoiArea();


    vtkDoubleArray* m_pointsArea;

    vtkPolyData* m_mesh;
};

#endif // POINTAREACOMPUTER_H
