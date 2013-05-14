/* ********************************************************************
 * Gradient Computer
 *
 * Copyright 2013 Joachim Giard
 * Université catholique de Louvain, Belgium
 * Apache License, Version 2.0
 *
 * Class containing various methods to compute gradients of a scalar
 * field on a mesh
 *
 * *******************************************************************/

#ifndef GRADIENTCOMPUTER_H
#define GRADIENTCOMPUTER_H

#include <vtkPolyData.h>

class GradientComputer
{
public:
    // Constructor.
    GradientComputer(vtkPolyData* mesh);

    // destructor
    ~GradientComputer();

    void WriteIntoFile(char* fileName);

    void ComputeGradient();
    void ComputeGradient(vtkDataArray* scalar);

private:

    void ProjectPointToPlane(double* pointToProject, double* inPlanePoint, double* normal, double* projectedPoint);
    void GetPointNeighbors(vtkIdType id, vtkIdList *neighbors);

    vtkDoubleArray* m_pointsArea;
    vtkPolyData* m_mesh;
    vtkDataArray* m_normals;
    vtkDoubleArray* m_gradientValueArray;
    vtkDoubleArray* m_steepestDescentArray;

};

#endif // GRADIENTCOMPUTER_H
