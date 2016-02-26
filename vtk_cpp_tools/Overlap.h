/* ********************************************************************
 * Overlap
 *
 * Copyright 2012 Joachim Giard
 * Université catholique de Louvain, Belgium
 * Apache License, Version 2.0
 *
 * Class containing various methods to compute overlaps
 * between pairs of labels
 *
 * *******************************************************************/

#ifndef OVERLAP_H
#define OVERLAP_H

#include "vtkPolyData.h"

class Overlap
{
public:
    // Constructor. mesh1 and mesh2 must be the same mesh
    // with a non-empty scalar field.
    Overlap(vtkPolyData* mesh1, vtkPolyData* mesh2);
    // destructor
    ~Overlap();

    //Compute Dice and Jaccard overlap score
    void ComputeOverlap();

    vtkDoubleArray* GetJaccard();
    vtkDoubleArray* GetDice();

    void WriteIntoFile(char* fileName);

private:

    vtkDoubleArray* m_pointsArea;
    vtkDoubleArray* m_labels1Areas;
    vtkDoubleArray* m_labels2Areas;
    vtkDoubleArray* m_labelsCommonAreas;

    vtkDoubleArray* m_jaccard;
    vtkDoubleArray* m_dice;

    vtkPolyData* m_mesh1;
    vtkPolyData* m_mesh2;

    int m_maxLabel;


};

#endif // OVERLAP_H
