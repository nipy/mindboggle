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

#include "GradientComputer.h"

#include "PointAreaComputer.h"

#include "MeshAnalyser.h"

#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkPolyDataNormals.h>

#include <stdio.h>

using namespace std;

GradientComputer::GradientComputer(vtkPolyData *mesh)
{
    m_mesh = mesh;

    vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
    pdn->SetInput(m_mesh);
    pdn->SetFeatureAngle(90);
    pdn->SplittingOff();
    pdn->Update();

    vtkPolyData* nor = pdn->GetOutput();
    nor->Update();
    m_normals = nor->GetPointData()->GetNormals();


    m_gradientValueArray = vtkDoubleArray::New();
    m_steepestDescentArray = vtkDoubleArray::New();
    m_steepestDescentArray->SetNumberOfComponents(3);
}


GradientComputer::~GradientComputer()
{
    m_gradientValueArray->Delete();
    m_steepestDescentArray->Delete();
}

void GradientComputer::WriteIntoFile(char *fileName)
{
    MeshAnalyser* ma = new MeshAnalyser(m_mesh);
    ma->WriteIntoFile(fileName,m_gradientValueArray);
    delete ma;
}

void GradientComputer::ComputeGradient()
{
    vtkDataArray* scalarField = m_mesh->GetPointData()->GetScalars();
    ComputeGradient(scalarField);
}

void GradientComputer::ComputeGradient(vtkDataArray *scalar)
{
    PointAreaComputer* pac = new PointAreaComputer(m_mesh);
    pac->ComputeArea();
    vtkDoubleArray* pointArea = pac->GetArea();
    double point[3];
    double neigPoint[3];
    double projNeigPoint[3];
    double normal[3];

    double scalarDiff;
    double gradientValue;

    double steepestDescent[3];
    double totalNeigArea;

    int neigId;

    vtkIdList* neig = vtkIdList::New();

    for(int i =0; i<m_mesh->GetNumberOfPoints();i++)
    {
        GetPointNeighbors(i,neig);
        m_mesh->GetPoint(i,point);
        m_normals->GetTuple(i,normal);

        totalNeigArea = 0;
        for(int k=0; k<3; k++)
        {
            steepestDescent[k] = 0;
        }

        for(int j=0; j<neig->GetNumberOfIds(); j++)
        {
            neigId = neig->GetId(j);
            m_mesh->GetPoint(neigId, neigPoint);
            ProjectPointToPlane(neigPoint,point,normal,projNeigPoint);
            totalNeigArea += pointArea->GetValue(neigId);
            scalarDiff = scalar->GetTuple1(neigId)-scalar->GetTuple1(i);

            for(int k=0; k<3; k++)
            {
                steepestDescent[k] += pointArea->GetValue(neigId)
                        * (projNeigPoint[k] - point[k])
                        * scalarDiff;
            }
        }
        for(int k=0; k<3; k++)
        {
            steepestDescent[k] /= totalNeigArea;
        }
        m_steepestDescentArray->InsertNextTuple(steepestDescent);
        gradientValue = vtkMath::Norm(steepestDescent);
        m_gradientValueArray->InsertNextValue(gradientValue);
    }

}

void GradientComputer::ProjectPointToPlane(double *pointToProject, double *inPlanePoint, double *normal, double *projectedPoint)
{
    double relativePoint[3];
    for(int k =0; k<3 ; k++)
    {
        relativePoint[k] = pointToProject[k] - inPlanePoint[k];
    }

    double dotProd = vtkMath::Dot(relativePoint, normal)/vtkMath::Norm(normal);

    for(int k =0; k<3 ; k++)
    {
        projectedPoint[k] = pointToProject[k] - dotProd * normal[k];
    }

}

void GradientComputer::GetPointNeighbors(vtkIdType id, vtkIdList* neighbors)
{
    vtkIdList *cellids = vtkIdList::New();
    m_mesh->GetPointCells(id,cellids);
    neighbors->Reset();
    for (int i = 0; i < cellids->GetNumberOfIds(); i++)
    {
        for(int j = 0; j < m_mesh->GetCell(cellids->GetId(i))->GetNumberOfPoints(); j++)
        {
            if(m_mesh->GetCell(cellids->GetId(i))->GetPointId(j)<m_mesh->GetNumberOfPoints())
            {
                neighbors->InsertUniqueId(m_mesh->GetCell(cellids->GetId(i))->GetPointId(j));
            }
        }
    }
    cellids->Delete();
}
