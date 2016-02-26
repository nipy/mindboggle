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

#include "Overlap.h"

#include "PointAreaComputer.h"

#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <stdio.h>

using namespace std;

Overlap::Overlap(vtkPolyData *mesh1, vtkPolyData *mesh2)
{
    m_pointsArea = vtkDoubleArray::New();
    m_labels1Areas = vtkDoubleArray::New();
    m_labels2Areas = vtkDoubleArray::New();
    m_labelsCommonAreas = vtkDoubleArray::New();

    m_jaccard = vtkDoubleArray::New();
    m_dice = vtkDoubleArray::New();

    m_mesh1 = mesh1;
    m_mesh2 = mesh2;

    m_maxLabel = 0;
}

Overlap::~Overlap()
{
    m_pointsArea->Delete();
    m_labels1Areas->Delete();
    m_labels2Areas->Delete();
    m_labelsCommonAreas->Delete();
}

void Overlap::ComputeOverlap()
{
//    ComputePointsArea();

    PointAreaComputer* areaComputer = new PointAreaComputer(m_mesh1);
    areaComputer->ComputeArea();
    m_pointsArea = areaComputer->GetArea();

    int cLabel1, cLabel2;
    double cArea1, cArea2, cAreaC;
    double iArea;
    double dice, jaccard;

    vtkDataArray* labels1 = m_mesh1->GetPointData()->GetScalars();
    vtkDataArray* labels2 = m_mesh2->GetPointData()->GetScalars();

    for(int i = 0; i<m_mesh1->GetNumberOfPoints(); i++)
    {
//        cout<<m_pointsArea->GetValue(i)<<" "<<labels1->GetTuple1(i)<<" "<<labels2->GetTuple1(i)<<endl;
        if(m_maxLabel < labels1->GetTuple1(i))
        {
            m_maxLabel = labels1->GetTuple1(i);
        }
        if(m_maxLabel < labels2->GetTuple1(i))
        {
            m_maxLabel = labels2->GetTuple1(i);
        }
    }

    for(int i = 0; i<m_maxLabel; i++)
    {
        m_labels1Areas->InsertNextValue(0);
        m_labels2Areas->InsertNextValue(0);
        m_labelsCommonAreas->InsertNextValue(0);
    }

    for(int i = 0; i<m_mesh1->GetNumberOfPoints(); i++)
    {
        cLabel1 = labels1->GetTuple1(i);
        cLabel2 = labels2->GetTuple1(i);

        cArea1 = m_labels1Areas->GetValue(cLabel1);
        cArea2 = m_labels2Areas->GetValue(cLabel2);

        iArea = m_pointsArea->GetValue(i);


        m_labels1Areas->SetValue(cLabel1, cArea1 + iArea);
        m_labels2Areas->SetValue(cLabel2, cArea2 + iArea);

        if(cLabel1 == cLabel2)
        {
            cAreaC = m_labelsCommonAreas->GetValue(cLabel1);
            m_labelsCommonAreas->SetValue(cLabel1, cAreaC + iArea);
        }
    }

    cout<<"Label  Dice  Jaccard Common Area1 Area2"<<endl;

    for(int i = 0; i<=m_maxLabel; i++)
    {
        cAreaC = m_labelsCommonAreas->GetValue(i);
        cArea1 = m_labels1Areas->GetValue(i);
        cArea2 = m_labels2Areas->GetValue(i);

        if(cAreaC > 0)
        {
            dice = cAreaC / ( (cArea1 + cArea2 ) / 2.0 );
            jaccard = cAreaC / (cArea1 + cArea2 - cAreaC);
        }
        else
        {
            dice = 0;
            jaccard = 0;
        }

        m_dice->InsertNextValue(dice);
        m_jaccard->InsertNextValue(jaccard);
        if(m_labels1Areas->GetValue(i) > 0 || m_labels2Areas->GetValue(i) > 0 )
        {
            cout<<i<<" "<<dice<<" "<<jaccard<<" "<< m_labelsCommonAreas->GetValue(i)<<" "<< m_labels1Areas->GetValue(i)<<" "<< m_labels2Areas->GetValue(i)<<endl;
        }
    }

}

void Overlap::WriteIntoFile(char* fileName)
{
    ofstream myfile;
    myfile.open (fileName, ios::out);


    myfile<<"Label  Dice  Jaccard Common Area1 Area2"<<endl;

    for(int i = 0; i<=m_maxLabel; i++)
    {
        if(m_labels1Areas->GetValue(i) > 0 || m_labels2Areas->GetValue(i) > 0 )
        {
            myfile <<i<<" "
                 <<m_dice->GetValue(i)<<" "
                 <<m_jaccard->GetValue(i)<<" "
                 << m_labelsCommonAreas->GetValue(i)<<" "
                 << m_labels1Areas->GetValue(i)<<" "
                 << m_labels2Areas->GetValue(i)<<endl;
        }
    }

    myfile.close();
}






vtkDoubleArray *Overlap::GetJaccard()
{
    return m_jaccard;
}

vtkDoubleArray *Overlap::GetDice()
{
    return m_dice;
}
