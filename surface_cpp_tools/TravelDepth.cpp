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

#include "TravelDepth.h"

#include <vtkPolyDataReader.h>
#include <vtkHull.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointLocator.h>
#include <math.h>

TravelDepth::TravelDepth(char *fileName)
{
    vtkPolyDataReader* reader=vtkPolyDataReader::New();
    reader->SetFileName(fileName);
    reader->Update();

    m_mesh=vtkPolyData::New();
    m_mesh->DeepCopy(reader->GetOutput());

    reader->Delete();

    Initialize();
}

void TravelDepth::ComputeDepth()
{
    ComputeConvexHull();

    m_hullLocator = vtkCellLocator::New();
    m_hullLocator->SetDataSet(m_hull);
    m_hullLocator->BuildLocator();

    AffectDepthToVisiblePoints();
    DoublePropagationLoop();


    m_mesh->GetPointData()->SetScalars(m_depth);
    m_mesh->Update();
}

vtkDoubleArray *TravelDepth::GetDepth()
{
    return m_depth;
}

void TravelDepth::WriteIntoFile(char *fileName)
{
    vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
    writer->SetFileName(fileName);
    writer->SetInput(m_mesh);
    writer->Update();
    writer->Write();
    writer->Delete();
}

TravelDepth::TravelDepth(vtkPolyData* mesh)
{
    m_mesh=vtkPolyData::New();
    m_mesh->DeepCopy(mesh);

    Initialize();
}

void TravelDepth::Initialize()
{
    m_meshAnalyser = new MeshAnalyser(m_mesh);
    m_depth = vtkDoubleArray::New();
    m_normalizedDepth = vtkDoubleArray::New();

    m_meshLocator = vtkCellLocator::New();
    m_meshLocator->SetDataSet(m_mesh);
    m_meshLocator->BuildLocator();

    m_meshPointLocator = vtkPointLocator::New();
    m_meshPointLocator->SetDataSet(m_mesh);
    m_meshPointLocator->BuildLocator();

    m_lowConfidenceIds = vtkIdList::New();
    m_confidence = vtkDoubleArray::New();
    m_referenceIds = vtkIdList::New();

    for(int i = 0; i<m_mesh->GetNumberOfPoints(); i++)
    {
        m_lowConfidenceIds->InsertNextId(i);
        m_confidence->InsertNextValue(0);
    }

    m_maxIt = 20;
    m_maxConfidence = 7;

    m_maxBound=500;
}

void TravelDepth::ComputeConvexHull()
{
    //convex hull construction
    int recPlan=3;

    vtkHull *hull = vtkHull::New();
    hull->SetInputConnection(m_mesh->GetProducerPort());
    hull->AddRecursiveSpherePlanes(recPlan);
    hull->Update();

    m_hull = vtkPolyData::New();
    m_hull->DeepCopy(hull->GetOutput());
    m_hull->Update();

    cout<<"Hull generated"<<endl;
}

void TravelDepth::AffectDepthToVisiblePoints()
{

    double point1[3], point2[3];
    int  subid;
    vtkIdType cellid;
    double ec=0;

    double threshold=0.2;

    for(int i = 0; i < m_mesh->GetNumberOfPoints(); i++)
    {
        m_mesh->GetPoint(i,point1);
        m_hullLocator->FindClosestPoint(point1,point2,cellid,subid, ec);
        ec=sqrt(ec);

        //If the point is very close
        if (ec<threshold)
        {
            m_depth->InsertNextValue(ec);
            ReachMaxConfidence(i);
        }
        else
        {
            //if there is no other point of the surface between the point and the convex hull
            if( !CheckIntersection(point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]) )
            {
                m_depth->InsertNextValue(ec);
                ReachMaxConfidence(i);
            }
            else
            {
                m_depth->InsertNextValue(m_maxBound);
            }
        }
    }
}

void TravelDepth::DoublePropagationLoop()
{
    int iteration = 0;

    while(m_lowConfidenceIds->GetNumberOfIds()>0 && iteration < m_maxIt)
    {
        iteration++;

        GeodesicPropagation();
//        EuclideanPushPropagation();
        EuclideanPullPropagation();

    }

}

void TravelDepth::GeodesicPropagation()
{
    vtkIdList* tempReferencePoints = vtkIdList::New();
    tempReferencePoints->DeepCopy(m_referenceIds);
    m_referenceIds->Reset();

    vtkIdList* newReferencePoints = vtkIdList::New();

    vtkIdType curId, localId;

    vtkIdList* neighbors = vtkIdList::New();

    double point1[3], point2[3];

    double ec;

    vtkIdList* changed = vtkIdList::New();
    vtkIdList* candidates = vtkIdList::New();


    while(tempReferencePoints->GetNumberOfIds()>0)
    {
        cout<<tempReferencePoints->GetNumberOfIds()<<endl;

        for(int i = 0; i<tempReferencePoints->GetNumberOfIds(); i++)
        {
            curId = tempReferencePoints->GetId(i);
            m_mesh->GetPoint(curId, point1);

            if(m_confidence->GetValue(curId)>0)
            {
                neighbors->Reset();
                GetPointNeighbors(curId,neighbors);

                for(int j=0;j<neighbors->GetNumberOfIds();j++)
                {
                    localId=neighbors->GetId(j);
                    m_mesh->GetPoint(localId, point2);

                    candidates->InsertUniqueId(localId);

                    ec=sqrt(vtkMath::Distance2BetweenPoints(point1,point2));

                    if(UpdateDepth(curId, localId, ec))
                    {
                        changed->InsertUniqueId(localId);
                        newReferencePoints->InsertUniqueId(localId);
                        m_confidence->SetValue(localId, m_confidence->GetValue(curId)-1);
                    }
                }
            }
            else
            {
//                cout<<"a point in the reference set with a low confidence..."<<endl;
            }
        }
        tempReferencePoints->Reset();
        tempReferencePoints->DeepCopy(newReferencePoints);
        newReferencePoints->Reset();
    }


    UpdateConfidence(candidates, changed);

    tempReferencePoints->Delete();
}

void TravelDepth::EuclideanPushPropagation()
{
    vtkIdList* tempReferencePoints = vtkIdList::New();
    tempReferencePoints->DeepCopy(m_referenceIds);
    m_referenceIds->Reset();

    vtkIdList* newReferencePoints = vtkIdList::New();

    vtkIdType curId, localId;

    vtkIdList* neighbors = vtkIdList::New();

    double point1[3], point2[3];

    double ec;

    vtkIdList* changed = vtkIdList::New();
    vtkIdList* candidates = vtkIdList::New();

    double radius = 2;


    while(tempReferencePoints->GetNumberOfIds()>0)
    {
        for(int i = 0; i<tempReferencePoints->GetNumberOfIds(); i++)
        {
            curId = tempReferencePoints->GetId(i);
            m_mesh->GetPoint(curId, point1);

            if(m_confidence->GetValue(curId)>0)
            {
                neighbors->Reset();

                m_meshPointLocator->FindPointsWithinRadius(radius, point1, neighbors);

                for(int j=0;j<neighbors->GetNumberOfIds();j++)
                {
                    localId=neighbors->GetId(j);
                    m_mesh->GetPoint(localId, point2);

                    candidates->InsertUniqueId(localId);

                    ec=sqrt(vtkMath::Distance2BetweenPoints(point1,point2));

                    if(UpdateDepth(curId, localId, ec))
                    {
                        changed->InsertUniqueId(localId);
                        newReferencePoints->InsertUniqueId(localId);
                    }
                }
            }
        }
        tempReferencePoints->Reset();
        tempReferencePoints->DeepCopy(newReferencePoints);
        newReferencePoints->Reset();
    }


    UpdateConfidence(candidates, changed);

    tempReferencePoints->Delete();
}

void TravelDepth::EuclideanPullPropagation()
{
    vtkIdList* tempCandidatePoints = vtkIdList::New();
    tempCandidatePoints->DeepCopy(m_lowConfidenceIds);

    vtkIdType localId;
    vtkIdType closestId;
    double point1[3], point2[3];

    double ec;

    vtkPoints* referencePoints = vtkPoints::New();

    vtkIdList* changed = vtkIdList::New();


//    cout<<m_referenceIds->GetNumberOfIds()<<endl;

//    if(m_referenceIds->GetNumberOfIds() < m_mesh->GetNumberOfPoints()*0.05)
//    {
//        RefillReferenceSet();
//    }

    for(int i = 0;i<m_referenceIds->GetNumberOfIds();i++)
    {
        referencePoints->InsertNextPoint(m_mesh->GetPoints()->GetPoint(m_referenceIds->GetId(i)));
    }

    vtkPolyData* referencePolyData = vtkPolyData::New();
    referencePolyData->SetPoints(referencePoints);
    referencePolyData->Update();

    vtkPointLocator* referenceLocator= vtkPointLocator::New();
    referenceLocator->SetDataSet(referencePolyData);
    referenceLocator->BuildLocator();
    referenceLocator->Update();

    for(int i = 0; i<m_lowConfidenceIds->GetNumberOfIds(); i++)
    {
        localId = m_lowConfidenceIds->GetId(i);

        m_mesh->GetPoint(localId,point1);
        closestId = referenceLocator->FindClosestPoint(point1);
        m_mesh->GetPoint(closestId,point2);

        //if there is no other point of the surface between the point and the convex hull
        if( !CheckIntersection(point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]) )
        {
            ec=sqrt(vtkMath::Distance2BetweenPoints(point1,point2));
            if(UpdateDepth(closestId,localId,ec))
            {
                changed->InsertUniqueId(localId);
            }
        }
    }

    DecreaseReferenceSet(changed);
    UpdateConfidence(m_lowConfidenceIds, changed);

}

void TravelDepth::CompleteGeodesicPropagation()
{
    vtkIdList* badIds = vtkIdList::New();
    vtkIdList* neighbors = vtkIdList::New();

    vtkIdType localId;
    double ec;
    double point1[3], point2[3];

    for(int i = 0; i<m_mesh->GetNumberOfPoints(); i++)
    {
        if(m_confidence->GetValue(i) <= 0)
        {
            badIds->InsertNextId(i)  ;
        }
    }

    while(badIds->GetNumberOfIds()>0)
    {
        for(int i = 0; i<m_mesh->GetNumberOfPoints(); i++)
        {
            neighbors->Reset();
            GetPointNeighbors(i,neighbors);

            for(int j=0;j<neighbors->GetNumberOfIds();j++)
            {
                localId=neighbors->GetId(j);
                m_mesh->GetPoint(localId, point2);

                ec=sqrt(vtkMath::Distance2BetweenPoints(point1,point2));

                if(UpdateDepth(i, localId, ec))
                {
                    badIds->DeleteId(localId);
                }
            }
        }
    }
}

void TravelDepth::DecreaseReferenceSet(vtkIdList* changed)
{
    vtkIdType currentId;

    for(int i = 0; i<m_referenceIds->GetNumberOfIds(); i++)
    {
        currentId = m_referenceIds->GetId(i);
        if(changed->IsId(currentId) < 0)
        {
            m_referenceIds->DeleteId(currentId);
        }
    }
}

void TravelDepth::RefillReferenceSet()
{
    int subNumber = floor(m_mesh->GetNumberOfPoints()/5.0);

    vtkIdType randomId;

    for(int i = 0; i<subNumber; i++)
    {
        randomId = rand() % m_mesh->GetNumberOfPoints();
        m_referenceIds->InsertUniqueId(randomId);
    }
}

void TravelDepth::ReachMaxConfidence(int i)
{
    m_confidence->SetValue(i, m_maxConfidence);
    m_lowConfidenceIds->DeleteId(i);
    m_referenceIds->InsertUniqueId(i);
}

void TravelDepth::GetPointNeighbors(vtkIdType id, vtkIdList* neighbors)
{
    vtkIdList *cellids = vtkIdList::New();
    m_mesh->GetPointCells(id,cellids);
    neighbors->Reset();
    for (int i = 0; i < cellids->GetNumberOfIds(); i++)
    {
        for(int j = 0; j < m_mesh->GetCell(cellids->GetId(i))->GetNumberOfPoints(); j++)
        {
            if(m_mesh->GetCell(cellids->GetId(i))->GetPointId(j)< m_mesh->GetNumberOfPoints())
            {
                neighbors->InsertUniqueId(m_mesh->GetCell(cellids->GetId(i))->GetPointId(j));
            }
        }
    }
    cellids->Delete();
}

bool TravelDepth::UpdateDepth(vtkIdType idParent, vtkIdType idChild, double valueDiff)
{
    double newValue = m_depth->GetValue(idParent) + valueDiff;

    bool changed = false;

    if(newValue < m_depth->GetValue(idChild) )
    {
        m_depth->SetValue(idChild, newValue);
        m_confidence->SetValue(idChild, m_confidence->GetValue(idParent)-1);
//        m_referenceIds->InsertUniqueId(idChild);
        changed = true;
    }
    return changed;
}

void TravelDepth::UpdateConfidence(vtkIdList *idsToCheck, vtkIdList *idsNotToUpdate)
{
    int curConf;
    vtkIdType localId;

    for(int i = 0; i<idsToCheck->GetNumberOfIds();i++)
    {
        localId = idsToCheck->GetId(i);

        if( idsNotToUpdate->IsId(localId) < 0)
        {
            curConf =m_confidence->GetValue(localId);
            if(curConf >= m_maxConfidence-1)
            {
                ReachMaxConfidence(localId);
            }
            else
            {
                m_confidence->SetValue(localId, curConf+1);
            }
        }
    }
}

bool TravelDepth::CheckIntersection(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double t, ptline[3], pcoords[3];
    int subId;

    double point1[3], point2[3];
    double shift[3];

    point1[0] = x1;
    point1[0] = y1;
    point1[0] = z1;

    point2[0] = x2;
    point2[0] = y2;
    point2[0] = z2;

    for(int k = 0; k<3; k++)
    {
        shift[k]=point2[k]-point1[k];
        point1[k]+=shift[k]*0.001;
        point2[k]+=shift[k]*0.001;
    }

    return (m_meshLocator->IntersectWithLine(point1, point2, 0.0001, t, ptline, pcoords, subId) ==1);
}

