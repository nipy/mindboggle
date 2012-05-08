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


#include "MedialSurfaceComputer.h"

#include "vtkCellArray.h"
#include "vtkCleanPolyData.h"
#include "vtkPointLocator.h"
#include "vtkPolyDataWriter.h"


MedialSurfaceComputer::MedialSurfaceComputer(vtkPolyData *originalMesh, vtkPolyData *fundi)
{
    m_originalFundi = fundi;
    m_originalMesh = originalMesh;

    m_originalMeshAnalyser = new MeshAnalyser(m_originalMesh);

    m_originalPointLocator = vtkPointLocator::New();
    m_originalPointLocator->SetDataSet(m_originalMesh);
    m_originalPointLocator->BuildLocator();
    m_originalPointLocator->Update();

    m_pitsNorm = vtkDoubleArray::New();
    m_pitsNorm->SetNumberOfComponents(3);

    m_nbLayers = 15;

    m_processedFundi = vtkPolyData::New();
    m_medialSurface = vtkPolyData::New();

    ProcessFundi();

    BuildMedialSurface();
}



MedialSurfaceComputer::~MedialSurfaceComputer()
{
    delete m_originalMeshAnalyser;
    m_pitsNorm->Delete();
}

void MedialSurfaceComputer::WriteIntoFile(char *fileName)
{
    vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
    writer->SetFileName(fileName);

    writer->SetInput(m_medialSurface);
    writer->Update();
    writer->Write();
    writer->Delete();
}

void MedialSurfaceComputer::ProcessFundi()
{
    vtkIdType npts;
    vtkIdType* ids;
    vtkIdType i0, i1;

    double p0[3], p1[3];

    int nbSub = 2;

    vtkPoints* layer0Points = vtkPoints::New();

    int cnt = 0;

    int l;

    double tempPoint[3];
    vtkCellArray* layer0Cells = vtkCellArray::New();

    vtkCellArray* lines = vtkCellArray::New();
    lines->DeepCopy(m_originalFundi->GetLines());
    int moreCells = lines->GetNextCell(npts,ids);

    while( moreCells != 0)
    {
        i0 = ids[0];
        m_originalFundi->GetPoint(i0,p0);
        layer0Points->InsertNextPoint(p0);
        cnt ++;

        i1 = ids[1];
        m_originalFundi->GetPoint(i1,p1);

        for(int j = 1 ; j<nbSub ; j++)
        {
            l = (double)j/(double)(nbSub-1);
            for(int k = 0; k < 3 ; k++)
            {
                tempPoint[k] = (1-l) * p0[k] +  l * p1[k];
            }
            layer0Points->InsertNextPoint(tempPoint);
            cnt ++;

            vtkIdList* cello1 = vtkIdList::New();
            cello1->InsertNextId(cnt-1);
            cello1->InsertNextId(cnt-2);

            layer0Cells->InsertNextCell(cello1);
        }
        moreCells = lines->GetNextCell(npts,ids);
    }

    vtkPolyData* layer0 = vtkPolyData::New();
    layer0->SetPoints(layer0Points);
    layer0->SetPolys(layer0Cells);
    layer0->Update();

    vtkCleanPolyData* clean = vtkCleanPolyData::New();
    clean->SetInput(layer0);
    clean->Update();

    m_processedFundi->DeepCopy(clean->GetOutput());

    m_layers.push_back(m_processedFundi);

    AffectNormalsToFundi();

    layer0Points->Delete();
    layer0Cells->Delete();
    lines->Delete();
    layer0->Delete();
    clean->Delete();

}

void MedialSurfaceComputer::BuildMedialSurface()
{
    int nbLayers = 15;

    for (int l = 1; l<nbLayers ; l++)
    {
        vtkPoints* layerPoints = vtkPoints::New();

        MeshAnalyser* mal0 = new MeshAnalyser(m_layers.at(l-1));

        for (int i = 0; i< m_processedFundi->GetNumberOfPoints() ; i++)
        {
            double pointn[3];
            SelectNextPoint(pointn, mal0 , i);
            layerPoints->InsertNextPoint(pointn);
        }

        vtkPolyData* lpd = vtkPolyData::New();
        lpd->SetPoints(layerPoints);
        lpd->SetLines(m_processedFundi->GetLines());
        lpd->Update();

        m_layers.push_back(lpd);

        delete mal0;
    }

    BuildMeshFromLayers();

}

void MedialSurfaceComputer::AffectNormalsToFundi()
{
    double point0[3];
    vtkIdType i0;
    double n[3];

    for(int i = 0; i < m_processedFundi->GetNumberOfPoints() ; i++)
    {
        m_processedFundi->GetPoint(i,point0);
        i0 = m_originalPointLocator->FindClosestPoint(point0);
        m_originalMeshAnalyser->GetNormals()->GetTuple(i0, n);
        m_pitsNorm->InsertNextTuple(n);
    }
}

void MedialSurfaceComputer::BuildMeshFromLayers()
{
    vtkPoints* medialPoints = vtkPoints::New();
    vtkCellArray* medialCells = vtkCellArray::New();

    vtkIdList* neib = vtkIdList::New();
    vtkIdType nbNeib;

    int nbPits = m_processedFundi->GetNumberOfPoints();

    for (int i = 0; i< nbPits ; i++)
    {
        medialPoints->InsertNextPoint(m_layers.at(0)->GetPoint(i));
    }

    for (int l = 1; l < m_nbLayers ; l++)
    {
        MeshAnalyser* mal0 = new MeshAnalyser(m_layers.at(l-1));

        for (int i = 0; i< nbPits ; i++)
        {
            medialPoints->InsertNextPoint(m_layers.at(l)->GetPoint(i));

            mal0->GetPointNeighbors(i,neib);
            nbNeib = neib->GetNumberOfIds();

            for(int j = 0; j<nbNeib ; j++)
            {
                if(neib->GetId(j) < i)
                {
                    vtkIdList* cello1 = vtkIdList::New();

                    cello1->InsertNextId(i + (l-1) * nbPits);
                    cello1->InsertNextId(i + l * nbPits);
                    cello1->InsertNextId(neib->GetId(j) + (l-1) * nbPits);

                    medialCells->InsertNextCell(cello1);

                }
                else if (neib->GetId(j) > i)
                {
                    vtkIdList* cello1 = vtkIdList::New();

                    cello1->InsertNextId(i + (l-1) * nbPits);
                    cello1->InsertNextId(i + l * nbPits);
                    cello1->InsertNextId(neib->GetId(j) + l * nbPits);

                    medialCells->InsertNextCell(cello1);
                }
            }
        }
        delete mal0;
    }

    m_medialSurface->SetPoints(medialPoints);
    m_medialSurface->SetPolys(medialCells);
    m_medialSurface->Update();

}

void MedialSurfaceComputer::SelectNextPoint(double point[3], MeshAnalyser* mal0, int i)
{
    double point0[3];
    vtkIdList* neib = vtkIdList::New();
    double p0[3], p1[3], p2[3];
    double vec[3];

    mal0->GetMesh()->GetPoint(i,point0);

    m_pitsNorm->GetTuple(i,vec);


    for(int k = 0 ; k<3 ; k++)
    {
        p0[k] = point0[k] + vec[k];
        point[k] = p0[k];
    }



//    //find neighbors of the parent point
//    mal0->GetPointNeighbors(i,neib);
//    int nbNeib = neib->GetNumberOfIds();
//    int cnt = 0;
//    for(int j = 0 ; j< nbNeib ; j++)
//    {
//        if(neib->GetId(j) != i)
//        {
//            if(cnt == 1)
//            {
//                mal0->GetMesh()->GetPoint(neib->GetId(j),p2);
//                cnt++;
//            }
//            else if (cnt == 0)
//            {
//                mal0->GetMesh()->GetPoint(neib->GetId(j),p1);
//                cnt++;
//            }
//        }
//    }

//    pl2->FindPointsWithinRadius(radius,point0,radiusPoints);

//    nbInRadius = radiusPoints->GetNumberOfIds();

//    minScore = 1000000;


//    for(int r1 = 0 ; r1<nbInRadius ; r1++)
//    {
//        midPoints->GetPoint(radiusPoints->GetId(r1), cp);

//        d0 = vtkMath::Distance2BetweenPoints(p0, cp);

//        if(cnt == 0)
//        {
//            score1 = d0;
//        }
//        if(cnt == 1)
//        {
//            d1 = vtkMath::Distance2BetweenPoints(p1, cp);
//            score1 = d0 - d1;
//        }
//        if(cnt == 2)
//        {
//            d1 = vtkMath::Distance2BetweenPoints(p1, cp);
//            d2 = vtkMath::Distance2BetweenPoints(p2, cp);

//            score1 = d0 - d1 - d2;
//        }
//        if(minScore > score1)
//        {
//            minScore = score1;
//            bi = radiusPoints->GetId(r1);
//        }
//    }


//    neib->Delete();

}
