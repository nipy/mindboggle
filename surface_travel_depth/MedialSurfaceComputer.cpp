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

#include "vtkDelaunay2D.h"


MedialSurfaceComputer::MedialSurfaceComputer(vtkPolyData *originalMesh, vtkPolyData *fundi)
{
    m_originalFundi = fundi;
    m_originalMesh = originalMesh;
    m_metric = originalMesh->GetPointData()->GetScalars();

    m_originalMeshAnalyser = new MeshAnalyser(m_originalMesh);

    m_originalPointLocator = vtkPointLocator::New();
    m_originalPointLocator->SetDataSet(m_originalMesh);
    m_originalPointLocator->BuildLocator();
    m_originalPointLocator->Update();

    m_pitsNorm = vtkDoubleArray::New();
    m_pitsNorm->SetNumberOfComponents(3);

    m_nbLayers = 10;

    m_processedFundi = vtkPolyData::New();
    m_medialSurface = vtkPolyData::New();

    ProcessFundi();
    FindCandidatePoints();

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

    m_fundiMeshAnalyser = new MeshAnalyser(m_processedFundi);

    AffectNormalsToFundi();

    layer0Points->Delete();
    layer0Cells->Delete();
    lines->Delete();
    layer0->Delete();
    clean->Delete();

}

void MedialSurfaceComputer::BuildMedialSurface()
{
    for (int l = 1; l<m_nbLayers ; l++)
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

        FilterPointPosition(lpd);

        m_layers.push_back(lpd);

        delete mal0;
    }

    BuildMeshFromLayers();
    AnchorAllPointsToCandidates();

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

void MedialSurfaceComputer::FindCandidatePoints()
{
    m_candidatePoints = vtkPoints::New();

    double maxScore;
    vtkIdList* radiusPoints = vtkIdList::New();
    int nbInRadius;
    vtkIdType gi1, gi2;
    double point0[3];
    double radius = 1;
    double score1, score2;
    double res, cs;
    double vec[3];
    double dp1, dp2;
    vtkIdType ti1, ti2;
    double pg1[3], pg2[3];
    double n1[3], n2[3];
    double thr = -0.7;

    for(int i = 0; i< m_originalMesh->GetNumberOfPoints() ; i++)
    {
        m_originalMesh->GetPoint(i, point0);

        m_originalPointLocator->FindPointsWithinRadius(radius,point0,radiusPoints);

        maxScore = -1000000;

        nbInRadius = radiusPoints->GetNumberOfIds();

        gi1 = -1;
        gi2 = -1;

        for(int r1 = 0 ; r1<nbInRadius ; r1++)
        {
            ti1 = radiusPoints->GetId(r1);
            m_originalMesh->GetPoint(ti1,pg1);

            score1 = - m_metric->GetTuple1(radiusPoints->GetId(r1));

            m_originalMeshAnalyser->GetNormals()->GetTuple(ti1,n1);

            for(int r2 = 0 ; r2<nbInRadius ; r2++)
            {
                ti2 = radiusPoints->GetId(r2);
                m_originalMesh->GetPoint(ti2,pg2);

                m_originalMeshAnalyser->GetNormals()->GetTuple(ti2,n2);

                res = vtkMath::Dot(n1,n2);
                cs = vtkMath::Distance2BetweenPoints(pg1,pg2);

                vec[0] = pg2[0] - pg1[0];
                vec[1] = pg2[1] - pg1[1];
                vec[2] = pg2[2] - pg1[2];

                dp1 = vtkMath::Dot(n1,vec) / vtkMath::Norm(vec);
                dp2 = vtkMath::Dot(n2,vec) / vtkMath::Norm(vec);

                if(res < thr)
                {

                    score2 = - m_metric->GetTuple1(radiusPoints->GetId(r2)) - cs + dp1 - dp2;

                    if(score1 + score2 > maxScore)
                    {
                        maxScore = score1 + score2;
                        gi1 = ti1;
                        gi2 = ti2;

                    }
                }
            }
        }

        if(gi1 > -1 && gi2 > -1)
        {
            m_originalMesh->GetPoint(gi1,pg1);
            m_originalMesh->GetPoint(gi2,pg2);

            for(int k =0; k<3;k++)
            {
                point0[k] = (pg1[k] + pg2[k])/2.0;
            }

            m_candidatePoints->InsertNextPoint(point0);
        }
    }


    m_candidatePolyData = vtkPolyData::New();
    m_candidatePolyData->SetPoints(m_candidatePoints);
    m_candidatePolyData->Update();

    m_candidateLocator = vtkPointLocator::New();
    m_candidateLocator->SetDataSet(m_candidatePolyData);
    m_candidateLocator->BuildLocator();
    m_candidateLocator->Update();

}


void MedialSurfaceComputer::FilterPointPosition(vtkPolyData *mesh)
{
    vtkIdList* neib = vtkIdList::New();
    int nbNeib;
    double pointn[3], pn1[3];

    MeshAnalyser* ma = new MeshAnalyser(mesh);


    for(int s = 0; s<5 ; s++)
    {
        for(int i = 0; i<m_processedFundi->GetNumberOfPoints() ; i++)
        {
            ma->GetPointNeighbors(i,neib);
            nbNeib = neib->GetNumberOfIds();

            pointn[0] = 0;
            pointn[1] = 0;
            pointn[2] = 0;

            for(int ne = 0; ne < nbNeib; ne++)
            {
                mesh->GetPoint(neib->GetId(ne),pn1);

                pointn[0] += pn1[0]/nbNeib;
                pointn[1] += pn1[1]/nbNeib;
                pointn[2] += pn1[2]/nbNeib;

            }
            mesh->GetPoints()->SetPoint(i, pointn);

        }
    }
}

void MedialSurfaceComputer::AnchorPointToCandidate(double point[3])
{
    int cpi = m_candidateLocator->FindClosestPoint(point);

    double pointc[3];
    m_candidatePoints->GetPoint(cpi, pointc);

    for(int k = 0 ; k<3 ; k++)
    {
        point[k] = pointc[k];
    }

}

void MedialSurfaceComputer::AnchorAllPointsToCandidates()
{
    double point[3];

    for(int i = 0; i<m_medialSurface->GetNumberOfPoints() ; i++)
    {
        m_medialSurface->GetPoint(i,point);
        AnchorPointToCandidate(point);
        m_medialSurface->GetPoints()->SetPoint(i,point);
    }


}

void MedialSurfaceComputer::SelectNextPoint(double point[3], MeshAnalyser* mal0, int i)
{
    double point0[3];
    vtkIdList* neib = vtkIdList::New();
    double p0[3], p1[3], p2[3];
    double cp[3];
    double vec[3];
    double radius = 3;
    vtkIdList* radiusPoints = vtkIdList::New();
    int nbInRadius;
    double minScore;
    double d0, d1, d2;
    double score1;
    vtkIdType bi;

    mal0->GetMesh()->GetPoint(i,point0);

    m_pitsNorm->GetTuple(i,vec);

    //find neighbors of the parent point
    mal0->GetPointNeighbors(i,neib);
    int nbNeib = neib->GetNumberOfIds();
    int cnt = 0;
    for(int j = 0 ; j< nbNeib ; j++)
    {
        if(neib->GetId(j) != i)
        {
            if(cnt == 1)
            {
                mal0->GetMesh()->GetPoint(neib->GetId(j),p2);
                cnt++;
            }
            else if (cnt == 0)
            {
                mal0->GetMesh()->GetPoint(neib->GetId(j),p1);
                cnt++;
            }
        }
    }

    for(int k = 0 ; k<3 ; k++)
    {
        p0[k] = point0[k] + 2 * vec[k];
    }

    bi = m_candidateLocator->FindClosestPoint(p0);

    m_candidatePoints->GetPoint(bi, point);
    m_processedFundi->GetPoint(i,cp);


    for(int k = 0 ; k<3 ; k++)
    {
        vec[k] = point[k] - cp[k];
    }

    double nn = vtkMath::Norm(vec);


    for(int k = 0 ; k<3 ; k++)
    {
        vec[k] /= nn;
    }

    m_pitsNorm->SetTuple(i,vec);


//    m_candidateLocator->FindPointsWithinRadius(radius,p0,radiusPoints);

//    nbInRadius = radiusPoints->GetNumberOfIds();

//    minScore = 1000000;

//    for(int r1 = 0 ; r1<nbInRadius ; r1++)
//    {
//        m_candidatePoints->GetPoint(radiusPoints->GetId(r1), cp);

//        d0 = vtkMath::Distance2BetweenPoints(p0, cp);

//        score1 = d0;

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


//    minScore = vtkMath::Distance2BetweenPoints(p0, p1);

//    if(minScore < radius)
//    {
//        for(int k = 0; k<3;k++)
//        {
//            point[k] = p1[k];
//        }
//    }
//    else
//    {
//        for(int k = 0; k<3;k++)
//        {
//            point[k] = point0[k];
//        }

//    }



    neib->Delete();
    radiusPoints->Delete();

}


