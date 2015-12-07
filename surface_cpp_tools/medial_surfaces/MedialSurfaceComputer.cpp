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

#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkPointLocator.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkMath.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkDelaunay3D.h>
#include <vtkDelaunay2D.h>

#include <vtkUnstructuredGrid.h>

MedialSurfaceComputer::MedialSurfaceComputer(vtkPolyData *originalMesh, vtkPolyData *fundi)
{
    m_originalFundi = fundi;
    m_originalMesh = originalMesh;
    m_metric = originalMesh->GetPointData()->GetScalars();

    m_originalPointLocator = vtkPointLocator::New();
    m_originalPointLocator->SetDataSet(m_originalMesh);
    m_originalPointLocator->BuildLocator();
    m_originalPointLocator->Update();

    m_originalCellLocator = vtkOBBTree::New();
    m_originalCellLocator->SetDataSet(m_originalMesh);
    m_originalCellLocator->BuildLocator();
    m_originalCellLocator->Update();

    m_pitsNormals = vtkDoubleArray::New();
    m_pitsNormals->SetNumberOfComponents(3);

    m_nbLayers = 10;

    m_processedFundi = vtkPolyData::New();
    m_medialSurface = vtkPolyData::New();


    vtkPolyDataNormals* pdn = vtkPolyDataNormals::New();
//  VTK6 Update: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
//  Old: pdn->SetInput(m_originalMesh);
    pdn->SetInputData(m_originalMesh);

    pdn->Update();
    m_normals = pdn->GetOutput()->GetPointData()->GetNormals();

    BuildMedialSurface();
}



MedialSurfaceComputer::~MedialSurfaceComputer()
{
    m_pitsNormals->Delete();
}

void MedialSurfaceComputer::WriteIntoFile(char *fileName)
{
    vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
    writer->SetFileName(fileName);

//  VTK6 Update: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
//  Old: writer->SetInput(m_medialSurface);
    writer->SetInputData(m_medialSurface);

//  Redundant?:
//  writer->Update();
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

//  VTK6 Update: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
//  Old: clean->SetInput(layer0);
    clean->SetInputData(layer0);

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
    ProcessFundi();
    FindCandidatePoints();
    FrontVoronoi();

//    vtkDelaunay3D* del = vtkDelaunay3D::New();
//    del->SetInputData(m_candidatePolyData);
//    del->SetAlpha(2);
//    del->Update();

//    vtkPolyData* delpd = vtkPolyData::New();
//    delpd->SetPoints(del->GetOutput()->GetPoints());
//    delpd->SetPolys(del->GetOutput()->GetCells());
//    delpd->Update();

//    vtkDelaunay2D* del2 = vtkDelaunay2D::New();
////    del2->SetInputData(m_candidatePolyData);
//    del2->SetAlpha(2);
//    del2->SetSource(delpd);
//    del2->Update();


//    m_medialSurface->DeepCopy(del->GetOutput());

//    m_medialSurface->SetPoints(del2->GetOutput()->GetPoints());
//    m_medialSurface->SetPolys(del2->GetOutput()->GetCells());
////    m_medialSurface->GetPointData()->SetScalars(this->test);
//    m_medialSurface->Update();



//    vtkSphereSource* sphere =  vtkSphereSource::New();
//    sphere->SetPhiResolution(6);
//    sphere->SetThetaResolution(6);
//    sphere->SetRadius(0.5);
//    sphere->Update();

//    vtkGlyph3D* glyph = vtkGlyph3D::New();
//    glyph->SetInputData(m_candidatePolyData);
//    glyph->SetSource(sphere->GetOutput());
//    glyph->Update();

//    m_medialSurface->DeepCopy(glyph->GetOutput());

    for (int l = 1; l<m_nbLayers ; l++)
    {
        vtkPoints* layerPoints = vtkPoints::New();

        for (int i = 0; i< m_processedFundi->GetNumberOfPoints() ; i++)
        {
            double pointn[3];
            SelectNextPoint(pointn, m_layers.at(l-1) , i);
            layerPoints->InsertNextPoint(pointn);
        }

        vtkPolyData* lpd = vtkPolyData::New();
        lpd->SetPoints(layerPoints);
        lpd->SetLines(m_processedFundi->GetLines());
        lpd->Update();

//        FilterPointPosition(lpd);

//        vtkIndent indent;
//        lpd->PrintSelf(cout, indent);

        m_layers.push_back(lpd);
    }

    BuildMeshFromLayers();
//    AnchorAllPointsToCandidates();

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
        m_normals->GetTuple(i0, n);
        m_pitsNormals->InsertNextTuple(n);
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
        for (int i = 0; i< nbPits ; i++)
        {
            GetPointNeighbors(m_layers.at(l-1), i, neib);
            nbNeib = neib->GetNumberOfIds();

            medialPoints->InsertNextPoint(m_layers.at(l)->GetPoint(i));

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
    }

    m_medialSurface->SetPoints(medialPoints);
    m_medialSurface->SetPolys(medialCells);
    m_medialSurface->Update();

}

void MedialSurfaceComputer::FindCandidatePoints()
{
    m_candidatePoints = vtkPoints::New();

    bool goodProject;
    double interPoint[3];
    double point[3];
    double candidatePoint[3];
    double radius = 3;
    vtkIdList* radiusPoints = vtkIdList::New();
    int nbInRadius;


    for(int i = 0; i< m_originalMesh->GetNumberOfPoints() ; i++)
    {
        goodProject = ProjectInNormalDirection(i, interPoint);
        m_originalMesh->GetPoint(i,point);


        if(goodProject)
        {
            for(int k = 0 ; k < 3 ; k++)
            {
                candidatePoint[k] = (point[k] + interPoint[k]) / 2;
            }

//            m_originalPointLocator->FindPointsWithinRadius(radius,candidatePoint,radiusPoints);

//            nbInRadius = radiusPoints->GetNumberOfIds();

//            for(int k = 0 ; k < 3 ; k++)
//            {
//                candidatePoint[k] = 0;
//            }

//            for(int j = 0 ; j<nbInRadius ; j++)
//            {
//                m_originalMesh->GetPoint(radiusPoints->GetId(j), point);
//                for(int k = 0 ; k < 3 ; k++)
//                {
//                    candidatePoint[k] += point[k]/nbInRadius;
//                }


//            }


            m_candidatePoints->InsertNextPoint(candidatePoint);
        }


    }




//    double maxScore;
//    vtkIdList* radiusPoints = vtkIdList::New();
//    int nbInRadius;
//    vtkIdType gi1, gi2;
//    double point0[3];
//    double radius = 3;
//    double score1, score2;
//    double res, cs;
//    double vec[3];
//    double dp1, dp2;
//    vtkIdType ti1, ti2;
//    double pg1[3], pg2[3];
//    double n1[3], n2[3];
//    double thr = -0.7;

//    for(int i = 0; i< m_originalMesh->GetNumberOfPoints() ; i++)
//    {
//        m_originalMesh->GetPoint(i, point0);

//        m_originalPointLocator->FindPointsWithinRadius(radius,point0,radiusPoints);

//        maxScore = -1000000;

//        nbInRadius = radiusPoints->GetNumberOfIds();

//        gi1 = -1;
//        gi2 = -1;

//        for(int r1 = 0 ; r1<nbInRadius ; r1++)
//        {
//            ti1 = radiusPoints->GetId(r1);
//            m_originalMesh->GetPoint(ti1,pg1);

//            score1 = - m_metric->GetTuple1(radiusPoints->GetId(r1));

//            m_normals->GetTuple(ti1,n1);

//            for(int r2 = 0 ; r2<nbInRadius ; r2++)
//            {
//                ti2 = radiusPoints->GetId(r2);
//                m_originalMesh->GetPoint(ti2,pg2);

//                m_normals->GetTuple(ti2,n2);

//                res = vtkMath::Dot(n1,n2);
//                cs = vtkMath::Distance2BetweenPoints(pg1,pg2);

//                vec[0] = pg2[0] - pg1[0];
//                vec[1] = pg2[1] - pg1[1];
//                vec[2] = pg2[2] - pg1[2];

//                dp1 = vtkMath::Dot(n1,vec) / vtkMath::Norm(vec);
//                dp2 = vtkMath::Dot(n2,vec) / vtkMath::Norm(vec);

//                if(res < thr)
//                {

//                    score2 = - m_metric->GetTuple1(radiusPoints->GetId(r2)) - cs + dp1 - dp2;

//                    if(score1 + score2 > maxScore)
//                    {
//                        maxScore = score1 + score2;
//                        gi1 = ti1;
//                        gi2 = ti2;
//                    }
//                }
//            }
//        }

//        if(gi1 > -1 && gi2 > -1)
//        {
//            m_originalMesh->GetPoint(gi1,pg1);
//            m_originalMesh->GetPoint(gi2,pg2);

//            for(int k =0; k<3;k++)
//            {
//                point0[k] = (pg1[k] + pg2[k])/2.0;
//            }

//            m_candidatePoints->InsertNextPoint(point0);
//        }
//    }


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

    for(int s = 0; s<10 ; s++)
    {
        for(int i = 0; i<m_processedFundi->GetNumberOfPoints() ; i++)
        {
            GetPointNeighbors(mesh, i, neib);
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

void MedialSurfaceComputer::SelectNextPoint(double point[3], vtkPolyData* prevLayer, int i)
{
    double point0[3];
    vtkIdList* neib = vtkIdList::New();
    double p0[3], p1[3], p2[3], op[3];
    double cp[3];
    double vec[3];
    double nvec[3];
    double vec1[3];
    double radius = 3;
    vtkIdList* radiusPoints = vtkIdList::New();
    int nbInRadius;
    double maxScore = -10000;
    double d0, d1, d2;
    double score1;
    vtkIdType bi;

    m_processedFundi->GetPoint(i,op);

    prevLayer->GetPoint(i,point0);

    m_pitsNormals->GetTuple(i,vec);

    //find neighbors of the parent point
    GetPointNeighbors(prevLayer, i, neib);
    int nbNeib = neib->GetNumberOfIds();
    int cnt = 0;
//    for(int j = 0 ; j< nbNeib ; j++)
//    {
//        if(neib->GetId(j) != i)
//        {
//            if(cnt == 1)
//            {
//                prevLayer->GetPoint(neib->GetId(j),p2);
//                cnt++;
//            }
//            else if (cnt == 0)
//            {
//                prevLayer->GetPoint(neib->GetId(j),p1);
//                cnt++;
//            }
//        }
//    }

    for(int k = 0 ; k<3 ; k++)
    {
        p0[k] = point0[k] + 2 * vec[k];
    }

    vtkIdList* cpl = vtkIdList::New();

    m_candidateLocator->FindClosestNPoints(20,point0,cpl);

    bool cf = false;

    for(int j = 0; j< cpl->GetNumberOfIds() ; j++)
    {
//        cout<<i<<" "<<m_bins->GetId(cpl->GetId(j))<<endl;

        if(i == m_bins->GetId(cpl->GetId(j)))
        {
            cf = true;

            bi = cpl->GetId(j);
            m_candidatePoints->GetPoint(bi,p1);

            for(int k = 0 ; k<3 ; k++)
            {
                vec1[k] = p1[k] - point0[k];
            }

            d0 = vtkMath::Norm(vec1);

            for(int k = 0 ; k<3 ; k++)
            {
                vec1[k] /= d0;
            }

            score1 = vtkMath::Dot(vec1,vec);

            if(score1 > maxScore)
            {
                maxScore = score1;
                for(int k = 0 ; k<3 ; k++)
                {
                    point[k] = p1[k];
                    nvec[k] = p1[k] - op[k];
                }

                d1 = vtkMath::Norm(nvec);

                for(int k = 0 ; k<3 ; k++)
                {
                    nvec[k] /= d1;
                }

            }
        }
    }

    if(!cf)
    {
        for(int k = 0 ; k<3 ; k++)
        {
            point[k] = p0[k];
        }
    }

//    m_pitsNormals->SetTuple(i,nvec);



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


void MedialSurfaceComputer::GetPointNeighbors(vtkPolyData *pd, vtkIdType id, vtkIdList *neighbors)
{
    vtkIdList *cellids = vtkIdList::New();
    pd->GetPointCells(id,cellids);
    neighbors->Reset();
    for (int i = 0; i < cellids->GetNumberOfIds(); i++)
    {
        for(int j = 0; j < pd->GetCell(cellids->GetId(i))->GetNumberOfPoints(); j++)
        {
            if(pd->GetCell(cellids->GetId(i))->GetPointId(j)<pd->GetNumberOfPoints())
            {
                neighbors->InsertUniqueId(pd->GetCell(cellids->GetId(i))->GetPointId(j));
            }
        }
    }
    cellids->Delete();
}


bool MedialSurfaceComputer::ProjectInNormalDirection(vtkIdType pointId, double interPoint[3])
{
    double step = 5;

    bool inter = false;

    vtkPoints *points = vtkPoints::New();
    vtkIdList *cellIds = vtkIdList::New();

    double point1[3], point2[3], n[3];
    m_originalMesh->GetPoint(pointId, point1);
    m_normals->GetTuple(pointId, n);

    for(int k = 0 ; k < 3 ; k++)
    {
        point2[k] = point1[k] + step*n[k];

//        point1[k] += 0.1*n[k];
    }


    double t, ptline[3], pcoords[3];
    int subId;

//    int isInter = m_originalCellLocator->IntersectWithLine(point1, point2, 0.0001, t, ptline, pcoords, subId);
     int isInter = m_originalCellLocator->IntersectWithLine(point1, point2, points, cellIds);

    int nbInter = points->GetNumberOfPoints();
    double  minT = 100000;
    double pointi[3];


    for(int i = 0 ; i<nbInter ; i++)
    {
        points->GetPoint(i,pointi);
        t = vtkMath::Distance2BetweenPoints(point1, pointi);

        if(minT > t && t > 0.1)
        {
            minT = t;

            for(int k = 0; k<3 ; k++)
            {
//                interPoint[k] = ptline[k];//point1[k] + t * (point2[k] - point1[k]);
                interPoint[k] = pointi[k];
    //            cout<<interPoint[k]<<" ";
            }
    //        cout<<endl;
            inter = true;
        }
    }



    return inter;
}

void MedialSurfaceComputer::FrontVoronoi()
{
    double point[3];
    m_bins = vtkIdList::New();

    vtkPointLocator* pl = vtkPointLocator::New();
    pl->SetDataSet(m_processedFundi);
    pl->BuildLocator();
    pl->Update();

    for(int i = 0 ; i< m_candidatePoints->GetNumberOfPoints(); i++)
    {
        m_candidatePoints->GetPoint(i,point);

        m_bins->InsertNextId(pl->FindClosestPoint(point));
    }

}



