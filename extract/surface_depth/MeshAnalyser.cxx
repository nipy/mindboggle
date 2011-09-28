/* ********************************************************************
 * MeshAnalyser
 * Author: Joachim Giard
 * Date (last update): 2 Nov 2009
 * 
 * Class containing various methods to apply on 3D meshes. The 
 * meshanalyser object takes a vtkPolyData object as input
 * 
 * *******************************************************************/


#include "MeshAnalyser.h"
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkHull.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkCellArray.h>
#include <vtkPointSet.h>
#include <vtkDataSet.h>
#include <vtkPoints.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkQuadricDecimation.h>
#include <map>
#include <vtkPolyDataNormals.h>
#include <vtkDataArray.h>
#include <vtkBox.h>
#include <vtkBoundingBox.h>
#include <vtkObjectFactory.h>
#include <vtkSphereSource.h>
#include <algorithm>

MeshAnalyser::MeshAnalyser(char* fileName)
{
	this->mesh=vtkPolyData::New();
	this->simpl=vtkPolyData::New();
	this->depth=vtkDoubleArray::New();
	this->derDepth=vtkDoubleArray::New();
	this->pointSurf=vtkDoubleArray::New();
	this->pointSurfSimple=vtkDoubleArray::New();
	this->geoDistRing=vtkDoubleArray::New();
	this->geoDistArray=vtkDoubleArray::New();
	this->geoDistRingSimple=vtkDoubleArray::New();
	this->curv=vtkDoubleArray::New();
	this->gCurv=vtkDoubleArray::New();
	this->rough=vtkDoubleArray::New();
	this->mask=vtkDoubleArray::New();
	this->mu=vtkDoubleArray::New();
	this->theta=vtkDoubleArray::New();
	this->thetaxy=vtkDoubleArray::New();
	this->L=vtkDoubleArray::New();
	this->C=vtkDoubleArray::New();
	this->S=vtkDoubleArray::New();
	this->lmsdm=vtkDoubleArray::New();
	this->msdm=vtkDoubleArray::New();
	this->totalSurface=0;
	this->totalSurfaceSimple=0;
	this->inRing=vtkIdList::New();
	this->inRingSimple=vtkIdList::New();
	this->predGeo=vtkIdList::New();
	this->test=vtkDoubleArray::New();
	this->close=vtkIdList::New();
	this->voronoiBin=vtkIdList::New();
	this->VoronoiPatchesColors2==vtkDoubleArray::New();
	
	vtkPolyDataReader* reader=vtkPolyDataReader::New();
	reader->SetFileName(fileName);
	reader->Update();

	this->mesh->DeepCopy(reader->GetOutput());
	this->nbPoints=reader->GetOutput()->GetNumberOfPoints();
	reader->Delete();
	
	ComputePointSurface();
	ComputeNormals();
}

MeshAnalyser::MeshAnalyser(vtkPolyData* mesh)
{
	this->mesh=vtkPolyData::New();
	this->simpl=vtkPolyData::New();
	this->depth=vtkDoubleArray::New();
	this->derDepth=vtkDoubleArray::New();
	this->pointSurf=vtkDoubleArray::New();
	this->pointSurfSimple=vtkDoubleArray::New();
	this->geoDistRing=vtkDoubleArray::New();
	this->geoDistArray=vtkDoubleArray::New();
	this->geoDistRingSimple=vtkDoubleArray::New();
	this->curv=vtkDoubleArray::New();
	this->rough=vtkDoubleArray::New();
	this->mask=vtkDoubleArray::New();
	this->mu=vtkDoubleArray::New();
	this->theta=vtkDoubleArray::New();
	this->thetaxy=vtkDoubleArray::New();
	this->L=vtkDoubleArray::New();
	this->C=vtkDoubleArray::New();
	this->S=vtkDoubleArray::New();
	this->lmsdm=vtkDoubleArray::New();
	this->msdm=vtkDoubleArray::New();
	this->gCurv=vtkDoubleArray::New();
	this->totalSurface=0;
	this->totalSurfaceSimple=0;
	this->inRing=vtkIdList::New();
	this->inRingSimple=vtkIdList::New();
	this->predGeo=vtkIdList::New();
	this->test=vtkDoubleArray::New();
	this->close=vtkIdList::New();
	this->voronoiBin=vtkIdList::New();
	this->VoronoiPatchesColors2==vtkDoubleArray::New();
	
	this->mesh->DeepCopy(mesh);
	this->nbPoints=mesh->GetNumberOfPoints();
	

	ComputePointSurface();
	ComputeNormals();
}


MeshAnalyser::~MeshAnalyser()
{
	this->mesh->Delete();
	this->simpl->Delete();
	this->depth->Delete();
	this->derDepth->Delete();
	this->pointSurf->Delete();
	this->pointSurfSimple->Delete();
	this->geoDistRing->Delete();
	this->geoDistArray->Delete();
	this->geoDistRingSimple->Delete();
	this->curv->Delete();
	this->gCurv->Delete();
	this->rough->Delete();
	this->mask->Delete();
	this->mu->Delete();
	this->theta->Delete();
	this->thetaxy->Delete();
	this->L->Delete();
	this->C->Delete();
	this->S->Delete();
	this->lmsdm->Delete();
	this->msdm->Delete();
	this->inRing->Delete();
	this->inRingSimple->Delete();
	this->predGeo->Delete();
	this->test->Delete();
	this->close->Delete();
	this->voronoiBin->Delete();
	this->VoronoiPatchesColors2->Delete();
}

void MeshAnalyser::SetMesh(char* fileName)
{
	vtkPolyDataReader* reader=vtkPolyDataReader::New();
	reader->SetFileName(fileName);
	reader->Update();

	this->mesh->DeepCopy(reader->GetOutput());
	this->nbPoints=reader->GetOutput()->GetNumberOfPoints();
	reader->Delete();
}

void MeshAnalyser::SetMesh(vtkPolyData* mesh)
{
	this->mesh->DeepCopy(mesh);
	this->nbPoints=mesh->GetNumberOfPoints();
}

void MeshAnalyser::GetPointNeighborsSimple(vtkIdType id, vtkIdList* neighbors)
{
	vtkIdList *cellids = vtkIdList::New();
	this->simpl->GetPointCells(id,cellids);
	neighbors->Reset();
	for (int i = 0; i < cellids->GetNumberOfIds(); i++)
	{
		for(int j = 0; j < this->simpl->GetCell(cellids->GetId(i))->GetNumberOfPoints(); j++)
		{
			if(this->simpl->GetCell(cellids->GetId(i))->GetPointId(j)<this->nbPointsSimple)
			{
				neighbors->InsertUniqueId(this->simpl->GetCell(cellids->GetId(i))->GetPointId(j));
			}

	}
	cellids->Delete();
}
}

void MeshAnalyser::GetPointNeighbors(vtkIdType id, vtkIdList* neighbors)
{
	vtkIdList *cellids = vtkIdList::New();
	this->mesh->GetPointCells(id,cellids);

	neighbors->Reset();
	for (int i = 0; i < cellids->GetNumberOfIds(); i++)
	{
		for(int j = 0; j < this->mesh->GetCell(cellids->GetId(i))->GetNumberOfPoints(); j++)
		{
			if(this->mesh->GetCell(cellids->GetId(i))->GetPointId(j)<this->nbPoints)
			{
				neighbors->InsertUniqueId(this->mesh->GetCell(cellids->GetId(i))->GetPointId(j));
			}
		}
	}
	cellids->Delete();
}

void MeshAnalyser::GetPointNeighbors2(vtkIdType id, vtkIdList* neighbors, vtkPolyData* Mesh)
{
	vtkIdList *cellids = vtkIdList::New();
	Mesh->GetPointCells(id,cellids);
	neighbors->Reset();
	for (int i = 0; i < cellids->GetNumberOfIds(); i++)
	{
		for(int j = 0; j < Mesh->GetCell(cellids->GetId(i))->GetNumberOfPoints(); j++)
		{
			if(Mesh->GetCell(cellids->GetId(i))->GetPointId(j)<this->nbPoints)
			{
				neighbors->InsertUniqueId(Mesh->GetCell(cellids->GetId(i))->GetPointId(j));
			}
		}
	}
	cellids->Delete();
}

/*
void MeshAnalyser::GeoDistRing(vtkIdType st, double maxDist, double approx)
{
	//number of closest points in this->simpl to relate to each point of this->mesh	
	int N=3;
	//minimal ray to effectuate front propagation without passing through this->simpl (to be upgraded)
	double ray=0;
	double point1[3],point2[3], pointt[3],ec;
	vtkIdType cellId;
	int subId;
	double dist2;
	vtkIdList* pointIds = vtkIdList::New();
	
	//Ids and dist to points from this->simpl close to st
	double closeIds[N];
	double closeDists[N];

	//compute a decimated mesh if not already done
	if(this->approxFactor!=approx)
	{
		Simplify(approx);
		this->approxFactor=approx;
		
		if(VERBOSE) cout<<endl<<"Decimation done"<<endl;
		
		//vtkPointLocator* pl=vtkPointLocator::New();
		vtkCellLocator* pl=vtkCellLocator::New();
		pl->SetDataSet(this->simpl);
		pl->BuildLocator();
		
		this->allDist.clear();
		
		//temporary list of N closest points
		vtkIdList* Nclo = vtkIdList::New();
			
		if(VERBOSE) cout<<"computing neighborhood"<<endl;
		for(int i=0;i<this->nbPoints;i++)
		{
			this->mesh->GetPoint(i,point1);
			//pl->FindClosestNPoints(N, point1, Nclo);
			pl->FindClosestPoint(point1, pointt, cellId, subId, dist2);
        	this->simpl->GetCellPoints(cellId,pointIds);	
        	//N=pointIds->GetNumberOfIds();
        	this->test->InsertNextValue(pointIds->GetId(2));
			for(int j=0;j<N;j++)
			{
				//this->close->InsertNextId(Nclo->GetId(j));
				this->close->InsertNextId(pointIds->GetId(j));
			}
		}

		if(VERBOSE) cout<<"close computed"<<endl;
	}
	else if(this->close->GetNumberOfIds()==0)
	{

		if(VERBOSE) cout<<"aquí no hauria d'entrar"<<endl;
		//vtkPointLocator* pl=vtkPointLocator::New();
		vtkCellLocator* pl=vtkCellLocator::New();
		pl->SetDataSet(this->simpl);
		pl->BuildLocator();
	
		this->allDist.clear();
		
		//temporary list of N closest points
		vtkIdList* Nclo = vtkIdList::New();
			
		if(VERBOSE) cout<<"computing neighborhood"<<endl;
		for(int i=0;i<this->nbPoints;i++)
		{
			this->mesh->GetPoint(i,point1);
			//pl->FindClosestNPoints(N, point1, Nclo);
			pl->FindClosestPoint(point1, pointt, cellId, subId, dist2);
        	this->simpl->GetCellPoints(cellId,pointIds);	
        	//N=pointIds->GetNumberOfIds();	
			this->test->InsertNextValue(cellId);
			for(int j=0;j<N;j++)
			{
				//this->close->InsertNextId(Nclo->GetId(j));
				this->close->InsertNextId(pointIds->GetId(j));
			}
		}
	
	}

	if(VERBOSE) cout<<"ha sortit"<<endl;

	this->mesh->GetPoint(st,point1);

	for(int j=0;j<N;j++)
	{
		this->simpl->GetPoint(this->close->GetId(st*N+j),point2);
		ec=vtkMath::Distance2BetweenPoints(point1,point2);
		ec=sqrt(ec);
		closeIds[j]=this->close->GetId(st*N+j);
		closeDists[j]=ec;
		if(ec>ray)ray=ec;	
	}

	if(VERBOSE) cout<<"for acabat"<<endl;
	double md=0;
	
	//Compute all distances on this->simpl
	if(this->allDist.empty())
	{
		//if(VERBOSE) cout<<"computing distances on this->simpl"<<endl;
		for(int i=0;i<this->nbPointsSimple;i++)
		{
			GeoDistRingSimple(i,maxDist);
			for(int j=0;j<this->nbPointsSimple;j++)
			{
				this->allDist.push_back(this->geoDistRingSimple->GetValue(j));
			}		
		}
	}

	if(VERBOSE) cout<<"if acabat"<<endl;
	//if(VERBOSE) cout<<"distances on this->simpl computed"<<endl;
	

	int nbInRing;
	float minDist,sd,ed,nd;
	vtkIdType sn,en,curId1,curId2;	

	this->geoDistRing->Reset();

	//set all distances to max
	for(int i=0;i<this->nbPoints;i++)
	{
		this->geoDistRing->InsertNextValue(maxDist);
	}

	//if(VERBOSE) cout<<"browse and update"<<endl;
	//browse points in this->simpl and update close points of this->mesh if necessary
	for(int i=0;i<this->nbPoints;i++)
	{
	
		this->mesh->GetPoint(i,point1);
	
		for(int j=0;j<N;j++)
		{
			curId1=this->close->GetId(i*N+j);
			this->simpl->GetPoint(curId1,point2);	//closest neighbor of the arrival point
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);					//distance between arrivel point and closest simple point
					
			for(int k=0;k<N;k++)
			{
				curId2=closeIds[k];					//closest neighbor of the starting point
				nd=this->allDist[curId1*this->nbPointsSimple+curId2];
				//nd=10;
				if(nd<maxDist)
				{
					sd=ec+nd+closeDists[k];
					//sd=ec;
					//sd=closeDists[k];
					//sd=nd;
					if(this->geoDistRing->GetValue(i)>sd)
					{
						this->geoDistRing->SetValue(i,sd);
					}						
				}				
			}
		}
	}

	GeoDistRing(st,ray);
	nbInRing=this->inRing->GetNumberOfIds();
	
	for(int j=0;j<nbInRing;j++)
	{
		curId1=this->inRing->GetId(j);
		this->geoDistRing->SetValue(curId1,this->geoDistRing->GetValue(curId1));
	}
	this->geoDistRing->SetValue(st,0);
	//if(VERBOSE) cout<<"distance computed"<<endl;

	if(VERBOSE) cout<<"GeoDistRing finished"<<endl;
}
*/

void MeshAnalyser::GeoDistRing(vtkIdType st, double maxDist, double approx)
{
	vtkDoubleArray* tempDist=vtkDoubleArray::New();

	//number of closest points in this->simpl to relate to each point of this->mesh
	int N=3;
	//minimal ray to effectuate front propagation without passing through this->simpl (to be upgraded)
	double ray=0;
	double point1[3],point2[3], pointt[3],ec;
	vtkIdType cellId;
	int subId;
	double dist2;
	vtkIdList* pointIds = vtkIdList::New();

	//Ids and dist to points from this->simpl close to st
	double *closeIds = (double*)malloc(sizeof(double)*N);
	double *closeDists = (double*)malloc(sizeof(double)*N);

	//compute a decimated mesh if not already done
	if(this->approxFactor!=approx)
	{
		Simplify(approx);
		this->approxFactor=approx;

		if(1) cout<<endl<<"Decimation done"<<endl;

		//vtkPointLocator* pl=vtkPointLocator::New();
		vtkCellLocator* pl=vtkCellLocator::New();
		pl->SetDataSet(this->simpl);
		pl->BuildLocator();

		this->allDist.clear();

		//temporary list of N closest points
		vtkIdList* Nclo = vtkIdList::New();

		if(1) cout<<"computing neighborhood"<<endl;
		for(int i=0;i<this->nbPoints;i++)
		{
			this->mesh->GetPoint(i,point1);
			//pl->FindClosestNPoints(N, point1, Nclo);
			pl->FindClosestPoint(point1, pointt, cellId, subId, dist2);
        	this->simpl->GetCellPoints(cellId,pointIds);
        	//N=pointIds->GetNumberOfIds();
        	this->test->InsertNextValue(pointIds->GetId(2));
			for(int j=0;j<N;j++)
			{
				//this->close->InsertNextId(Nclo->GetId(j));
				this->close->InsertNextId(pointIds->GetId(j));
			}
		}
	}
	else if(this->close->GetNumberOfIds()==0)
	{
		//vtkPointLocator* pl=vtkPointLocator::New();
		vtkCellLocator* pl=vtkCellLocator::New();
		pl->SetDataSet(this->simpl);
		pl->BuildLocator();

		this->allDist.clear();

		//temporary list of N closest points
		vtkIdList* Nclo = vtkIdList::New();

		if(1) cout<<"computing neighborhood"<<endl;
		for(int i=0;i<this->nbPoints;i++)
		{
			this->mesh->GetPoint(i,point1);
			//pl->FindClosestNPoints(N, point1, Nclo);
			pl->FindClosestPoint(point1, pointt, cellId, subId, dist2);
        	this->simpl->GetCellPoints(cellId,pointIds);
        	//N=pointIds->GetNumberOfIds();
			this->test->InsertNextValue(cellId);
			for(int j=0;j<N;j++)
			{
				//this->close->InsertNextId(Nclo->GetId(j));
				this->close->InsertNextId(pointIds->GetId(j));
			}
		}

	}

	this->mesh->GetPoint(st,point1);

	for(int j=0;j<N;j++)
	{
		this->simpl->GetPoint(this->close->GetId(st*N+j),point2);
		ec=vtkMath::Distance2BetweenPoints(point1,point2);
		ec=sqrt(ec);
		closeIds[j]=this->close->GetId(st*N+j);
		closeDists[j]=ec;
		if(ec>ray)ray=ec;
	}

	double md=0;

	//Compute all distances on this->simpl
	if(this->allDist.empty())
	{
		//if(VERBOSE) cout<<"computing distances on this->simpl"<<endl;
		for(int i=0;i<this->nbPointsSimple;i++)
		{
			GeoDistRingSimple(i,maxDist);
			for(int j=0;j<this->nbPointsSimple;j++)
			{
				this->allDist.push_back(this->geoDistRingSimple->GetValue(j));
			}
		}
	}

	//if(VERBOSE) cout<<"distances on this->simpl computed"<<endl;


	int nbInRing;
	float minDist,sd,ed,nd;
	vtkIdType sn,en,curId1,curId2;

	this->geoDistRing->Reset();

	//set all distances to max
	for(int i=0;i<this->nbPoints;i++)
	{
		tempDist->InsertNextValue(maxDist);
	}

	//if(VERBOSE) cout<<"browse and update"<<endl;
	//browse points in this->simpl and update close points of this->mesh if necessary
	for(int i=0;i<this->nbPoints;i++)
	{

		this->mesh->GetPoint(i,point1);

		for(int j=0;j<N;j++)
		{
			curId1=this->close->GetId(i*N+j);
			this->simpl->GetPoint(curId1,point2);	//closest neighbor of the arrival point
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);					//distance between arrivel point and closest simple point

			for(int k=0;k<N;k++)
			{
				curId2=closeIds[k];					//closest neighbor of the starting point
				nd=this->allDist[curId1*this->nbPointsSimple+curId2];
				//nd=10;
				if(nd<maxDist)
				{
					sd=ec+nd+closeDists[k];
					//sd=ec;
					//sd=closeDists[k];
					//sd=nd;
					if(tempDist->GetValue(i)>sd)
					{
						tempDist->SetValue(i,sd);
					}
				}
			}
		}
	}

	GeoDistRing(st,ray);

	vtkPointLocator* pll=vtkPointLocator::New();
	    pll->SetDataSet(this->mesh);
	    pll->BuildLocator();

	    vtkIdList* inRay=vtkIdList::New();

	    pll->FindPointsWithinRadius(ray,point1,inRay);

	    nbInRing=inRay->GetNumberOfIds();

	    for(int j=0;j<nbInRing;j++)
	    {
	        this->mesh->GetPoint(inRay->GetId(j),point2);

	        ec=vtkMath::Distance2BetweenPoints(point1,point2);
	        ec=sqrt(ec);

	        if(ec>this->geoDistRing->GetValue(j)/2)this->geoDistRing->SetValue(j,ec);
	    }

	    pll->Delete();
	    inRay->Delete();


/*
	nbInRing=this->inRing->GetNumberOfIds();

	for(int j=0;j<this->nbPoints;j++)
	{
		if(this->geoDistRing->GetValue(j)>=ray)
		{
			this->geoDistRing->SetValue(j,tempDist->GetValue(j));
		}
	}
	this->geoDistRing->SetValue(st,0);
	tempDist->Delete();
	//if(VERBOSE) cout<<"distance computed"<<endl;

	 */
}

void MeshAnalyser::GeoDistRingSimple(vtkIdType stPoint, double maxDist)
{		
	multimap<float,vtkIdType> frontier;
		
	//temporary indexes
	vtkIdType curId;
	vtkIdType curFrontierId;
		
	vtkIdList* curNeib=vtkIdList::New();
	int nbNeib;

	this->geoDistRingSimple->Reset();
	this->geoDistRingSimple->Squeeze();
	
	this->inRingSimple->Reset();
	this->inRingSimple->Squeeze();

	//fill the distance vector with maximum value
	for(int i=0;i<this->nbPointsSimple;i++)
	{
		this->geoDistRingSimple->InsertNextValue(maxDist);
	}

	//in the front set, there is only the starting point
	frontier.insert(pair<float,vtkIdType>(0,stPoint));
	this->inRingSimple->InsertNextId(stPoint);
	int nbFrontierPoints=1;
	this->geoDistRingSimple->SetValue(stPoint,0);

	double point1[3], point2[3];
	double ec=0;
	double upgDist;

	multimap<float,vtkIdType>::iterator iter;

	//iteration screening the neighbors of the front point with the smallest value
	while(!frontier.empty())
	{
		iter=frontier.begin();
		curFrontierId=iter->second;
		
		this->simpl->GetPoint(curFrontierId,point1);
		curNeib->Reset();
		curNeib->Squeeze();
		GetPointNeighborsSimple(curFrontierId,curNeib);
		nbNeib=curNeib->GetNumberOfIds();
		//neighbors of the front point 
		for(int j=0;j<nbNeib;j++)
		{	
			curId=curNeib->GetId(j);		

			this->simpl->GetPoint(curId,point2);
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);
			upgDist=this->geoDistRingSimple->GetValue(curFrontierId)+ec;
			//distance is upgraded if necessary
			if(upgDist<this->geoDistRingSimple->GetValue(curId))
			{
				this->geoDistRingSimple->SetValue(curId,upgDist);
				//if point is upgraded, it enters in the front
				frontier.insert(pair<float,vtkIdType>(upgDist,curId));
				this->inRingSimple->InsertUniqueId(curId);
			}
		}
		frontier.erase(iter);
	}
	curNeib->Delete();
}


void MeshAnalyser::GeoDistRing(vtkIdType stPoint, double maxDist)
{

	this->geoCenter=stPoint;
		
	multimap<float,vtkIdType> frontier;
		
	//temporary indexes
	vtkIdType curId;
	vtkIdType curFrontierId;
		
	vtkIdList* curNeib=vtkIdList::New();
	int nbNeib;

	this->geoDistRing->Reset();
	this->geoDistRing->Squeeze();
	
	this->inRing->Reset();
	//this->inRing->Squeeze();
	
	this->predGeo->Reset();

	//fill the distance vector with maximum value
	for(int i=0;i<this->nbPoints;i++)
	{
		this->geoDistRing->InsertNextValue(maxDist);
		this->predGeo->InsertNextId(-1);
	}

	//in the front set, there is only the starting point
	frontier.insert(pair<float,vtkIdType>(0,stPoint));
	this->inRing->InsertNextId(stPoint);
	int nbFrontierPoints=1;
	this->geoDistRing->SetValue(stPoint,0);

	double point1[3], point2[3];
	double ec=0;
	double upgDist;

	multimap<float,vtkIdType>::iterator iter;

	//iteration screening the neighbors of the front point with the smallest value
	while(!frontier.empty())
	{
		iter=frontier.begin();
		curFrontierId=iter->second;
		
		this->mesh->GetPoint(curFrontierId,point1);
		curNeib->Reset();
		curNeib->Squeeze();
		GetPointNeighbors(curFrontierId,curNeib);
		nbNeib=curNeib->GetNumberOfIds();
		//neighbors of the front point 
		for(int j=0;j<nbNeib;j++)
		{	
			curId=curNeib->GetId(j);		

			this->mesh->GetPoint(curId,point2);
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);
			upgDist=this->geoDistRing->GetValue(curFrontierId)+ec;
			//distance is upgraded if necessary
			if(upgDist<this->geoDistRing->GetValue(curId))
			{
				this->geoDistRing->SetValue(curId,upgDist);
				//if point is upgraded, it enters in the front
				frontier.insert(pair<float,vtkIdType>(upgDist,curId));
				this->inRing->InsertUniqueId(curId);
				this->predGeo->SetId(curId,curFrontierId);
			}
		}
		frontier.erase(iter);
	}
	curNeib->Delete();
}

void MeshAnalyser::WriteIntoFile(char* fileName, char* prop)
{
	vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
	writer->SetFileName(fileName);
	
	if(strcmp("geoDist",prop)==0) this->mesh->GetPointData()->SetScalars(this->geoDistRing);
	//else if(strcmp("geoDistApp",prop)==0) this->mesh->GetPointData()->SetScalars(this->geoDistRing);
	else if(strcmp("depth",prop)==0) this->mesh->GetPointData()->SetScalars(this->depth);
	else if(strcmp("derDepth",prop)==0) this->mesh->GetPointData()->SetScalars(this->derDepth);
	else if(strcmp("curv",prop)==0) this->mesh->GetPointData()->SetScalars(this->curv);
	else if(strcmp("rough",prop)==0) this->mesh->GetPointData()->SetScalars(this->rough);
	else if(strcmp("gCurv",prop)==0) this->mesh->GetPointData()->SetScalars(this->gCurv);
	else if(strcmp("test",prop)==0) this->mesh->GetPointData()->SetScalars(this->test);
	else if(strcmp("surf",prop)==0) this->mesh->GetPointData()->SetScalars(this->pointSurf);
	else if(strcmp("mask",prop)==0) this->mesh->GetPointData()->SetScalars(this->mask);
	else if(strcmp("L",prop)==0) this->mesh->GetPointData()->SetScalars(this->L);
	else if(strcmp("C",prop)==0) this->mesh->GetPointData()->SetScalars(this->C);
	else if(strcmp("S",prop)==0) this->mesh->GetPointData()->SetScalars(this->S);
	else if(strcmp("lmsdm",prop)==0) this->mesh->GetPointData()->SetScalars(this->lmsdm);
	else if(strcmp("msdm",prop)==0) this->mesh->GetPointData()->SetScalars(this->msdm);
	else if(strcmp("geoDistArray",prop)==0) this->mesh->GetPointData()->SetScalars(this->geoDistArray);
	//convert the voronoi bin indices from Id to double
	else if(strcmp("voronoi",prop)==0) 
	{
		vtkDoubleArray* value=vtkDoubleArray::New();
		
		for(int i=0;i<this->nbPoints;i++)
		{
			value->InsertNextValue(this->voronoiBin->GetId(i));
		}
		this->mesh->GetPointData()->SetScalars(value);
		value->Delete();

	}
	//If no valid code is used, the index of the point is the scalar
	else 
	{
		vtkDoubleArray* noValue=vtkDoubleArray::New();
		
		for(int i=0;i<this->nbPoints;i++)
		{
			noValue->InsertNextValue(i);
		}
		this->mesh->GetPointData()->SetScalars(noValue);
		//noValue->Delete();
	}
	
	
	this->mesh->Update();
	if(strcmp("simple",prop)==0) writer->SetInput(this->simpl);
	else if(strcmp("geoDistSimple",prop)==0)
	{
		Simplify(500);
		GeoDistRingSimple(250,1000);
		this->simpl->GetPointData()->SetScalars(this->geoDistRingSimple);
		this->simpl->Update();
		writer->SetInput(this->simpl);
	}
	else writer->SetInput(this->mesh);
	writer->Update();
	writer->Write();
	writer->Delete();

}

void MeshAnalyser::WriteIntoFile(char* fileName)
{
	vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
	writer->SetFileName(fileName);
	
	writer->SetInput(this->mesh);
	writer->Update();
	writer->Write();
	writer->Delete();
}


void MeshAnalyser::WriteIntoFile(char* fileName, vtkDataArray* propExt)
{
	vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
	writer->SetFileName(fileName);
	
	this->mesh->GetPointData()->SetScalars(propExt);	
	this->mesh->Update();
	
	writer->SetInput(this->mesh);
	writer->Update();
	writer->Write();
	writer->Delete();

}

void MeshAnalyser::ComputePointSurface()
{	
	//initialisation
	vtkCellArray* cells=this->mesh->GetPolys();
	int nbPolys = cells->GetNumberOfCells();
	this->totalSurface=0;
	this->pointSurf->Reset();

	for(int i=0;i<this->nbPoints;i++)
	{
		this->pointSurf->InsertNextValue(0);
	}
	
	int cellIds[3];

	double pos1[3];
	double pos2[3];
	double pos3[3];
	
	float a;
	float b;
	float c;
	float p;
	
	float area;

	for (int i=0;i<nbPolys;i++)
	{
		cellIds[0]=this->mesh->GetCell(i)->GetPointId(0);
		cellIds[1]=this->mesh->GetCell(i)->GetPointId(1);
		cellIds[2]=this->mesh->GetCell(i)->GetPointId(2);
		
		if(cellIds[0]>this->nbPoints||cellIds[1]>this->nbPoints||cellIds[2]>this->nbPoints)
		{
			if(1) cout<<"point surface cell id problem: "<<i<<" "<<cellIds[0]<<" "<<cellIds[1]<<" "<<cellIds[2]<<endl;
		}
		else
		{
			this->mesh->GetPoint(cellIds[0],pos1);
			this->mesh->GetPoint(cellIds[1],pos2);
			this->mesh->GetPoint(cellIds[2],pos3);
			
			//distances between points of each triangle
			a=sqrt(vtkMath::Distance2BetweenPoints(pos1,pos2));
			b=sqrt(vtkMath::Distance2BetweenPoints(pos1,pos3));
			c=sqrt(vtkMath::Distance2BetweenPoints(pos2,pos3));
		
			//Herons's formula
			area=0.25*sqrt((a+b+c)*(b+c-a)*(a-b+c)*(a+b-c));
			
			this->totalSurface+=area;
			
			//add a third the face area to each point of the triangle	
			this->pointSurf->SetValue(cellIds[0],this->pointSurf->GetValue(cellIds[0])+area/3);
			this->pointSurf->SetValue(cellIds[1],this->pointSurf->GetValue(cellIds[1])+area/3);
			this->pointSurf->SetValue(cellIds[2],this->pointSurf->GetValue(cellIds[2])+area/3);
		}
	}
}

void MeshAnalyser::ComputePointSurfaceSimple()
{
	//initialisation
	vtkCellArray* cells=this->simpl->GetPolys();
	int nbPolys = cells->GetNumberOfCells();
	this->totalSurfaceSimple=0;
	this->pointSurfSimple->Reset();
	
	for(int i=0;i<this->nbPointsSimple;i++)
	{
		this->pointSurfSimple->InsertNextValue(0);
	}
	
	int cellIds[3];

	double pos1[3];
	double pos2[3];
	double pos3[3];
	
	float a;
	float b;
	float c;
	
	float area;

	for (int i=0;i<nbPolys;i++)
	{
		cellIds[0]=this->simpl->GetCell(i)->GetPointId(0);
		cellIds[1]=this->simpl->GetCell(i)->GetPointId(1);
		cellIds[2]=this->simpl->GetCell(i)->GetPointId(2);
		
		this->simpl->GetPoint(cellIds[0],pos1);
		this->simpl->GetPoint(cellIds[1],pos2);
		this->simpl->GetPoint(cellIds[2],pos3);
		
		//distances between points of each triangle
		a=sqrt(vtkMath::Distance2BetweenPoints(pos1,pos2));
		b=sqrt(vtkMath::Distance2BetweenPoints(pos1,pos3));
		c=sqrt(vtkMath::Distance2BetweenPoints(pos2,pos3));
		
		//Herons's formula	
		area=0.25*sqrt((a+b+c)*(b+c-a)*(a-b+c)*(a+b-c));
		
		this->totalSurfaceSimple+=area;
		
		//add a third the face area to each point of the triangle	
		this->pointSurfSimple->SetValue(cellIds[0],this->pointSurfSimple->GetValue(cellIds[0])+area/3);
		this->pointSurfSimple->SetValue(cellIds[1],this->pointSurfSimple->GetValue(cellIds[1])+area/3);
		this->pointSurfSimple->SetValue(cellIds[2],this->pointSurfSimple->GetValue(cellIds[2])+area/3);
	}
}


void MeshAnalyser::ComputeTravelDepth(bool norm)
{
	int maxBound=5000;

	//convex hull construction
	int recPlan=3;
	
	vtkHull *hull = vtkHull::New();
	hull->SetInput(this->mesh);
	hull->AddRecursiveSpherePlanes(recPlan);
	vtkPolyData *pq =hull->GetOutput();
	pq->Update();

	//Point locator for the closest point of the convex hull
	vtkCellLocator *pl = vtkCellLocator::New();
	pl->SetDataSet(pq);
	pl->BuildLocator();

	vtkDoubleArray* depths=vtkDoubleArray::New();

	vtkIntArray* label=vtkIntArray::New();
	vtkIntArray* labelT=vtkIntArray::New();
	
	for (int i=0;i<this->nbPoints;i++)
	{
		label->InsertNextValue(-1);
		labelT->InsertNextValue(-1);
	}

	vtkIdList* frontier=vtkIdList::New();
	double dist, point1[3], point2[3], epsilon = 1;
	int cellid, subid;

//-----------------------------------------------//
//direct allocation of the distance with the     //
//convex hull for the mesh                       //
//-----------------------------------------------// 

	double ec=0;
	vtkIdType closestPointId;

	hull->Delete();

//-----------------------------------------------//
//allocation of the distance to the              //
//mesh if the point is visible                   //
//-----------------------------------------------// 

	double threshold=0.30;
	double minDist;	
	double vector[3];

	vtkCellLocator* cellLocator=vtkCellLocator::New();
	cellLocator->SetDataSet(this->mesh);
	cellLocator->BuildLocator();
	
	double t, ptline[3], pcoords[3];
  	int subId;	
	int isInter;
	vtkPoints* nHPoints=vtkPoints::New();
	vtkIdList* nhList=vtkIdList::New();
	int doneLabel=7;
	int goodLabel=6;
	double tempSurf=0;
	
	for(int i = 0; i < this->nbPoints; i++)
	{
		this->mesh->GetPoint(i,point1);
		pl->FindClosestPoint(point1,point2,cellid,subid, ec);
		ec=sqrt(ec);
			
		//If the point is very close
		if (ec<threshold)
		{
			depths->InsertNextValue(ec);
			label->SetValue(i,doneLabel);
			nHPoints->InsertNextPoint(point1);
			nhList->InsertNextId(i);
			tempSurf+=this->pointSurf->GetValue(i);
		}
		else
		{
			//small shift, because the intersection is the starting point otherwise
			vector[0]=point2[0]-point1[0];
			vector[1]=point2[1]-point1[1];
			vector[2]=point2[2]-point1[2];

			point1[0]+=vector[0]*0.001;
			point1[1]+=vector[1]*0.001;
			point1[2]+=vector[2]*0.001;
			
			point2[0]+=vector[0]*0.1;
			point2[1]+=vector[1]*0.1;
			point2[2]+=vector[2]*0.1;

			t=1;
			
			isInter=cellLocator->IntersectWithLine(point1, point2, 0.0001, t, ptline, pcoords, subId);
			
			//if there is no other point a the surface between the surface and the coarse surface
			//if (t < 10e-5| t > 1-10e-5)
			if(isInter==0)
			{
				depths->InsertNextValue(ec);
				label->SetValue(i,doneLabel);
				nHPoints->InsertNextPoint(point1);
				nhList->InsertNextId(i);
				tempSurf+=this->pointSurf->GetValue(i);
			}
			else
			{
				frontier->InsertUniqueId(i);
				depths->InsertNextValue(maxBound);
			}
		
		}
	}

	pl->Delete();

	vtkPolyData* nHPointSet=vtkPolyData::New();
	nHPointSet->SetPoints(nHPoints);
		
	int nbNH=nHPoints->GetNumberOfPoints();
	
	int nbFrontierPoints=frontier->GetNumberOfIds();
	vtkIdList* curNeib=vtkIdList::New();
	int nbNeib;
	vtkIdType curId;
	vtkIdType curFrontierId;
	
	vtkIdList* oldFront=vtkIdList::New();	
	int nbOld;	
	double minNeibDist;	
	double curDist;	
	double n1[3],n2[3];
	
	vtkIdList* nhListr=vtkIdList::New();
	vtkIdList *result=vtkIdList::New();
	
//-----------------------------------------------//
//allocation of the distance to the              //
//mesh if the point is  not visible              //
//-----------------------------------------------// 
	
	int nbLeftPoints;	
	vtkIdList* frontierT=vtkIdList::New();	
	int curLab;
	int iteration=0;
	
	while(nbFrontierPoints>0&&nbNH>1)
	{
		iteration++;
		
		frontierT->DeepCopy(frontier);
		labelT->DeepCopy(label);
	
		nbLeftPoints=nbFrontierPoints;
	
		while(nbLeftPoints>0)
		{
			for(int i=0;i<nbLeftPoints;i++)
			{
				curFrontierId=frontierT->GetId(i);
				this->mesh->GetPoint(curFrontierId,point2);
				
				curNeib->Reset();
				GetPointNeighbors(curFrontierId,curNeib);
				nbNeib=curNeib->GetNumberOfIds();
				minNeibDist=depths->GetValue(curFrontierId);
				curLab=-1;
				for(int j=0;j<nbNeib;j++)
				{	
					curId=curNeib->GetId(j);
					
					if (depths->GetValue(curId)<maxBound)
					{
						this->mesh->GetPoint(curId,point1);
						ec=vtkMath::Distance2BetweenPoints(point1,point2);	
						ec=sqrt(ec);				
						curDist=depths->GetValue(curId)+ec;

						if(curDist<minNeibDist && label->GetValue(curId)-1>0)
						{ 	
							minNeibDist=curDist;
							curLab=label->GetValue(curId)-1;	
						}			
					}	
				}
				if (curLab<0)
				{
					oldFront->InsertUniqueId(curFrontierId);
					labelT->SetValue(curFrontierId,labelT->GetValue(curFrontierId)+1);
					if(labelT->GetValue(curFrontierId)>goodLabel)oldFront->InsertUniqueId(curFrontierId);
				}
				if(minNeibDist<maxBound)
				{
					depths->SetValue(curFrontierId,minNeibDist);
					label->SetValue(curFrontierId,curLab);
					labelT->SetValue(curFrontierId,curLab);
					oldFront->InsertUniqueId(curFrontierId);
				}
			}

			nbOld= oldFront->GetNumberOfIds();
			
			for(int i=0;i<nbOld;i++)
			{
				frontierT->DeleteId(oldFront->GetId(i));
			}
			
			oldFront->Reset();
			labelT->Reset();
			nbLeftPoints=frontierT->GetNumberOfIds();
		}
		
	
		vtkPointLocator* tcl=vtkPointLocator::New();
		tcl->SetDataSet(nHPointSet);
		tcl->BuildLocator();
		
		vtkPoints* tnh=vtkPoints::New();

		nhListr->Reset();

		tempSurf=0;
		for(int i=0;i<nbFrontierPoints;i++)
		{
			curFrontierId=frontier->GetId(i);
			this->mesh->GetPoint(curFrontierId,point1);
			
			tcl->FindClosestNPoints (2, point1,result);
			
			if(nhList->GetId(result->GetId(0))!=curFrontierId)
			{
				curId=nhList->GetId(result->GetId(0));
			}
			else curId=nhList->GetId(result->GetId(1));
			
			this->mesh->GetPoint(curId,point2);
					
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);
			
			this->normals->GetTuple(curFrontierId,n1);
			this->normals->GetTuple(curId,n2);
			
			point1[0]+=n1[0]*0.1*ec;
			point1[1]+=n1[1]*0.1*ec;
			point1[2]+=n1[2]*0.1*ec;
			
			point2[0]+=n2[0]*0.1*ec;
			point2[1]+=n2[1]*0.1*ec;
			point2[2]+=n2[2]*0.1*ec;
			
			isInter=cellLocator->IntersectWithLine(point1, point2, 0.0001, t, ptline, pcoords, subId);
			
			if(isInter==0)
			{
				
				if(ec+depths->GetValue(curId)<depths->GetValue(curFrontierId))
				{
					label->SetValue(curFrontierId,label->GetValue(curId)-1);
					depths->SetValue(curFrontierId,ec+depths->GetValue(curId));
				}
				else if(label->GetValue(curFrontierId)>0)\
				label->SetValue(curFrontierId,label->GetValue(curFrontierId)+1);
				if(label->GetValue(curFrontierId)>goodLabel)
				{
					//frontier->DeleteId(curFrontierId);
					tnh->InsertNextPoint(point1);
					nhListr->InsertNextId(curFrontierId);
					tempSurf+=this->pointSurf->GetValue(curFrontierId);
				}
			}
			else 
			{	
				if(label->GetValue(curFrontierId)>0)label->SetValue(curFrontierId,label->GetValue(curFrontierId)+1);
				if(label->GetValue(curFrontierId)>goodLabel)
				{
					tnh->InsertNextPoint(point1);
					nhListr->InsertNextId(curFrontierId);
					tempSurf+=this->pointSurf->GetValue(curFrontierId);
				}
			}
		}
		
		tcl->Delete();
		nHPointSet->Reset();
		nHPointSet->SetPoints(tnh);
		nbFrontierPoints=frontier->GetNumberOfIds();
		nbNH=tnh->GetNumberOfPoints();
		nhList->Reset();
		nhList->DeepCopy(nhListr);
		tnh->Delete();		
	}
	
	labelT->Delete();

	frontier->Reset();
	int nbNotDone=0;	
	
	for(int i=0;i<this->nbPoints;i++)
	{
		if(label->GetValue(i)<doneLabel)
		{
			frontier->InsertNextId(i);
			label->SetValue(i,-1);
		}
		if(depths->GetValue(i)==maxBound)nbNotDone++;
	}
	
	nbLeftPoints=frontier->GetNumberOfIds();
	
	while(nbLeftPoints>0)
	{
		for(int i=0;i<nbLeftPoints;i++)
		{
			curFrontierId=frontier->GetId(i);
			this->mesh->GetPoint(curFrontierId,point2);
			
			curNeib->Reset();
			GetPointNeighbors(curFrontierId,curNeib);
			nbNeib=curNeib->GetNumberOfIds();
			minNeibDist=depths->GetValue(curFrontierId);
			
			for(int j=0;j<nbNeib;j++)
			{	
				curId=curNeib->GetId(j);
				
				if (depths->GetValue(curId)<maxBound)
				{
					this->mesh->GetPoint(curId,point1);
					ec=vtkMath::Distance2BetweenPoints(point1,point2);
					
					ec=sqrt(ec);
					
					curDist=depths->GetValue(curId)+ec;

					if(curDist<minNeibDist)
					{ 	
						minNeibDist=curDist;	
					}
				}
			}
			
			if(minNeibDist<maxBound)
			{
				if(minNeibDist==depths->GetValue(curFrontierId))
				{	
					label->SetValue(curFrontierId,label->GetValue(curFrontierId)+1);
				}
				else
				{
					depths->SetValue(curFrontierId,minNeibDist);
					label->SetValue(curFrontierId,label->GetValue(curId)-1);
				}
			}
			
			if(label->GetValue(curFrontierId)>goodLabel)
			{
				oldFront->InsertUniqueId(curFrontierId);
			}
		}

		nbOld= oldFront->GetNumberOfIds();
		
		for(int i=0;i<nbOld;i++)
		{
			frontier->DeleteId(oldFront->GetId(i));
		}
		
		oldFront->Reset();

		nbLeftPoints=frontier->GetNumberOfIds();
		//if(VERBOSE) cout<<"nb left points: "<<nbLeftPoints<<endl;
	}

	label->Delete();
	frontier->Delete();

	double MIN=1000;
	double MAX=-1;

	if(norm==true)
	{
		for(int i = 0; i < this->nbPoints; i++)
		{
			if(MAX<depths->GetValue(i)&depths->GetValue(i)<maxBound)MAX=depths->GetValue(i);
			if(MIN>depths->GetValue(i))MIN=depths->GetValue(i);
		}
	}

	double cc;
	double tot=0;
	for(int i = 0; i < this->nbPoints; i++)
	{
		if(norm==true)cc = (depths->GetValue(i)-MIN)/(MAX-MIN);
		else cc=depths->GetValue(i);
		this->depth->InsertNextValue(cc);
		tot+=cc*this->pointSurf->GetValue(i);
	}
	
	depths->Delete();
}

void MeshAnalyser::ComputeTravelDepth(bool norm, bool sphere)
{

	double pt1[3];
	double pt2[3];
	double xMin=0.0,yMin=0.0,zMin=0.0,xMax=0.0,yMax=0.0,zMax=0.0;
	double diag = 0.0;

	for(int i=0;i<nbPoints;i++)
	{
			this->mesh->GetPoint(i,pt1);

			if(xMin>pt1[0])xMin=pt1[0];
			if(yMin>pt1[1])yMin=pt1[1];
			if(zMin>pt1[2])zMin=pt1[2];
			if(xMax<pt1[0])xMax=pt1[0];
			if(yMax<pt1[1])yMax=pt1[1];
			if(zMax<pt1[2])zMax=pt1[2];

	}

	this->diag = sqrt( (xMax-xMin)*(xMax-xMin)+(yMax-yMin)*(yMax-yMin)+(zMax-zMin)*(zMax-zMin) );

	if(1) cout<<"diag: "<<this->diag<<endl;





    int maxBound=5000000;

    //convex hull construction
    int recPlan=3;
    vtkPolyData *pq;

    double center[3];
    double dist, point1[3], point2[3], epsilon = 1;
    this->mesh->GetCenter(center);
    double ec=0, maxDist=-1;

    if(sphere)
    {
        for(int i=0;i<this->nbPoints;i++)
        {
            this->mesh->GetPoint(i,point1);
            ec=vtkMath::Distance2BetweenPoints(center,point1);
            ec=sqrt(ec);
            if(maxDist<ec)maxDist=ec;
        }

        vtkIdList* nhListr=vtkIdList::New();

        vtkSphereSource* ss=vtkSphereSource::New();
        ss->SetCenter(center);
        ss->SetRadius(maxDist);
        ss->SetThetaResolution(16);
        ss->SetPhiResolution(16);
        ss->Update();

        pq=ss->GetOutput();
        pq->Update();

        vtkPolyDataWriter* w=vtkPolyDataWriter::New();
        w->SetInput(ss->GetOutput());
        w->SetFileName((char*)"sphere.vtk");
        w->Write();
        w->Update();
        w->Delete();
    }

    else
    {
        vtkHull *hull = vtkHull::New();
        hull->SetInput(this->mesh);
        hull->AddRecursiveSpherePlanes(recPlan);
        hull->Update();

        pq=hull->GetOutput();
        pq->Update();
        hull->Delete();

    }


    //Point locator for the closest point of the convex hull
    vtkCellLocator *pl = vtkCellLocator::New();
    pl->SetDataSet(pq);
    pl->BuildLocator();

   vtkDoubleArray* depths=vtkDoubleArray::New();

	vtkIntArray* label=vtkIntArray::New();
	vtkIntArray* labelT=vtkIntArray::New();



	for (int i=0;i<this->nbPoints;i++)
	{
		label->InsertNextValue(-1);
		labelT->InsertNextValue(-1);
	}

	vtkIdList* frontier=vtkIdList::New();

	int cellid, subid;

//-----------------------------------------------//
//direct allocation of the distance with the     //
//convex hull for the mesh                       //
//-----------------------------------------------//


	vtkIdType closestPointId;




//-----------------------------------------------//
//allocation of the distance to the              //
//mesh if the point is visible                   //
//-----------------------------------------------//

	double threshold=0.30;
	double minDist;
	double vector[3];

	vtkCellLocator* cellLocator=vtkCellLocator::New();
	cellLocator->SetDataSet(this->mesh);
	cellLocator->BuildLocator();

	double t, ptline[3], pcoords[3];
  	int subId;
	int isInter;
	vtkPoints* nHPoints=vtkPoints::New();
	vtkIdList* nhList=vtkIdList::New();
	int doneLabel=7;
	int goodLabel=6;
	double tempSurf=0;

	for(int i = 0; i < this->nbPoints; i++)
	{
		this->mesh->GetPoint(i,point1);
		pl->FindClosestPoint(point1,point2,cellid,subid, ec);
		ec=sqrt(ec);

		//If the point is very close
		if (ec<threshold)
		{
			depths->InsertNextValue(ec);
			label->SetValue(i,doneLabel);
			nHPoints->InsertNextPoint(point1);
			nhList->InsertNextId(i);
			tempSurf+=this->pointSurf->GetValue(i);
		}
		else
		{
			//small shift, because the intersection is the starting point otherwise
			vector[0]=point2[0]-point1[0];
			vector[1]=point2[1]-point1[1];
			vector[2]=point2[2]-point1[2];

			point1[0]+=vector[0]*0.001;
			point1[1]+=vector[1]*0.001;
			point1[2]+=vector[2]*0.001;

			point2[0]+=vector[0]*0.1;
			point2[1]+=vector[1]*0.1;
			point2[2]+=vector[2]*0.1;

			t=1;

			isInter=cellLocator->IntersectWithLine(point1, point2, 0.0001, t, ptline, pcoords, subId);

			//if there is no other point a the surface between the surface and the coarse surface
			//if (t < 10e-5| t > 1-10e-5)
			if(isInter==0)
			{
				depths->InsertNextValue(ec);
				label->SetValue(i,doneLabel);
				nHPoints->InsertNextPoint(point1);
				nhList->InsertNextId(i);
				tempSurf+=this->pointSurf->GetValue(i);
			}
			else
			{
				frontier->InsertUniqueId(i);
				depths->InsertNextValue(maxBound);
			}

		}
	}

	pl->Delete();

	vtkPolyData* nHPointSet=vtkPolyData::New();
	nHPointSet->SetPoints(nHPoints);

	int nbNH=nHPoints->GetNumberOfPoints();

	int nbFrontierPoints=frontier->GetNumberOfIds();
	vtkIdList* curNeib=vtkIdList::New();
	int nbNeib;
	vtkIdType curId;
	vtkIdType curFrontierId;

	vtkIdList* oldFront=vtkIdList::New();
	int nbOld;
	double minNeibDist;
	double curDist;
	double n1[3],n2[3];

	vtkIdList* nhListr=vtkIdList::New();
	vtkIdList *result=vtkIdList::New();

//-----------------------------------------------//
//allocation of the distance to the              //
//mesh if the point is  not visible              //
//-----------------------------------------------//

	int nbLeftPoints;
	vtkIdList* frontierT=vtkIdList::New();
	int curLab;
	int iteration=0;

	while(nbFrontierPoints>0&&nbNH>1)
	{
		iteration++;

		frontierT->DeepCopy(frontier);
		labelT->DeepCopy(label);

		nbLeftPoints=nbFrontierPoints;

		while(nbLeftPoints>0)
		{
			for(int i=0;i<nbLeftPoints;i++)
			{
				curFrontierId=frontierT->GetId(i);
				this->mesh->GetPoint(curFrontierId,point2);

				curNeib->Reset();
				GetPointNeighbors(curFrontierId,curNeib);
				nbNeib=curNeib->GetNumberOfIds();
				minNeibDist=depths->GetValue(curFrontierId);
				curLab=-1;
				for(int j=0;j<nbNeib;j++)
				{
					curId=curNeib->GetId(j);

					if (depths->GetValue(curId)<maxBound)
					{
						this->mesh->GetPoint(curId,point1);
						ec=vtkMath::Distance2BetweenPoints(point1,point2);
						ec=sqrt(ec);
						curDist=depths->GetValue(curId)+ec;

						if(curDist<minNeibDist && label->GetValue(curId)-1>0)
						{
							minNeibDist=curDist;
							curLab=label->GetValue(curId)-1;
						}
					}
				}
				if (curLab<0)
				{
					oldFront->InsertUniqueId(curFrontierId);
					labelT->SetValue(curFrontierId,labelT->GetValue(curFrontierId)+1);
					if(labelT->GetValue(curFrontierId)>goodLabel)oldFront->InsertUniqueId(curFrontierId);
				}
				if(minNeibDist<maxBound)
				{
					depths->SetValue(curFrontierId,minNeibDist);
					label->SetValue(curFrontierId,curLab);
					labelT->SetValue(curFrontierId,curLab);
					oldFront->InsertUniqueId(curFrontierId);
				}
			}

			nbOld= oldFront->GetNumberOfIds();

			for(int i=0;i<nbOld;i++)
			{
				frontierT->DeleteId(oldFront->GetId(i));
			}

			oldFront->Reset();
			labelT->Reset();
			nbLeftPoints=frontierT->GetNumberOfIds();
		}


		vtkPointLocator* tcl=vtkPointLocator::New();
		tcl->SetDataSet(nHPointSet);
		tcl->BuildLocator();

		vtkPoints* tnh=vtkPoints::New();

		nhListr->Reset();

		tempSurf=0;
		for(int i=0;i<nbFrontierPoints;i++)
		{
			curFrontierId=frontier->GetId(i);
			this->mesh->GetPoint(curFrontierId,point1);

			tcl->FindClosestNPoints (2, point1,result);

			if(nhList->GetId(result->GetId(0))!=curFrontierId)
			{
				curId=nhList->GetId(result->GetId(0));
			}
			else curId=nhList->GetId(result->GetId(1));

			this->mesh->GetPoint(curId,point2);

			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);

			this->normals->GetTuple(curFrontierId,n1);
			this->normals->GetTuple(curId,n2);

			point1[0]+=n1[0]*0.1*ec;
			point1[1]+=n1[1]*0.1*ec;
			point1[2]+=n1[2]*0.1*ec;

			point2[0]+=n2[0]*0.1*ec;
			point2[1]+=n2[1]*0.1*ec;
			point2[2]+=n2[2]*0.1*ec;

			isInter=cellLocator->IntersectWithLine(point1, point2, 0.0001, t, ptline, pcoords, subId);

			if(isInter==0)
			{

				if(ec+depths->GetValue(curId)<depths->GetValue(curFrontierId))
				{
					label->SetValue(curFrontierId,label->GetValue(curId)-1);
					depths->SetValue(curFrontierId,ec+depths->GetValue(curId));
				}
				else if(label->GetValue(curFrontierId)>0)\
				label->SetValue(curFrontierId,label->GetValue(curFrontierId)+1);
				if(label->GetValue(curFrontierId)>goodLabel)
				{
					//frontier->DeleteId(curFrontierId);
					tnh->InsertNextPoint(point1);
					nhListr->InsertNextId(curFrontierId);
					tempSurf+=this->pointSurf->GetValue(curFrontierId);
				}
			}
			else
			{
				if(label->GetValue(curFrontierId)>0)label->SetValue(curFrontierId,label->GetValue(curFrontierId)+1);
				if(label->GetValue(curFrontierId)>goodLabel)
				{
					tnh->InsertNextPoint(point1);
					nhListr->InsertNextId(curFrontierId);
					tempSurf+=this->pointSurf->GetValue(curFrontierId);
				}
			}
		}

		tcl->Delete();
		nHPointSet->Reset();
		nHPointSet->SetPoints(tnh);
		nbFrontierPoints=frontier->GetNumberOfIds();
		nbNH=tnh->GetNumberOfPoints();
		nhList->Reset();
		nhList->DeepCopy(nhListr);
		tnh->Delete();
	}

	labelT->Delete();

	frontier->Reset();
	int nbNotDone=0;

	for(int i=0;i<this->nbPoints;i++)
	{
		if(label->GetValue(i)<doneLabel)
		{
			frontier->InsertNextId(i);
			label->SetValue(i,-1);
		}
		if(depths->GetValue(i)==maxBound)nbNotDone++;
	}

	nbLeftPoints=frontier->GetNumberOfIds();

	while(nbLeftPoints>0)
	{
		for(int i=0;i<nbLeftPoints;i++)
		{
			curFrontierId=frontier->GetId(i);
			this->mesh->GetPoint(curFrontierId,point2);

			curNeib->Reset();
			GetPointNeighbors(curFrontierId,curNeib);
			nbNeib=curNeib->GetNumberOfIds();
			minNeibDist=depths->GetValue(curFrontierId);

			for(int j=0;j<nbNeib;j++)
			{
				curId=curNeib->GetId(j);

				if (depths->GetValue(curId)<maxBound)
				{
					this->mesh->GetPoint(curId,point1);
					ec=vtkMath::Distance2BetweenPoints(point1,point2);

					ec=sqrt(ec);

					curDist=depths->GetValue(curId)+ec;

					if(curDist<minNeibDist)
					{
						minNeibDist=curDist;
					}
				}
			}

			if(minNeibDist<maxBound)
			{
				if(minNeibDist==depths->GetValue(curFrontierId))
				{
					label->SetValue(curFrontierId,label->GetValue(curFrontierId)+1);
				}
				else
				{
					depths->SetValue(curFrontierId,minNeibDist);
					label->SetValue(curFrontierId,label->GetValue(curId)-1);
				}
			}

			if(label->GetValue(curFrontierId)>goodLabel)
			{
				oldFront->InsertUniqueId(curFrontierId);
			}
		}

		nbOld= oldFront->GetNumberOfIds();

		for(int i=0;i<nbOld;i++)
		{
			frontier->DeleteId(oldFront->GetId(i));
		}

		oldFront->Reset();

		nbLeftPoints=frontier->GetNumberOfIds();
		//if(VERBOSE) cout<<"nb left points: "<<nbLeftPoints<<endl;
	}

	label->Delete();
	frontier->Delete();

	double MIN=1000;
	double MAX=-1;

	if(norm==true)
	{
		for(int i = 0; i < this->nbPoints; i++)
		{
			if(MAX<depths->GetValue(i)&depths->GetValue(i)<maxBound)MAX=depths->GetValue(i);
			if(MIN>depths->GetValue(i))MIN=depths->GetValue(i);
		}
	}

	double cc;
	double tot=0;
	for(int i = 0; i < this->nbPoints; i++)
	{
		if(norm==true)cc = (depths->GetValue(i)-MIN)/(MAX-MIN);
		else cc=depths->GetValue(i);
		this->depth->InsertNextValue(cc);
		tot+=cc*this->pointSurf->GetValue(i);
	}

	depths->Delete();




}



void MeshAnalyser::Simplify(double factor)
{
	//if the last simplification was not exactly the same
	if(this->approxFactor!=factor)
	{
		this->nbPointsSimple=10000000;
		this->simpl->DeepCopy(this->mesh);
		
		if(factor<=1)
		{
			vtkQuadricDecimation* dec=vtkQuadricDecimation::New();
			
			dec->SetInput(this->simpl);
			dec->SetTargetReduction(factor);
			dec->Update();

			this->simpl=dec->GetOutput();
			this->nbPointsSimple=this->simpl->GetNumberOfPoints();
		}

		else
		{
			//if the decimation factor is a number of points, 
			//the decimation is done until this number of points is reached
			while(this->nbPointsSimple>factor)
			{
				vtkQuadricDecimation* dec=vtkQuadricDecimation::New();
				dec->SetInput(this->simpl);
				dec->SetTargetReduction(0.5);
				dec->Update();

				this->simpl=dec->GetOutput();
				if(this->simpl->GetNumberOfPoints()==this->nbPointsSimple)break;
				this->nbPointsSimple=this->simpl->GetNumberOfPoints();
				if(1) cout<<this->nbPointsSimple<<" ";
			}
		}
		this->approxFactor=factor;
	}
}

void MeshAnalyser::ComputeNormals()
{
	vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
	pdn->SetInput(this->mesh);
	pdn->SetFeatureAngle(90);
	pdn->Update();

	vtkPolyData* no=pdn->GetOutput();
	no->Update();
	this->normals=no->GetPointData()->GetNormals();
}

void MeshAnalyser::ComputeBothCurvatures(double ray)
{
	//initialisation
	vtkPolyData* upNorm=vtkPolyData::New();
	vtkPoints* upNormPoints=vtkPoints::New();
	
	double point1[3], point2[3], norm[3];
	double eps=0.01;
	
	//computation of the mesh modification by projection of the points
	//in the direction of the normal.
	for(int i=0;i<this->nbPoints;i++)
	{
		this->mesh->GetPoint(i,point1);
		this->normals->GetTuple(i,norm);
		
		point2[0]=point1[0]+norm[0]*eps;
		point2[1]=point1[1]+norm[1]*eps;
		point2[2]=point1[2]+norm[2]*eps;
		
		upNormPoints->InsertNextPoint(point2);
	}	
	
	upNorm->DeepCopy(this->mesh);
	upNorm->SetPoints(upNormPoints);
	upNorm->Update();
	
	upNormPoints->Delete();
	
	//computation of the laplacian smoothed mesh
	vtkSmoothPolyDataFilter* smoothed = vtkSmoothPolyDataFilter::New();
	smoothed->SetInput(this->mesh);
	smoothed->SetRelaxationFactor(0.7);
	smoothed->SetNumberOfIterations(200);
	smoothed->FeatureEdgeSmoothingOff();
	smoothed->Update();
		
	int nbInRing;

	//using meshanalyser to compute the areas affected to each point
	MeshAnalyser* mat=new MeshAnalyser(upNorm);
	
	vtkDoubleArray* upPs=mat->GetPointSurface();
	
	upNorm->Delete();
	
	MeshAnalyser* mas=new MeshAnalyser(smoothed->GetOutput());
	
	vtkDoubleArray* sPs=mas->GetPointSurface();
	
	smoothed->Delete();
	
	double surf,normSurf, siSurf;
	
	//point locator for the smoothing of the curvature field
	vtkPointLocator* pl=vtkPointLocator::New();
	pl->SetDataSet(this->mesh);
	pl->BuildLocator();
	vtkIdList* Nclo = vtkIdList::New();	
		
	vtkIdType curId;
	
	double maxCurv=-10;
	double minCurv=1000;
	double maxgCurv=-10;
	double mingCurv=1000;
	
	for(int i=0;i<this->nbPoints;i++)
	{
		this->mesh->GetPoint(i,point1);
		pl->FindPointsWithinRadius(ray,point1,Nclo);
		//GeoDistRing(i,ray);
		//nbInRing=this->inRing->GetNumberOfIds();
		nbInRing=Nclo->GetNumberOfIds();
		
		surf=0;
		normSurf=0;
		siSurf=0;
		
		for(int j=0;j<nbInRing;j++)
		{
			//curId=this->inRing->GetId(j);
			curId=Nclo->GetId(j);
			
			surf+=this->pointSurf->GetValue(curId);
			normSurf+=upPs->GetValue(curId);
			siSurf+=sPs->GetValue(curId);
		}
		if(normSurf/surf<minCurv)minCurv=normSurf/surf;
		if(normSurf/surf>maxCurv)maxCurv=normSurf/surf;
		if(siSurf/surf<mingCurv)mingCurv=siSurf/surf;
		if(siSurf/surf>maxgCurv)maxgCurv=siSurf/surf;	
		
		this->curv->InsertNextValue(normSurf/surf);
		this->gCurv->InsertNextValue(siSurf/surf);
		
	}
	
	delete mat;
	delete mas;
	
	double curCurv;
	
	for(int i=0;i<this->nbPoints;i++)
	{	
		curCurv=this->curv->GetValue(i);
		this->curv->SetValue(i,(maxCurv-curCurv)/(maxCurv-minCurv)*(-2)+1);
		curCurv=this->gCurv->GetValue(i);
		this->gCurv->SetValue(i,(maxgCurv-curCurv)/(maxgCurv-mingCurv)*2-1);
	}
	
}


void MeshAnalyser::ComputeCurvature(double res)
{
	//if(1) cout<<"Dins ComputeCurvature"<<endl;

	double pt1[3],pt2[3],N[3]={0,0,0};

	double j;

	vtkSmoothPolyDataFilter* smoothed = vtkSmoothPolyDataFilter::New();
	smoothed->SetInput(this->mesh);
	smoothed->SetRelaxationFactor(res);
	smoothed->SetNumberOfIterations(200);
	smoothed->FeatureEdgeSmoothingOff();
	smoothed->Update();

	this->curv->Reset();

	double pt1pt2[3];
	for( int i = 0; i < this->nbPoints; i++)
	{
		this->mesh->GetPoint(i,pt1);
		this->normals->GetTuple(i,N);
		smoothed->GetOutput()->GetPoint(i,pt2);
		for(int w=0; w<3;w++)
		{
			pt1pt2[w]=pt2[w]-pt1[w];
		}

		if(vtkMath::Norm(pt1pt2)==0) this->curv->InsertNextTuple1(0);

	//	j=((1+vtkMath::Dot(N,pt1pt2)/(vtkMath::Norm(N)*vtkMath::Norm(pt1pt2)))/2.0);

		//if(//isnan(j)) if(VERBOSE) cout<<"j is nan "<<vtkMath::Dot(N,pt1pt2)<<" "<<vtkMath::Norm(N)<<" "<<vtkMath::Norm(pt1pt2)<<endl;

		//else this->curv->InsertNextTuple1((1+vtkMath::Dot(N,pt1pt2)/(vtkMath::Norm(N)*vtkMath::Norm(pt1pt2)))/2.0);
		else this->curv->InsertNextTuple1((1+vtkMath::Dot(N,pt1pt2)/(vtkMath::Norm(pt1pt2)))/2.0);
		//else this->curv->InsertNextTuple1((1+vtkMath::Dot(N,pt1pt2)/(vtkMath::Norm(N)))/2.0);
	}
	
	smoothed->Delete();

	//if(1) cout<<"ComputeCurvature done"<<endl;

}

void MeshAnalyser::ComputeRoughness()
{
	//if(1) cout<<"Dins ComputeRoughness"<<endl;

    double maxDist=0,ec=0;
    double  nmaxDist=0;
    //GeoDistRing(0,1000);
    double point1[3],point2[3];

    int nbInRing;
    double std=0,ringSurf=0,mean=0;
    vtkIdType crp;
    this->mesh->Update();
    this->nbPoints = this->mesh->GetNumberOfPoints();
    for(int i=0;i<this->nbPoints;i++)
    {
        vtkIdList* curNeib=vtkIdList::New();
    	GetPointNeighbors(i,curNeib);
        nbInRing=curNeib->GetNumberOfIds();
        this->mesh->GetPoint(i,point1);
        for(int j=0;j<nbInRing;j++)
        {
            crp=curNeib->GetId(j);
            this->mesh->GetPoint(crp,point2);
            ec=vtkMath::Distance2BetweenPoints(point1,point2);
            ec=sqrt(ec);
            if(maxDist<ec)maxDist=ec;
        }
        //this->rough->InsertNextValue(0);
        curNeib->Delete();

    }
    //if(1) cout<<"max Dist:"<<maxDist<<endl;

    //if(1) cout<<"Entremig"<<endl;

  //  getc(stdin);

    nmaxDist=maxDist/2.0;

    //if(1) cout<<"n max Dist:"<<nmaxDist<<endl;

/*    ComputeCurvature(0.01);
    for(int i=0;i<this->nbPoints;i++)
    {
        std=0;
        mean=0;
        ringSurf=0;
        curNeib->Reset();
        GetPointNeighbors(i,curNeib);
        nbInRing=curNeib->GetNumberOfIds();
        for(int j=0;j<nbInRing;j++)
        {
            crp=curNeib->GetId(j);
            std+=this->pointSurf->GetValue(crp)*pow(this->curv->GetValue(crp)-0.5,2);
            //if(VERBOSE) cout<<this->curv->GetValue(crp)<<endl;
            mean+=this->pointSurf->GetValue(crp)*(this->curv->GetValue(crp)-0.5);
            ringSurf+=this->pointSurf->GetValue(crp);
        }
        std=sqrt((std/ringSurf-pow(mean,2)));
        this->rough->InsertNextValue(std);
    }
*/

    vtkPointLocator* pl=vtkPointLocator::New();
    pl->SetDataSet(this->mesh);
    pl->BuildLocator();
    //vtkIdList* Nclo = vtkIdList::New();

    double maxC,minC;
    ComputeCurvature(0.2);
    for(int i=0;i<this->nbPoints;i++)
    {
        //maxC=-10;
        //minC=10;
        std=0;
        mean=0;
        ringSurf=0;
        vtkIdList* Nclo = vtkIdList::New();
        Nclo->Reset();
        this->mesh->GetPoint(i,point1);
        pl->FindPointsWithinRadius(nmaxDist,point1,Nclo);
        nbInRing=Nclo->GetNumberOfIds();

        if (i==0){if(1) cout<<"Numero de punts dins Ring per la roughness:"<<nbInRing<<endl;}

        //GeoDistRing(i,maxDist);
        //nbInRing=this->inRing->GetNumberOfIds();

        for(int j=0;j<nbInRing;j++)
        {
        	//dist euclidienne
            crp=Nclo->GetId(j);
            //dist geodesique
            //crp=this->inRing->GetId(j);
            std+=this->pointSurf->GetValue(crp)*pow(this->curv->GetValue(crp)-0.5,2);
            //if(VERBOSE) cout<<this->curv->GetValue(crp)<<endl;
            mean+=this->pointSurf->GetValue(crp)*(this->curv->GetValue(crp)-0.5);
            ringSurf+=this->pointSurf->GetValue(crp);
            //if (maxC<this->curv->GetValue(crp))maxC=this->curv->GetValue(crp);
            //if (minC>this->curv->GetValue(crp))minC=this->curv->GetValue(crp);
       }

        if(ringSurf==0)std=0;
        else std=sqrt(fabs(std/ringSurf-(mean/ringSurf)*(mean/ringSurf)));

        //if(isnan(std))if(VERBOSE) cout<<"Ringsurf"<<ringSurf<<endl;
        this->rough->InsertNextValue(std);
        //this->rough->InsertNextValue(maxC-minC);
    }

    ComputeCurvature(0.5);
    double maxRough=-100, minRough=100,curRough;
    for(int i=0;i<this->nbPoints;i++)
    {
        //curRough=(this->visibility->GetValue(i)+1)*(this->rough->GetValue(i)+this->depth->GetValue(i)/2+fabs(this->curv->GetValue(i)-0.5));
        //curRough=this->rough->GetValue(i)+this->depth->GetValue(i)/2+fabs(this->curv->GetValue(i)-0.5);
        //curRough=this->rough->GetValue(i)+fabs(this->curv->GetValue(i)-0.5);
        curRough=this->rough->GetValue(i);

        this->rough->SetValue(i,curRough);
        if(maxRough<curRough)maxRough=curRough;
        if(minRough>curRough)minRough=curRough;
    }
    for(int i=0;i<this->nbPoints;i++)
    {

        this->rough->SetValue(i,(this->rough->GetValue(i)-minRough)/(maxRough-minRough));
    }

    //if(1) cout<<"ComputeRoughness done"<<endl;
}

void MeshAnalyser::ComputeMask()
{
	if(1) cout<<"Dins ComputeMask"<<endl;

	this->mask->Reset();

	this->ComputeRoughness();
	this->ComputeCurvature(0.7);


	int i=0;

	if(1) cout<<"ComputeCurvature and Roughness done"<<endl;
	//getc(stdin);
	for(i=0;i<nbPoints;i++)
	{
		//if(isnan(this->curv->GetValue(i))) if(VERBOSE) cout<<"computecurv:"<<this->curv->GetValue(i)<<endl;
	}
/*	if(VERBOSE) cout<<"computecurvature:"<<this->curv->GetValue(1)<<endl;
	if(VERBOSE) cout<<"computecurvature:"<<this->curv->GetValue(2)<<endl;
	if(VERBOSE) cout<<"computecurvature:"<<this->curv->GetValue(3)<<endl;
	if(VERBOSE) cout<<"computecurvature:"<<this->curv->GetValue(4)<<endl;
	if(VERBOSE) cout<<"computecurvature:"<<this->curv->GetValue(5)<<endl;
*/
	//getc(stdin);

	for(i=0;i<nbPoints;i++)
	{
		//if(isnan(this->rough->GetValue(i))) if(VERBOSE) cout<<"computerough:"<<this->rough->GetValue(i)<<endl;
	}
/*	if(VERBOSE) cout<<"computeroughness:"<<this->rough->GetValue(1)<<endl;
	if(VERBOSE) cout<<"computeroughness:"<<this->rough->GetValue(2)<<endl;
	if(VERBOSE) cout<<"computeroughness:"<<this->rough->GetValue(3)<<endl;
	if(VERBOSE) cout<<"computeroughness:"<<this->rough->GetValue(4)<<endl;
	if(VERBOSE) cout<<"computeroughness:"<<this->rough->GetValue(5)<<endl;
*/
	for(int i=0;i<this->nbPoints;i++)
	{
		//if(VERBOSE) cout<<"iteració:"<<i<<endl;
		this->mask->InsertNextValue((this->curv->GetValue(i) * this->rough->GetValue(i)));
		//this->mask->SetValue(i,((this->curv->GetValue(i))*(this->rough->GetValue(i))));
	}

	for(i=0;i<nbPoints;i++)
	{
		//if(isnan(this->mask->GetValue(i))) if(VERBOSE) cout<<"computemask:"<<this->mask->GetValue(i)<<endl;
	}


/*	if(VERBOSE) cout<<"ComputeMask done. Mask(v):"<<this->mask->GetValue(1)<<endl;
	if(VERBOSE) cout<<"ComputeMask done. Mask(v):"<<this->mask->GetValue(2)<<endl;
	if(VERBOSE) cout<<"ComputeMask done. Mask(v):"<<this->mask->GetValue(3)<<endl;
	if(VERBOSE) cout<<"ComputeMask done. Mask(v):"<<this->mask->GetValue(4)<<endl;
	if(VERBOSE) cout<<"ComputeMask done. Mask(v):"<<this->mask->GetValue(5)<<endl;
*/
	//getc(stdin);
}

void MeshAnalyser::ComputeParamCorsini()
{

	this->ComputeCurvature(0.7);

    double maxDist=0,ec,minDist=0;
    //GeoDistRing(0,1000);
    double point1[3],point2[3];

    int nbInRing;

    double std=0,ringSurf=0,mean=0,sum=0,sum2=0;
    vtkIdType crp;

    vtkIdList* curNeib=vtkIdList::New();

    //this->IdClosNeib->Reset();

    for(int i=0;i<this->nbPoints;i++)
    {
        GetPointNeighbors(i,curNeib);
        nbInRing=curNeib->GetNumberOfIds();
        this->mesh->GetPoint(i,point1);
        for(int j=0;j<nbInRing;j++)
        {
            crp=curNeib->GetId(j);
            this->mesh->GetPoint(crp,point2);
            ec=vtkMath::Distance2BetweenPoints(point1,point2);
            ec=sqrt(ec);
            if(maxDist<ec)maxDist=ec;
        }
    }



    vtkPointLocator* pl=vtkPointLocator::New();
    pl->SetDataSet(this->mesh);
    pl->BuildLocator();
    vtkIdList* Nclo = vtkIdList::New();

    double maxC,minC;
    //ComputeCurvature(0.2);
    for(int i=0;i<this->nbPoints;i++)
    {
        //maxC=-10;
        //minC=10;
        std=0;
        mean=0;
        ringSurf=0;
        Nclo->Reset();
        this->mesh->GetPoint(i,point1);
        pl->FindPointsWithinRadius(maxDist/2.0,point1,Nclo);
        nbInRing=Nclo->GetNumberOfIds();

        //if(VERBOSE) cout<<"NbInRing:"<<nbInRing<<endl;

       // getc(stdin);

        //bucle per calcular mu

        for(int j=0;j<nbInRing;j++)
        {
        	//if(VERBOSE) cout<<"Dins for per calcular mu"<<endl;
        	//dist euclidienne
            crp=Nclo->GetId(j);
            //dist geodesique
            //crp=this->inRing->GetId(j);
            mean+=(this->curv->GetValue(crp))/nbInRing;
            //if(VERBOSE) cout<<"mean:"<<mean<<endl;
        }

        this->mu->InsertNextValue(mean);

        //if(VERBOSE) cout<<"mu:"<<mu->GetValue(i)<<endl;

        //bucle per calcular sigma

        sum=0;

        for(int j=0;j<nbInRing;j++)
        {
        	//dist euclidienne
            crp=Nclo->GetId(j);
            //dist geodesique
            //crp=this->inRing->GetId(j);
            sum+=pow(this->curv->GetValue(crp)-this->mu->GetValue(i),2);

            if(i==1)
            {
				if(1) cout<<"curv "<<this->curv->GetValue(crp)<<endl;
	    		if(1) cout<<"crp "<<crp<<endl;
            }

        }

		if(nbInRing==0)nbInRing=1;
        std=sqrt(sum/nbInRing);

        this->theta->InsertNextValue(std);
        if(fabs(this->theta->GetValue(i))<1.0e-15)
        {
        	this->theta->SetValue(i,0.0);
        	//if(1) cout<<"fabs std <1.0e-100"<<endl;
        }

        if(i==1) if(1) cout<<"nbInring(i=1):"<<nbInRing<<endl;

        //if(VERBOSE) cout<<"theta:"<<theta->GetValue(i)<<endl;
}
		//if(1) cout<<"mu "<<this->mu->GetValue(1)<<endl;
		//if(1) cout<<"theta(i=1)"<<this->theta->GetValue(1)<<endl;
		//    if(VERBOSE) cout<<"mu 7 dins compute params"<<this->mu->GetValue(7)<<endl;
}

void MeshAnalyser::ComputeSigmaCorsini2(MeshAnalyser* ma2)
{
	double thetaxyx=0,thetaxyy=0;
	double p1[3];
	double ec=0,sum2=0,sum3=0;
	double point1[3],point2[3];
	int nbInRingm1=0, nbInRingm2=0 ,nbInRing1=0,nbInRing2=0;
	vtkIdList* curNeib=vtkIdList::New();
	vtkIdType crp;

	int i,j;
	double maxDist;

	vtkPolyData* mesh2 = ma2->GetMesh();

	vtkPointLocator *pointLocatora = vtkPointLocator::New();
	pointLocatora->SetDataSet(mesh2);
	pointLocatora->BuildLocator();

	vtkPointLocator *pointLocatorb = vtkPointLocator::New();
	pointLocatorb->SetDataSet(this->mesh);
	pointLocatorb->BuildLocator();

	vtkIdList* Nclo = vtkIdList::New();

	// Find closest point

	for(i=0;i<this->nbPoints;i++)
	{
		ma2->GetPointNeighbors(i,curNeib);
		nbInRingm1=curNeib->GetNumberOfIds();
		mesh2->GetPoint(i,point1);

		for(int j=0;j<nbInRingm1;j++)
		{
			crp=curNeib->GetId(j);
			mesh2->GetPoint(crp,point2);
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);
			if(maxDist<ec)maxDist=ec;
		}
		//if(i==1)if(VERBOSE) cout<<"nbInRingm1: "<<nbInRingm1<<endl;
	}

	for(i=0;i<this->nbPoints;i++)
	{
		GetPointNeighbors(i,curNeib);
		nbInRingm2=curNeib->GetNumberOfIds();
		//if(i==1)if(VERBOSE) cout<<"nbInRingm2: "<<nbInRingm2<<endl;
	}

	for(i=0;i<this->nbPoints;i++)
	{
	    Nclo->Reset();
	    mesh2->GetPoint(i,point1);
	    pointLocatora->FindPointsWithinRadius(maxDist/2.0,point1,Nclo);
	    nbInRing1=Nclo->GetNumberOfIds();

	    //if(VERBOSE) cout<<"NbInRing1 dins ComputeSigmaCorsini:"<<nbInRing1<<endl;

	//bucle per calcular thetaxyx

	    sum2=0;

	    for(int j=0;j<nbInRing1;j++)
	    {
	    	crp=Nclo->GetId(j);
	    	//dist euclidienne
	    	sum2+=(ma2->curv->GetValue(crp)-ma2->mu->GetValue(i))*(this->curv->GetValue(crp)-this->mu->GetValue(i));
	    }
	    if(nbInRing1==0)
	    {
	    	nbInRing1=1;
	    	thetaxyx=0;
	    }
	    else thetaxyx=sum2/nbInRing1;
	//    if(i==7) if(VERBOSE) cout<<"thetaxyx 7 "<<thetaxyx<<endl;

	 Nclo->Reset();
	 this->mesh->GetPoint(i,point1);
	 pointLocatorb->FindPointsWithinRadius(maxDist/2.0,point1,Nclo);
	 nbInRing2=Nclo->GetNumberOfIds();
	 //if(VERBOSE) cout<<"NbInRing2 dins ComputeSigmaCorsini:"<<nbInRing2<<endl;

	 sum3=0;

	 for(int j=0;j<nbInRing2;j++)
	 {
	  	crp=Nclo->GetId(j);
	   	//dist euclidienne
	   	sum3+=(this->curv->GetValue(crp)-this->mu->GetValue(i))*(ma2->curv->GetValue(crp)-ma2->mu->GetValue(i));
	 }
	 if(nbInRing2==0)
	 {
		 nbInRing2=1;
		 thetaxyy=0;
	 }
	 else thetaxyy=sum3/nbInRing2;

    this->thetaxy->InsertNextValue((thetaxyx+thetaxyy)/2);
    //if(VERBOSE) cout<<"thetaxy:"<<thetaxy->GetValue(i)<<endl;

    if(fabs(this->thetaxy->GetValue(i))<1.0e-15)
    {
    	this->thetaxy->SetValue(i,0.0);
    	//if(1) cout<<"fabs thetaxy <1.0e-15"<<endl;
    }

    if(i==1) if(1) cout<<"nbInRing1(i=1):"<<nbInRing1<<endl;
    if(i==1) if(1) cout<<"nbInRing2(i=1):"<<nbInRing2<<endl;
}
//	if(VERBOSE) cout<<"thetaxy 7 "<<this->thetaxy->GetValue(7);
	if(1) cout<<"thetaxy(i=1)"<<this->thetaxy->GetValue(1)<<endl;
}

void MeshAnalyser::ComputeSigmaCorsini(vtkPolyData* wMesh,vtkPolyData* Mesh, MeshAnalyser* ma2, MeshAnalyser* ma3)
{
	double thetaxyx=0,thetaxyy=0;
	double p1[3];
	double ec=0,sum2=0,sum3=0;
	double point1[3],point2[3];
	vtkIdType ptId;
	int nbInRingm1=0, nbInRingm2=0 ,nbInRing1=0,nbInRing2=0;
	vtkIdList* curNeib=vtkIdList::New();
	vtkIdType crp;

	vtkPointLocator* p3=vtkPointLocator::New();
	p3->SetDataSet(this->mesh);
	p3->BuildLocator();
	vtkIdList* Nclo = vtkIdList::New();


	vtkPointLocator *pointLocatora = vtkPointLocator::New();
	pointLocatora->SetDataSet(Mesh);
	pointLocatora->BuildLocator();

	vtkPointLocator *pointLocatorb = vtkPointLocator::New();
	pointLocatorb->SetDataSet(wMesh);
	pointLocatorb->BuildLocator();

	int i,j;
	double maxDist1;
	double maxDist2;



	// Find closest point

	for(i=0;i<this->nbPoints;i++)
	{


		GetPointNeighbors2(i,curNeib, ma2->GetMesh());
		nbInRingm1=curNeib->GetNumberOfIds();
		ma2->GetMesh()->GetPoint(i,point1);

		for(int j=0;j<nbInRingm1;j++)
		{
			crp=curNeib->GetId(j);
			ma2->GetMesh()->GetPoint(crp,point2);
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);
			if(maxDist1<ec)maxDist1=ec;
		}
		if(i==1)if(1) cout<<"nbInRingm1: "<<nbInRingm1<<endl;
	}

	for(i=0;i<this->nbPoints;i++)
	{
		GetPointNeighbors(i,curNeib);
		nbInRingm2=curNeib->GetNumberOfIds();
		this->mesh->GetPoint(i,point1);

		for(int j=0;j<nbInRingm2;j++)
		{
			crp=curNeib->GetId(j);
			this->mesh->GetPoint(crp,point2);
			ec=vtkMath::Distance2BetweenPoints(point1,point2);
			ec=sqrt(ec);
			if(maxDist2<ec)maxDist2=ec;
		}
		if(i==1)if(1) cout<<"nbInRingm2: "<<nbInRingm2<<endl;
	}




	for(i=0;i<this->nbPoints;i++)
	{
	    Nclo->Reset();
	    ma2->GetMesh()->GetPoint(i,point1);
	    pointLocatora->FindPointsWithinRadius(maxDist1/2.0,point1,Nclo);
	    nbInRing1=Nclo->GetNumberOfIds();

	    //if(VERBOSE) cout<<"NbInRing1 dins ComputeSigmaCorsini:"<<nbInRing1<<endl;

	//bucle per calcular thetaxy x-based

	    for(int j=0;j<nbInRing1;j++)
	    {
	    	crp=Nclo->GetId(j);
	    	ma2->GetMesh()->GetPoint(crp,p1);
	    	ptId = pointLocatorb->FindClosestPoint(p1);
	    	//dist euclidienne

	    	sum2+=(ma2->curv->GetValue(crp)-ma2->mu->GetValue(i))*(ma3->curv->GetValue(ptId)-ma3->mu->GetValue(i));
	    	if(i==1)
	    	{
	    		if(1) cout<<"ma2 curv inRing1 "<<ma2->curv->GetValue(crp)<<endl;
	    		if(1) cout<<"ma2 mu inRing1 "<<ma2->mu->GetValue(i)<<endl;
	    		if(1) cout<<"ma3 curv inRing1 "<<ma3->curv->GetValue(ptId)<<endl;
	    		if(1) cout<<"ma3 mu inRing1 "<<ma3->mu->GetValue(i)<<endl;
	    		if(1) cout<<"crp "<<crp<<"ha de ser igual a ptId "<<ptId<<endl;

	    	}
	    	//if(crp!=ptId) if(VERBOSE) cout<<"Iteració "<<i<<"crp: "<<crp<<" ptId: "<<ptId<<endl;

	    }

	    if(nbInRing1==0) nbInRing1=1;
	    else thetaxyx=sum2/sqrt((double)nbInRing1);
	//    if(i==7) if(VERBOSE) cout<<"thetaxyx 7 "<<thetaxyx<<endl;



	    Nclo->Reset();
	    this->mesh->GetPoint(i,point1);
	    pointLocatorb->FindPointsWithinRadius(maxDist2/2.0,point1,Nclo);
	    nbInRing2=Nclo->GetNumberOfIds();

	    //if(VERBOSE) cout<<"NbInRing2 dins ComputeSigmaCorsini:"<<nbInRing2<<endl;

	    for(int j=0;j<nbInRing2;j++)
	    {
	    	crp=Nclo->GetId(j);
	    	wMesh->GetPoint(crp,p1);
	    	ptId = pointLocatora->FindClosestPoint(p1);
	    	//dist euclidienne

	    	sum3+=(ma3->curv->GetValue(crp)-ma3->mu->GetValue(i))*(ma2->curv->GetValue(ptId)-ma2->mu->GetValue(i));
	    	if(i==1)
	    	{
	    		if(1) cout<<"ma2 curv inRing 2 "<<ma2->curv->GetValue(ptId)<<endl;
	    		if(1) cout<<"ma2 mu inRing 2 "<<ma2->mu->GetValue(i)<<endl;
	    		if(1) cout<<"ma3 curv inRing 2 "<<ma3->curv->GetValue(crp)<<endl;
	    		if(1) cout<<"ma3 mu inRing 2 "<<ma3->mu->GetValue(i)<<endl;
	    		if(1) cout<<"crp "<<crp<<"ha de ser igual a ptId "<<ptId<<endl;

	    	}
	    	//if(crp!=ptId) if(VERBOSE) cout<<"Iteració "<<i<<"crp: "<<crp<<" ptId: "<<ptId<<endl;
	    }

	    if(nbInRing2==0) nbInRing2=1;
	    thetaxyy=sum2/sqrt((double)nbInRing2);
//	    if(i==7) if(VERBOSE) cout<<"thetaxyy 7 "<<thetaxyy<<endl;


	//bucle per calcular thetaxy y-based

    this->thetaxy->InsertNextValue(((thetaxyx)+(thetaxyy))/2);
    //if(VERBOSE) cout<<"thetaxy:"<<thetaxy->GetValue(i)<<endl;

    if(fabs(this->thetaxy->GetValue(i))<1.0e-15)
    {
    	this->thetaxy->SetValue(i,0.0);
    	//if(1) cout<<"fabs thetaxy <1.0e-15"<<endl;
    }

    if(i==1) if(1) cout<<"nbInRing1(i=1):"<<nbInRing1<<endl;
    if(i==1) if(1) cout<<"nbInRing2(i=1):"<<nbInRing2<<endl;
}
//	if(VERBOSE) cout<<"thetaxy 7 "<<this->thetaxy->GetValue(7);
	if(1) cout<<"thetaxy(i=1)"<<this->thetaxy->GetValue(1)<<endl;
}


void MeshAnalyser::ComputeCorsini(MeshAnalyser* ma2)
{
	if(1) cout<<"Dins ComputeCorsini"<<endl;

	double msdm, sumlmsdm =0.0;
	double stand_dev_cross=0.0;
	double a=3.0;
	double alph=0.4;
	double beta=0.4;
	double gamm=0.2;
	double maxu=0.0;
	double maxtheta=0.0;

	//L(curvature comparison) calculation

	double maxL=0.0;
	double imaxL=0.0;

	for(int i=0;i<this->nbPoints;i++)
	{
		if(this->mu->GetValue(i)>ma2->GetMu()->GetValue(i))maxu=this->mu->GetValue(i);
		else maxu=ma2->GetMu()->GetValue(i);
		if(maxu==0.0)
		{
			this->L->InsertNextValue(0);
			//if(1) cout<<"L és 0"<<endl;
		}
		else this->L->InsertNextValue(fabs((this->mu->GetValue(i))-(ma2->GetMu()->GetValue(i)))/maxu);
		if(maxL<this->L->GetValue(i))
		{
			maxL=this->L->GetValue(i);
			imaxL=i;
		}
	}

	//if(1) cout<<"L(1)"<<this->L->GetValue(1)<<endl;
	//if(1) cout<<"L(2)"<<this->L->GetValue(2)<<endl;
	//if(1) cout<<"L(3)"<<this->L->GetValue(3)<<endl;
	//if(1) cout<<"L(4)"<<this->L->GetValue(4)<<endl;
	//if(1) cout<<"L(1000)"<<this->L->GetValue(1000)<<endl;
	if(1) cout<<"L maxim"<<maxL<<endl;
	//if(1) cout<<"id L maxim"<<imaxL<<endl;
	if(1) cout<<"L calculada"<<endl;

	//C(contrast comparison) calculation

	double maxC=0.0;
	double imaxC=0.0;

	for(int i=0;i<this->nbPoints;i++)
	{
		//if(VERBOSE) cout<<"dins bucle C iteració"<<i<<endl;
		if(this->theta->GetValue(i)>ma2->GetTheta()->GetValue(i))maxtheta=this->theta->GetValue(i);
		else maxtheta=ma2->GetTheta()->GetValue(i);
		if(maxtheta==0.0)
		{
			this->C->InsertNextValue(0);
			//if(1) cout<<"C és 0"<<endl;
		}
		else this->C->InsertNextValue(fabs(this->theta->GetValue(i)-ma2->GetTheta()->GetValue(i))/maxtheta);
		//if(VERBOSE) cout<<"C ("<<i<<")"<<this->C->GetValue(i)<<endl;
		if(maxC<this->C->GetValue(i))
		{
			maxC=this->C->GetValue(i);
			imaxC=i;
		}
		if(i==415) if(1) cout<<"ma3 theta 415 "<<this->theta->GetValue(415)<<" ma2 theta "<<ma2->GetTheta()->GetValue(415)<<endl;
	}



/*	if(1) cout<<"C(1)"<<this->C->GetValue(1)<<endl;
	if(1) cout<<"C(2)"<<this->C->GetValue(2)<<endl;
	if(1) cout<<"C(3)"<<this->C->GetValue(3)<<endl;
	if(1) cout<<"C(4)"<<this->C->GetValue(4)<<endl;
	if(1) cout<<"C(1000)"<<this->C->GetValue(1000)<<endl;
*/	if(1) cout<<"C maxim"<<maxC<<endl;
//	if(1) cout<<"id C maxim"<<imaxC<<endl;
	if(1) cout<<"C calculada"<<endl;
	//S(structure comparison) calculation

	double maxS=0.0;
	double imaxS=0.0;

	for(int i=0;i<this->nbPoints;i++)
	{
		if (this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i)==0.0)
		{
			this->S->InsertNextValue(0);
			//if(1) cout<<"S és 0"<<endl;
		}
		else this->S->InsertNextValue(fabs(this->thetaxy->GetValue(i))/this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i));
		if(maxS<this->S->GetValue(i))
		{
			maxS=this->S->GetValue(i);
			imaxS=i;
		}

	}


/*	for(int i=0;i<this->nbPoints;i++)
	{
		if (this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i)==0.0)
		{
			this->S->InsertNextValue(0);
			if(VERBOSE) cout<<"S és 0"<<endl;
		}
		else this->S->InsertNextValue((fabs((this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i))-this->thetaxy->GetValue(i)))/(this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i)));
		//if (i<=10) if(VERBOSE) cout<<"C ("<<i<<")"<<this->C->getValue(i)<<endl;
		if(maxS<this->S->GetValue(i))
		{
			maxS=this->S->GetValue(i);
			imaxS=i;


		}

		/*if (this->S->GetValue(i) == std::numeric_limits<double>::infinity())
		{
			this->S->SetValue(i,0);
		}

		if((this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i))<this->S->GetValue(i))
		{
			maxS=this->S->GetValue(i);
			imaxS=i;
		}
		if(VERBOSE) cout<<"S maxim dins bucle 2n for"<<maxS<<endl;
		if(VERBOSE) cout<<"id S maxim dins bucle 2n for"<<imaxS<<endl;
		*/
/*
}

*/
/*	if(1) cout<<"S(0)"<<this->S->GetValue(0)<<endl;
	if(1) cout<<"S(1)"<<this->S->GetValue(1)<<endl;
	if(1) cout<<"S(2)"<<this->S->GetValue(2)<<endl;
	if(1) cout<<"S(3)"<<this->S->GetValue(3)<<endl;
	if(1) cout<<"S(4)"<<this->S->GetValue(4)<<endl;
	if(1) cout<<"S(1000)"<<this->S->GetValue(1000)<<endl;
	if(1) cout<<"S(8356)"<<this->S->GetValue(8356)<<endl;
	//(this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i))-this->thetaxy->GetValue(i)))/(this->theta->GetValue(i)*ma2->GetTheta()->GetValue(i)));
	if(1) cout<<"this->theta->GetValue(0): "<<this->theta->GetValue(0)<<endl;
	if(1) cout<<"ma2->GetTheta()->GetValue(0): "<<ma2->GetTheta()->GetValue(0)<<endl;
	if(1) cout<<"this->thetaxy->GetValue(0): "<<this->thetaxy->GetValue(0)<<endl;
	if(VERBOSE) cout<<"thetaxy(1)"<<this->thetaxy->GetValue(1)<<endl;
	if(VERBOSE) cout<<"thetaxy(2)"<<this->thetaxy->GetValue(2)<<endl;
	if(VERBOSE) cout<<"thetaxy(3)"<<this->thetaxy->GetValue(3)<<endl;
	if(VERBOSE) cout<<"thetaxy(4)"<<this->thetaxy->GetValue(4)<<endl;
	if(VERBOSE) cout<<"thetaxy(1000)"<<this->thetaxy->GetValue(1000)<<endl;
*/	if(1) cout<<"S maxim fora bucle"<<maxS<<endl;
//	if(1) cout<<"id S maxim fora bucle"<<imaxS<<endl;
	if(1) cout<<"S calculada"<<endl;

	//LMSDM Calculation

	this->lmsdm->Reset();

	cout<<"after lmsdm reset"<<endl;
	double maxlsmdm=0.0;
	double imaxlsmdm=0.0;

	for(int i=0;i<this->nbPoints;i++)
	{
		//if(1) cout<<"Dins el bucle lsmdm iteració"<<i<<endl;
		if(((alph*pow(this->L->GetValue(i),a))+(beta*pow(this->C->GetValue(i),a))+(gamm*pow(this->S->GetValue(i),a)))<0.0) if(1) cout<<"msdm negative"<<endl;
		this->lmsdm->InsertNextValue(pow(((alph*pow(this->L->GetValue(i),a))+(beta*pow(this->C->GetValue(i),a))+(gamm*pow(this->S->GetValue(i),a))),1/a));
		if(maxlsmdm<this->lmsdm->GetValue(i))
		{
			maxlsmdm=this->lmsdm->GetValue(i);
			imaxlsmdm=i;
		}

	}


	//if(1) cout<<"Lmsdm calculada"<<this->lmsdm->GetValue(1)<<endl;
	//if(1) cout<<"Lmsdm calculada"<<this->lmsdm->GetValue(2)<<endl;
	//if(1) cout<<"Lmsdm calculada"<<this->lmsdm->GetValue(3)<<endl;
	//if(1) cout<<"Lmsdm calculada"<<this->lmsdm->GetValue(4)<<endl;
	//if(1) cout<<"Lmsdm calculada"<<this->lmsdm->GetValue(1000)<<endl;
	if(1) cout<<"Lmsdm maxim"<<maxlsmdm<<endl;
	//if(1) cout<<"id lmsdm maxim"<<imaxlsmdm<<endl;
	//MSDM Calculation
	for(int i=0;i<this->nbPoints;i++)
	{
		if(this->lmsdm->GetValue(i)<0) if(1) cout<<"msdm negative"<<endl;
		sumlmsdm+= ((pow(this->lmsdm->GetValue(i),a))/(double)nbPoints);
	}

	if(sumlmsdm<=0) if(1) cout<<"sumlmsdm és negatiu o 0"<<endl;
	msdm = pow(sumlmsdm,(1/a));
	if(1) cout<<"msdm:"<<msdm<<endl;
}



vtkDoubleArray* MeshAnalyser::ComputeVoronoi(vtkIdList* centroidsList)
{
	//initialisation
	int nbCentroids	=centroidsList->GetNumberOfIds();
	vtkDoubleArray* shortestDists=vtkDoubleArray::New();
	vtkDoubleArray* VoronoiColors =vtkDoubleArray::New();
	
	double maxDist=1000;
	
	for(int i=0; i<this->nbPoints;i++)
	{
		shortestDists->InsertNextValue(maxDist);
		this->voronoiBin->InsertNextId(-1);
		VoronoiColors->InsertNextValue(0.0);
	}
	
	vtkIdType curCentroid;
	
	//compute the geodesic distance from centroids to all points
	for(int i=0;i<nbCentroids;i++)
	{
		curCentroid=centroidsList->GetId(i);
		GeoDistRing(curCentroid,maxDist);
		for(int j=0;j<this->nbPoints;j++)
		{
			if(shortestDists->GetValue(j)>this->geoDistRing->GetValue(j))
			{
					shortestDists->SetValue(j,this->geoDistRing->GetValue(j));
					this->voronoiBin->SetId(j,curCentroid);
					VoronoiColors->SetValue(j,(double) i/ (double) nbCentroids);
			}		
		}	
	}

	shortestDists->Delete();
	if(1) cout<<"Voronoi computed"<<endl;
	return VoronoiColors;
}


void MeshAnalyser::ProngsDetection()
{

	if(1) cout<<"dins prongsdetection"<<endl;

	double pt1[3];
	double pt2[3];
	double xMin=0.0,yMin=0.0,zMin=0.0,xMax=0.0,yMax=0.0,zMax=0.0;
	double diag = 0.0;

	for(int i=0;i<nbPoints;i++)
	{
		this->mesh->GetPoint(i,pt1);

		if(xMin>pt1[0])xMin=pt1[0];
		if(yMin>pt1[1])yMin=pt1[1];
		if(zMin>pt1[2])zMin=pt1[2];
		if(xMax<pt1[0])xMax=pt1[0];
		if(yMax<pt1[1])yMax=pt1[1];
		if(zMax<pt1[2])zMax=pt1[2];

	}

	diag = sqrt( (xMax-xMin)*(xMax-xMin)+(yMax-yMin)*(yMax-yMin)+(zMax-zMin)*(zMax-zMin) );

	if(1) cout<<"diag: "<<diag<<endl;

	//if(VERBOSE) cout<<"bounds: xmin: "<<this->BBox->GetBounds(1)<<endl;

	if(1) cout <<"abans intro"<<endl;
	 

	this->geoDistArray->Reset();
	double dist=0.0;
	if(1) cout <<"abans 2n intro"<<endl;
	 

	if(1) cout<<"reset done"<<endl;

/*	//code pour le calcul de la distance géodesique
	for(int i = 0; i<this->nbPoints; i++)
	{
		if(VERBOSE) cout<<"dins for iteració "<<i<<endl;

		dist=0.0;
		this->geoDistRing->Reset();
		//this->mesh->GetPoint(i,pt1);
		GeoDistRing(i, diag,400);

		if(VERBOSE) cout<<"dist geod: "<<this->GetGeoDistRing()->GetValue(i)<<endl;

		for(int j = 0; j < this->nbPoints; j++)
		{
			//if(VERBOSE) cout<<"dins 2n for iteració "<<j<<endl;
			if(i != j) dist+= this->GetGeoDistRing()->GetValue(j);
		}



		this->geoDistArray->InsertNextValue(dist);
	}
*/
	//code pour la distance euclidienne
	for	(int i = 0; i<this->nbPoints; i++)
	{
		//if(VERBOSE) cout<<"dins for distance euclidienne"<<endl;

		dist=0.0;

		this->geoDistRing->Reset();
		this->mesh->GetPoint(i,pt1);
		for(int j=0;j<this->nbPoints;j++)
		{
			this->mesh->GetPoint(j,pt2);
			//double *pt3= pt2;
			if(i!=j) dist+=sqrt((pt1[0]-pt2[0])*(pt1[0]-pt2[0])+(pt1[1]-pt2[1])*(pt1[1]-pt2[1])+(pt1[2]-pt2[2])*(pt1[2]-pt2[2]));
		}
		this->geoDistArray->InsertNextValue(dist);
	}

	if(1) cout<<"done"<<endl;
	if(1) cout<<"geodesic distance point 0"<<this->geoDistArray->GetValue(0)<<endl;
	if(1) cout<<"geodesic distance point 1"<<this->geoDistArray->GetValue(1)<<endl;
	if(1) cout<<"geodesic distance point 2"<<this->geoDistArray->GetValue(2)<<endl;
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////

//CENTER DATA

void MeshAnalyser::CenterData(vtkPolyData *p)
{
	int n = p->GetNumberOfPoints();
	double m[3] = {0,0,0};
	for(int i = 0; i < n ; i++){
		double a[3] = {0,0,0};
		p->GetPoint(i,a);
		for(int w =0; w < 3; w++)
			m[w] += a[w]/((double) n);
		}
	for(int i = 0; i < n ; i++){
		double a[3] = {0,0,0};
		p->GetPoint(i,a);
		for(int w =0; w < 3; w++)
			a[w] -= m[w];
		p->GetPoints()->SetPoint(i,a);
		}
	p->Update();

}


