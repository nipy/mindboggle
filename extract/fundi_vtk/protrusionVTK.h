//
// File: protrusionVTK.h
// Created by:  Patrice Rondao Alface
// Created on: Tue Dec  2 11:37:19 2003
//

#ifndef _PROTRUSION_H_
#define _PROTRUSION_H_


//------------------------------------------------------------------------------
// OmbilicRobustness Declaration
//------------------------------------------------------------------------------
#include "geodesic_distanceVTK.cxx"
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTriangle.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPolyDataWriter.h>
#include <vtkMath.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkDecimatePro.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkDelaunay3D.h>
#include <map>
#include <list>
#include <set>
#include <vtkIdListCollection.h>

#include <math.h>
#include "PolyDataUtilities.cxx"
//#define INFINITE 1e15





class Protrusion
{
	//typedef double COST_TYPE;
  public:
		Protrusion();
		~Protrusion();
		
		vtkPolyData* GetOutput() {return polyData;} //retourne le polydata avec
		//la protrusion en dataarray. Appeler execute() pour que la classe
		//fasse qqcho.
		double* GetValues(){return values;}//rend juste le DataArray avec
		//la protrusion. Appeler Execute() avant.
		vtkPolyData* GetProngs(){return prongs;} //donne les prongs comme un
		//polydata. Appeler Execute() avant
		vtkIdList* GetProngsArray(){return prongsarray;}
		vtkDoubleArray* GetPatches(){return patches;}
		vtkIdList* GetPrec(){return prec;}
		vtkIdList* GetPatchesId(){return patchesids;}
		vtkIdListCollection* GetProngsIdList(){return ProngsIdList;}
		
		void SetInput(vtkPolyData * p){ polyData->DeepCopy(p); }
		void SetSmoothingFactor(double d){ smoothingFactor = d;}//permet de lisser
		//le maillage pour eviter trop de prongs (choisis d=1 pour pas de lissage
		//ou d=100 pour un filtrage raisonnable.
		void SetNeighBorhoodExtension(int s){ extension = s;}//etendue du voisinage
		//pour chercher un prong. s=1 c'est le One-Ring, si s=2 c'est le 2-ring, etc.
		
		void SetPrecision(double d){decimFactor = d;}//c'est le facteur de decimation.
		void SetVisualizationOn() { VISUALIZE = true; }//pour visualiser le polydata
		//avant de le sauver dans un fichier
		void SetGlyphFactor(double d){ sg = d;}	//pour la taille des boules pour les
		//prongs. Depend du mesh.
		void SetMoment(double m){moment = m;}//exposant des distances geodesiques dans
		//le calcul de protusion: m=2->prot_i=sum_j(geodesic_dist(i,j)^2*area_pj)
		void Execute();//Obligatoire de l'appeler pour executer le programme.
		void print_3D_balls(vtkPolyData *p, vtkPolyData *q, double ff);
		void print_3D_balls2(vtkPolyData *p, double ff);
		void print_3D_balls3(vtkPolyData *p, double ff);
		bool IsOnBorder(vtkPolyData *p, int i);
		void getPointNeighbours(vtkPolyData * polyData,int pt, vtkIdList * Neighbours);
		void ProngsDet(vtkPolyData *polyData);
		//void SetInput(vtkPolyData * p){ polyData = p; }
		vtkPolyData * find_max(vtkDoubleArray* prot, vtkPolyData* polydata);
		void print_3D(vtkPolyData *p);
		void ComputeNeighborhoodExtension(vtkPolyData *p, int num, int pt, int extension, int numpatch);
		vtkDoubleArray* ComputeFlaggedNeighborhoodExtension(vtkPolyData *p, int num, int pt, int extension, int numpatch, vtkDoubleArray*da, double flag, vtkDoubleArray *colorring);

		set<int>* ComputeNeighborhood(vtkPolyData *p, int i);

		
		vtkDoubleArray* patches;
		vtkIdList *prec;
		vtkIdList *patchesids;
		vtkIdListCollection* ProngsIdList;

	protected:
		bool is_max( int i, float f);
		bool is_max( int i, float f, vtkPolyData *p);
		bool is_min( int i, float f);

		set<int>* ComputeNeighborhood(int i);
		void ComputeProngs();
		void ComputeProtrusion(vtkPolyData *p, float ray, char* message);	
		
	private:
		
		int n;
		int extension;
		bool VISUALIZE;
		vtkPolyData* polyData;
		vtkPolyData* Out;
		vtkPolyData* prongs;

		double smoothingFactor;
		double decimFactor;
		double moment;
		double * values;
		double sg;
		vtkIdList* prongsarray;
};


//------------------------------------------------------------------------------
// Protrusion Implementation
//------------------------------------------------------------------------------


Protrusion::Protrusion()
{
		n = 0;
		VISUALIZE = false;
		polyData = vtkPolyData::New();
		prongs = vtkPolyData::New();
		prongsarray = vtkIdList::New();
		Out = vtkPolyData::New();
		extension = 1;
		smoothingFactor = 1;
		decimFactor = .9;
		sg = .1;
		moment = 1;
		patches=vtkDoubleArray::New();
		patchesids=vtkIdList::New();
		ProngsIdList=vtkIdListCollection::New();
}


Protrusion::~Protrusion()
{
	// TODO: put destructor code here
	polyData->Delete();
	Out->Delete();
	prongs->Delete();
	delete [] values;
	patches->Delete();
	patchesids->Delete();
	ProngsIdList->Delete();
}


#endif	//_PROTRUSION_H_
