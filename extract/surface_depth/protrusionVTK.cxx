
#include "protrusionVTK.h"
#include <vtkDecimatePro.h>
#include <vtkSplineWidget.h>
#include "vtkPolyDataReader.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"


#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkSphereSource.h"
#include "vtkGlyph3D.h"
#include "vtkTriangle.h"

#include "vtkClipPolyData.h"
#include "vtkSphere.h"
//#include "PolySect.cxx"
#include <vtkPolyDataNormals.h>
//#include "vtkDecimate.h"
#include "vtkPointLocator.h"
//#include <protrusionVTK.cxx>
#include "vtkHull.h"
#include "vtkPolyDataWriter.h"
#include "vtkHull.h"
#include "vtkCellLocator.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "vtkDoubleArray.h"
#include "vtkGlyph3D.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkPlaneSource.h"
#include "vtkPlane.h"
#include "vtkProperty.h"
#include "vtkIdListCollection.h"
#include "vtkKochanekSpline.h"
#include "vtkParametricSpline.h"
#include "vtkSplineWidget.h"
#include "vtkGenericClip.h"
#include "vtkDataSetMapper.h"

#include <set>

//#include "geodesic_distanceVTK.cxx"

#define VERBOSE 1


using namespace std;







/*void print_list(list<double*> a, char *s){
	
	if(VERBOSE) cout<<"\n"<<s<<": ";
	for(list<double*>::iterator it = a.begin(); it!= a.end(); it++){
		if(VERBOSE) cout<<"\n"<<*it<<" "<<(*it)[0]<<", "<<(*it)[1]<<", "<<(*it)[2];
		
		
		}
		getc(stdin);
	
	}
*/
/*void getPointNeighbours(vtkPolyData * polyData,int pt, vtkIdList * Neighbours){


  vtkIdList *NC=vtkIdList::New();
  

  Neighbours->InsertNextId(1);
  Neighbours->Reset();

  polyData->GetPointCells(pt,NC);
  for(int i=0;i<NC->GetNumberOfIds();i++){

    vtkIdList *CP=vtkIdList::New();
    polyData->GetCellPoints(NC->GetId(i),CP);
    for(int j=0;j<CP->GetNumberOfIds();j++){

      if(Neighbours->IsId(CP->GetId(j))==-1 && CP->GetId(j)!=pt )
		Neighbours->InsertNextId(CP->GetId(j));
    }
    CP->Delete();
  }
  NC->Delete();
}
*/

//pt=id du point
set<int>* Protrusion::ComputeNeighborhood(int pt){
	
	set<int> *s = new set<int>  [1];
	int count = this->extension;
	
	
	vtkIdList *temp=vtkIdList::New();
	vtkIdList *prec=vtkIdList::New();
	vtkIdList *currentcrown=vtkIdList::New();
	prec->InsertNextId(pt);
	while( count  > 0){
		
		for(int k = 0; k < prec->GetNumberOfIds(); k++){
			vtkIdList *NC=vtkIdList::New();
			this->polyData->GetPointCells(prec->GetId(k),NC);
  	
			for(int i=0;i<NC->GetNumberOfIds();i++){

    			vtkIdList *CP=vtkIdList::New();
    			polyData->GetCellPoints(NC->GetId(i),CP);
    			for(int j=0;j<CP->GetNumberOfIds();j++){

      				if(temp->IsId(CP->GetId(j))==-1 && CP->GetId(j)!=pt )
						temp->InsertNextId(CP->GetId(j));
    				}
    			CP->Delete();
  				}
  			
				for(int i = 0; i < temp->GetNumberOfIds();i++){
	 				int vid = temp->GetId(i);
					if(s->find(vid)==s->end()){
						s->insert(vid);
						currentcrown->InsertNextId(vid);
					}
				}
				temp->Reset();
 				NC->Delete();
  				
			}
			prec->Reset();
			for(int i = 0; i < currentcrown->GetNumberOfIds();i++)
				prec->InsertNextId(currentcrown->GetId(i));
			currentcrown->Reset();
			
	--count;
	}//end while
	
	temp->Delete();
	prec->Delete();
  	return s;
 }
//extension = taille du ring extension =1 on a le 1-ring
set<int>* Protrusion::ComputeNeighborhood(vtkPolyData *p, int pt){

	set<int> *s = new set<int>  [1];
	//int count = this->extension;
	int count = 4;

	vtkIdList *temp=vtkIdList::New();
	vtkIdList *prec=vtkIdList::New();
	vtkIdList *currentcrown=vtkIdList::New();
	prec->InsertNextId(pt);
	while( count  > 0){

		for(int k = 0; k < prec->GetNumberOfIds(); k++){
			vtkIdList *NC=vtkIdList::New();
			p->GetPointCells(prec->GetId(k),NC);

			for(int i=0;i<NC->GetNumberOfIds();i++){

    			vtkIdList *CP=vtkIdList::New();
    			p->GetCellPoints(NC->GetId(i),CP);
    			for(int j=0;j<CP->GetNumberOfIds();j++){

      				if(temp->IsId(CP->GetId(j))==-1 && CP->GetId(j)!=pt )
						temp->InsertNextId(CP->GetId(j));
    				}
    			CP->Delete();
  				}

				for(int i = 0; i < temp->GetNumberOfIds();i++){
	 				int vid = temp->GetId(i);
					if(s->find(vid)==s->end()){
						s->insert(vid);
						currentcrown->InsertNextId(vid);
					}
				}
				temp->Reset();
 				NC->Delete();

			}
			prec->Reset();
			for(int i = 0; i < currentcrown->GetNumberOfIds();i++)
				prec->InsertNextId(currentcrown->GetId(i));
			currentcrown->Reset();

	--count;




	}//end while

	temp->Delete();
	prec->Delete();
  	return s;
 }

//-------------------------------------------------------------------------------

//extension = taille du ring extension =1 on a le 1-ring
void Protrusion::ComputeNeighborhoodExtension(vtkPolyData *p, int num, int pt, int extension, int numpatch){

	if(VERBOSE) cout<<"computeneighborhoodextension start"<<endl;
	set<int> *s = new set<int>  [1];
	int count = extension;

	vtkIdList *temp=vtkIdList::New();
	vtkIdList *prec=vtkIdList::New();
	vtkIdList *currentcrown=vtkIdList::New();
	vtkIdList *temppatchids=vtkIdList::New();

	this->patchesids->Reset();
	prec->InsertNextId(pt);

	if(VERBOSE) cout<<"computeneighborhoodextension abans while"<<endl;
	while( count  > 0){

		for(int k = 0; k < prec->GetNumberOfIds(); k++){
			vtkIdList *NC=vtkIdList::New();
			p->GetPointCells(prec->GetId(k),NC);

			for(int i=0;i<NC->GetNumberOfIds();i++)
			{
    			vtkIdList *CP=vtkIdList::New();
    			p->GetCellPoints(NC->GetId(i),CP);
    			for(int j=0;j<CP->GetNumberOfIds();j++)
    			{
      				if(temp->IsId(CP->GetId(j))==-1 && CP->GetId(j)!=pt )
      					temp->InsertNextId(CP->GetId(j));
   				}

    			CP->Delete();
  			}

			for(int i = 0; i < temp->GetNumberOfIds();i++)
			{
	 			int vid = temp->GetId(i);
	 			//if(VERBOSE) cout<<"iteració: "<<i<<endl;
				if(s->find(vid)==s->end())
				{
					if(VERBOSE) cout<<"entra dins l'if"<<endl;
					s->insert(vid);
					currentcrown->InsertNextId(vid);
					//if(VERBOSE) cout<<"abans patches->setvalue"<<endl;
					//this->patches->SetValue(temp->GetId(i),(double) numpatch/ (double) num);
					//if(VERBOSE) cout<<"patches->setvalue"<<endl;
					if(VERBOSE) cout<<"abans pacthes id insertnextid"<<endl;
					this->patchesids->InsertNextId(vid);
					if(VERBOSE) cout<<"patchesids insertnextid"<<endl;
					temppatchids->InsertNextId(vid);

				}

			}
			if(VERBOSE) cout<<"abans reset temp"<<endl;
			temp->Reset();
			if(VERBOSE) cout<<"després reset temp"<<endl;
 			NC->Delete();

		}
		prec->Reset();

		//if(VERBOSE) cout<<"anem a ficar un item"<<endl;
//		this->ProngsIdList->AddItem(this->patchesids);
//		this->ProngsIdList->AddItem(temppatchids);
		//if(VERBOSE) cout<<"no hi ha cap item ficat"<<endl;
		for(int i = 0; i < currentcrown->GetNumberOfIds();i++)
		{
			prec->InsertNextId(currentcrown->GetId(i));

		}
		currentcrown->Reset();

	--count;
	}//end while

	temp->Delete();
	prec->Delete();
	if(VERBOSE) cout<<"computeneighborhoodextension fin"<<endl;
  	//return VoronoiColors;
 }


//---------------------------------------------------------------------------------------------------------------------
	
bool Protrusion::is_max(int i, float f){

	bool out = true;
	set<int>*s = this->ComputeNeighborhood(i);
	set<int>::iterator sit = s->begin();
	double lambda = f * this->values[i];
	for(; sit != s->end(); ++sit){
		if(lambda<=values[*sit]){
			out = false;
			break;
			}
		}
	s->clear();
	return out;

	}

bool Protrusion::is_max(int i, float f, vtkPolyData *p){
/*	if(IsOnBorder(p,i))
	{
		return 0;
	}

	else
	{
*/
		bool out = true;
		set<int>*s = this->ComputeNeighborhood(p,i);
		set<int>::iterator sit = s->begin();
		double lambda = f * this->values[i];
		for(; sit != s->end(); ++sit){
			if(lambda<=values[*sit]){
				out = false;
				break;
				}
			}
		s->clear();
		return out;
//	}
}

//-------------------------------------------------------------------------------
vtkPolyData * Protrusion::find_max(vtkDoubleArray *prot, vtkPolyData *polydata){
	vtkPolyData *out = vtkPolyData::New();
	int n = polydata->GetNumberOfPoints();
	vtkPoints *pts = vtkPoints::New();
	vtkIdList * voisins = vtkIdList::New();
	double pt[3];

	this->prongsarray->Reset();

	for(int i = 0; i < n ; i++){
		this->getPointNeighbours(polydata, i, voisins);
		bool valid = true;
		for(int j = 0; j < voisins->GetNumberOfIds(); j++){
			if(IsOnBorder(polydata,i)){
				valid = false; break;
			}
			if(prot->GetValue(i) <= prot->GetValue(voisins->GetId(j))){
				valid = false; break;
			}
		}
		if(valid){
			polydata->GetPoint(i,pt);
			pts->InsertNextPoint(pt);
			this->prongsarray->InsertNextId(i);
		}
	}
	out->SetPoints(pts);
	out->Update();
 	return out;
}

//----------------------------------------------------------------------------------------------------------------------------
bool Protrusion::is_min(int i, float f ){
	
	bool out = true;
	set<int>*s = this->ComputeNeighborhood(i);	
	set<int>::iterator sit = s->begin();
	double lambda = f * this->values[i]; 
	for(; sit != s->end(); ++sit){
		if(lambda>values[*sit]){
			out = false;
			break;
			}
		}
	s->clear();	
	return out;	
	
	}


bool Protrusion::IsOnBorder(vtkPolyData *p, int i){
	vtkIdList *l = vtkIdList::New();
	vtkIdList *l2 = vtkIdList::New();
	p->GetPointCells(i,l);
	getPointNeighbours(p,i,l2);
	return (l->GetNumberOfIds()!=l2->GetNumberOfIds());
}

//---------------------------------------------------------------------------------------------------------------------

void Protrusion::getPointNeighbours(vtkPolyData * polyData,int pt, vtkIdList * Neighbours){

  vtkIdList *NC=vtkIdList::New();
  vtkIdList *CP=vtkIdList::New();
  Neighbours->InsertNextId(1);
  Neighbours->Reset();

  polyData->GetPointCells(pt,NC);
  for(int i=0;i<NC->GetNumberOfIds();i++){
    polyData->GetCellPoints(NC->GetId(i),CP);
    for(int j=0;j<CP->GetNumberOfIds();j++){

      if(Neighbours->IsId(CP->GetId(j))==-1 && CP->GetId(j)!=pt )
			Neighbours->InsertNextId(CP->GetId(j));
    }
    CP->Reset();
  }
  CP->Delete();
  NC->Delete();
}


void Protrusion::print_3D_balls(vtkPolyData *p, vtkPolyData *q, double ff){
	if(VERBOSE) cout<<"number points mesh q"<<q->GetNumberOfPoints()<<endl;
	if(VERBOSE) cout<<"number points mesh p"<<p->GetNumberOfPoints()<<endl;
	//create sphere geometry
	vtkSphereSource *s = vtkSphereSource::New();
	s->SetRadius(1.0);
	s->SetPhiResolution(6);
	s->SetThetaResolution(6);

	//vtkGlyph3D is a filter that copies a geometric representation (called a glyph) to every point in
	//the input dataset.
	vtkGlyph3D *gl = vtkGlyph3D::New();
//	vtkDoubleArray *da = vtkDoubleArray::New();
//	for(int i = 0; i < p->GetNumberOfPoints(); i++)
//		da->InsertNextTuple1(1);
//	p->GetPointData()
	gl->SetInput(p);
	gl->SetSource(s->GetOutput());
	gl->SetScaleFactor(ff);
	vtkPolyDataMapper *mm = vtkPolyDataMapper::New();
	mm->SetInput(gl->GetOutput());
	vtkActor *aa = vtkActor::New();
	aa->SetMapper(mm);

	//aa->GetProperty()->SetColor(1.0,1.0,1.0);

	vtkPolyDataMapper *m = vtkPolyDataMapper::New();
	vtkActor *a = vtkActor::New();
	vtkRenderer *r = vtkRenderer::New();
	vtkRenderWindow* rw = vtkRenderWindow::New();
	vtkRenderWindowInteractor *rwi = vtkRenderWindowInteractor::New();
	m->SetInput(q);
	a->SetMapper(m);
	//a->GetProperty()->SetColor(1,1,1);
	r->AddActor(a);r->AddActor(aa);
	r->SetBackground(1,1,1);
	rw->AddRenderer(r);


	rwi->SetRenderWindow(rw);

	rw->Render();
	rwi->Start();

	m->Delete();
	a->Delete();
	aa->Delete();
	r->Delete();
	rw->Delete();
	rwi->Delete();

}

void Protrusion::print_3D_balls2(vtkPolyData *p, double ff){

	if(VERBOSE) cout<<"number points mesh p"<<p->GetNumberOfPoints()<<endl;
	//create sphere geometry

	vtkCubeSource * c = vtkCubeSource::New();
	c->SetXLength(3.5);//200
	//c->SetYLength(0.01);
	c->SetYLength(0.5);//40
	c->SetZLength(0.5);//40

	vtkRenderer *r = vtkRenderer::New();
	vtkRenderWindow* rw = vtkRenderWindow::New();
	vtkRenderWindowInteractor *rwi = vtkRenderWindowInteractor::New();



//		vtkPlane *implicitPlane = vtkPlane::New();
//			          implicitPlane->SetOrigin(0.5, 0, 0); // (0.5, 0, 0);
//			          implicitPlane->SetNormal(1, 1, 1); // (1, 1, 1);


		vtkPlaneSource *planeSource= vtkPlaneSource::New();
		planeSource->SetPoint1(4,0,4);
		planeSource->SetPoint2(-4,0,-4);
		planeSource->SetXResolution(20);
		planeSource->SetYResolution(20);

		planeSource->Update();

		vtkPolyDataMapper* planeSourceMapper = vtkPolyDataMapper::New();
	    planeSourceMapper->SetInput(planeSource->GetOutput());
	    vtkActor* planeSourceActor = vtkActor::New();
	       planeSourceActor->SetMapper(planeSourceMapper);


	          // Create the filter



/*	     vtkSplineWidget* spline = vtkSplineWidget::New();
	       spline->SetInteractor(rwi);
	       spline->SetInput(planeSource->GetOutput());
	       spline->SetPriority(1.0);
	       spline->KeyPressActivationOff();
	       spline->PlaceWidget();
	       spline->ProjectToPlaneOn();
	       spline->SetProjectionNormal(0);
	       spline->SetProjectionPosition(102.4);  //initial plane oriented position
	       spline->SetProjectionNormal(3); //allow arbitrary oblique orientations
	       spline->SetPlaneSource(planeSource);




	     // Specify the type of spline (change from default vtkCardinalSpline)
	     vtkKochanekSpline* xspline = vtkKochanekSpline::New();
	     vtkKochanekSpline* yspline = vtkKochanekSpline::New();
	     vtkKochanekSpline* zspline = vtkKochanekSpline::New();

	     vtkParametricSpline* para = spline->GetParametricSpline();

	     para->SetXSpline(xspline);
	     para->SetYSpline(yspline);
	     para->SetZSpline(zspline);

	     vtkPolyData* poly = vtkPolyData::New();
	       spline->GetPolyData(poly);


	     //vtkSplineWidgetCallback* swcb = vtkSplineWidget::New();
	     //  swcb->Poly = poly;

	     //spline->AddObserver(vtkCommand::InteractionEvent,swcb);

*/



	//vtkSphereSource *s = vtkSphereSource::New();
	//s->SetRadius(1.0);
	//s->SetPhiResolution(6);
	//s->SetThetaResolution(6);

	//vtkGlyph3D is a filter that copies a geometric representation (called a glyph) to every point in
	//the input dataset.
	vtkGlyph3D *gl = vtkGlyph3D::New();
//	vtkDoubleArray *da = vtkDoubleArray::New();
//	for(int i = 0; i < p->GetNumberOfPoints(); i++)
//		da->InsertNextTuple1(1);
//	p->GetPointData()
	gl->SetInput(p);
	gl->SetSource(c->GetOutput());
	gl->SetScaleFactor(ff);
	vtkPolyDataMapper *mm = vtkPolyDataMapper::New();
	mm->SetInput(gl->GetOutput());
	vtkActor *aa = vtkActor::New();
	aa->SetMapper(mm);

	//aa->GetProperty()->SetColor(1.0,1.0,1.0);


//	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	//vtkDataSetMapper *mapper = vtkDataSetMapper::New();
//	mapper->SetInput(implicitPlane);
//	vtkActor *actor = vtkActor::New();
//	actor->SetMapper(mapper);




	//a->GetProperty()->SetColor(1,1,1);
	//vtkRenderer *r = vtkRenderer::New();
	//vtkRenderWindow* rw = vtkRenderWindow::New();
	//vtkRenderWindowInteractor *rwi = vtkRenderWindowInteractor::New();
	r->AddActor(planeSourceActor);
	r->AddActor(aa);
	//r->AddActor(actor);
	r->SetBackground(1,1,1);
	rw->AddRenderer(r);


	rwi->SetRenderWindow(rw);

	rw->Render();
	rwi->Start();


	aa->Delete();
	r->Delete();
	rw->Delete();
	rwi->Delete();

}

void Protrusion::print_3D_balls3(vtkPolyData *p, double ff){

	if(VERBOSE) cout<<"number points mesh p"<<p->GetNumberOfPoints()<<endl;
	//create sphere geometry
	vtkSphereSource *s = vtkSphereSource::New();
	s->SetRadius(1.0);
	s->SetPhiResolution(6);
	s->SetThetaResolution(6);

	//vtkGlyph3D is a filter that copies a geometric representation (called a glyph) to every point in
	//the input dataset.
	vtkGlyph3D *gl = vtkGlyph3D::New();
//	vtkDoubleArray *da = vtkDoubleArray::New();
//	for(int i = 0; i < p->GetNumberOfPoints(); i++)
//		da->InsertNextTuple1(1);
//	p->GetPointData()
	gl->SetInput(p);
	gl->SetSource(s->GetOutput());
	gl->SetScaleFactor(ff);
	vtkPolyDataMapper *mm = vtkPolyDataMapper::New();
	mm->SetInput(gl->GetOutput());
	vtkActor *aa = vtkActor::New();
	aa->SetMapper(mm);

	//aa->GetProperty()->SetColor(1.0,1.0,1.0);

	//vtkPolyDataMapper *m = vtkPolyDataMapper::New();
	//vtkActor *a = vtkActor::New();
	vtkRenderer *r = vtkRenderer::New();
	vtkRenderWindow* rw = vtkRenderWindow::New();
	vtkRenderWindowInteractor *rwi = vtkRenderWindowInteractor::New();
	//m->SetInput(q);
	//a->SetMapper(m);
	//a->GetProperty()->SetColor(1,1,1);
	//r->AddActor(a);
	r->AddActor(aa);
	r->SetBackground(0,0,0);
	rw->AddRenderer(r);


	rwi->SetRenderWindow(rw);

	rw->Render();
	rwi->Start();

	aa->Delete();
	r->Delete();
	rw->Delete();
	rwi->Delete();

}


void Protrusion::print_3D(vtkPolyData *p){

	vtkPolyDataMapper *m = vtkPolyDataMapper::New();
	vtkActor *a = vtkActor::New();
	vtkRenderer *r = vtkRenderer::New();
	vtkRenderWindow* rw = vtkRenderWindow::New();
	vtkRenderWindowInteractor *rwi = vtkRenderWindowInteractor::New();
	m->SetInput(p);
	a->SetMapper(m);
	r->AddActor(a);
	r->SetBackground(1,1,1);
	rw->AddRenderer(r);
	rwi->SetRenderWindow(rw);
	rw->Render();
	rwi->Start();

}


void Protrusion::ProngsDet(vtkPolyData *polyData)
{

	vtkDoubleArray *da = vtkDoubleArray::New();
	vtkDoubleArray *db = vtkDoubleArray::New();
	//vtkPoints*qts = vtkPoints::New();
	vtkPoints *qts2 = vtkPoints::New();
	//this->polyData->DeepCopy(polyData);
	//polyData->Update();
	//this->polyData->Update();
	this->n = polyData->GetNumberOfPoints();
	this->values = new double [n];
	//double val =0.0;

	da = (vtkDoubleArray*) polyData->GetPointData()->GetScalars();

	//qts = polyData->GetPoints();


	//vtkPolyData *q = vtkPolyData::New();

	if(VERBOSE) cout<<"n: "<<n<<endl;
	if(VERBOSE) cout<<"getpoint1:"<<polyData->GetPoint(1)<<endl;
	if(VERBOSE) cout<<"polydata scalar 0"<<da->GetValue(0)<<endl;
	if(VERBOSE) cout<<"polydata scalar 1"<<da->GetValue(1)<<endl;
	if(VERBOSE) cout<<"polydata scalar 2"<<da->GetValue(2)<<endl;
	if(VERBOSE) cout<<"abans for"<<endl;
	 
	int count = 0;
	polyData->BuildLinks();

	this->prongsarray->Reset();

	for(int i = 0; i < n; i++)
	{
		this->values[i] = da->GetValue(i);
	}
	for(int i = 0; i < n; i++){
		//if(VERBOSE) cout<<"iteració "<<i<<endl;
		//this->values[i] = da->GetValue(i);
		//if(VERBOSE) cout<<"value "<<i<<" :"<<this->values[i]<<endl;
		int p = int (10*this->values[i]);

	if((this->is_max( i, this->smoothingFactor, polyData))&&(!IsOnBorder(polyData,i))){
	//if((!IsOnBorder(polyData,i))){
	//if((this->is_max( i, this->smoothingFactor, polyData))){
		if(VERBOSE) cout<<"i: "<<i<<endl;
		if(VERBOSE) cout<<"count: "<<count<<endl;
			double a[3];
			polyData->GetPoint(i,a);
			this->prongsarray->InsertNextId(i);
			if(VERBOSE) cout<<"point a: "<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
			qts2->InsertNextPoint(a);
			if(VERBOSE) cout<<"qts2: "<<qts2->GetPoint(count)<<endl;
			db->InsertNextTuple1(this->values[i]);
			if(VERBOSE) cout<<"db ("<<count<<") es "<<db->GetValue(count)<<endl;
			count++;
			}
		}
	if(VERBOSE) cout<<"després for"<<endl;

	this->prongs->SetPoints(qts2);

	this->prongs->GetPointData()->SetScalars(db);
	this->prongs->Update();
	qts2->Delete();
	db->Delete();
	if (this->VISUALIZE){

		if(VERBOSE) cout<<"visualize es true"<<endl;

		vtkGlyph3D *gl = vtkGlyph3D::New();
		vtkSphereSource *sp = vtkSphereSource::New();

		//boules prennent moins de place en la memoire GPU
		sp->SetThetaResolution(4);
		sp->SetPhiResolution(4);

		gl->SetInput(this->prongs);
		gl->SetSource(sp->GetOutput());
		gl->SetScaleFactor(this->sg);
		//gl->SetColorModeToDataOff();

		vtkPolyDataMapper *map1 =  vtkPolyDataMapper::New();
		map1->SetInput(gl->GetOutput());
		vtkPolyDataMapper *map2 =  vtkPolyDataMapper::New();
		map2->SetInput(polyData);

		vtkActor *act1 = vtkActor::New();
		//act1->GetProperty()->SetColor(1,1,1);
		vtkActor *act2 = vtkActor::New();

		act1->SetMapper(map1);
		act2->SetMapper(map2);

		vtkRenderer *ren = vtkRenderer::New();

		ren->AddActor(act1);
		ren->AddActor(act2);
		vtkRenderWindow *rw = vtkRenderWindow::New();
		rw->AddRenderer(ren);

		vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
		iren->SetRenderWindow(rw);

		iren->Start();

	}
}


void Protrusion::Execute(){
	
	vtkDecimatePro * decim = vtkDecimatePro::New();
	vtkPolyData *polydata2 = vtkPolyData::New();
	polydata2->DeepCopy(this->polyData);
	polydata2->Update();
	decim->SetInput(polydata2);
	decim->SetTargetReduction(this->decimFactor);
	vtkPolyData * BasePoints = decim->GetOutput();
	BasePoints->Update();
	

	int nb = BasePoints->GetNumberOfPoints();
	this->n = this->polyData->GetNumberOfPoints();
	
	this->values = new double [n];
	for (int i = 0; i < n; i++) this->values[i] = 0.0;
	double *area_bi = new double [nb];
	
	printf("\ncompute area weights: nb = %d \n",nb);
	
	for(int j = 0; j < nb ; j++){
		
		vtkIdList *liste = vtkIdList::New();
		area_bi[j] = 0;
		BasePoints->GetPointCells(j,liste);
		for(int k = 0; k < liste->GetNumberOfIds(); k++){
			int id0 = BasePoints->GetCell(liste->GetId(k))->GetPointId(0);
			int id1 = BasePoints->GetCell(liste->GetId(k))->GetPointId(1);
			int id2 = BasePoints->GetCell(liste->GetId(k))->GetPointId(2);
			double A[3],B[3],C[3];
			BasePoints->GetPoint(id0,A);
			BasePoints->GetPoint(id1,B);
			BasePoints->GetPoint(id2,C);
			area_bi[j]+=vtkTriangle::TriangleArea( A, B, C); 			
		}
		area_bi[j] = area_bi[j]/3.0;
	//	printf("\n area_bi[%d] = %f\n",j,area_bi[j]);
		liste->Delete();
		}
		
	int pid = 0;
	vtkPointLocator *pl = vtkPointLocator::New(); 
	pl->SetDataSet(this->polyData);
	pl->BuildLocator();

	for(int j=0;j<nb;j++){
		cout<<"STOP1, j:"<<j<<endl;
		double A[3],B[3],C[3];
		BasePoints->GetPoint(j,A);
		int id = pl->FindClosestPoint(A);
		this->polyData->GetPoint(id,B);
		cout<<"STOP2"<<endl;
		
		Geodesic geo = Geodesic();
		this->polyData->GetPoint(0,C);
		cout<<"STOP3"<<endl;
		geo.set_polyhedron(this->polyData);
		geo.set_start_point(B);
		geo.set_end_point(C);
		geo.set_compute_all();
		//if(j%100==0) cout<<"\nexecution: "<<j<<" on "<<nb<<"\n"<<endl;
		cout<<"Abans execute"<<endl;
		geo.execute();
		cout<<"\nexecution: "<<j;
		//if(VERBOSE) cout<<"\nexecution: "<<j;
		//if(VERBOSE) cout<<"\nexecuvtkPolyData * ption";getc(stdin);
		double* dist = geo.get_all_distances();
		//if(VERBOSE) cout<<"\nend";getc(stdin);
		cout<<"Abans bucle"<<endl;
			for(int i = 0; i < n ; i++){
				this->values[i] += pow(dist[i],this->moment) * area_bi[j];
			}
		//delete [] dist;
		cout<<"Després bucle"<<endl;
		}
	pl->Delete();

	double min = 1e10;
	double max = 0.0;
	
	for(int i = 0; i < this->n; ++i) if (this->values[i]<min) min = this->values[i];
	for(int i = 0; i < this->n; ++i) if (this->values[i]>max) max = this->values[i];
	for(int i = 0; i < this->n; ++i) this->values[i] = (this->values[i]-min)/(max-min);	
	
	vtkDoubleArray *da = vtkDoubleArray::New();
	vtkDoubleArray *db = vtkDoubleArray::New();
	for(int i = 0; i < this->n; ++i) 
		da->InsertNextTuple1(this->values[i]);
	this->polyData->GetPointData()->SetScalars(da);
	da->Delete();
	//vtkPolyData *q = vtkPolyData::New();
	vtkPoints*qts = vtkPoints::New();
	
	//float sf[10] = {0,.9,.9,.95,.95,.99,.99,.99,1,1};
	
	for(int i = 0; i < n; ++i){ 
		int p = int (10*this->values[i]);
		if((this->is_max( i, this->smoothingFactor))&&(!IsOnBorder(this->polyData,i))){
			
			double a[3];
			polyData->GetPoint(i,a);
			qts->InsertNextPoint(a);
			db->InsertNextTuple1(this->values[i]);
			}
		}

	this->prongs->SetPoints(qts);
	this->prongs->GetPointData()->SetScalars(db);
	this->prongs->Update();	
	qts->Delete();
	db->Delete();

	if (this->VISUALIZE){
		
		vtkGlyph3D *gl = vtkGlyph3D::New();
		vtkSphereSource *sp = vtkSphereSource::New();
	
		//boules prennent moins de place en la memoire GPU
		sp->SetThetaResolution(4);
		sp->SetPhiResolution(4);

		gl->SetInput(this->prongs);
		gl->SetSource(sp->GetOutput());
		gl->SetScaleFactor(this->sg);
		//gl->SetColorModeToDataOff();
	
		vtkPolyDataMapper *map1 =  vtkPolyDataMapper::New();
		map1->SetInput(gl->GetOutput());
		vtkPolyDataMapper *map2 =  vtkPolyDataMapper::New();
		map2->SetInput(polyData);
	
		vtkActor *act1 = vtkActor::New();
		//act1->GetProperty()->SetColor(1,1,1);
		vtkActor *act2 = vtkActor::New();
	
		act1->SetMapper(map1);
		act2->SetMapper(map2);
	
		vtkRenderer *ren = vtkRenderer::New();
	
		ren->AddActor(act1);
		ren->AddActor(act2);
		vtkRenderWindow *rw = vtkRenderWindow::New();
		rw->AddRenderer(ren);

		vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
		iren->SetRenderWindow(rw);
	
		iren->Start();
		

		}
	delete [] area_bi;	
	decim->Delete();
	
}

//---------------------------------------------------------------------------------------------------------------------


vtkDoubleArray * Protrusion::ComputeFlaggedNeighborhoodExtension(vtkPolyData *p, int num, int pt, int extension, int numpatch, vtkDoubleArray *da, double flag, vtkDoubleArray *colorring){

	if(VERBOSE) cout<<"computeneighborhoodextension start"<<endl;
	set<int> *s = new set<int>  [1];
	int count = extension;

	vtkIdList *temp=vtkIdList::New();
	vtkIdList *prec=vtkIdList::New();
	vtkIdList *currentcrown=vtkIdList::New();
	vtkIdList *temppatchids=vtkIdList::New();

	this->patchesids->Reset();
	prec->InsertNextId(pt);

	if(VERBOSE) cout<<"computeneighborhoodextension abans while"<<endl;
	while( count  > 0){

		for(int k = 0; k < prec->GetNumberOfIds(); k++){
			vtkIdList *NC=vtkIdList::New();
			p->GetPointCells(prec->GetId(k),NC);

			for(int i=0;i<NC->GetNumberOfIds();i++)
			{
    			vtkIdList *CP=vtkIdList::New();
    			p->GetCellPoints(NC->GetId(i),CP);
    			for(int j=0;j<CP->GetNumberOfIds();j++)
    			{
      				if(temp->IsId(CP->GetId(j))==-1 && CP->GetId(j)!=pt )
      					temp->InsertNextId(CP->GetId(j));
   				}

    			CP->Delete();
  			}



			for(int i = 0; i < temp->GetNumberOfIds();i++)
			{


				int vid = temp->GetId(i);
	 			if(VERBOSE) cout<<"iteració: "<<i<<endl;
				if(s->find(vid)==s->end())
				{
					if(VERBOSE) cout<<"entra dins l'if"<<endl;
//here
					if(da->GetValue(vid)==flag){
						s->insert(vid);
						currentcrown->InsertNextId(vid);
						

						int coloret = ((int)(extension-count)) % 9;
						if(VERBOSE) cout<<"Val: "<<da->GetValue(vid)<<" prong_id "<<flag<<" numpointsinset "<<this->patchesids->GetNumberOfIds()<<endl; getc(stdin);
					//this->patches->SetValue(temp->GetId(i),(double) numpatch/ (double) num);
						colorring->SetValue(vid, (double) coloret / 9.0);
					if(VERBOSE) cout<<"patches->setvalue"<<endl;
					if(VERBOSE) cout<<"abans pacthes id insertnextid"<<endl;
						this->patchesids->InsertNextId(vid);
					if(VERBOSE) cout<<"patchesids insertnextid"<<endl;
						temppatchids->InsertNextId(vid);
					}

					else if(this->patchesids->GetNumberOfIds()>=64){
						//cout<<"we are touching a Voronoi boundary"<<endl;
						s->insert(vid);
						currentcrown->InsertNextId(vid);


						//int coloret = ((int)(extension-count)) % 9;

						//colorring->SetValue(vid, (double) coloret / 9.0);

						//this->patchesids->InsertNextId(vid);

						temppatchids->InsertNextId(vid);
					//	flag=-1;
						//cout<<"flag a -1"<<endl;
					}

				}

			}
			//if(VERBOSE) cout<<"abans reset temp"<<endl;
			temp->Reset();
			//if(VERBOSE) cout<<"després reset temp"<<endl;
 			NC->Delete();

		}
		prec->Reset();

		//if(VERBOSE) cout<<"anem a ficar un item"<<endl;
//		this->ProngsIdList->AddItem(this->patchesids);
//		this->ProngsIdList->AddItem(temppatchids);
		//if(VERBOSE) cout<<"no hi ha cap item ficat"<<endl;
		for(int i = 0; i < currentcrown->GetNumberOfIds();i++)
		{
			prec->InsertNextId(currentcrown->GetId(i));

		}
		currentcrown->Reset();

	--count;

    p->GetPointData()->SetScalars(colorring);

//    print_3D(p);





	}//end while

	temp->Delete();
	prec->Delete();
	if(VERBOSE) cout<<"computeneighborhoodextension fin"<<endl;
  	return colorring;
 }

//vtkDoubleArray * Protrusion::ComputeScalarsFromMesh(vtkPolyData *p, int num, int pt, int extension, int numpatch, vtkDoubleArray *da, double flag, vtkDoubleArray *colorring){

