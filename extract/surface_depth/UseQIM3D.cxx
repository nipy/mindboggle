#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataCollection.h>
#include <vtkPointLocator.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkMath.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkTransform.h>
#include <vtkCleanPolyData.h>
#include <vtkTubeFilter.h>
#include <vtkSplineWidget.h>
#include <vtkPolyLine.h>

#include <vtkSmartPointer.h>

#include <vtkBoostPrimMinimumSpanningTree.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkTree.h>

#include <vtkGraphLayoutView.h>
#include <vtkGraphLayoutStrategy.h>
#include <vtkGraphWriter.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSimple2DLayoutStrategy.h>

#include "dlib/image_io.h"
#include "dlib/gui_widgets.h"
#include "dlib/image_io.h"
#include "dlib/image_keypoint.h"
#include <fstream>

//#include "QIM3D3.cxx"
//#include "QIM3D.cxx"
#include "MeshAnalyser.cxx"
#include "meshAttack.cxx"
#include "protrusionVTK.cxx"

#include <set>

using namespace std;
using namespace dlib;

#define VERBOSE 0


int main(int argc,char** argv)
{

       char FileNameIn[50];
       sprintf(FileNameIn,"%s.vtk",argv[1]);

//---------//
//read mesh//
//---------//

       cout<<"Before reading mesh"<<endl;
       //getc(stdin);

       vtkPolyDataReader* meshReader=vtkPolyDataReader::New();
       meshReader->SetFileName(FileNameIn);
       meshReader->Update();

       vtkPolyDataConnectivityFilter *connectivity = vtkPolyDataConnectivityFilter::New();
       connectivity->SetInput(meshReader->GetOutput());
       connectivity->SetExtractionModeToLargestRegion();
       vtkPolyData *polydata = connectivity->GetOutput();
       polydata->Update();

       //render clourpoints
       vtkPolyData* wMesh1;
       wMesh1 = meshReader->GetOutput();

       Protrusion* prot=new Protrusion();
       //prot->print_3D(wMesh1);
       //prot->print_3D_balls2(wMesh1, 0.00001);

       //getc(stdin);


//-------------------------------//
//Hollows detection. Travel depth//
//-------------------------------//
#define TRAVELDEPTH 1

#if TRAVELDEPTH
       vtkPolyData* MeshDepth;
       MeshDepth = meshReader->GetOutput();

       MeshAnalyser* madepth = new MeshAnalyser(MeshDepth);

       if(VERBOSE) cout<<"abans ComputeTravelDepth"<<endl;

       //madepth->ComputeTravelDepth(false,true);
       madepth->ComputeTravelDepth(false);

       if(VERBOSE) cout<<"després ComputeTravelDepth"<<endl;

       double curDepth, tDepth, difDepth, meanDepth;
       vtkIdList* inRing=vtkIdList::New();
       int nbInRing;
       vtkIdType curId;
//       vtkDoubleArray* derDepth=vtkDoubleArray::New();

       int nbp = madepth->GetNumberOfPoints();

       if(VERBOSE) cout<<"abans moyenage"<<endl;

       int radi=madepth->GetDiag();

       if(VERBOSE) cout<<"diag: "<<radi<<endl;

       //getc(stdin);

       madepth->derDepth->Reset();

       if(VERBOSE) cout<<"després moyenage"<<endl;
       if(VERBOSE) cout<<"derDepth(1)"<<madepth->derDepth->GetValue(1)<<endl;

       char FileNameOutdepth[50];
       char FileNameOutdepth2[50];
       //sprintf(FileNameOutdepth,"%s-derDepth.vtk",argv[1]);
       sprintf(FileNameOutdepth2,"%s-depth.vtk",argv[1]);

       //madepth->WriteIntoFile(FileNameOutdepth,(char*)"derDepth");
       madepth->WriteIntoFile(FileNameOutdepth2,(char*)"depth");

       MeshDepth->GetPointData()->SetScalars(madepth->depth);


       Protrusion* print=new Protrusion();
       vtkPolyData *pq = vtkPolyData::New();
       //OPCIO 1: amb DerDepth
       //pq=print->find_max(madepth->GetDerDepth(),MeshDepth);

       //OPCIO 2: directament amb travel depth
       pq=print->find_max(madepth->GetTravelDepth(),MeshDepth);

       //OPCIO 1
       //MeshDepth->GetPointData()->SetScalars(madepth->GetDerDepth());
       print->print_3D_balls(pq,MeshDepth,0.9);

       int hollows = pq->GetNumberOfPoints();

       cout<<"hollows: "<<hollows<<endl;

       vtkDoubleArray* colorhol = vtkDoubleArray::New();

       for(int i=0;i<pq->GetNumberOfPoints();i++)
       {
    	   colorhol->InsertNextValue(1.0);
       }
       pq->GetPointData()->SetScalars(colorhol);

       prot->print_3D_balls3(pq,0.6);

       //OPCIO 2
       MeshDepth->GetPointData()->SetScalars(madepth->GetTravelDepth());

       //print->print_3D_balls(pq,MeshDepth,0.9);

       for(int i=0;i<print->GetProngsArray()->GetNumberOfIds();i++){

    	   cout<<"hollow num: "<<i<<" id: "<<print->GetProngsArray()->GetId(i)<<endl;
    	   cout<<"depth value:"<<madepth->depth->GetValue(print->GetProngsArray()->GetId(i))<<endl;
    	   int a = print->GetProngsArray()->GetId(i);
    }
       cout<<"això va"<<endl;

       //delete madepth;
       //delete print;

#endif

       cout<<"arriba fins aquí"<<endl;
       return 0;

}
