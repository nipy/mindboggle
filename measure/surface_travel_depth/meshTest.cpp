#include "MeshAnalyser.h"
#include "BrainMeshing.h"

#include "vtkMetaImageReader.h"

#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkUnstructuredGrid.h>

#include <vtkDelaunay3D.h>
#include <vtkDecimatePro.h>

#include <vtkImageDifference.h>

#include <vtkXMLPolyDataReader.h>

#include <vtkTransform.h>
#include <vtkImageReader.h>
#include "VtkFileEditor.h"

#include "FsSurfaceReader.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    FsSurfaceReader* fsr1 = new FsSurfaceReader(argv[1]);

    // the reader for the vtp files
/*    vtkXMLPolyDataReader* reader = vtkXMLPolyDataReader::New();
    reader->SetFileName(argv[1]);
    reader->Update();*/

    MeshAnalyser* ma = new MeshAnalyser(fsr1->GetVTKData());

//    MeshAnalyser* ma1 = new MeshAnalyser(fsr1->GetVTKData());
//    ma1->Simplify(20000);


//    MeshAnalyser* ma = new MeshAnalyser(argv[1]);


    ma->ComputeTravelDepth(false);

//    MeshAnalyser* ma = new MeshAnalyser(argv[1]);

//    ma->ComputeClosedMesh();

    //may be replaced by ComputeTravelDepthFromClosed
    //which is slower but does not allow to cross the surface
//    ma->ComputeEuclideanDepth(true);

//    ma->ComputeCurvature(0.7);

//    ma->WriteIntoFile((char*)"closed.vtk",(char*)"closed");

//    ma->ComputeMedialSurfaces();

    ma->WriteIntoFile(argv[2], (char*)"depth");

//    VtkFileEditor* vfe = new VtkFileEditor(argv[2]);
//    vfe->CreateField((char*)"DEPTH", ma->GetEuclideanDepth());
//    vfe->CreateField((char*)"CURVATURE", ma->GetCurvature());
//    delete vfe;

    delete ma;
//    delete fsr1;


    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;



}



