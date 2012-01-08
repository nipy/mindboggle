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

    ma->ComputeClosedMesh();

    //may be replaced by ComputeTravelDepthFromInflated
    //which is slower but does not allow to cross the surface
    ma->ComputeEuclideanDepthFromClosed(true);    // Added by Forrest 2012-01-05, after getting Joachim's email

//    ma->ComputeTravelDepthFromClosed(true);  // Added by Forrest 2012-01-05, for testing

    ma->ComputeCurvature(0.7);

//    ma->WriteIntoFile((char*)"closed.vtk",(char*)"closed");

//    ma->WriteIntoFile(argv[2],(char*)"depth");  // Commented by Forrest 2012-01-05
//    ma->WriteIntoFile(argv[2],(char*)"Travel_Depth_and_Curv");  // Added by Forrest 2012-01-05, for testing
    ma->WriteIntoFile(argv[2],(char*)"Euclidean_Depth_and_Curv");  // Added by Forrest 2012-01-05, for testing

/*  commented by Forrest  2012-01-07 
    VtkFileEditor* vfe = new VtkFileEditor(argv[2]);
    vfe->CreateField((char*)"EuclideanDepth", ma->GetEuclideanDepth());  
    vfe->CreateField((char*)"curvature", ma->GetCurvature());

    delete vfe;
*/
    delete ma;
    delete fsr1;

    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;

}



