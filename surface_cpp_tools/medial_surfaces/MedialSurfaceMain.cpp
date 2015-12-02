#include "MedialSurfaceComputer.h"

#include <vtkPolyDataReader.h>

#include <iostream>
#include <fstream>
#include <stdio.h>

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    vtkPolyDataReader* travelReader = vtkPolyDataReader::New();
    travelReader->SetFileName(argv[1]);
//  VTK6 Update: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Removal_of_Update
//  ???
    travelReader->Update();

    vtkPolyDataReader* fundiReader = vtkPolyDataReader::New();
    fundiReader->SetFileName(argv[2]);
//  VTK6 Update: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Removal_of_Update
//  ???
    fundiReader->Update();

    MedialSurfaceComputer* msc = new MedialSurfaceComputer( travelReader->GetOutput(), fundiReader->GetOutput() );

    msc->WriteIntoFile(argv[3]);

    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



