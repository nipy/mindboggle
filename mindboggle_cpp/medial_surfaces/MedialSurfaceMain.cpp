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
    travelReader->Update();

    vtkPolyDataReader* fundiReader = vtkPolyDataReader::New();
    fundiReader->SetFileName(argv[2]);
    fundiReader->Update();

    MedialSurfaceComputer* msc = new MedialSurfaceComputer( travelReader->GetOutput(), fundiReader->GetOutput() );

    msc->WriteIntoFile(argv[3]);

    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



