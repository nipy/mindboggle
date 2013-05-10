#include "../TravelDepth.h"
#include "../MeshAnalyser.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

int main(int argc, char** argv)
{
    time_t start= time(NULL);

//    TravelDepth* depthComputer = new TravelDepth(argv[1]);
//    depthComputer->ComputeDepth();
//    depthComputer->WriteIntoFile(argv[2]);

    MeshAnalyser* depthComputer = new MeshAnalyser(argv[1]);
    depthComputer->ComputeGeodesicDepthFromClosed(true);
    depthComputer->WriteIntoFile(argv[2],(char*)"geoDepth");

    cout<<"Elapsed time (geodesic depth): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



