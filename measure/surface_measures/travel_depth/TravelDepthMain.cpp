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
    depthComputer->ComputeTravelDepth(true);
    depthComputer->WriteIntoFile(argv[2],(char*)"depth");



    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



