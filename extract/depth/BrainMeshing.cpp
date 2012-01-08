/* ********************************************************************
 * BrainMeshing
 *
 * Copyright 2009 Joachim Giard
 * Université catholique de Louvain, Belgium
 * Apache License, Version 2.0
 *
 * Class containing methods to construct a 3D meshes from
 * a 3D image with various parameters.
 *
 * *******************************************************************/

#include <vtkImageContinuousErode3D.h>
#include <vtkImageContinuousDilate3D.h>
//#include <vtkPolyDataWriter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkMarchingContourFilter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkTriangleFilter.h>

#include <vtkImageFFT.h>

#include <vtkImageShrink3D.h>

#include <vtkImageIdealLowPass.h>
#include <vtkImageRFFT.h>

#include <vtkTransformPolyDataFilter.h>

#include "BrainMeshing.h"

using namespace std;

#define vround(a) ((int) (a+0.5))
#define lround(a) ((int) (a-0.5))
#define max(a,b) (a>b?a:b)


BrainMeshing::BrainMeshing(vtkImageData *inputImage)
{
    m_inputImage = inputImage;
    m_shrinkFactor = 3;
    m_probeRadius = 5;
    m_threshold = 70;

    m_modifiedImage = vtkImageData::New();

    m_transform = vtkTransform::New();
}

//    vtkImageFFT* fft=vtkImageFFT::New();
//    fft->SetInputConnection(is3->GetOutputPort());
//    fft->Update();

//    cout<<"fft calculated"<<endl;

//    //vtkImageButterworthLowPass* lowPass=vtkImageButterworthLowPass::New();
//    vtkImageIdealLowPass* lowPass=vtkImageIdealLowPass::New();
//    //lowPass->SetInputConnection(highPass->GetOutputPort());
//    lowPass->SetInputConnection(fft->GetOutputPort());
//    //lowPass->SetInputConnection(im->GetOutputPort());
//    lowPass->SetCutOff(0.04);
//    lowPass->ReleaseDataFlagOff();
//    lowPass->Update();



//    vtkImageRFFT* rfft=vtkImageRFFT::New();
//    rfft->SetInputConnection(lowPass->GetOutputPort());
//    rfft->Update();

vtkPolyData* BrainMeshing::computeInflatedPolyData(){


    vtkImageShrink3D* is3 = vtkImageShrink3D::New();
    is3->SetInputConnection(m_inputImage->GetProducerPort());
    is3->SetShrinkFactors(m_shrinkFactor,m_shrinkFactor,m_shrinkFactor);
    is3->Update();


    if(m_probeRadius == 0){
        m_modifiedImage = is3->GetOutput();
    }
    else{

        //Dilatation
              vtkImageContinuousDilate3D *dilate = vtkImageContinuousDilate3D::New();
       //3DImage is the image of interest
    //          dilate->SetInputConnection(m_inputImage->GetProducerPort());
              dilate->SetInputConnection(is3->GetOutputPort());
       //sk is the size of the dilation kernel in each direction
              double sk = m_probeRadius;

              dilate->SetKernelSize(sk,sk,sk);
              dilate->UpdateWholeExtent();


              //Erosion
              vtkImageContinuousErode3D *erode = vtkImageContinuousErode3D::New();
              erode->SetInputConnection(dilate->GetOutputPort());
              erode->SetKernelSize(sk,sk,sk);
              erode->UpdateWholeExtent();

              m_modifiedImage ->DeepCopy( erode->GetOutput());

    }


          //Marching cubes : 3D image -> 3D mesh

 //Tg is the value for the isosurface
          double Tg = m_threshold;

          vtkMarchingContourFilter* imc3=vtkMarchingContourFilter::New();
//          imc3->SetInputConnection(m_modifiedImage->GetProducerPort());
          imc3->SetInputConnection(m_modifiedImage->GetProducerPort());
          imc3->SetValue(0,Tg);
          imc3->UpdateWholeExtent();

//          dilate->Delete();
//          erode->Delete();

          //Connectivity filters (in order to clean small residues and keep only 1 connected surface)

          vtkPolyDataConnectivityFilter* con3 =vtkPolyDataConnectivityFilter::New();
          con3->SetInput(imc3->GetOutput());
          con3->SetExtractionModeToLargestRegion();
          con3->UpdateWholeExtent();
          //con3->Update();

          imc3->Delete();

          //Cleaners (remove double points and connections)

          vtkCleanPolyData* cleaner3 =vtkCleanPolyData::New();
          cleaner3->SetInput(con3->GetOutput());
          cleaner3->UpdateWholeExtent();

          con3->Delete();

          vtkTriangleFilter* tFilter = vtkTriangleFilter::New();
          tFilter->SetInput(cleaner3->GetOutput());
          tFilter->Update();

          cleaner3->Delete();

          vtkTransformPolyDataFilter* transf = vtkTransformPolyDataFilter::New();
          transf->SetInput(tFilter->GetOutput());
          transf->SetTransform(m_transform);
          transf->Update();

          //Polydatawriter (no need of it if there is an output to the function)

//          sprintf(FileNameOut3,"%s-SES.vtk",argv[1]);     // argv[1] = filename without extension ->vtk

   // may need to be recentered. The method for proteins may not be applied in this case.
//          vtkPolyDataWriter * writer3 = vtkPolyDataWriter::New();
//          writer3->SetFileName(FileNameOut3);

//          writer3->SetInput(cleaner3->GetOutput());
//          writer3->UpdateWholeExtent();
//          writer3->Write();

//          writer3->Delete();

          return transf->GetOutput();

}

void BrainMeshing::SetProbeRadius(double pr)
{
    m_probeRadius = pr;
}

void BrainMeshing::SetShrinkFactor(double sf)
{
    m_shrinkFactor = sf;
}

void BrainMeshing::SetThreshold(double th)
{
    m_threshold = th;
}

vtkImageData* BrainMeshing::GetModifiedImage()
{
    return m_modifiedImage;
}

void BrainMeshing::SetTransform(vtkTransform *transform)
{
    m_transform = transform;
}
