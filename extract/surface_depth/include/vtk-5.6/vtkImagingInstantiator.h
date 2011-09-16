#ifndef __vtkImagingInstantiator_h
#define __vtkImagingInstantiator_h

#include "vtkInstantiator.h"



class VTK_IMAGING_EXPORT vtkImagingInstantiator
{
  public:
  vtkImagingInstantiator();
  ~vtkImagingInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkImagingInstantiator vtkImagingInstantiatorInitializer;

#endif
