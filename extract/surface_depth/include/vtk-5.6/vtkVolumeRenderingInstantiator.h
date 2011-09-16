#ifndef __vtkVolumeRenderingInstantiator_h
#define __vtkVolumeRenderingInstantiator_h

#include "vtkInstantiator.h"



class VTK_VOLUMERENDERING_EXPORT vtkVolumeRenderingInstantiator
{
  public:
  vtkVolumeRenderingInstantiator();
  ~vtkVolumeRenderingInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkVolumeRenderingInstantiator vtkVolumeRenderingInstantiatorInitializer;

#endif
