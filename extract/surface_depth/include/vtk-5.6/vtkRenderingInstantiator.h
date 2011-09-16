#ifndef __vtkRenderingInstantiator_h
#define __vtkRenderingInstantiator_h

#include "vtkInstantiator.h"



class VTK_RENDERING_EXPORT vtkRenderingInstantiator
{
  public:
  vtkRenderingInstantiator();
  ~vtkRenderingInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkRenderingInstantiator vtkRenderingInstantiatorInitializer;

#endif
