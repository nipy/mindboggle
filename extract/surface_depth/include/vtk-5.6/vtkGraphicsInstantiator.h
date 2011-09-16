#ifndef __vtkGraphicsInstantiator_h
#define __vtkGraphicsInstantiator_h

#include "vtkInstantiator.h"



class VTK_GRAPHICS_EXPORT vtkGraphicsInstantiator
{
  public:
  vtkGraphicsInstantiator();
  ~vtkGraphicsInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkGraphicsInstantiator vtkGraphicsInstantiatorInitializer;

#endif
