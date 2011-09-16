#ifndef __vtkViewsInstantiator_h
#define __vtkViewsInstantiator_h

#include "vtkInstantiator.h"



class VTK_VIEWS_EXPORT vtkViewsInstantiator
{
  public:
  vtkViewsInstantiator();
  ~vtkViewsInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkViewsInstantiator vtkViewsInstantiatorInitializer;

#endif
