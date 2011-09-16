#ifndef __vtkWidgetsInstantiator_h
#define __vtkWidgetsInstantiator_h

#include "vtkInstantiator.h"



class VTK_WIDGETS_EXPORT vtkWidgetsInstantiator
{
  public:
  vtkWidgetsInstantiator();
  ~vtkWidgetsInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkWidgetsInstantiator vtkWidgetsInstantiatorInitializer;

#endif
