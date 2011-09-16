#ifndef __vtkCommonInstantiator_h
#define __vtkCommonInstantiator_h

#include "vtkInstantiator.h"



class VTK_COMMON_EXPORT vtkCommonInstantiator
{
  public:
  vtkCommonInstantiator();
  ~vtkCommonInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkCommonInstantiator vtkCommonInstantiatorInitializer;

#endif
