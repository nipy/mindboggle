#ifndef __vtkGenericFilteringInstantiator_h
#define __vtkGenericFilteringInstantiator_h

#include "vtkInstantiator.h"



class VTK_GENERIC_FILTERING_EXPORT vtkGenericFilteringInstantiator
{
  public:
  vtkGenericFilteringInstantiator();
  ~vtkGenericFilteringInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkGenericFilteringInstantiator vtkGenericFilteringInstantiatorInitializer;

#endif
