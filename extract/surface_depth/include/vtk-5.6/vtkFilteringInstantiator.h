#ifndef __vtkFilteringInstantiator_h
#define __vtkFilteringInstantiator_h

#include "vtkInstantiator.h"



class VTK_FILTERING_EXPORT vtkFilteringInstantiator
{
  public:
  vtkFilteringInstantiator();
  ~vtkFilteringInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkFilteringInstantiator vtkFilteringInstantiatorInitializer;

#endif
