#ifndef __vtkHybridInstantiator_h
#define __vtkHybridInstantiator_h

#include "vtkInstantiator.h"



class VTK_HYBRID_EXPORT vtkHybridInstantiator
{
  public:
  vtkHybridInstantiator();
  ~vtkHybridInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkHybridInstantiator vtkHybridInstantiatorInitializer;

#endif
