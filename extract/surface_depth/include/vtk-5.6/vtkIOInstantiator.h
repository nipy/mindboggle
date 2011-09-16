#ifndef __vtkIOInstantiator_h
#define __vtkIOInstantiator_h

#include "vtkInstantiator.h"



class VTK_IO_EXPORT vtkIOInstantiator
{
  public:
  vtkIOInstantiator();
  ~vtkIOInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkIOInstantiator vtkIOInstantiatorInitializer;

#endif
