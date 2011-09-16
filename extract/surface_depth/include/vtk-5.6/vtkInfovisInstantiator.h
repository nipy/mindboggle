#ifndef __vtkInfovisInstantiator_h
#define __vtkInfovisInstantiator_h

#include "vtkInstantiator.h"



class VTK_INFOVIS_EXPORT vtkInfovisInstantiator
{
  public:
  vtkInfovisInstantiator();
  ~vtkInfovisInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkInfovisInstantiator vtkInfovisInstantiatorInitializer;

#endif
