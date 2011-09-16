#ifndef __vtkGeovisInstantiator_h
#define __vtkGeovisInstantiator_h

#include "vtkInstantiator.h"



class VTK_GEOVIS_EXPORT vtkGeovisInstantiator
{
  public:
  vtkGeovisInstantiator();
  ~vtkGeovisInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkGeovisInstantiator vtkGeovisInstantiatorInitializer;

#endif
