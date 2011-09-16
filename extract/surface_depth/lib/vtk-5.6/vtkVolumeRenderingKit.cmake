# Directory containing class headers.
SET(VTK_VOLUMERENDERING_HEADER_DIR "${VTK_INSTALL_PREFIX}/include/vtk-5.6")

# Classes in vtkVolumeRendering.
SET(VTK_VOLUMERENDERING_CLASSES
  "vtkDirectionEncoder"
  "vtkEncodedGradientEstimator"
  "vtkEncodedGradientShader"
  "vtkFiniteDifferenceGradientEstimator"
  "vtkFixedPointRayCastImage"
  "vtkFixedPointVolumeRayCastCompositeGOHelper"
  "vtkFixedPointVolumeRayCastCompositeGOShadeHelper"
  "vtkFixedPointVolumeRayCastCompositeHelper"
  "vtkFixedPointVolumeRayCastCompositeShadeHelper"
  "vtkFixedPointVolumeRayCastHelper"
  "vtkFixedPointVolumeRayCastMIPHelper"
  "vtkFixedPointVolumeRayCastMapper"
  "vtkGPUVolumeRayCastMapper"
  "vtkHAVSVolumeMapper"
  "vtkProjectedTetrahedraMapper"
  "vtkRayCastImageDisplayHelper"
  "vtkRecursiveSphereDirectionEncoder"
  "vtkSphericalDirectionEncoder"
  "vtkVolumeMapper"
  "vtkVolumeOutlineSource"
  "vtkVolumePicker"
  "vtkVolumeProMapper"
  "vtkVolumeRayCastCompositeFunction"
  "vtkVolumeRayCastFunction"
  "vtkVolumeRayCastIsosurfaceFunction"
  "vtkVolumeRayCastMIPFunction"
  "vtkVolumeRayCastMapper"
  "vtkVolumeRenderingFactory"
  "vtkVolumeTextureMapper"
  "vtkVolumeTextureMapper2D"
  "vtkVolumeTextureMapper3D"
  "vtkUnstructuredGridBunykRayCastFunction"
  "vtkUnstructuredGridHomogeneousRayIntegrator"
  "vtkUnstructuredGridLinearRayIntegrator"
  "vtkUnstructuredGridPartialPreIntegration"
  "vtkUnstructuredGridPreIntegration"
  "vtkUnstructuredGridVolumeMapper"
  "vtkUnstructuredGridVolumeRayCastFunction"
  "vtkUnstructuredGridVolumeRayCastIterator"
  "vtkUnstructuredGridVolumeRayIntegrator"
  "vtkUnstructuredGridVolumeRayCastMapper"
  "vtkUnstructuredGridVolumeZSweepMapper"
  "vtkOpenGLGPUVolumeRayCastMapper"
  "vtkOpenGLHAVSVolumeMapper"
  "vtkOpenGLProjectedTetrahedraMapper"
  "vtkOpenGLRayCastImageDisplayHelper"
  "vtkOpenGLVolumeTextureMapper2D"
  "vtkOpenGLVolumeTextureMapper3D")

# Abstract classes in vtkVolumeRendering.
SET(VTK_VOLUMERENDERING_CLASSES_ABSTRACT
  "vtkDirectionEncoder"
  "vtkEncodedGradientEstimator"
  "vtkFixedPointVolumeRayCastHelper"
  "vtkRayCastImageDisplayHelper"
  "vtkVolumeMapper"
  "vtkVolumeRayCastFunction"
  "vtkVolumeTextureMapper"
  "vtkUnstructuredGridVolumeMapper"
  "vtkUnstructuredGridVolumeRayCastFunction"
  "vtkUnstructuredGridVolumeRayCastIterator"
  "vtkUnstructuredGridVolumeRayIntegrator")

# Wrap-exclude classes in vtkVolumeRendering.
SET(VTK_VOLUMERENDERING_CLASSES_WRAP_EXCLUDE)

# Set convenient variables to test each class.
FOREACH(class ${VTK_VOLUMERENDERING_CLASSES})
  SET(VTK_CLASS_EXISTS_${class} 1)
ENDFOREACH(class)
FOREACH(class ${VTK_VOLUMERENDERING_CLASSES_ABSTRACT})
  SET(VTK_CLASS_ABSTRACT_${class} 1)
ENDFOREACH(class)
FOREACH(class ${VTK_VOLUMERENDERING_CLASSES_WRAP_EXCLUDE})
  SET(VTK_CLASS_WRAP_EXCLUDE_${class} 1)
ENDFOREACH(class)
