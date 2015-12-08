# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# This project is designed to be built outside the Insight source tree.
PROJECT(travel_depth)

# Find VTK.
FIND_PACKAGE(VTK REQUIRED NO_MODULE)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(
TravelDepthMain TravelDepthMain.cpp
)
INSTALL(TARGETS TravelDepthMain
    RUNTIME DESTINATION /bin)

TARGET_LINK_LIBRARIES(
TravelDepthMain
TravelDepth
MeshAnalyser
${VTK_LIBRARIES}

# VTK6 migration: www.vtk.org/Wiki/VTK/Build_System_Migration
# OLD:
#vtkCommon
#vtkRendering
#vtkHybrid
# NEW:
#vtkCommonCore
#vtkCommonSystem
#vtkCommonColor
#vtkCommonDataModel
#vtkCommonMath
#vtkCommonMisc
#vtkCommonTransforms
#vtkCommonExecutionModel
#vtkCommonComputationalGeometry
#vtkRenderingCore
#vtkRenderingContext2D
#vtkRenderingFreeType
#vtkRenderingFreeTypeOpenGL
#vtkRenderingHybridOpenGL
#vtkRenderingImage
#vtkRenderingLabel
#vtkRenderingLOD
#vtkRenderingOpenGL
#vtkRenderingAnnotation
#vtkRenderingVolume
#vtkRenderingVolumeAMR
#vtkRenderingVolumeOpenGL
#vtkRenderingGL2PS
#vtkImagingHybrid
#vtkFiltersHybrid
)
