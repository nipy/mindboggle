# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
PROJECT(medial_surface)

# Find VTK.
FIND_PACKAGE(VTK REQUIRED NO_MODULE)
INCLUDE(${VTK_USE_FILE})

ADD_LIBRARY(
MedialSurfaceComputer MedialSurfaceComputer.cpp
)

ADD_EXECUTABLE(
MedialSurfaceMain MedialSurfaceMain.cpp
)
INSTALL(TARGETS MedialSurfaceMain
    RUNTIME DESTINATION /bin)

TARGET_LINK_LIBRARIES(
MedialSurfaceMain
MedialSurfaceComputer
${VTK_LIBRARIES}
)

# VTK6 migration: www.vtk.org/Wiki/VTK/Build_System_Migration
# Not in VTK6:
#vtkCommon
#vtkRendering
#vtkHybrid
