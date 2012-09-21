# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
PROJECT(surface_overlap)

# Find VTK.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(
SurfaceOverlapMain SurfaceOverlapMain.cpp
)
INSTALL(TARGETS SurfaceOverlapMain
    RUNTIME DESTINATION /bin)

TARGET_LINK_LIBRARIES(
SurfaceOverlapMain
Overlap
vtkCommon
vtkRendering
vtkHybrid
)
