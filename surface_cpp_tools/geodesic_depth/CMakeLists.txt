# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
PROJECT(geodesic_depth)

# Find VTK.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(
GeodesicDepthMain GeodesicDepthMain.cpp
)
INSTALL(TARGETS GeodesicDepthMain
    RUNTIME DESTINATION /bin)

TARGET_LINK_LIBRARIES(
GeodesicDepthMain
MeshAnalyser
vtkCommon
vtkRendering
vtkHybrid
)
