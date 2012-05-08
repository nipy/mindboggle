# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(brain)

# Find VTK.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})


ADD_SUBDIRECTORY(
travel_depth
medial_surface
)

ADD_LIBRARY(
MeshAnalyser MeshAnalyser.cpp
VtkFileEditor VtkFileEditor.cpp
FsSurfaceReader FsSurfaceReader.cpp
MedialSurfaceComputer MedialSurfaceComputer.cpp
)

