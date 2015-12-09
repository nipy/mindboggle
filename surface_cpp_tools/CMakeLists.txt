# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

#https://forums.developer.apple.com/thread/17334
SET(CMAKE_OSX_DEPLOYMENT_TARGET "")
SET(CMAKE_OSX_SYSROOT /)
# https://cmake.org/Bug/view.php?id=14695
#SET(MACOSX_DEPLOYMENT_TARGET ${DARWIN_MAJOR_SDK_VERSION}.${DARWIN_MINOR_SDK_VERSION})
#SET(CMAKE_CXX_FLAGS "--sysroot ${CMAKE_OSX_SYSROOT} ${CMAKE_CXX_FLAGS}")
#MESSAGE("Setting MACOSX_DEPLOYMENT_TARGET to '${MACOSX_DEPLOYMENT_TARGET}'.")

PROJECT(brain)

# Find VTK.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(
travel_depth
)

ADD_SUBDIRECTORY(
curvature
)

ADD_SUBDIRECTORY(
surface_overlap
)

ADD_SUBDIRECTORY(
area
)

ADD_SUBDIRECTORY(
gradient
)

ADD_SUBDIRECTORY(
geodesic_depth
)

#ADD_SUBDIRECTORY(
#medial_surfaces
#)

ADD_LIBRARY(
MeshAnalyser MeshAnalyser.cpp
)
INSTALL(TARGETS MeshAnalyser
        LIBRARY DESTINATION /lib
        ARCHIVE DESTINATION /lib)

ADD_LIBRARY(
FsSurfaceReader FsSurfaceReader.cpp
)

ADD_LIBRARY(
VtkFileEditor VtkFileEditor.cpp
)

ADD_LIBRARY(
TravelDepth TravelDepth.cpp
)

ADD_LIBRARY(
Overlap Overlap.cpp
)

ADD_LIBRARY(
GradientComputer GradientComputer.cpp
)

ADD_LIBRARY(
PointAreaComputer PointAreaComputer.cpp
)

TARGET_LINK_LIBRARIES(
MeshAnalyser
)

TARGET_LINK_LIBRARIES(
TravelDepth
)

TARGET_LINK_LIBRARIES(
GradientComputer
)

TARGET_LINK_LIBRARIES(
Overlap
PointAreaComputer
)
