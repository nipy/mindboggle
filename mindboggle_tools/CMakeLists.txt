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

<<<<<<< HEAD
ADD_SUBDIRECTORY(
gradient
)

ADD_SUBDIRECTORY(
geodesic_depth
)

ADD_SUBDIRECTORY(
medial_surfaces
)

=======
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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
<<<<<<< HEAD
GradientComputer GradientComputer.cpp
)

ADD_LIBRARY(
=======
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
PointAreaComputer PointAreaComputer.cpp
)

TARGET_LINK_LIBRARIES(
MeshAnalyser
)

TARGET_LINK_LIBRARIES(
TravelDepth
)

TARGET_LINK_LIBRARIES(
<<<<<<< HEAD
GradientComputer
)

TARGET_LINK_LIBRARIES(
=======
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
Overlap
PointAreaComputer
)
