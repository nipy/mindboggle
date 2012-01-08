# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
PROJECT(brain)

# Find VTK.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(meshTest
                meshTest.cpp
                MeshAnalyser.h
                MeshAnalyser.cpp
                BrainMeshing.h
                BrainMeshing.cpp
                FsSurfaceReader.h
                FsSurfaceReader.cpp
VtkFileEditor.h
VtkFileEditor.cpp)

TARGET_LINK_LIBRARIES(meshTest vtkCommon vtkRendering vtkHybrid)
