cmake_minimum_required(VERSION 3.0)

# Find VTK:
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# This project is designed to be built outside the Insight source tree.
project(curvature)

add_executable(
CurvatureMain CurvatureMain.cpp
)
install(TARGETS CurvatureMain
    RUNTIME DESTINATION /bin)

target_link_libraries(
CurvatureMain
MeshAnalyser
${VTK_LIBRARIES}
)
