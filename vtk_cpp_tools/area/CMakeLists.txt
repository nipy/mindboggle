cmake_minimum_required(VERSION 3.0)

# Find VTK:
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# This project is designed to be built outside the Insight source tree:
project(point_area)

add_executable(
PointAreaMain PointAreaMain.cpp
)
install(TARGETS PointAreaMain
    RUNTIME DESTINATION /bin)

target_link_libraries(
PointAreaMain
MeshAnalyser
PointAreaComputer
${VTK_LIBRARIES}
)
