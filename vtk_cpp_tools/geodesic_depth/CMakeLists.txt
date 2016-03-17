cmake_minimum_required(VERSION 3.0)

# Find VTK:
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# This project is designed to be built outside the Insight source tree:
project(geodesic_depth)

add_executable(
GeodesicDepthMain GeodesicDepthMain.cpp
)
install(TARGETS GeodesicDepthMain
    RUNTIME DESTINATION /bin)

target_link_libraries(
GeodesicDepthMain
MeshAnalyser
${VTK_LIBRARIES}
)
