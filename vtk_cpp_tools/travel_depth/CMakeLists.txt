cmake_minimum_required(VERSION 3.0)

# Find VTK:
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# This project is designed to be built outside the Insight source tree:
project(travel_depth)

add_executable(
TravelDepthMain TravelDepthMain.cpp
)
install(TARGETS TravelDepthMain
    RUNTIME DESTINATION /bin)

target_link_libraries(
TravelDepthMain
TravelDepth
MeshAnalyser
${VTK_LIBRARIES}
)
