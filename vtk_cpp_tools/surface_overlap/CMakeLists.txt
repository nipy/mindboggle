cmake_minimum_required(VERSION 3.0)

# Find VTK:
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# This project is designed to be built outside the Insight source tree:
project(surface_overlap)

add_executable(
SurfaceOverlapMain SurfaceOverlapMain.cpp
)
install(TARGETS SurfaceOverlapMain
    RUNTIME DESTINATION /bin)

target_link_libraries(
SurfaceOverlapMain
Overlap
${VTK_LIBRARIES}
)
