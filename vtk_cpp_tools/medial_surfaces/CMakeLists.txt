cmake_minimum_required(VERSION 3.0)

# Find VTK:
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# This project is designed to be built outside the Insight source tree:
project(medial_surface)

add_library(
MedialSurfaceComputer MedialSurfaceComputer.cpp
)

add_executable(
MedialSurfaceMain MedialSurfaceMain.cpp
)
install(TARGETS MedialSurfaceMain
    RUNTIME DESTINATION /bin)

target_link_libraries(
MedialSurfaceMain
MedialSurfaceComputer
${VTK_LIBRARIES}
)
