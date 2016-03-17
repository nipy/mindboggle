cmake_minimum_required(VERSION 3.0)

# Find VTK:
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# This project is designed to be built outside the Insight source tree:
project(gradient)

add_executable(
GradientMain GradientMain.cpp
)
install(TARGETS GradientMain
    RUNTIME DESTINATION /bin)

target_link_libraries(
GradientMain
GradientComputer
MeshAnalyser
PointAreaComputer
${VTK_LIBRARIES}
)
