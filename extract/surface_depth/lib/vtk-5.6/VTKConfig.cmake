#-----------------------------------------------------------------------------
#
# VTKConfig.cmake - VTK CMake configuration file for external projects.
#
# This file is configured by VTK and used by the UseVTK.cmake module
# to load VTK's settings for an external project.

# Compute the installation prefix from VTK_DIR.
SET(VTK_INSTALL_PREFIX "${VTK_DIR}")
GET_FILENAME_COMPONENT(VTK_INSTALL_PREFIX "${VTK_INSTALL_PREFIX}" PATH)
GET_FILENAME_COMPONENT(VTK_INSTALL_PREFIX "${VTK_INSTALL_PREFIX}" PATH)


# The VTK include file directories.
SET(VTK_INCLUDE_DIRS "${VTK_INSTALL_PREFIX}/include/vtk-5.6;C:/msys/1.0/boost_1_47_0/boost_1_47_0")

# The VTK library directories.
SET(VTK_LIBRARY_DIRS "${VTK_INSTALL_PREFIX}/lib/vtk-5.6")

# The VTK binary executable directories.  Note that if
# VTK_CONFIGURATION_TYPES is set (see below) then these directories
# will be the parent directories under which there will be a directory
# of runtime binaries for each configuration type.
SET(VTK_EXECUTABLE_DIRS "${VTK_INSTALL_PREFIX}/bin")

# The VTK runtime library directories.  Note that if
# VTK_CONFIGURATION_TYPES is set (see below) then these directories
# will be the parent directories under which there will be a directory
# of runtime libraries for each configuration type.
SET(VTK_RUNTIME_LIBRARY_DIRS "${VTK_INSTALL_PREFIX}/bin")

# The runtime library path variable name e.g. LD_LIBRARY_PATH,
# this environment variable should be set to VTK_RUNTIME_LIBRARY_DIRS
SET(VTK_RUNTIME_PATH_VAR_NAME "")

# The C and C++ flags added by VTK to the cmake-configured flags.
SET(VTK_REQUIRED_C_FLAGS " -mwin32 -mthreads")
SET(VTK_REQUIRED_CXX_FLAGS " -Wno-deprecated -mwin32 -mthreads")
SET(VTK_REQUIRED_EXE_LINKER_FLAGS " -mthreads")
SET(VTK_REQUIRED_SHARED_LINKER_FLAGS " -mthreads")
SET(VTK_REQUIRED_MODULE_LINKER_FLAGS " -mthreads")

# The VTK version number
SET(VTK_MAJOR_VERSION "5")
SET(VTK_MINOR_VERSION "6")
SET(VTK_BUILD_VERSION "1")

# The location of the UseVTK.cmake file.
SET(VTK_USE_FILE "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/UseVTK.cmake")

# The build settings file.
SET(VTK_BUILD_SETTINGS_FILE "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/VTKBuildSettings.cmake")

# The directory containing class list files for each kit.
SET(VTK_KITS_DIR "${VTK_INSTALL_PREFIX}/lib/vtk-5.6")

# The wrapping hints file.
SET(VTK_WRAP_HINTS "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/hints")

# CMake extension module directory and macro file.
SET(VTK_LOAD_CMAKE_EXTENSIONS_MACRO "C:/VTK/vtk-5.6.1-src/CMake/vtkLoadCMakeExtensions.cmake")
SET(VTK_CMAKE_DIR "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/CMake")
SET(VTK_CMAKE_EXTENSIONS_DIR "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/CMake")

# The list of available kits.
SET(VTK_KITS "COMMON;FILTERING;IO;GRAPHICS;GENERIC_FILTERING;IMAGING;RENDERING;VOLUMERENDERING;HYBRID;WIDGETS;INFOVIS;GEOVIS;VIEWS;CHARTS")

# The list of available languages.
SET(VTK_LANGUAGES "")

# VTK Configuration options.
SET(VTK_BUILD_SHARED_LIBS "OFF")
SET(VTK_DEBUG_LEAKS "OFF")
SET(VTK_HAVE_VP1000 "")
SET(VTK_LEGACY_REMOVE "OFF")
SET(VTK_LEGACY_SILENT "OFF")
SET(VTK_MPIRUN_EXE "")
SET(VTK_MPI_CLIENT_POSTFLAGS "")
SET(VTK_MPI_CLIENT_PREFLAGS "")
SET(VTK_MPI_MAX_NUMPROCS "")
SET(VTK_MPI_PRENUMPROC_FLAGS "")
SET(VTK_MPI_NUMPROC_FLAG "")
SET(VTK_MPI_POSTFLAGS "")
SET(VTK_MPI_PREFLAGS "")
SET(VTK_MPI_SERVER_POSTFLAGS "")
SET(VTK_MPI_SERVER_PREFLAGS "")
SET(VTK_OPENGL_HAS_OSMESA "OFF")
SET(VTK_USE_64BIT_IDS "OFF")
SET(VTK_USE_ANSI_STDLIB "ON")
SET(VTK_USE_BOOST "ON")
SET(VTK_USE_QT "OFF")
SET(VTK_USE_CARBON "OFF")
SET(VTK_USE_CG_SHADERS "OFF")
SET(VTK_USE_CHARTS "ON")
SET(VTK_USE_COCOA "OFF")
SET(VTK_USE_GEOVIS "ON")
SET(VTK_USE_GL2PS "OFF")
SET(VTK_USE_GLSL_SHADERS "ON")
SET(VTK_USE_GNU_R "OFF")
SET(VTK_USE_GUISUPPORT "OFF")
SET(VTK_USE_INFOVIS "ON")
SET(VTK_USE_MANGLED_MESA "OFF")
SET(VTK_USE_MATLAB_MEX "OFF")
SET(VTK_USE_MATROX_IMAGING "OFF")
SET(VTK_USE_METAIO "ON")
SET(VTK_USE_MFC "")
SET(VTK_USE_MPI "OFF")
SET(VTK_USE_MYSQL "OFF")
SET(VTK_USE_N_WAY_ARRAYS "ON")
SET(VTK_USE_ODBC "OFF")
SET(VTK_USE_PARALLEL "OFF")
SET(VTK_USE_PARALLEL_BGL "OFF")
SET(VTK_USE_POSTGRES "OFF") 
SET(VTK_USE_QVTK "OFF")
SET(VTK_USE_RENDERING "ON")
SET(VTK_USE_TDX "OFF")
SET(VTK_USE_TEXT_ANALYSIS "OFF")
SET(VTK_USE_TK "OFF")
SET(VTK_USE_TK "OFF")
SET(VTK_USE_VIDEO_FOR_WINDOWS "OFF")
SET(VTK_USE_VIEWS "ON")
SET(VTK_USE_VOLUMEPRO_1000 "OFF")
SET(VTK_USE_X "0")
SET(VTK_WRAP_JAVA "OFF")
SET(VTK_WRAP_PYTHON "OFF")
SET(VTK_WRAP_TCL "OFF")

# The Hybrid and VolumeRendering kits are now switched with Rendering.
SET(VTK_USE_HYBRID "ON")
SET(VTK_USE_VOLUMERENDERING "ON")

# The Tcl/Tk configuration.
SET(VTK_TCL_TK_STATIC "")
SET(VTK_TCL_TK_COPY_SUPPORT_LIBRARY "")
SET(VTK_TCL_SUPPORT_LIBRARY_PATH "")
SET(VTK_TK_SUPPORT_LIBRARY_PATH "")
SET(VTK_TCL_TK_MACROS_MODULE "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/CMake/vtkTclTkMacros.cmake")
SET(VTK_TCL_HOME "")
SET(VTK_WRAP_TCL_EXE "")
SET(VTK_WRAP_TCL_INIT_EXE "")
SET(VTK_TK_INTERNAL_DIR "")
SET(VTK_TK_RESOURCES_DIR "")

# The Java configuration.
SET(VTK_JAVA_JAR "")
SET(VTK_PARSE_JAVA_EXE "")
SET(VTK_WRAP_JAVA_EXE "")

# The Matlab configuration.
SET(VTK_MATLAB_ROOT_DIR "")
SET(VTK_MATLAB_INCLUDE_DIR "")
SET(VTK_MATLAB_LIB_DIR "")

# The Python configuration.
# If VTK_CONFIGURATION_TYPES is set (see below) then the VTK_PYTHONPATH_DIRS
# will have subdirectories for each configuration type.
SET(VTK_PYTHONPATH_DIRS "")
SET(VTK_WRAP_PYTHON_EXE "")
SET(VTK_WRAP_PYTHON_INIT_EXE "")

# Other executables
SET(VTK_ENCODESTRING_EXE "${VTK_INSTALL_PREFIX}/bin/vtkEncodeString.exe")

# The Doxygen configuration.
SET(VTK_DOXYGEN_HOME "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/doxygen")

# The VTK test script locations.
SET(VTK_HEADER_TESTING_PY "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/testing/HeaderTesting.py")
SET(VTK_FIND_STRING_TCL "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/testing/FindString.tcl")
SET(VTK_PRINT_SELF_CHECK_TCL "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/testing/PrintSelfCheck.tcl")
SET(VTK_RT_IMAGE_TEST_TCL "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/testing/rtImageTest.tcl")
SET(VTK_PRT_IMAGE_TEST_TCL "")

# The names of utility libraries used by VTK.
SET(VTK_PNG_LIBRARIES      "vtkpng")
SET(VTK_ZLIB_LIBRARIES     "vtkzlib")
SET(VTK_JPEG_LIBRARIES     "vtkjpeg")
SET(VTK_TIFF_LIBRARIES     "vtktiff")
SET(VTK_EXPAT_LIBRARIES    "vtkexpat")
SET(VTK_FREETYPE_LIBRARIES "vtkfreetype")
SET(VTK_LIBXML2_LIBRARIES  "vtklibxml2")
SET(VTK_LIBPROJ4_LIBRARIES  "vtkproj4")

# The VTK Qt configuration.
IF(VTK_USE_QVTK)
  INCLUDE(${VTK_DIR}/VTKConfigQt.cmake)
ENDIF(VTK_USE_QVTK)

# Relative install paths in the VTK install tree
SET(VTK_INSTALL_BIN_DIR "/bin")
SET(VTK_INSTALL_INCLUDE_DIR "/include/vtk-5.6")
SET(VTK_INSTALL_LIB_DIR "/lib/vtk-5.6")
SET(VTK_INSTALL_PACKAGE_DIR "/lib/vtk-5.6")

# A VTK install tree always provides one build configuration.  A VTK
# build tree may provide either one or multiple build configurations
# depending on the CMake generator used.  Since VTK can be used either
# from a build tree or an install tree it is useful for outside
# projects to know the configurations available.  If this
# VTKConfig.cmake is in a VTK install tree VTK_CONFIGURATION_TYPES
# will be empty and VTK_BUILD_TYPE will be set to the value of
# CMAKE_BUILD_TYPE used to build VTK.  If VTKConfig.cmake is in a VTK
# build tree then VTK_CONFIGURATION_TYPES and VTK_BUILD_TYPE will have
# values matching CMAKE_CONFIGURATION_TYPES and CMAKE_BUILD_TYPE for
# that build tree (only one will ever be set).
SET(VTK_CONFIGURATION_TYPES )
SET(VTK_BUILD_TYPE )

# The VTK library dependencies.
IF(NOT VTK_NO_LIBRARY_DEPENDS AND EXISTS "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/VTKLibraryDepends.cmake")
  INCLUDE("${VTK_INSTALL_PREFIX}/lib/vtk-5.6/VTKLibraryDepends.cmake")
ENDIF(NOT VTK_NO_LIBRARY_DEPENDS AND EXISTS "${VTK_INSTALL_PREFIX}/lib/vtk-5.6/VTKLibraryDepends.cmake")

# The old, less clear name for VTK_RUNTIME_LIBRARY_DIRS.  Kept here
# for compatibility.
SET(VTK_RUNTIME_DIRS ${VTK_RUNTIME_LIBRARY_DIRS})

# The name of the encoders used by VTK
SET( VTK_USE_FFMPEG_ENCODER "OFF" )
SET( VTK_USE_OGGTHEORA_ENCODER "OFF")


