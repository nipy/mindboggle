/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkToolkits.h.in

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef __vtkToolkits_h
#define __vtkToolkits_h

/* This header is configured by VTK's build process.  */

/*--------------------------------------------------------------------------*/
/* Selected VTK Toolkits                                                    */

#define VTK_USE_VIEWS
#define VTK_USE_GEOVIS
#define VTK_USE_INFOVIS
/* #undef VTK_USE_TEXT_ANALYSIS */
#define VTK_USE_N_WAY_ARRAYS
/* #undef VTK_USE_PARALLEL */
#define VTK_USE_RENDERING
/* #undef VTK_USE_GL2PS */

/* The Hybrid and VolumeRendering kits are now switched with Rendering.  */
#ifdef VTK_USE_RENDERING
# define VTK_USE_HYBRID
# define VTK_USE_VOLUMERENDERING
#endif

/*--------------------------------------------------------------------------*/
/* Rendering Configuration                                                  */

/* #undef VTK_USE_X */
/* #undef VTK_USE_MANGLED_MESA */
#define VTK_USE_OPENGL_LIBRARY
/* #undef VTK_OPENGL_HAS_OSMESA */
/* #undef VTK_USE_OFFSCREEN */

/* #undef VTK_USE_CG_SHADERS */
#define VTK_USE_GLSL_SHADERS

#define VTK_MATERIALS_DIRS
#ifdef VTK_MATERIALS_DIRS
#  undef VTK_MATERIALS_DIRS
#  define VTK_MATERIALS_DIRS "C:/VTK/cmakebins2/Utilities/MaterialLibrary/Repository"
#endif

/* #undef VTK_USE_VOLUMEPRO_1000 */

/*--------------------------------------------------------------------------*/
/* Wrapping Configuration                                                   */

/* #undef VTK_WRAP_TCL */
/* #undef VTK_WRAP_PYTHON */
/* #undef VTK_WRAP_JAVA */

/*--------------------------------------------------------------------------*/
/* Other Configuration Options                                              */

/* Whether we are building MPI support.  */
/* #undef VTK_USE_MPI */

/* Should VTK use the display?  */
#define VTK_USE_DISPLAY

/* Is VTK_DATA_ROOT defined? */
/* #undef VTK_DATA_ROOT */
#ifdef VTK_DATA_ROOT
#  undef VTK_DATA_ROOT
#  define VTK_DATA_ROOT "VTK_DATA_ROOT-NOTFOUND"
#endif

/* Should VTK use PostgreSQL?  */
/* #undef VTK_USE_POSTGRES */

/* Should VTK use MySQL?  */
/* #undef VTK_USE_MYSQL */

/* Should VTK use ODBC?  */
/* #undef VTK_USE_ODBC */
/* #undef VTK_ODBC_TEST_DSN */
#ifdef VTK_ODBC_TEST_DSN
#  undef VTK_ODBC_TEST_DSN
#  define VTK_ODBC_TEST_DSN ""
#endif

/* Is a test PostgreSQL database URL defined? */
/* #undef VTK_PSQL_TEST_URL */
#ifdef VTK_PSQL_TEST_URL
#  undef VTK_PSQL_TEST_URL
#  define VTK_PSQL_TEST_URL ""
#endif

/* Is a test MySQL database URL defined? */
/* #undef VTK_MYSQL_TEST_URL */
#ifdef VTK_MYSQL_TEST_URL
#  undef VTK_MYSQL_TEST_URL
#  define VTK_MYSQL_TEST_URL ""
#endif

/* Debug leaks support.  */
/* #undef VTK_DEBUG_LEAKS */

/* Whether VTK is using its own utility libraries.  */
/* #undef VTK_USE_SYSTEM_PNG */
/* #undef VTK_USE_SYSTEM_ZLIB */
/* #undef VTK_USE_SYSTEM_JPEG */
/* #undef VTK_USE_SYSTEM_TIFF */
/* #undef VTK_USE_SYSTEM_EXPAT */
/* #undef VTK_USE_SYSTEM_FREETYPE */
/* #undef VTK_USE_SYSTEM_LIBXML2 */
/* #undef VTK_USE_SYSTEM_LIBPROJ4 */
/* #undef VTK_USE_SYSTEM_OGGTHEORA */

/* Whether VTK is using vfw32 and if it supports video capture */
/* #undef VTK_USE_VIDEO_FOR_WINDOWS */
/* #undef VTK_VFW_SUPPORTS_CAPTURE */

/* Whether FFMPEG is found or not  */
/* #undef VTK_USE_FFMPEG_ENCODER */

/* Whether the user has linked in the MPEG2 library or not  */
/* #undef VTK_USE_MPEG2_ENCODER */

/* Whether the user wants to use the Ogg/Theora encoder or not */
/* #undef VTK_USE_OGGTHEORA_ENCODER */

/* Whether the user has built the MetaIO library or not  */
#define VTK_USE_METAIO

/* Options for GPUInfo */
/* #undef VTK_USE_DIRECTX */
/* #undef VTK_USE_CORE_GRAPHICS */
/* #undef VTK_USE_NVCONTROL */

/*--------------------------------------------------------------------------*/
/* Setup VTK based on platform features and configuration.                  */

/* OGLR */
/* OGLR is for GLX. It can be on Unix, Mac or Windows */
#if ((defined(VTK_USE_OPENGL_LIBRARY) && defined(VTK_USE_X) && \
      !defined(_WIN32)) ||                                     \
     (defined(VTK_USE_X) && defined(_WIN32)) || \
     (defined(VTK_USE_MANGLED_MESA) && !defined(_WIN32))) && \
     !(defined (VTK_USE_CARBON) || defined(VTK_USE_COCOA))
# define VTK_USE_OGLR
#endif

#if (defined(VTK_OPENGL_HAS_OSMESA) && !defined(VTK_USE_OGLR) &&  \
    !defined(_WIN32) && !defined(VTK_USE_CARBON) && !defined(VTK_USE_COCOA))
# define VTK_USE_OSMESA
#endif

#endif
