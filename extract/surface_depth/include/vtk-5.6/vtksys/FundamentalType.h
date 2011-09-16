/*============================================================================
  KWSys - Kitware System Library
  Copyright 2000-2009 Kitware, Inc., Insight Software Consortium

  Distributed under the OSI-approved BSD License (the "License");
  see accompanying file Copyright.txt for details.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the License for more information.
============================================================================*/
#ifndef vtksys_FundamentalType_h
#define vtksys_FundamentalType_h

#include <vtksys/Configure.h>

/* Redefine all public interface symbol names to be in the proper
   namespace.  These macros are used internally to kwsys only, and are
   not visible to user code.  Use kwsysHeaderDump.pl to reproduce
   these macros after making changes to the interface.  */
#if !defined(KWSYS_NAMESPACE)
# define kwsys_ns(x) vtksys##x
# define kwsysEXPORT vtksys_EXPORT
#endif

#if !vtksys_NAME_IS_KWSYS
# define kwsysFundamentalType        kwsys_ns(FundamentalType)
# define kwsysFundamentalType_Int8   kwsys_ns(FundamentalType_Int8)
# define kwsysFundamentalType_UInt8  kwsys_ns(FundamentalType_UInt8)
# define kwsysFundamentalType_Int16  kwsys_ns(FundamentalType_Int16)
# define kwsysFundamentalType_UInt16 kwsys_ns(FundamentalType_UInt16)
# define kwsysFundamentalType_Int32  kwsys_ns(FundamentalType_Int32)
# define kwsysFundamentalType_UInt32 kwsys_ns(FundamentalType_UInt32)
# define kwsysFundamentalType_Int64  kwsys_ns(FundamentalType_Int64)
# define kwsysFundamentalType_UInt64 kwsys_ns(FundamentalType_UInt64)
#endif

/* The size of fundamental types.  Types that do not exist have size 0.  */
#define vtksys_SIZEOF_CHAR 1
#define vtksys_SIZEOF_SHORT 2
#define vtksys_SIZEOF_INT 4
#define vtksys_SIZEOF_LONG 4
#define vtksys_SIZEOF_LONG_LONG 8
#define vtksys_SIZEOF___INT64 8

/* Whether types "long long" and "__int64" are enabled.  If a type is
   enabled then it is a unique fundamental type.  */
#define vtksys_USE_LONG_LONG 1
#define vtksys_USE___INT64 1

/* Whether type "char" is signed (it may be signed or unsigned).  */
#define vtksys_CHAR_IS_SIGNED 1

#if defined(__cplusplus)
extern "C"
{
#endif

/* Select an 8-bit integer type.  */
#if vtksys_SIZEOF_CHAR == 1
typedef signed char kwsysFundamentalType_Int8;
typedef unsigned char kwsysFundamentalType_UInt8;
#else
# error "No native data type can represent an 8-bit integer."
#endif

/* Select a 16-bit integer type.  */
#if vtksys_SIZEOF_SHORT == 2
typedef short kwsysFundamentalType_Int16;
typedef unsigned short kwsysFundamentalType_UInt16;
#elif vtksys_SIZEOF_INT == 2
typedef int kwsysFundamentalType_Int16;
typedef unsigned int kwsysFundamentalType_UInt16;
#else
# error "No native data type can represent a 16-bit integer."
#endif

/* Select a 32-bit integer type.  */
#if vtksys_SIZEOF_INT == 4
typedef int kwsysFundamentalType_Int32;
typedef unsigned int kwsysFundamentalType_UInt32;
#elif vtksys_SIZEOF_LONG == 4
typedef long kwsysFundamentalType_Int32;
typedef unsigned long kwsysFundamentalType_UInt32;
#else
# error "No native data type can represent a 32-bit integer."
#endif

/* Select a 64-bit integer type.  */
#if vtksys_SIZEOF_LONG == 8
typedef signed long   kwsysFundamentalType_Int64;
typedef unsigned long kwsysFundamentalType_UInt64;
/* Whether UInt64 can be converted to double.  */
# define vtksys_CAN_CONVERT_UI64_TO_DOUBLE 1
#elif vtksys_USE_LONG_LONG && vtksys_SIZEOF_LONG_LONG == 8
typedef signed long long   kwsysFundamentalType_Int64;
typedef unsigned long long kwsysFundamentalType_UInt64;
/* Whether UInt64 can be converted to double.  */
# define vtksys_CAN_CONVERT_UI64_TO_DOUBLE 1
#elif vtksys_USE___INT64 && vtksys_SIZEOF___INT64 == 8
typedef signed __int64   kwsysFundamentalType_Int64;
typedef unsigned __int64 kwsysFundamentalType_UInt64;
/* Whether UInt64 can be converted to double.  */
# define vtksys_CAN_CONVERT_UI64_TO_DOUBLE 0
#else
# error "No native data type can represent a 64-bit integer."
#endif

#if defined(__cplusplus)
} /* extern "C" */
#endif

/* If we are building a kwsys .c or .cxx file, let it use these macros.
   Otherwise, undefine them to keep the namespace clean.  */
#if !defined(KWSYS_NAMESPACE)
# undef kwsys_ns
# undef kwsysEXPORT
# if !defined(KWSYS_NAMESPACE) && !vtksys_NAME_IS_KWSYS
#  undef kwsysFundamentalType
#  undef kwsysFundamentalType_Int8
#  undef kwsysFundamentalType_UInt8
#  undef kwsysFundamentalType_Int16
#  undef kwsysFundamentalType_UInt16
#  undef kwsysFundamentalType_Int32
#  undef kwsysFundamentalType_UInt32
#  undef kwsysFundamentalType_Int64
#  undef kwsysFundamentalType_UInt64
# endif
#endif

/* If building a C or C++ file in kwsys itself, give the source file
   access to the configured macros without a configured namespace.  */
#if defined(KWSYS_NAMESPACE)
# define KWSYS_SIZEOF_CHAR vtksys_SIZEOF_CHAR
# define KWSYS_SIZEOF_SHORT vtksys_SIZEOF_SHORT
# define KWSYS_SIZEOF_INT vtksys_SIZEOF_INT
# define KWSYS_SIZEOF_LONG vtksys_SIZEOF_LONG
# define KWSYS_SIZEOF_LONG_LONG vtksys_SIZEOF_LONG_LONG
# define KWSYS_SIZEOF___INT64 vtksys_SIZEOF___INT64
# define KWSYS_USE_LONG_LONG vtksys_USE_LONG_LONG
# define KWSYS_USE___INT64 vtksys_USE___INT64
# define KWSYS_CHAR_IS_SIGNED vtksys_CHAR_IS_SIGNED
# define KWSYS_CAN_CONVERT_UI64_TO_DOUBLE vtksys_CAN_CONVERT_UI64_TO_DOUBLE
#endif

#endif
