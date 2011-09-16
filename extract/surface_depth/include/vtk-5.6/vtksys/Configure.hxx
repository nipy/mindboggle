/*============================================================================
  KWSys - Kitware System Library
  Copyright 2000-2009 Kitware, Inc., Insight Software Consortium

  Distributed under the OSI-approved BSD License (the "License");
  see accompanying file Copyright.txt for details.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the License for more information.
============================================================================*/
#ifndef vtksys_Configure_hxx
#define vtksys_Configure_hxx

/* Include C configuration.  */
#include <vtksys/Configure.h>

/* Disable cleanup of putenv memory for issues with GCOV */
#if 0
#define KWSYS_DO_NOT_CLEAN_PUTENV 
#else
#undef KWSYS_DO_NOT_CLEAN_PUTENV
#endif



/* Whether ANSI C++ stream headers are to be used.  */
#define vtksys_IOS_USE_ANSI 1

/* Whether ANSI C++ streams are in std namespace.  */
#define vtksys_IOS_HAVE_STD 1

/* Whether ANSI C++ <sstream> header is to be used.  */
#define vtksys_IOS_USE_SSTREAM 1

/* Whether old C++ <strstream.h> header is to be used.  */
#define vtksys_IOS_USE_STRSTREAM_H 0

/* Whether old C++ <strstrea.h> header is to be used.  */
#define vtksys_IOS_USE_STRSTREA_H 0

/* Whether C++ streams support the ios::binary openmode.  */
#define vtksys_IOS_HAVE_BINARY 1

/* Whether STL is in std namespace.  */
#define vtksys_STL_HAVE_STD 1

/* Whether the STL string has operator<< for ostream.  */
#define vtksys_STL_STRING_HAVE_OSTREAM 1

/* Whether the STL string has operator>> for istream.  */
#define vtksys_STL_STRING_HAVE_ISTREAM 1

/* Whether the STL string has operator!= for char*.  */
#define vtksys_STL_STRING_HAVE_NEQ_CHAR 1

/* Define the stl namespace macro.  */
#if vtksys_STL_HAVE_STD
# define vtksys_stl std
#else
# define vtksys_stl
#endif

/* Define the ios namespace macro.  */
#if vtksys_IOS_HAVE_STD
# define vtksys_ios_namespace std
#else
# define vtksys_ios_namespace
#endif
#if vtksys_IOS_USE_SSTREAM
# define vtksys_ios vtksys_ios_namespace
#else
# define vtksys_ios vtksys_ios
#endif

/* Define the ios::binary openmode macro.  */
#if vtksys_IOS_HAVE_BINARY
# define vtksys_ios_binary vtksys_ios::ios::binary
#else
# define vtksys_ios_binary 0
#endif

/* Whether the cstddef header is available.  */
#define vtksys_CXX_HAS_CSTDDEF 1

/* Whether the compiler supports null template arguments.  */
#define vtksys_CXX_HAS_NULL_TEMPLATE_ARGS 1

/* Define the null template arguments macro.  */
#if vtksys_CXX_HAS_NULL_TEMPLATE_ARGS
# define vtksys_CXX_NULL_TEMPLATE_ARGS <>
#else
# define vtksys_CXX_NULL_TEMPLATE_ARGS
#endif

/* Whether the compiler supports member templates.  */
#define vtksys_CXX_HAS_MEMBER_TEMPLATES 1

/* Whether the compiler supports argument dependent lookup.  */
#define vtksys_CXX_HAS_ARGUMENT_DEPENDENT_LOOKUP 1

/* Whether the compiler supports standard full specialization syntax.  */
#define vtksys_CXX_HAS_FULL_SPECIALIZATION 1

/* Define the specialization definition macro.  */
#if vtksys_CXX_HAS_FULL_SPECIALIZATION
# define vtksys_CXX_DEFINE_SPECIALIZATION template <>
#else
# define vtksys_CXX_DEFINE_SPECIALIZATION
#endif

/* Define typename keyword macro for use in declarations.  */
#if defined(_MSC_VER) && _MSC_VER < 1300
# define vtksys_CXX_DECL_TYPENAME
#else
# define vtksys_CXX_DECL_TYPENAME typename
#endif

/* Whether the stl has iterator_traits.  */
#define vtksys_STL_HAS_ITERATOR_TRAITS 1

/* Whether the stl has iterator_category.  */
#define vtksys_STL_HAS_ITERATOR_CATEGORY 0

/* Whether the stl has __iterator_category.  */
#define vtksys_STL_HAS___ITERATOR_CATEGORY 0

/* Whether the stl allocator is the standard template.  */
#define vtksys_STL_HAS_ALLOCATOR_TEMPLATE 1

/* Whether the stl allocator is not a template.  */
#define vtksys_STL_HAS_ALLOCATOR_NONTEMPLATE 0

/* Whether the stl allocator has rebind.  */
#define vtksys_STL_HAS_ALLOCATOR_REBIND 1

/* Whether the stl allocator has a size argument for max_size.  */
#define vtksys_STL_HAS_ALLOCATOR_MAX_SIZE_ARGUMENT 0

/* Whether the stl containers support allocator objects.  */
#define vtksys_STL_HAS_ALLOCATOR_OBJECTS 1

/* Whether struct stat has the st_mtim member for high resolution times.  */
#define vtksys_STAT_HAS_ST_MTIM 

/* If building a C++ file in kwsys itself, give the source file
   access to the macros without a configured namespace.  */
#if defined(KWSYS_NAMESPACE)
# if !vtksys_NAME_IS_KWSYS
#  define kwsys_stl vtksys_stl
#  define kwsys_ios vtksys_ios
#  define kwsys     vtksys
#  define kwsys_ios_binary vtksys_ios_binary
# endif
# define KWSYS_NAME_IS_KWSYS            vtksys_NAME_IS_KWSYS
# define KWSYS_STL_HAVE_STD             vtksys_STL_HAVE_STD
# define KWSYS_IOS_HAVE_STD             vtksys_IOS_HAVE_STD
# define KWSYS_IOS_USE_ANSI             vtksys_IOS_USE_ANSI
# define KWSYS_IOS_USE_SSTREAM          vtksys_IOS_USE_SSTREAM
# define KWSYS_IOS_USE_STRSTREAM_H      vtksys_IOS_USE_STRSTREAM_H
# define KWSYS_IOS_USE_STRSTREA_H       vtksys_IOS_USE_STRSTREA_H
# define KWSYS_IOS_HAVE_BINARY          vtksys_IOS_HAVE_BINARY
# define KWSYS_STAT_HAS_ST_MTIM         vtksys_STAT_HAS_ST_MTIM
# define KWSYS_CXX_HAS_CSTDDEF          vtksys_CXX_HAS_CSTDDEF
# define KWSYS_STL_STRING_HAVE_OSTREAM  vtksys_STL_STRING_HAVE_OSTREAM
# define KWSYS_STL_STRING_HAVE_ISTREAM  vtksys_STL_STRING_HAVE_ISTREAM
# define KWSYS_STL_STRING_HAVE_NEQ_CHAR vtksys_STL_STRING_HAVE_NEQ_CHAR
# define KWSYS_CXX_NULL_TEMPLATE_ARGS   vtksys_CXX_NULL_TEMPLATE_ARGS
# define KWSYS_CXX_HAS_MEMBER_TEMPLATES vtksys_CXX_HAS_MEMBER_TEMPLATES
# define KWSYS_CXX_HAS_FULL_SPECIALIZATION vtksys_CXX_HAS_FULL_SPECIALIZATION
# define KWSYS_CXX_DEFINE_SPECIALIZATION vtksys_CXX_DEFINE_SPECIALIZATION
# define KWSYS_CXX_DECL_TYPENAME        vtksys_CXX_DECL_TYPENAME
# define KWSYS_STL_HAS_ALLOCATOR_REBIND vtksys_STL_HAS_ALLOCATOR_REBIND
# define KWSYS_STL_HAS_ALLOCATOR_MAX_SIZE_ARGUMENT vtksys_STL_HAS_ALLOCATOR_MAX_SIZE_ARGUMENT
# define KWSYS_CXX_HAS_ARGUMENT_DEPENDENT_LOOKUP vtksys_CXX_HAS_ARGUMENT_DEPENDENT_LOOKUP
# define KWSYS_STL_HAS_ITERATOR_TRAITS vtksys_STL_HAS_ITERATOR_TRAITS
# define KWSYS_STL_HAS_ITERATOR_CATEGORY vtksys_STL_HAS_ITERATOR_CATEGORY
# define KWSYS_STL_HAS___ITERATOR_CATEGORY vtksys_STL_HAS___ITERATOR_CATEGORY
# define KWSYS_STL_HAS_ALLOCATOR_TEMPLATE vtksys_STL_HAS_ALLOCATOR_TEMPLATE
# define KWSYS_STL_HAS_ALLOCATOR_NONTEMPLATE vtksys_STL_HAS_ALLOCATOR_NONTEMPLATE
# define KWSYS_STL_HAS_ALLOCATOR_OBJECTS vtksys_STL_HAS_ALLOCATOR_OBJECTS
#endif

#endif
