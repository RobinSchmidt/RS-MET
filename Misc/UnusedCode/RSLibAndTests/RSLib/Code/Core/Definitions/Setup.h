#ifndef RS_SETUP_H
#define RS_SETUP_H

// identify operating system:
#if defined(_WIN32) || defined(__WIN32__)
  #define RS_SYSTEM_WIN32
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN // affects the declarations when including <windows.h>
  #endif
  #ifndef NOMINMAX
    #define NOMINMAX            // prevents MSVC system headers from defining MIN and MAX
  #endif

#elif defined(linux) || defined(__linux)
  #define RS_SYSTEM_LINUX

#elif defined(__APPLE__) || defined(MACOSX) || defined(macintosh) || defined(Macintosh)
  #define RS_SYSTEM_OSX

#else
  #error operating system not supported by RSLib

#endif

// identify compiler:
#if defined(_MSC_VER)
 #define RS_COMPILER_MICROSOFT  // rename to RS_COMPILER_MS
 #if defined(_M_X64)
  #define RS_COMPILER_MS_X64
 #else
  #define RS_COMPILER_MS_X86
 #endif
#elif defined(__GNUC__)
 #define RS_COMPILER_GCC        // maybe, we need to distinguish X86 and X64 versions here too?
#else
 #warning(Unknown Compiler)
#endif



// define portable debug macro:
#if !defined(NDEBUG)  // NDEBUG is added to the list of defined macros for Release
                      // builds when the MSVC project is created.
  #define RS_DEBUG
#endif

// macros for declaring that classes and functions should be exported to the dll when the library
// is built and imported from a dll when the client is built:
#if defined(_MSC_VER)

  #pragma warning(disable : 4241) // appears on dll-exported template instances
  #pragma warning(disable : 4251) // appears on (template) classes without RSLib_API prefix
  #ifndef _CRT_SECURE_NO_WARNINGS 
    #define _CRT_SECURE_NO_WARNINGS // MS compiler annoys with recommending secure functions
  #endif

  #if RS_BUILD_DLL
    // compiling the RSLib project as dll
    #define RSLib_API __declspec(dllexport)
    #define RSLib_TEMPLATE_INSTANCE

  #elif RS_USE_DLL
    // compiling a client project, that uses the RSLib as dll
    #define RSLib_API __declspec(dllimport)
    #define RSLib_TEMPLATE_INSTANCE extern

  #else
    #define RSLib_API
    // neither building as nor linking to a dll - nothing to do - empty macro
  #endif
  
#else
  #define RSLib_API // empty definition - nothing to do on compilers other than the MS

#endif

#endif
