#ifndef RAPT_MACRODEFINITIONS_H_INCLUDED
#define RAPT_MACRODEFINITIONS_H_INCLUDED

// ToDo: 
// -Try to minimize the visibility of these macros to client code by
//  -organzing the #inclusions in such a way, that not every client file that includes rapt.h
//   is forced to see them
//  -using #undef in the main rapt.h include file to undo the definitions (maybe make a special
//   file MacroUndefinitions.h that can be included at the end of rapt.h)
//  -inclusion of the undef-file will make the macros unavailable in rosic. ...unless we include 
//   this file here from rosic.h again
//  -we may also need to re-include the def-file it in rapt.cpp

//#ifdef _MSC_VER
//  #define RS_INLINE __forceinline
//  #define ASM(x) __asm {x}
//#else
//  #define RS_INLINE inline
//  //#define RS_INLINE __attribute__((always_inline))
//  #define ASM(x) asm("x");
//#endif
#define RS_INLINE inline

#if defined(DEBUG) || defined(_DEBUG)
  #define RS_DEBUG
  #define RS_ASSERT(expression) { if (! (expression)) RS_DEBUG_BREAK }
  #ifdef _MSC_VER
    #pragma intrinsic (__debugbreak)
    #define RS_DEBUG_BREAK __debugbreak();
  #else
    #define RS_DEBUG_BREAK __builtin_trap();  // preliminary - gcc only
  #endif
#else
  #define RS_DEBUG_BREAK { }
  #define RS_ASSERT(expression) { }
#endif
#define RS_ASSERT_FALSE RS_ASSERT(false) 


// Compiler identification:
// preliminary - later wrap that into #ifdef conditons
#define RS_COMPILER_MSC 1        // microsoft compiler
//#define RS_COMPILER_GCC 1     // gnu compiler collection
//#define RS_COMPILER_CLANG 1


// Deprecation:
#define RS_WARN_DEPRECATED 1 // comment to disable deprecation warnings
#ifdef DOXYGEN
  #define RS_DEPRECATED(functionDef)
  #define RS_DEPRECATED_WITH_BODY(functionDef, body)
#elif RS_COMPILER_MSC && RS_WARN_DEPRECATED
  #define RS_DEPRECATED_ATTRIBUTE                  __declspec(deprecated)
  #define RS_DEPRECATED(funcDef)                   RS_DEPRECATED_ATTRIBUTE funcDef
  #define RS_DEPRECATED_WITH_BODY(funcDef, body)   RS_DEPRECATED_ATTRIBUTE funcDef body
#elif (RS_COMPILER_GCC || RS_COMPILER_CLANG) && RS_WARN_DEPRECATED
  #define RS_DEPRECATED_ATTRIBUTE                  __attribute__ ((deprecated))
  #define RS_DEPRECATED(funcDef)                   funcDef RS_DEPRECATED_ATTRIBUTE
  #define RS_DEPRECATED_WITH_BODY(funcDef, body)   funcDef RS_DEPRECATED_ATTRIBUTE body
#else
  #define RS_DEPRECATED_ATTRIBUTE
  #define RS_DEPRECATED(funcDef)                   funcDef
  #define RS_DEPRECATED_WITH_BODY(funcDef, body)   funcDef body
#endif
// see juce_PlatformDefs.h, juce has further macros: JUCE_DEPRECATED_STATIC, 
// JUCE_DECLARE_DEPRECATED_STATIC, JUCE_CATCH_DEPRECATED_CODE_MISUSE
// -> figure out what these do and maybe mimick them, too

// Compiler hinting:
#if defined(__GNUC__) && __GNUC__ >= 4
#define RS_LIKELY(x)   (__builtin_expect((x), 1))
#define RS_UNLIKELY(x) (__builtin_expect((x), 0))
#else
#define RS_LIKELY(x)   (x)
#define RS_UNLIKELY(x) (x)
#endif
// this can be used to hint the compiler that some condition is likely to evaluate to true (or 
// false, for example:
//
// if( sampleCount >= delayLineLength )
//   wrapAround(sampleCount, delayLineLength)
//
// becomes:
//
// if( RS_UNLIKELY(sampleCount >= delayLineLength) )
//   wrapAround(sampleCount, delayLineLength)
//
// ...not really a practical example because a bitmask for wraparound is even better but anyway

// see https://www.youtube.com/watch?v=vrfYLlR8X8k at 1:11:00


// todo: figure out, what to do for the microsoft compiler


// bit twiddling:
/*
// bit-rotataion left and right - todo: implement this for other compilers, for the moment, we have a preliminary mapping to 0
#ifdef _MSVC_VER
#define ROTL(x, k) ( _rotl( (x), (k) ) )
#define ROTR(x, k) ( _rotr( (x), (k) ) )
#else
#define ROTL(x, k) (0)
#define ROTR(x, k) (0)
#endif
*/
//_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON)

#endif
