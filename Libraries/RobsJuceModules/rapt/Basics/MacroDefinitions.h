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


// Build configuration macros:
#define RS_USE_SSE
#define RS_USE_SSE2  // if you define this, RS_USE_SSE should also be defined
// These macros are actually supposed to be set by client code to determine the build config, so 
// they should perhaps reside in another file. Maybe a file BuildConfig.h that resides next to
// rapt.h and is included right before this one. Or maybe they should be placed into the .jucer 
// files? But no: this would imply a proliferation of exporters in Projucer, which are a
// pain to maintain anyway. Also, we don't really want to depend on Projucer for all builds. 
// We'll see... For the time being, they are here.

// Identify compiler:
#ifdef _MSC_VER  
  #define RS_COMPILER_MSC 1        // microsoft compiler
#elif __GNUC__   
  #define RS_COMPILER_GCC 1        // gnu compiler collection
#elif __clang__  
  #define RS_COMPILER_CLANG 1      // clang compiler
#endif

// Identify build target architecture:
#ifdef RS_COMPILER_MSC
  #ifdef _M_X64
    #define RS_ARCHITECTURE_X64
  #endif
#elif defined(RS_COMPILER_GCC) || defined(RS_COMPILER_CLANG)
  #ifdef __x86_64__
    #define RS_ARCHITECTURE_X64
  #endif
#endif
// todo: 
// -Verify the gcc/clang branch and add macro definition RS_ARCHITECTURE_ARM64 for ARM processors
// -Maybe we should also check the "-arch" settings of the compiler? Currently, we just assume that
//  the build machine is the same architecture as the target machine, but that doesn't need to be 
//  the case.

// Identify operating system:
// ...


// Identify instruction set architecture (ISA) to build for. This depends on the build target 
// architecture and some configuration macros set by client code.
#if defined(RS_ARCHITECTURE_X64) && defined(RS_USE_SSE)
  #define RS_INSTRUCTION_SET_SSE
#endif
#if defined(RS_ARCHITECTURE_X64) && defined(RS_USE_SSE2)
  #define RS_INSTRUCTION_SET_SSE2
#endif
// maybe we can use SSE also fo X86 in 32bit builds?

// If this macros is defined, it indicates that none of SIMD instruction sets should be used, i.e.
// the scalar fallback implementations should be used:
#if !defined(RS_USE_SSE)  // later, use "or", i.e. || with other simd instruction sets
  #define RS_NO_SIMD_FLOAT32X4
#endif
#if !defined(RS_USE_SSE2)  // later, use "or", i.e. || with other simd instruction sets
  #define RS_NO_SIMD_FLOAT64X2
#endif



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



/** Deprecation macros. To depracte a function declaration (in this case rsIFFT) use:

  template<class T> // replacement: rsLinearTransforms::fourierInvRadix2DIF
  RS_DEPRECATED(void rsIFFT(std::complex<T> *buffer, int N));

In this case, the function (rsIFFT) is defined elsewhere, and just calls the new function that 
acts as the new replacement for the old, deprecated function:

  template<class T>
  void rsIFFT(std::complex<T> *a, int N) { rsLinearTransforms::fourierInvRadix2DIF(a, N); }

i.e. the old, deprecated name will remain available as alias for the new function name, but using
it will issue deprecation warnings during compilation. We can also do the following:

  template<class T> RS_DEPRECATED_WITH_BODY(
    void rsRadix2FFT(std::complex<T>* x, int N),
    { rsLinearTransforms::fourierRadix2DIF(x, N); })  // use that directly instead!

to deprecate a function and call the new function directly inside the old in just a single macro. 
When using these macros, take care of the placement of the commas and semicolons. It's easy to 
make mistakes and then it will not compile. If the deprecated functions are not templatized, then
the template<class T> thing can just be removed and the T will be replaced by an actual type. It's
good practice to always hint in a comment, what the replacement is, i.e. how client code should be
updated to get rid of the deprecation warnings. */

//#define RS_WARN_DEPRECATED 1 // comment to disable deprecation warnings
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
// Will these macros work also for classes, such that we can declare a class-name deprecated but at
// the same time define the old name as alias for the new one? Maybe the old should inherit from 
// the new, also inheriting constructors? maybe define a macor RS_DEPRECTE_CLASS
// maybe try it with rsPolynomial -> rename it to rsPolynom ...but maybe not - that's not widely 
// used as an english word



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

// See:
// https://blog.kowalczyk.info/article/j/guide-to-predefined-macros-in-c-compilers-gcc-clang-msvc-etc..html
// https://sourceforge.net/p/predef/wiki/Home/
// https://sourceforge.net/p/predef/wiki/Compilers/
// https://sourceforge.net/p/predef/wiki/Architectures/

// https://abseil.io/docs/cpp/platforms/macros


// https://stackoverflow.com/questions/23934862/what-predefined-macro-can-i-use-to-detect-the-target-architecture-in-clang/41666292

#endif
