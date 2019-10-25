#ifndef RAPT_MACRODEFINITIONS_H_INCLUDED
#define RAPT_MACRODEFINITIONS_H_INCLUDED

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

// compiler hinting:

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
