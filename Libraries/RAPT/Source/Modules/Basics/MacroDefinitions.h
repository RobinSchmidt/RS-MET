#ifndef RAPT_MACRODEFINITIONS_H_INCLUDED
#define RAPT_MACRODEFINITIONS_H_INCLUDED

#ifdef _MSC_VER
  #define RS_INLINE __forceinline
  #define ASM(x) __asm {x}
#else
  #define RS_INLINE inline
  //#define RS_INLINE __attribute__((always_inline))
  #define ASM(x) asm("x");
#endif

#ifdef RS_DEBUG
  #ifdef _MSC_VER
    #pragma intrinsic (__debugbreak)
    #define RS_DEBUG_BREAK __debugbreak();
  #else
    #define RS_DEBUG_BREAK __builtin_trap();  // preliminary - gcc only
  #endif
  //#define rsAssert(expression) { if (! (expression)) RS_DEBUG_BREAK }
#else
  #define RS_DEBUG_BREAK { }
  #define rsAssert(expression) { }
#endif
#define rsAssertFalse rsAssert(false) 

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
