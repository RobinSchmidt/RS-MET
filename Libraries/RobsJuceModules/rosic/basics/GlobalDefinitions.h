#ifndef GlobalDefinitions_h
#define GlobalDefinitions_h

namespace rosic
{

/** This file contains a bunch of useful macros which are not wrapped into the
rosic namespace to facilitate their global use. 

many macros should be either transfered to RAPT or maybe removed altogether - only those that 
clearly have a datatype (are non-generic) should remain here
*/

#ifdef _MSC_VER
#define INLINE __forceinline
#else
#define INLINE inline
//#define INLINE __attribute__((always_inline))
#endif


#define MSC_X86_ASM 0  // define this for certain asm number manipulations on X86 with MS compiler

// something better to do here...

//_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON)

/*
#ifdef LINUX
  #define BEGIN_ASM asm(
  #define END_ASM );
  #define ASM(x) "x"
#else
  #define BEGIN_ASM __asm{
  #define END_ASM }
  #define ASM(x) (x)
#endif
*/

#ifdef LINUX
#define ASM(x) asm("(x)");
#else
#define ASM(x) __asm x
#endif

//-------------------------------------------------------------------------------------------------
// mathematical constants:

// get rid of them...

//#define PI 3.1415926535897932384626433832795
//#define TWO_PI 6.283185307179586476925286766559
//#define EULER 2.7182818284590452353602874713527
//#define SQRT2 1.4142135623730950488016887242097
//#define SQRT2_INV 0.70710678118654752440084436210485
//#define LN10 2.3025850929940456840179914546844
//#define LN10_INV 0.43429448190325182765112891891661
//#define LN2 0.69314718055994530941723212145818
//#define LN2_INV 1.4426950408889634073599246810019
//#define SEMITONE_FACTOR 1.0594630943592952645618252949463

//#define MILLI 0.001
//#define MICRO 0.000001
//#define NANO 0.000000001  // used by romos - get rid of usage there - define a constexpr there


//-------------------------------------------------------------------------------------------------
// type definitions:

// UnaryFunctionPointer is a pointer to a function that takes a double and returns a double:
typedef double(*UnaryFunctionPointer) (double);

/** Pointer to a function that takes a double and returns nothing (void). */
typedef void(*FunctionPointerDoubleToVoid) (double);

/** Pointer to a function that takes an int and a double and returns nothing (void). */
typedef void(*FunctionPointerIntDoubleToVoid) (int, double);

/** Pointer to a function that takes two ints and a double and returns nothing (void). */
typedef void(*FunctionPointerIntIntDoubleToVoid) (int, int, double);


// doubles, aligned at 64-bit (8 byte) boundaries:
#ifdef _MSC_VER
typedef __declspec(align(8)) double doubleA;
#else
typedef double doubleA; // something to do here...
#endif

// doubles, aligned at 128-bit (16 byte) boundaries:
#ifdef _MSC_VER
typedef __declspec(align(16)) double doubleA16;
#else
typedef double doubleA16;
#endif

// integers, aligned at 64-bit (8 byte) boundaries:
#ifdef _MSC_VER
typedef __declspec(align(8)) int intA;
#else
typedef int intA;
#endif

// unsigned 64 bit integers:
#ifdef _MSC_VER
typedef unsigned __int64 UINT64;
#else
typedef unsigned long long UINT64;
#endif

// signed 64 bit integers:
#ifdef _MSC_VER
typedef signed __int64 INT64;
#else
typedef signed long long INT64;
#endif

// unsigned 32 bit integers:
#ifdef _MSC_VER
typedef unsigned __int32 UINT32;
#else
//typedef unsigned long UINT32;
#endif

// 32 bit integers:
//typedef signed int INT32;

// ...constants for numerical precision issues, denorm, etc.:
#define TINY FLT_MIN
#define EPS DBL_EPSILON

// define infinity values:
//#define INF (1.0 / sqrt(0.0))       // this cheats the compiler to generate infinity
//#define NEG_INF (-1.0 / sqrt(0.0))

//UINT64 INF_ULL (0x7FF0000000000000ULL);
//UINT64 NEG_INF_ULL (0xFFF0000000000000ULL);
//#define INF (*(reinterpret_cast<double *>(&INF_ULL)))
//#define NEG_INF (*(reinterpret_cast<double *>(&NEG_INF_ULL)))

//// old:
//#ifdef _MSC_VER
//  INLINE double dummyFunction(double x) { return x; }
//  #define INF (1.0/dummyFunction(0.0))
//  #define NEG_INF (-1.0/dummyFunction(0.0))  // get rid of that - use -INF in the code
//  //#ifndef NAN
//  //  #define NAN (dummyFunction(0.0)/dummyFunction(0.0))
//  //#endif
//  #define INDEF (dummyFunction(0.0) / dummyFunction(0.0))
//#else
//  #define INF (1.0/0.0)
//  #define NEG_INF (-1.0/0.0)
//  #define INDEF (0.0/0.0)
//#endif

// new - move to rapt ...or get rid entirely:
static const double INF     =  std::numeric_limits<double>::infinity();
static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double INDEF   =  std::numeric_limits<double>::quiet_NaN();

//-------------------------------------------------------------------------------------------------
// debug stuff:

// this will try to break the debugger if one is currently hosting this app:
#ifdef _DEBUG

#ifdef _MSC_VER
#pragma intrinsic (__debugbreak)
#define DEBUG_BREAK __debugbreak();
#else
#define DEBUG_BREAK {}  // preliminary
#endif

#define DEBUG_HOOK { int debugDummy = 0; } // hook that can be placed for setting breakpoints

#else

// evaluate macros to no-op in release builds:
#define DEBUG_HOOK  {}
#define DEBUG_BREAK {}

#endif


// an replacement of the ASSERT macro
#define rassert(expression)         { if (! (expression)) DEBUG_BREAK }

//-------------------------------------------------------------------------------------------------
// bit twiddling:

//extract the exponent from a IEEE 754 floating point number (single and double precision):
#define EXPOFFLT(value) (((*((reinterpret_cast<UINT32 *>(&value)))&0x7FFFFFFF)>>23)-127)
#define EXPOFDBL(value) (((*((reinterpret_cast<UINT64 *>(&value)))&0x7FFFFFFFFFFFFFFFULL)>>52)-1023)
  // ULL indicates an unsigned long long literal constant

// bit-rotataion left and right - todo: implement this for other compilers to make key-genereation
// and -validation work, for the moment, we have a preliminary mapping to 0
#ifdef _MSVC_VER
#define ROTL(x, k) ( _rotl( (x), (k) ) )
#define ROTR(x, k) ( _rotr( (x), (k) ) )
#else
#define ROTL(x, k) (0)
#define ROTR(x, k) (0)
#endif

//-------------------------------------------------------------------------------------------------
// functions macros (some are still to be optimized algebraically):

//some useful macros for VST-Plug-In-Development:
//calculates the logarithm to an arbitrary base:
#define LOGB(x, B) (log((double)x)/log((double)B))

//calculates the logarithm to base 2:
#define LOG2(x) (1.4426950408889634073599246810019 * log((double)x))

//conversion from a MIDI-Note Number to a frequency in Hz:
#define PITCH2FREQ(pitch) (440*( pow(2, ((((double)pitch)-69)/12))))

//conversion from a frequency in Hz to a MIDI-Note Number:
//#define FREQ2PITCH(freq) (12*(LOG2(((double)freq)/440))+69) //macro is not tested yet

//conversion from a pitch offset value in semitones to a frequency scale factor:
#define PITCHOFFSET2FREQFACTOR(pitchOffset) (pow(2, (((double)pitchOffset)/12)))

//conversion from a frequency scale factor to a pitch offset in semitones:
#define FREQFACTOR2PITCHOFFSET(freqFactor) (12*LOG2((double)freqFactor))

//conversion from an amplitude-value to a level value in dB (amplitude 1 corresponds to 0 dB):
#define AMP2DB(amp) (20*log10((double)amp))

//conversion from a level value in dB to an amplitude-value:
#define DB2AMP(dB) (pow(10, (((double)dB)/20)))

//clip the signal at a specified threshold:
#define CLIP(signal, thresh) {if(signal > thresh) signal=thresh; else if(signal < (-thresh)) signal=-thresh;}

}

#endif
