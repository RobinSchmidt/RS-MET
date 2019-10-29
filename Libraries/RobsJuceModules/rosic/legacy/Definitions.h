#ifndef Definitions_h
#define Definitions_h

#include <intrin.h> 
#include <xmmintrin.h>
#include <float.h>

enum sampleFileStates
{
 EXISTENT_AND_VALID = 1,
 EXISTENT_BUT_INVALID,
 NONEXISTENT,
 NOT_ENOUGH_MEMORY
};

#define INLINE __forceinline

namespace MoreMath
{


//#define INLINE inline

//#define TINY 0.0000000001
#define TINY FLT_MIN




/*
#define _CRT_SECURE_NO_DEPRECATE // to prevent the compiler from complaining 
                                 // about using non-secure functions like 
                                 // strcpy
*/

//type definitions:
typedef __declspec( align(8) ) double doubleA; // doubles, aligned at 64-bit
                                               // (8 byte) boundaries

typedef __declspec( align(8) ) int intA; // integers, aligned at 64-bit
                                         // (8 byte) boundaries

/*
//floating point types:
typedef float flt32;
typedef double flt64;
typedef long double flt80;
typedef __declspec( align(8) ) double flt64A; //aligned at 64-bit (8 byte) boundaries
typedef __declspec( align(4) ) double flt32A; //aligned at 32-bit (4 byte) boundaries

typedef double sample;
typedef __declspec( align(8) ) double sampleA; //aligned at 64-bit boundaries

//unsigned integer types:
typedef unsigned __int64 uint64;
//typedef unsigned __int32 uint32; // clashes with the similar definition in JUCE -> use a namespace
typedef unsigned __int16 uint16;
typedef unsigned __int8 uint8;
typedef __declspec( align(8) ) unsigned __int64 uint64A; //for aligned declarations
typedef __declspec( align(4) ) unsigned __int32 uint32A;
typedef __declspec( align(2) ) unsigned __int16 uint16A;
typedef __declspec( align(1) ) unsigned __int8 uint8A;

//signed integer types:
typedef signed __int64 int64;
typedef signed __int32 int32;
typedef signed __int16 int16;
typedef signed __int8 int8;
typedef __declspec( align(8) ) signed __int64 int64A;
typedef __declspec( align(4) ) signed __int32 int32A;
typedef __declspec( align(2) ) signed __int16 int16A;
typedef __declspec( align(1) ) signed __int8 int8A;

//signed integer types:
#define INT64 signed __int64
#define INT32 signed __int32
#define INT16 signed __int16
#define INT8 signed __int8

*/

//some mathematical contants:
#define PI 3.1415926535897932384626433832795
#define SQRT2 1.4142135623730950488016887242097
#define SQRT2_INV 0.70710678118654752440084436210485

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

//extract the exponent from a IEEE 754 floating point number (single and double precision):
#define EXPOFFLT(value) (((*((reinterpret_cast<unsigned __int32 *>(&value)))&0x7FFFFFFF)>>23)-127)
#define EXPOFDBL(value) (((*((reinterpret_cast<unsigned __int64 *>(&value)))&0x7FFFFFFFFFFFFFFF)>>52)-1023)

} // end namespace



#endif // Definitions_h