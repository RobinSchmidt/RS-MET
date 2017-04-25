#include "math.h"
#include "string.h"

typedef unsigned __int64 uint64;
typedef unsigned __int32 uint32;
typedef unsigned __int16 uint16;
typedef unsigned __int8 uint8;

typedef signed __int64 int64;
typedef signed __int32 int32;
typedef signed __int16 int16;
typedef signed __int8 int8;

typedef float flt32;
typedef double flt64;
typedef long double flt80;

typedef double sample;


#pragma inline_depth(10)

//some mathematical contants:
#define PI 3.1415926535897932384626433832795

//some useful macros for VST-Plug-In-Development:

//calculates the logarithm to an arbitrary base:
#define LOGB(x, B) (log((sample)x)/log((sample)B))

//calculates the logarithm to base 2:
#define LOG2(x) (1.4426950408889634073599246810019 * log((sample)x))

//conversion from a MIDI-Note Number to a frequency in Hz:
#define PITCH2FREQ(pitch) (440*( pow(2, ((((sample)pitch)-69)/12))))

//conversion from a frequency in Hz to a MIDI-Note Number:
//#define FREQ2PITCH(freq) (12*(LOG2(((sample)freq)/440))+69) //macro is not tested yet

//conversion from a pitch offset value in semitones to a frequency scale factor:
#define PITCHOFFSET2FREQFACTOR(pitchOffset) (pow(2, (((sample)pitchOffset)/12)))

//conversion from a frequency scale factor to a pitch offset in semitones:
#define FREQFACTOR2PITCHOFFSET(freqFactor) (12*LOG2((sample)freqFactor))

//conversion from an amplitude-value to a level value in dB (amplitude 1 corresponds to 0 dB):
#define AMP2DB(amp) (20*log10((sample)amp))

//conversion from a level value in dB to an amplitude-value:
#define DB2AMP(dB) (pow(10, (((sample)dB)/20)))

//clip the signal at a specified threshold:
#define CLIP(signal, thresh) {if(signal > thresh) signal=thresh; else if(signal < (-thresh)) signal=-thresh;}

//extract the exponent from a IEEE 754 floating point number (single and double precision):
#define EXPOFFLT(value) (((*((reinterpret_cast<uint32 *>(&value)))&0x7FFFFFFF)>>23)-127)
#define EXPOFDBL(value) (((*((reinterpret_cast<uint64 *>(&value)))&0x7FFFFFFFFFFFFFFF)>>52)-1023)

//some useful functions for VST-Plug-In-Development:
extern "C"
{

sample pitch2frequency(sample pitch);
//conversion from a MIDI-Note Number to a
//frequency in Hz

sample frequency2pitch(sample frequency);
//conversion from a frequency in Hz to a
//MIDI-Note Number

sample pitchOffset2freqFactor(sample pitchOffset);
//conversion from a pitch offset value in semitones
//to a frequency scale factor

sample freqFactor2pitchOffset(sample freqFactor);
//conversion from a frequency scale factor to a
//pitch offset in semitones

sample amp2dB(sample amp);
//conversion from an amplitude-value (0-1)
//to a level value in dB

sample dB2amp(sample db);
//conversion from a level value in dB
//to an amplitude-value (0-1)

sample logB(sample basis, sample x);
//calculates the logarithm to an arbitrary base

sample log2(sample x);
//calculates the logarithm to base 2

sample absDbl(sample x);
//calculates the absolute value

sample sign(sample x);
//returns the sign of the input (+1 or -1)

sample clip(sample signal, sample threshold);
//clips the signal at the specified threshold

inline short getExponentFromFloat(float valueFlt);
//extracts the exponent from a IEEE 754 
//single precision floating point number 

inline short getExponentFromDouble(double valueDbl);
//extracts the exponent from a IEEE 754 
//double precision floating point number 

void pitch2string(long pitch, char *notestring);
//converts a MIDI-Note number to a string (4 charachters + ASCII-NULL)

}


