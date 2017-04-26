/*

*/

#include "string.h"
#include "math.h"
//#include "VSTools.h" //leads to error when functions are declared as extern "C" in VSTools.h

typedef double sample;

//-----------------------------------------------------------------------------------------------------------
//calculates the logarithm to an arbitrary basis B
// Parameters: 1: basis B    (sample-format)
//	            2: argument x (sample-format)
// Return:     logarithm of x to basis B (sample-format)
sample logB(sample B, sample x)
{
 return log(x)/log(B);
}

//-----------------------------------------------------------------------------------------------------------
//calculates the logarithm to basis 2
// Parameters: argument x (sample-format)
// Return:     logarithm of x to basis 2 (sample-format)
sample log2(sample x)
{
 //very ineficcient: return logB(2, x);
 return 1.4426950408889634073599246810019 * log(x); //1.44... = 1/log(2) 
}

//-----------------------------------------------------------------------------------------------------------
//calculates the absolute value of its parameter   
//	Paramter:  any number (sample format)
// Return:    absolute value (sample-format)
sample absDbl(sample value)
{
	sample abs;

	if(value>=0)
		abs = value;
	else
		abs = -value;

	return abs;	
}

//-----------------------------------------------------------------------------------------------------------
//returns the sign of its parameter
//	Paramter:  any number (sample format)
// Return:    sign (+1 or -1) (sample format)
sample sign(sample x)
{
	if(x>=0)
		return 1;
	else
		return -1;
}

//-----------------------------------------------------------------------------------------------------------
//Conversion from an amplitude value to decibels (dB) where an amplitude of 1 corresponds to 0 dB    
//	Paramter:  amplitude value (sample format)
// Return:    level in dB (sample-format)
sample amp2dB(sample amp)
{
	return 20*log10(amp);
}

//-----------------------------------------------------------------------------------------------------------
//Conversion from a dB-value ton an amplitude value (where 0 dB corresponds to 1) 
//	Paramter:  level in dB (sample format)
// Return:    amplitude value (sample-format)
sample dB2amp(sample dB)
{
	return pow(10, (dB/20));
}

//-----------------------------------------------------------------------------------------------------------
//Conversion from a frequency scale factor to a pitch offset in semitones 
//	Paramter: frequency scale factor (sample format)
// Return:   corresponding pitch offset in semitones (sample format)
sample freqFactor2pitchOffset(sample freqFactor)
{
	return 12*log2(freqFactor);
}

//-----------------------------------------------------------------------------------------------------------
//Conversion from a frequency in Hz to a MIDI-Note number   
//	Paramter:  frequency in Hz (sample format)
// Return:    MIDI-Note number (sample-format)
sample frequency2pitch(sample frequency)
{
 return 12*(log2(frequency/440))+69;
}

//-----------------------------------------------------------------------------------------------------------
//Conversion from a MIDI-Note Number to a frequency in Hz
//	Paramter: MIDI-Note number (sample-format) 
// Return:   corresponding frequency in Hz (sample format) 
sample pitch2frequency(sample pitch)
{
 return 440*( pow(2, ((pitch-69)/12)));
}

//-----------------------------------------------------------------------------------------------------------
//Conversion from a pitch offset in semitones to a frequency scale factor
//	Paramter: pitch offset in semitones (sample format) 
// Return:   corresponding frequency scale factor (sample format)
sample pitchOffset2freqFactor(sample pitchOffset)
{
	return pow(2, (pitchOffset/12));
}

//-----------------------------------------------------------------------------------------------------------
//clip the signal at the specified threshold (seems not to work yet)
// parameter: signal to be clipped and clip-treshold
// return:    clipped signal
sample clip(sample signal, sample threshold)
{
	if(signal > threshold)
		return threshold;
	else if(signal < (-threshold))
		return -threshold;
	else
		return signal;
}


short getExponentFromFloat(float ValueFlt)
{
 const unsigned long bitMask = 0x7FFFFFFF;  //unsigned 32 bit integer 
                                             //leftmost bit is 0, all others are 1

 static unsigned long exponent = 0;  //this is the variable which should finally hold
                                     //the exponent

 static unsigned long *ptrInt; //pointer to a 32-bit integer
 static void          *ptr;    //pointer to to anything

 ptr      = &ValueFlt;               //determine the adress of the float-variable
 ptrInt   = (unsigned long*) ptr;    //store this adress in the int-pointer

 exponent = *ptrInt;                 //get the value and interpret is as int

 //mask out the sign-bit via bitwise and:
 exponent = exponent & bitMask;

 //rightshift by 23 to shift the mantissa away:
 exponent =  exponent >> 23;

 //subtract the bias:
 exponent -= 127;

 return (short) exponent;
}

short getExponentFromDouble(double ValueDbl)
{
 const unsigned __int64 bitMask  = 0x7FFFFFFFFFFFFFFF;  //unsigned 64 bit integer 
                                                         //leftmost bit is 0, all others are 1

 static unsigned __int64 exponent = 0;  //this is the variable which should finally hold
                                        //the exponent

 static unsigned __int64 *ptrInt; //pointer to a 64-bit integer
 static void             *ptr;    //pointer to to anything


 //reinterpret the bit pattern of the double as an unsigned 64 bit integer by
 //meas of pointer-casts:
 ptr      = &ValueDbl;               //determine the adress of the double-variable
 ptrInt   = (unsigned __int64*) ptr; //store this adress in the int-pointer
 exponent = *ptrInt;                 //get the value and interpret is as int

 //mask out the sign-bit via bitwise and:
 exponent = exponent & bitMask;

 //rightshift by 53 to shift the mantissa away:
 exponent =  exponent >> 52;

 //subtract the bias:
 exponent -= 1023;

 return (short) exponent;
}


//-----------------------------------------------------------------------------------------------------------
//converts an integer midi-note-value into a corresponding note-string
// Parameters: 1: MIDI-Note number (long-format)
//	            2: Pointer to char for writing the string (writes always 4 characters plus ASCII-NULL)
// Return:     no return value, the return-string is stored in the second parameter 
void pitch2string(long pitch, char *notestring)
{
	switch(pitch)
	{
 	case   0: strcpy(notestring, "C-1 ");break;
		case   1: strcpy(notestring, "C#-1");break;
		case   2: strcpy(notestring, "D-1 ");break;
		case   3: strcpy(notestring, "D#-1");break;
	 case   4: strcpy(notestring, "E-1 ");break;
		case   5: strcpy(notestring, "F-1 ");break;
		case   6: strcpy(notestring, "F#-1");break;
		case   7: strcpy(notestring, "G-1 ");break;
		case   8: strcpy(notestring, "G#-1");break;
		case   9: strcpy(notestring, "A-1 ");break;
		case  10: strcpy(notestring, "A#-1");break;
	 case  11: strcpy(notestring, "B-1 ");break;

		case  12: strcpy(notestring, "C0  ");break;
		case  13: strcpy(notestring, "C#0 ");break;
		case  14: strcpy(notestring, "D0  ");break;
		case  15: strcpy(notestring, "D#0 ");break;
	 case  16: strcpy(notestring, "E0  ");break;
		case  17: strcpy(notestring, "F0  ");break;
		case  18: strcpy(notestring, "F#0 ");break;
		case  19: strcpy(notestring, "G0  ");break;
		case  20: strcpy(notestring, "G#0 ");break;
		case  21: strcpy(notestring, "A0  ");break;
		case  22: strcpy(notestring, "A#0 ");break;
	 case  23: strcpy(notestring, "B0  ");break;

		case  24: strcpy(notestring, "C1  ");break;
		case  25: strcpy(notestring, "C#1 ");break;
		case  26: strcpy(notestring, "D1  ");break;
		case  27: strcpy(notestring, "D#1 ");break;
	 case  28: strcpy(notestring, "E1  ");break;
		case  29: strcpy(notestring, "F1  ");break;
		case  30: strcpy(notestring, "F#1 ");break;
		case  31: strcpy(notestring, "G1  ");break;
		case  32: strcpy(notestring, "G#1 ");break;
		case  33: strcpy(notestring, "A1  ");break;
		case  34: strcpy(notestring, "A#1 ");break;
	 case  35: strcpy(notestring, "B1  ");break;

		case  36: strcpy(notestring, "C2  ");break;
		case  37: strcpy(notestring, "C#2 ");break;
		case  38: strcpy(notestring, "D2  ");break;
		case  39: strcpy(notestring, "D#2 ");break;
	 case  40: strcpy(notestring, "E2  ");break;
		case  41: strcpy(notestring, "F2  ");break;
		case  42: strcpy(notestring, "F#2 ");break;
		case  43: strcpy(notestring, "G2  ");break;
		case  44: strcpy(notestring, "G#2 ");break;
		case  45: strcpy(notestring, "A2  ");break;
		case  46: strcpy(notestring, "A#2 ");break;
	 case  47: strcpy(notestring, "B2  ");break;

		case  48: strcpy(notestring, "C3  ");break;
		case  49: strcpy(notestring, "C#3 ");break;
		case  50: strcpy(notestring, "D3  ");break;
		case  51: strcpy(notestring, "D#3 ");break;
	 case  52: strcpy(notestring, "E3  ");break;
		case  53: strcpy(notestring, "F3  ");break;
		case  54: strcpy(notestring, "F#3 ");break;
		case  55: strcpy(notestring, "G3  ");break;
		case  56: strcpy(notestring, "G#3 ");break;
		case  57: strcpy(notestring, "A3  ");break;
		case  58: strcpy(notestring, "A#3 ");break;
	 case  59: strcpy(notestring, "B3  ");break;

		case  60: strcpy(notestring, "C4  ");break;
		case  61: strcpy(notestring, "C#4 ");break;
		case  62: strcpy(notestring, "D4  ");break;
		case  63: strcpy(notestring, "D#4 ");break;
	 case  64: strcpy(notestring, "E4  ");break;
		case  65: strcpy(notestring, "F4  ");break;
		case  66: strcpy(notestring, "F#4 ");break;
		case  67: strcpy(notestring, "G4  ");break;
		case  68: strcpy(notestring, "G#4 ");break;
		case  69: strcpy(notestring, "A4  ");break;
		case  70: strcpy(notestring, "A#4 ");break;
	 case  71: strcpy(notestring, "B4  ");break;

		case  72: strcpy(notestring, "C5  ");break;
		case  73: strcpy(notestring, "C#5 ");break;
		case  74: strcpy(notestring, "D5  ");break;
		case  75: strcpy(notestring, "D#5 ");break;
	 case  76: strcpy(notestring, "E5  ");break;
		case  77: strcpy(notestring, "F5  ");break;
		case  78: strcpy(notestring, "F#5 ");break;
		case  79: strcpy(notestring, "G5  ");break;
		case  80: strcpy(notestring, "G#5 ");break;
		case  81: strcpy(notestring, "A5  ");break;
		case  82: strcpy(notestring, "A#5 ");break;
	 case  83: strcpy(notestring, "B5  ");break;

		case  84: strcpy(notestring, "C6  ");break;
		case  85: strcpy(notestring, "C#6 ");break;
		case  86: strcpy(notestring, "D6  ");break;
		case  87: strcpy(notestring, "D#6 ");break;
	 case  88: strcpy(notestring, "E6  ");break;
		case  89: strcpy(notestring, "F6  ");break;
		case  90: strcpy(notestring, "F#6 ");break;
		case  91: strcpy(notestring, "G6  ");break;
		case  92: strcpy(notestring, "G#6 ");break;
		case  93: strcpy(notestring, "A6  ");break;
		case  94: strcpy(notestring, "A#6 ");break;
	 case  95: strcpy(notestring, "B6  ");break;

		case  96: strcpy(notestring, "C7  ");break;
		case  97: strcpy(notestring, "C#7 ");break;
		case  98: strcpy(notestring, "D7  ");break;
		case  99: strcpy(notestring, "D#7 ");break;
	 case 100: strcpy(notestring, "E7  ");break;
		case 101: strcpy(notestring, "F7  ");break;
		case 102: strcpy(notestring, "F#7 ");break;
		case 103: strcpy(notestring, "G7  ");break;
		case 104: strcpy(notestring, "G#7 ");break;
		case 105: strcpy(notestring, "A7  ");break;
		case 106: strcpy(notestring, "A#7 ");break;
	 case 107: strcpy(notestring, "B7  ");break;

		case 108: strcpy(notestring, "C8  ");break;
		case 109: strcpy(notestring, "C#8 ");break;
		case 110: strcpy(notestring, "D8  ");break;
		case 111: strcpy(notestring, "D#8 ");break;
	 case 112: strcpy(notestring, "E8  ");break;
		case 113: strcpy(notestring, "F8  ");break;
		case 114: strcpy(notestring, "F#8 ");break;
		case 115: strcpy(notestring, "G8  ");break;
		case 116: strcpy(notestring, "G#8 ");break;
		case 117: strcpy(notestring, "A8  ");break;
		case 118: strcpy(notestring, "A#8 ");break;
	 case 119: strcpy(notestring, "B8  ");break;

		case 120: strcpy(notestring, "C9  ");break;
		case 121: strcpy(notestring, "C#9 ");break;
		case 122: strcpy(notestring, "D9  ");break;
		case 123: strcpy(notestring, "D#9 ");break;
	 case 124: strcpy(notestring, "E9  ");break;
		case 125: strcpy(notestring, "F9  ");break;
		case 126: strcpy(notestring, "F#9 ");break;
		case 127: strcpy(notestring, "G9  ");break;

		default: strcpy(notestring, "ERR "); //value represents no valid MIDI-Note
	}
	return;
}
