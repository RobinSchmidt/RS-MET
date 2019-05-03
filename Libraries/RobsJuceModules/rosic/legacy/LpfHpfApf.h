#ifndef LpfHpfApf_h
#define LpfHpfApf_h

#include "AudioModule.h"

/**

This class combines a first order lowpass-, first order highpass- and first
order allpass-filter into a single object. The 3 filters are in series.
This class is makes it more convenient to use such a combination of filters,
when desired - for example in feedback loops for delay-lines to make their
feedback frequency dependent and/or dispersive.

*/

class LpfHpfApf
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 LpfHpfApf();   ///< Constructor.
 ~LpfHpfApf();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 void setLpfCutoff (double newLpfCutoff);
 ///< Sets the cutoff frequency of the lowpass-filter.

 void setHpfCutoff (double newHpfCutoff);
 ///< Sets the cutoff frequency of the highpass-filter.

 void setApfCutoff (double newApfCutoff);
 ///< Sets the characteristic frequency of the allpass-filter.

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE double getSample(double in);
 ///< Calculates one output sample at a time.

 //---------------------------------------------------------------------------
 //others:

 void resetBuffers();
 ///< Resets the internal buffers of the filter to zero.

 //===========================================================================

protected:

 // buffering:
 doubleA out;
 doubleA x1, x2, x3;  // past input values
 doubleA y1, y2, y3;  // past output values

 // overall 3-pole filter-coefficients:
 doubleA b0, b1, b2, b3; // feedforward coeffs
 doubleA     a1, a2, a3; // feedback coeffs

 // filter parameters:
 doubleA lpfCutoff;
 doubleA hpfCutoff;
 doubleA apfCutoff;

 doubleA sampleRate;
 doubleA sampleRateRec;  // reciprocal of the sampleRate
 
 // private member functions:
 void calcCoeffs(); 
  // calculates the filter coefficients from filter parameters. The 
  // design-equations for the first order lowpass- and highpass-filters come
  // from "The Scientist and Engineers Guide to Digital Signal Processing" 
  // (www.dspguide.com) and the allpass-design is from the DAFX-book
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE double LpfHpfApf::getSample(double in)
{
 // calculate output-sample:
 out = b0*in + b1*x1 + b2*x2 + b3*x3
             + a1*y1 + a2*y2 + a3*y3;

 // update buffer-variables:
 x3 = x2;
 x2 = x1;
 x1 = in;
 y3 = y2;
 y2 = y1;
 y1 = out;

 return out;
}

#endif // LpfHpfApf_h
