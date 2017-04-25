#ifndef LowpassHighpass_h
#define LowpassHighpass_h

#include "AudioModule.h"

/**

This class combines a first order lowpass-, and a first order highpass-filter
into a single object. The 2 filters, which are originally in a series 
connection, are implemented as a single biquad-stage.

*/

class LowpassHighpass
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 LowpassHighpass();   ///< Constructor.
 ~LowpassHighpass();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 void setLpfCutoff (double newLpfCutoff);
 ///< Sets the cutoff frequency of the lowpass-filter.

 void setHpfCutoff (double newHpfCutoff);
 ///< Sets the cutoff frequency of the highpass-filter.

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE double getSample(double in);
 ///< Calculates one output sample at a time.

 //---------------------------------------------------------------------------
 // others:

 void resetBuffers();
 ///< Resets the internal buffers of the filter to zero.

 //===========================================================================

protected:

 // buffering:
 doubleA out;
 doubleA x1, x2;           // past input values
 doubleA y1, y2;           // past output values

 // biquad-coefficients:
 doubleA b0, b1, b2;       // feedforward coeffs
 doubleA     a1, a2;       // feedback coeffs

 // filter parameters:
 doubleA lpfCutoff;
 doubleA hpfCutoff;

 doubleA sampleRate;
 doubleA sampleRateRec;    // reciprocal of the sampleRate
 
 // private member functions:
 void calcCoeffs();     
  // calculates the filter coefficients from filter parameters. The 
  // design-equations for the first order filters come from "The Scientist and
  // Engineers Guide to Digital Signal Processing" (www.dspguide.com)
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE double LowpassHighpass::getSample(double in)
{
 // calculate output-sample:
 out = b0*in + b1*x1 + b2*x2
             + a1*y1 + a2*y2 
             + TINY;          // to avoid denorm problems

 // update buffer-variables:
 x2 = x1;
 x1 = in;
 y2 = y1;
 y1 = out;  

 return out;
}

#endif // LowpassHighpass_h
