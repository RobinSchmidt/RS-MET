#ifndef Upsampler_h
#define Upsampler_h

#include "AudioModule.h"

/**

This module peforms upsampling of an input signal by means of linear 
interpolation. Simply call setInputSample() for each sample on the old 
sample-rate and for each setInputSample()-call you call n times the 
getSample() function where n is the oversampling factor. This makes n samples 
from 1 sample.

*/

class Upsampler : public AudioModule
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          Upsampler();  ///< Constructor.
 virtual ~Upsampler();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setOversamplingFactor(int newOversamplingFactor);
 ///< Sets the factor of oversampling.

 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE void setInputSample(double in);
 /**< Puts one sample into the upsampler. For each call to this function you 
      should call n times the getSample()-function where n is the oversampling
      factor. */

 virtual INLINE double getSample();
 /**< Calculates an linearly interpolated output sample. */


 //---------------------------------------------------------------------------
 // others:
 virtual void reset();

 //===========================================================================

protected:

 intA oversamplingFactor;
 doubleA oversamplingFactorRec;

 doubleA previousInput;
 //doubleA currentInput;
 doubleA slope;

 doubleA position;

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):


INLINE void Upsampler::setInputSample(double in)
{
 //calculate the interpolation slope:
 slope = (in-previousInput)*oversamplingFactorRec;

 previousInput = in;
 position      = 0.0;
}

INLINE double Upsampler::getSample()
{
 static double out;

 out = previousInput + slope*position;

 position += 1.0;

 return out;
}

#endif // Upsampler_h
