#ifndef NoiseGenerator_h
#define NoiseGenerator_h

#include <stdlib.h>
#include "AudioModule.h"
#include "IirFilter.h"

/**

This is a noise generator capapble of producing various kinds of noise.

*/

class NoiseGenerator : public AudioModule
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          NoiseGenerator();  ///< Constructor.
 virtual ~NoiseGenerator();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual void setLowpassOrder(int newLowpassOrder);
 /**< sets the order of the lowpass-filter. */

 virtual void setLowpassCutoff(double newLowpassCutoff);
 /**< sets the cutoff-frequency of the lowpass-filter. */

 virtual void setHighpassOrder(int newHighpassOrder);
 /**< sets the order of the highpass-filter. */

 virtual void setHighpassCutoff(double newHighpassCutoff);
 /**< sets the cutoff-frequency of the highpass-filter. */


 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE double getSampleWhite();
 /**< Calculates a sample of white noise in the range -1.0...+1.0. */

 virtual INLINE double getSampleLowpassWhite();
 /**< Calculates a sample of lowpass-filtered white noise. A gain factor is 
      applied to compensate for the loss of energy in the filter. */

 virtual INLINE double getSampleHighpassWhite();
 /**< Calculates a sample of highpass-filtered white noise. A gain factor is 
      applied to compensate for the loss of energy in the filter. */

 virtual INLINE double getSampleBandpassWhite();
 /**< Calculates a sample of bandpass-filtered white noise. A gain factor is 
      applied to compensate for the loss of energy in the filter. */

 virtual INLINE double getSamplePink();
 /**< Calculates a sample of pink noise in the range, A gain factor is 
      applied to compensate for the loss of energy in the pinking-filter. */

 virtual INLINE double getSampleLowpassPink();
 /**< Calculates a sample of lowpass-filtered pink noise. A gain factor is 
      applied to compensate for the loss of energy in the filter. */

 virtual INLINE double getSampleHighpassPink();
 /**< Calculates a sample of highpass-filtered pink noise. A gain factor is 
      applied to compensate for the loss of energy in the filter. */

 virtual INLINE double getSampleBandpassPink();
 /**< Calculates a sample of bandpass-filtered pink noise. A gain factor is 
      applied to compensate for the loss of energy in the filter. */


 //---------------------------------------------------------------------------
 // others:

protected:

 // embedded audio-modules:
 IirFilter lowpassFilter, highpassFilter;

 double scale, offset;

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double NoiseGenerator::getSampleWhite()
{
 return scale * ((double) rand()) + offset;
}

INLINE double NoiseGenerator::getSampleLowpassWhite()
{
 return lowpassFilter.getSample(getSampleWhite());
}

INLINE double NoiseGenerator::getSampleHighpassWhite()
{
 return highpassFilter.getSample(getSampleWhite());
}

INLINE double NoiseGenerator::getSampleBandpassWhite()
{
 return highpassFilter.getSample(getSampleLowpassWhite());
}




INLINE double NoiseGenerator::getSamplePink()
{
 return scale * ((double) rand()) + offset;
}

INLINE double NoiseGenerator::getSampleLowpassPink()
{
 return lowpassFilter.getSample(getSamplePink());
}

INLINE double NoiseGenerator::getSampleHighpassPink()
{
 return highpassFilter.getSample(getSamplePink());
}

INLINE double NoiseGenerator::getSampleBandpassPink()
{
 return highpassFilter.getSample(getSampleLowpassPink());
}



#endif // NoiseGenerator_h
