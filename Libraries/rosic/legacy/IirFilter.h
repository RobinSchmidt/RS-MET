#ifndef IirFilter_h
#define IirFilter_h

//#include "AudioModule.h"
#include "BiquadCascade.h"
#include "IirDesigner.h"

/**

This class implements an infinite impulse response filter by utilizing objects
of class "IirDesigner" and "BiquadCascade".

*/

class IirFilter : public AudioModule
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          IirFilter();  ///< Constructor.     
 virtual ~IirFilter();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual void setMode(int newMode);       
 ///< Chooses one of the modes (LP, HP, BP, BR) as enumerated above.

 virtual void setSlope(int newSlope);      
 /** Determines the order of the filter - 1 results in 6 dB/oct slope, 
     2->12dB/oct 3->18, 4->24 and so on. The actual order of the filter
     depends on its characteristic: in LPF's and HPF's the order is
     identical to the value of "slope", in BPFs and BRFs the order is
     twice as high. */

 virtual void setFreq1(double newFreq1);      
 /**< Sets the cutoff frequency for lowpass- and highpass-filters, or the 
      lower cutoff frequency for bandpass- and bandreject-filters. */

 virtual void setFreq2(double newFreq2);
 /**< Sets the upper cutoff frequency for bandpass- and bandreject-filters. */

 virtual void setEqGain(double newEqGain);
 /**< Sets the gain for the shelving and peaking modes. The value is expected 
      as raw multiplicative value - not in decibels. */

 virtual void setGlobalGain(double newGlobalGain);
 /**< Sets an overall gain factor for the signal. */


 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE double getSample(double in); 
 ///< Calculates a single filtered output-sample.

 //---------------------------------------------------------------------------
 // others:

 virtual void resetBuffers (); 
 /**< Sets the buffers for the previous input and output samples of all biquad
      stages to zero. */

 //===========================================================================

protected:

 // embedded audio-modules:
 BiquadCascade filter;
 IirDesigner   designer;

 // maximum number of biquad-stages:
 static const int maxNumStages = 12;

 // filter parameters:
 intA    mode, slope;
 doubleA freq1, freq2, eqGain;

 // filter coefficients:
 doubleA globalGain;       // a global gain-factor
 doubleA a0[maxNumStages]; // a0-coefficients for the individual stages
 doubleA a1[maxNumStages]; // a1-coefficients for the individual stages
 doubleA a2[maxNumStages]; // a2-coefficients for the individual stages
 doubleA b0[maxNumStages]; // b0-coefficients for the individual stages
 doubleA b1[maxNumStages]; // b1-coefficients for the individual stages
 doubleA b2[maxNumStages]; // b2-coefficients for the individual stages

 void resetBiquadCoeffs();
 void updateBiquadCoeffs();
};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double IirFilter::getSample(double in)
{
 return filter.getSampleDirect1(in);
}

#endif // IirFilter_h
