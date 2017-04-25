#ifndef MagicCarpetEqualizer_h
#define MagicCarpetEqualizer_h

#include "BiquadCascade.h"
#include "BiquadDesigner.h"
#include "IirDesigner.h"

/**

This class implements a 3-band parameteric equalizer with additional first order
lowpass- and highpass-filters.

*/

class MagicCarpetEqualizer
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 MagicCarpetEqualizer();  ///< Constructor.
 ~MagicCarpetEqualizer();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 void setEqFreq(double newFreq, int stage);
 /**< Sets the frequency for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 void setEqGain(double newGain, int stage);
 /**< Sets the gain for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 void setEqQ(double newQ, int stage);
 /**< Sets the Q-factor for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double* inL, 
                                  double* inR,
                                  double* outL,
                                  double* outR);
 ///< Calculates one streo output frame at a time.

 //---------------------------------------------------------------------------
 // others:


 //---------------------------------------------------------------------------
 // embedded objects (in the order of the signal flow):


 //===========================================================================

protected:

 // number of equalizer-stages:
 static const int numStages  = 3;

 // parameter variables:
 doubleA eqFreq[numStages]; // frequencies for the eq-bands (in Hz).
 doubleA eqGain[numStages]; // gain-values for the eq-bands (in dB).
 doubleA eqQ[numStages];    // Q-factors for the eq-bands.

 // embedded objects:
 BiquadDesigner eqDesigners[numStages];

 BiquadCascade  eqChainL, eqChainR;

 // biquad-coefficients for the 3 equalizer-stages:
 double b0[numStages];
 double b1[numStages];
 double b2[numStages];
 double a1[numStages];
 double a2[numStages];

 bool   isOff;

 void resetEqCoeffs();

 double sampleRate;
};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE void MagicCarpetEqualizer::getSampleFrameStereo(double* inL, 
                                                       double* inR,
                                                       double* outL,
                                                       double* outR)
{
 // just copy input into output when we are in bypass-mode:
 if( isOff == true )
 {
  *outL = *inL;
  *outR = *inR;
  return;
 }

 static doubleA tmpL, tmpR; // these variables will contain the left and right
                            // output at the end

 // write the inputs into the temporary (accumulating) variables and add some
 // anti-denorm signal to it (has to be different on both channels to work
 // also in M/S-mode):
 tmpL = *inL; // + getAntiDenormSaw(8, 1.0e-035);
 tmpR = *inR; // + getAntiDenormSaw(8, 0.5e-035);

 // apply the equalizer-chain:
 tmpL = eqChainL.getSampleDirect1(tmpL);
 tmpR = eqChainR.getSampleDirect1(tmpR);

 // store the output-signals in the slots:
 *outL = tmpL;
 *outR = tmpR;

 return;
}

#endif // MagicCarpetEqualizer_h