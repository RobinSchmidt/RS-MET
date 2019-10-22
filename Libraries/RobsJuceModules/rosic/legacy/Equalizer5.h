#ifndef Equalizer5_h
#define Equalizer5_h

#include "AudioModule.h"
#include "BiquadCascade.h"
#include "BiquadDesigner.h"
#include "IirDesigner.h"
//#include "DirectFormIir.h"

/**

This class implements a 5-band parameteric equalizer with additional 
butterworth lowpass- and highpass-filters of variable order. The low- and 
high-frequency equalizer-band can be switched between first and second order
shelving and peak/bell characteristics.

*/

class Equalizer5 : public AudioModule 
{

public:

 /** This is an enumeration of the available filter modes for the equalizer 
     stages. */
 enum equalizerModes
 {
  BYPASS = 0,           // bypass
  PEAK,                 // peak/bell filter
  ONE_POLE_LOW_SHELV,   // first order low-shelving filter
  TWO_POLE_LOW_SHELV,   // second order low-shelving filter
  ONE_POLE_HIGH_SHELV,  // first order high-shelving filter
  TWO_POLE_HIGH_SHELV,  // second order high-shelving filter
  //BANDREJECT,           // second order band reject
  NOTCH                 // second order notch filter
 };

 /** This is an enumeration of the available filter modes for the lowpass and 
     highpass stages.  */
 enum filterModes
 {
  //BYPASS = 0,        // bypass
  BUTTERWORTH = 1,     // Butterworth characteristic
  LINKWITZ_RILEY,      // Linkwitz/Riley characteristic
 };

 enum stereoModes
 {
  MONO_5 = 0,
  STEREO_LINKED,
  STEREO_LR,
  STEREO_MS,
  MONO_10
 };

 //---------------------------------------------------------------------------
 // construction/destruction:

          Equalizer5();  ///< Constructor.
 virtual ~Equalizer5();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(flt64 newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual void setStereoMode(int64 newStereoMode);
 /**< Switches between the 4 available stereo-modes: mono, stereo-linked, 
      stereo and mid/side */

 virtual void setLpfOrder(int64 newLpfOrder, int64 channel);
 /**< Sets the order of the Butterworth-lowpass. */

 virtual void setLpfCutoff(flt64 newLpfCutoff, int64 channel);
 /**< Sets the cutoff-frequency of the Butterworth-lowpass. */

 virtual void setHpfOrder(int64 newHpfOrder, int64 channel);
 /**< Sets the order of the Butterworth-highpass. */

 virtual void setHpfCutoff(flt64 newHpfCutoff, int64 channel);
 /**< Sets the cutoff-frequency of the Butterworth-highpass. */

 virtual void setEqMode(int64 newMode, int64 stage, int64 channel);
 /**< Sets the mode for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side). See enum above for the available modes. */

 virtual void setEqFreq(flt64 newFreq, int64 stage, int64 channel);
 /**< Sets the frequency for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 virtual void setEqGain(flt64 newGain, int64 stage, int64 channel);
 /**< Sets the gain for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 virtual void setEqQ(flt64 newQ, int64 stage, int64 channel);
 /**< Sets the Q-factor for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 virtual void setLevel(flt64 newLevel, int64 channel);
 /**< Adjusts the global output-volume for each channel. Values are expected in
      decibels. */

 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE flt64 getSample(flt64 in);
 ///< Calculates one output sample at a time.

 virtual INLINE void getSampleFrameStereo(flt64* inL, 
                                          flt64* inR,
                                          flt64* outL,
                                          flt64* outR);
 ///< Calculates one streo output frame at a time.

 //---------------------------------------------------------------------------
 // others:

 virtual void getMagnitudeResponse(float *omegas, float *magnitudes, 
                                   int numBins, int channel);
 /** Calculates the magnitudes of the frequency-response of the "channel"-th
     channel at the normalized radian frequencies given in the array "omegas"
     and stores them in the array "magnitudes". Both arrays are assumed to be
     "numBins" long. */

 //---------------------------------------------------------------------------
 // embedded objects (in the order of the signal flow):


 //===========================================================================

protected:

 // number of equalizer-stages:
 static const int64 numStages  = 5;

 // number of biquad-stages in the highpass- and lowpass-filters:
 //static const int64 numFltStages = 4;

 // parameter variables:

 int64A stereoMode;           // the stereo mode (mono, stereo-linked, stereo,
                              // mid/side
 int64A eqMode[numStages][2]; // modes for the eq-bands. 1st index: stage(1-5)
                              // 2nd index: channel 0:left/mid 1:right/side
 flt64A eqFreq[numStages][2]; // frequencies for the eq-bands (in Hz).
 flt64A eqGain[numStages][2]; // gain-values for the eq-bands (in dB).
 flt64A eqQ[numStages][2];    // Q-factors for the eq-bands.

 int64A hpfOrder[2];          // orders for the HPFs. index indicates channel
 int64A hpfMode[2];           // modes of the HPFs 
 flt64A hpfCutoff[2];         // cutoff-frequencies of the HPFs

 int64A lpfOrder[2];          // orders for the LPFs. index indicates channel
 int64A lpfMode[2];           // modes of the LPFs 
 flt64A lpfCutoff[2];         // cutoff-frequencies of the LPFs

 flt64  volL, volR;           // volume-factors for left and right (or mid and
                              // side) channel

 // embedded objects:
 BiquadDesigner eqDesigners[numStages][2];
 IirDesigner    hpfDesigners[2];
 IirDesigner    lpfDesigners[2];

 BiquadCascade  eqChainL, eqChainR;
 BiquadCascade  hpfL, hpfR;
 BiquadCascade  lpfL, lpfR;

 /*
 DirectFormIir  hpfL, hpfR;
 DirectFormIir  lpfL, lpfR;
 */

 // biquad-coefficients for highpass-filters (2 channels):
 flt64 b0hpL[numStages];
 flt64 b1hpL[numStages];
 flt64 b2hpL[numStages];
 flt64 a1hpL[numStages];
 flt64 a2hpL[numStages];

 flt64A b0hpR[numStages];
 flt64A b1hpR[numStages];
 flt64A b2hpR[numStages];
 flt64A a1hpR[numStages];
 flt64A a2hpR[numStages];

 // biquad-coefficients for the 5 equalizer-stages (2 channels):
 flt64 b0eqL[numStages];
 flt64 b1eqL[numStages];
 flt64 b2eqL[numStages];
 flt64 a1eqL[numStages];
 flt64 a2eqL[numStages];

 flt64A b0eqR[numStages];
 flt64A b1eqR[numStages];
 flt64A b2eqR[numStages];
 flt64A a1eqR[numStages];
 flt64A a2eqR[numStages];

 // biquad-coefficients for lowpass-filters (2 channels):
 flt64 b0lpL[numStages];
 flt64 b1lpL[numStages];
 flt64 b2lpL[numStages];
 flt64 a1lpL[numStages];
 flt64 a2lpL[numStages];

 flt64A b0lpR[numStages];
 flt64A b1lpR[numStages];
 flt64A b2lpR[numStages];
 flt64A a1lpR[numStages];
 flt64A a2lpR[numStages];

 void resetHpfCoeffsL();
 void resetHpfCoeffsR();
 void resetEqCoeffsL();
 void resetEqCoeffsR();
 void resetLpfCoeffsL();
 void resetLpfCoeffsR();
};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE flt64 Equalizer5::getSample(flt64 in)
{
 return in;
}

INLINE void Equalizer5::getSampleFrameStereo(flt64* inL, 
                                             flt64* inR,
                                             flt64* outL,
                                             flt64* outR)
{
 static flt64A tmpL, tmpR; // these variables will contain the left and right
                           // output at the end
 static flt64A tmpM, tmpS; // these variables are for the mid- and side-signal 
                           // in the mid-side mode

 // write the inputs into the temporary (accumulating) variables and add some
 // anti-denorm signal to it (has to be different on both channels to work
 // also in M/S-mode):
 tmpL = *inL + getAntiDenormSaw(8, 1.0e-035);
 tmpR = *inR + getAntiDenormSaw(8, 0.5e-035);

 switch( stereoMode )
 {
  case MONO_5:
  {
   // apply the highpass:
   tmpL = hpfL.getSample(tmpL);

   // apply the equalizer-chain:
   tmpL = eqChainL.getSample(tmpL);

   // apply the lowpass and global volume:
   tmpL = volL * lpfL.getSample(tmpL);

   // copy the result into the right channel:
   tmpR = tmpL;
  }
  break;

  case STEREO_LINKED:
  {
   // apply the highpass:
   tmpL = hpfL.getSample(tmpL);
   tmpR = hpfR.getSample(tmpR);

   // apply the equalizer-chain:
   tmpL = eqChainL.getSample(tmpL);
   tmpR = eqChainR.getSample(tmpR);

   // apply the lowpass and global volume:
   tmpL = volL * lpfL.getSample(tmpL);
   tmpR = volL * lpfR.getSample(tmpR);

   // not implemented yet
   //tmpL = 0;
   //tmpR = 0;
  }
  break;

  case STEREO_LR:
  {
   // apply the highpass:
   tmpL = hpfL.getSample(tmpL);
   tmpR = hpfR.getSample(tmpR);

   // apply the equalizer-chain:
   tmpL = eqChainL.getSample(tmpL);
   tmpR = eqChainR.getSample(tmpR);

   // apply the lowpass and global volume:
   tmpL = volL * lpfL.getSample(tmpL);
   tmpR = volR * lpfR.getSample(tmpR);
  }
  break;

  case STEREO_MS:
  {
   // calculate the mid- and side-signal from the left- and right-signal:
   tmpM = SQRT2_INV * (tmpL + tmpR);
   tmpS = SQRT2_INV * (tmpL - tmpR);

   // apply the highpass:
   tmpM = hpfL.getSample(tmpM);
   tmpS = hpfR.getSample(tmpS);

   // apply the equalizer-chain:
   tmpM = eqChainL.getSample(tmpM);
   tmpS = eqChainR.getSample(tmpS);

   // apply the lowpass and global volume:
   tmpM = volL * lpfL.getSample(tmpM);
   tmpS = volR * lpfR.getSample(tmpS);

   // calculate the left- and right-signal from the mid- and side-signal:
   tmpL = SQRT2_INV * (tmpM + tmpS);
   tmpR = SQRT2_INV * (tmpM - tmpS);
  }
  break;

  case MONO_10:
  {
   // apply the highpasses:
   tmpL = hpfL.getSample(tmpL);
   tmpL = hpfR.getSample(tmpL);

   // apply the equalizer-chains:
   tmpL = eqChainL.getSample(tmpL);
   tmpL = eqChainR.getSample(tmpL);

   // apply the lowpasses and global volumes:
   tmpL = volL * lpfL.getSample(tmpL);
   tmpL = volR * lpfR.getSample(tmpL);

   // copy the result into the right channel:
   tmpR = tmpL;
  }
  break;

 } // end of switch(stereoMode)

 // store the output-signals in the slots:
 *outL = tmpL;
 *outR = tmpR;

 return;
}

#endif // Equalizer5_h