#ifndef HighOrderEqualizer_h
#define HighOrderEqualizer_h

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

class HighOrderEqualizer : public AudioModule 
{

public:

 /** This is an enumeration of the available filter modes for the equalizer 
     stages. */
 enum equalizerModes
 {
  BYPASS = 0,           // bypass
  LOW_SHELV,
  PEAK,
  NOTCH,
  HIGH_SHELV,
 };

 enum channels
 {
  LEFT  = 0,
  RIGHT
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

          HighOrderEqualizer();  ///< Constructor.
 virtual ~HighOrderEqualizer();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual void setStereoMode(int newStereoMode);
 /**< Switches between the 4 available stereo-modes: mono, stereo-linked, 
      stereo and mid/side */

 virtual void setLpfOrder(int newLpfOrder, int channel);
 /**< Sets the order of the Butterworth-lowpass. */

 virtual void setLpfCutoff(double newLpfCutoff, int channel);
 /**< Sets the cutoff-frequency of the Butterworth-lowpass. */

 virtual void setHpfOrder(int newHpfOrder, int channel);
 /**< Sets the order of the Butterworth-highpass. */

 virtual void setHpfCutoff(double newHpfCutoff, int channel);
 /**< Sets the cutoff-frequency of the Butterworth-highpass. */

 virtual void setEqMode(int newMode, int stage, int channel);
 /**< Sets the mode for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side). See enum above for the available modes. */

 virtual void setEqOrder(int newOrder, int stage, int channel);
 /**< Sets the order for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side). See enum above for the available modes. */

 virtual void setEqFreq(double newFreq, int stage, int channel);
 /**< Sets the frequency for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 virtual void setEqGain(double newGain, int stage, int channel);
 /**< Sets the gain for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 virtual void setEqQ(double newQ, int stage, int channel);
 /**< Sets the Q-factor for one equalizer-stage. "stage" indicates the index
      of the stage (0...4) and channel indicates the signal channel (0: left or 
      mid, 1: right or side) */

 virtual void setLevel(double newLevel, int channel);
 /**< Adjusts the global output-volume for each channel. Values are expected in
      decibels. */

 virtual void setBypassAll(bool newBypassAll);
 /**< Switches the whole equalizer on and off. */

 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE double getSample(double in);
 ///< Calculates one output sample at a time.

 virtual INLINE void getSampleFrameStereo(double* inL, 
                                          double* inR,
                                          double* outL,
                                          double* outR);
 ///< Calculates one streo output frame at a time.

 //---------------------------------------------------------------------------
 // others:

 //virtual void getMagnitudeResponse(float *omegas, float *magnitudes, 
   //                                int numBins, int channel);
 /** Calculates the magnitudes of the frequency-response of the "channel"-th
     channel at the normalized radian frequencies given in the array "omegas"
     and stores them in the array "magnitudes". Both arrays are assumed to be
     "numBins" long. */

 virtual void getMagnitudeResponse(int channel, int* arrayLengths, 
                                   float** freqArray, float** magArray);
 /** Returns the magnitude response for one of the channels. The number of
     bins will be returned in *arrayLenghts, the frequencies (in Hz) will be 
     returned in the array *freqArray and the magnitudes (in dB) will be 
     returned in the array *magArray. */

 //---------------------------------------------------------------------------
 // embedded objects (in the order of the signal flow):


 //===========================================================================

protected:

 static const int maxBiquadsInFilter = 4;
  // maximum number of biquad-stages in the highpass- and lowpass-filters

 static const int maxBiquadsInEq = 8;
  // maximum number of biquad-stages in the equalizer-filters

 static const int numEqStages = 5;
  // number of equalizer-stages in the equalizer-chain

 static const int numChannels = 2;
  // number of channels:

 static const int numBins = 1648;
  // the number of bins for the frequency responses

 // parameter variables:
 intA stereoMode; // the stereo mode (mono5, stereo-linked, 
                  // stereo L/R, stereo M/S, mono10

 intA eqMode[numEqStages][numChannels];  // modes for the eq-bands. 
  // 1st index: stage(0-4), 2nd index: channel 0:left/mid 1:right/side

 intA eqOrder[numEqStages][numChannels];   // orders for the eq-bands. 
  // 1st index: stage(0-4), 2nd index: channel 0:left/mid 1:right/side

 doubleA eqFreq[numEqStages][numChannels]; 
  // center-/corner- frequencies for the eq-bands (in Hz)

 doubleA eqGain[numEqStages][numChannels]; 
  // gain-values for the eq-bands (in dB)

 doubleA eqQ[numEqStages][numChannels];    
  // Q-factors for the eq-bands.

 intA    hpfOrder[numChannels];   // orders for the HPFs
 doubleA hpfCutoff[numChannels];  // cutoff-frequencies of the HPFs

 intA    lpfOrder[numChannels];   // orders for the LPFs
 doubleA lpfCutoff[numChannels];  // cutoff-frequencies of the LPFs

 doubleA vol[numChannels];  // volume-factors for left and right (or mid and
                            // side) channel

 bool bypassAll; // switch for a global bypass

 // embedded objects:
 IirDesigner designer; // this object will do all the filter-design

 // the biquad-cascade-filters for the channels:
 BiquadCascade hpf[numChannels]; 
 BiquadCascade eq[numEqStages][numChannels]; 
 BiquadCascade lpf[numChannels];

 // biquad-coefficients for the highpass-filters (2 channels):
 doubleA b0hp[numChannels][maxBiquadsInFilter];
 doubleA b1hp[numChannels][maxBiquadsInFilter];
 doubleA b2hp[numChannels][maxBiquadsInFilter];
 doubleA a1hp[numChannels][maxBiquadsInFilter];
 doubleA a2hp[numChannels][maxBiquadsInFilter];

 // biquad-coefficients for the equalizer-filters (2 channels):
 doubleA b0eq[numEqStages][numChannels][maxBiquadsInEq];
 doubleA b1eq[numEqStages][numChannels][maxBiquadsInEq];
 doubleA b2eq[numEqStages][numChannels][maxBiquadsInEq];
 doubleA a1eq[numEqStages][numChannels][maxBiquadsInEq];
 doubleA a2eq[numEqStages][numChannels][maxBiquadsInEq];

 // biquad-coefficients for the lowpass-filters (2 channels):
 doubleA b0lp[numChannels][maxBiquadsInFilter];
 doubleA b1lp[numChannels][maxBiquadsInFilter];
 doubleA b2lp[numChannels][maxBiquadsInFilter];
 doubleA a1lp[numChannels][maxBiquadsInFilter];
 doubleA a2lp[numChannels][maxBiquadsInFilter];

 // member variables which hold the frequency-values for the magnitude 
 // response and some related quantities:
 doubleA freqs[numBins];      // bin frequencies in Hz
 doubleA omegas[numBins];     // normalized radian bin frequencies
 doubleA cosOmegas[numBins];  // cosines of the omegas[]
 doubleA cos2Omegas[numBins]; // cosines of the 2*omegas[]

 // member variables which hold the magnitudes-responses of the individual
 // biquad-cascade-filters:
 doubleA hpfResponses[numChannels][numBins];
 doubleA eqResponses[numEqStages][numChannels][numBins];
 doubleA lpfResponses[numChannels][numBins];

 // member-variables which hold the magnitude response of the whole filter
 // chain (for each channel, as floats for display in the JUCE-GUI):
 float floatFrequencies[numBins];
 float floatMagnitudes[numChannels][numBins];

 // functions for resetting the coefficients to neutral values:
 void resetHpfCoeffs(int channel);
 void resetEqCoeffs (int stage, int channel);
 void resetLpfCoeffs(int channel);

 // functions, that calculate the new filter coefficients and pass them to
 // the appropriate BiquadCascade-objects:
 void updateHpfCoeffs(int channel);
 void updateEqCoeffs (int stage, int channel);
 void updateLpfCoeffs(int channel);

 // functions, that calculate the new filter magnitude responses:
 void updateHpfFreqResponse(int channel);
 void updateEqFreqResponse (int stage, int channel);
 void updateLpfFreqResponse(int channel);
 void updateOverallFreqResponse(int channel);

 // we use some multi-dimensional arrays in this class. the hierarchy of the 
 // indices is as follows: 
 //  -1st index: equalizer-stage (for the eq-array)
 //  -2nd index: the channel
 //  -3rd index: the biquad-stage or the frequency bin
 // not all of these hierarchy-levels need to be present - if some level is 
 // absent, the others are shifted to higher levels accordingly

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double HighOrderEqualizer::getSample(double in)
{
 return in;
}

INLINE void HighOrderEqualizer::getSampleFrameStereo(double* inL, 
                                                     double* inR,
                                                     double* outL,
                                                     double* outR)
{
 static doubleA tmpL, tmpR; // these variables will contain the left and right
                            // output at the end
 static doubleA tmpM, tmpS; // these variables are for the mid- and side-signal 
                            // in the mid-side mode
 static intA    s;

 // just copy the inputs into the outputs and return in bypass-mode:
 if( bypassAll )
 {
  *outL = *inL;
  *outR = *inR;
  return;
 }

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
   tmpL = hpf[0].getSampleDirect1(tmpL);

   // apply the equalizer-chain:
   for(s=0; s<numEqStages; s++)
    tmpL = eq[s][0].getSampleDirect1(tmpL); 

   // apply the lowpass and global volume:
   tmpL = vol[0] * lpf[0].getSampleDirect1(tmpL);

   // copy the result into the right channel:
   tmpR = tmpL;
  }
  break;

  case STEREO_LINKED:
  {
   // apply the highpass:
   tmpL = hpf[0].getSampleDirect1(tmpL);
   tmpR = hpf[1].getSampleDirect1(tmpR);

   // apply the equalizer-chain:
   for(s=0; s<numEqStages; s++)
   {
    tmpL = eq[s][0].getSampleDirect1(tmpL); 
    tmpR = eq[s][1].getSampleDirect1(tmpR); 
   }

   // apply the lowpass and global volume:
   tmpL = vol[0] * lpf[0].getSampleDirect1(tmpL);
   tmpR = vol[0] * lpf[1].getSampleDirect1(tmpR);
  }
  break;

  case STEREO_LR:
  {
   // apply the highpass:
   tmpL = hpf[0].getSampleDirect1(tmpL);
   tmpR = hpf[1].getSampleDirect1(tmpR);

   // apply the equalizer-chain:
   for(s=0; s<numEqStages; s++)
   {
    tmpL = eq[s][0].getSampleDirect1(tmpL); 
    tmpR = eq[s][1].getSampleDirect1(tmpR); 
   }

   // apply the lowpass and global volume:
   tmpL = vol[0] * lpf[0].getSampleDirect1(tmpL);
   tmpR = vol[1] * lpf[1].getSampleDirect1(tmpR);
  }
  break;

  case STEREO_MS:
  {
   // calculate the mid- and side-signal from the left- and right-signal:
   tmpM = SQRT2_INV * (tmpL + tmpR);
   tmpS = SQRT2_INV * (tmpL - tmpR);

   // apply the highpass:
   tmpM = hpf[0].getSampleDirect1(tmpM);
   tmpS = hpf[1].getSampleDirect1(tmpS);

   // apply the equalizer-chain:
   for(s=0; s<numEqStages; s++)
   {
    tmpM = eq[s][0].getSampleDirect1(tmpM); 
    tmpS = eq[s][1].getSampleDirect1(tmpS); 
   }

   // apply the lowpass and global volume:
   tmpM = vol[0] * lpf[0].getSampleDirect1(tmpM);
   tmpS = vol[1] * lpf[1].getSampleDirect1(tmpS);

   // calculate the left- and right-signal from the mid- and side-signal:
   tmpL = SQRT2_INV * (tmpM + tmpS);
   tmpR = SQRT2_INV * (tmpM - tmpS);
  }
  break;

  case MONO_10:
  {
   // apply the highpasses:
   tmpL = hpf[0].getSampleDirect1(tmpL);
   tmpL = hpf[1].getSampleDirect1(tmpL);

   // apply the equalizer-chains:
   for(s=0; s<numEqStages; s++)
   {
    tmpL = eq[s][0].getSampleDirect1(tmpL); 
    tmpL = eq[s][1].getSampleDirect1(tmpL); 
   }

   // apply the lowpasses and global volumes:
   tmpL = vol[0] * lpf[0].getSampleDirect1(tmpL);
   tmpL = vol[1] * lpf[1].getSampleDirect1(tmpL);

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

#endif // HighOrderEqualizer_h