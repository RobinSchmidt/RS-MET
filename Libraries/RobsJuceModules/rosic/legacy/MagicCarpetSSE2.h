#ifndef MagicCarpetSSE2_h
#define MagicCarpetSSE2_h

#include "Definitions.h"
#include "FourierTransformer.h"
#include "MoreMath.h"
#include "AudioModule.h"
//#include "LowpassHighpassStereoSSE2.h"
//#include "MagicCarpetEqualizerSSE2.h"
#include "MagicCarpetModulator.h"
//#include "MagicCarpetVoiceSSE2.h"

/**

MagicCarpetSSE2 is a synthesizer based on a mix of special audio-samples. blah...

*/


class MagicCarpetSSE2
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 MagicCarpetSSE2();   ///< Constructor.
 ~MagicCarpetSSE2();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Sets the sample-rate.

 void setNumVoices(int newNumVoices); 
 ///< Sets the maximum number of voices (polyphony).

 void setWaveTable(int     destinationSlot,
                   double* newWaveTable, 
                   int     newTableLength);

 void setSpectrum(int     destinationSlot,
                  double* newSpectrum, 
                  int     newSpectrumLength);

 // switches for sample-playback options:
 void setSlotLoopMode       (int whichSlot, bool newLoopModeSwitch);
 void setSlotSingleCycleMode(int whichSlot, bool newSingleCycleSwitch);
 void setSlotStereoMode     (int whichSlot, bool newStereoSwitch);
 void setSlotMute           (int whichSlot, bool newMuteSwitch);

 // numeric parameters for the sample-playback:
 void setSlotDetune    (int whichSlot, double newDetune);
 void setSlotTuneByKey (int whichSlot, double newTuneByKey);
 void setSlotLevel     (int whichSlot, double newLevel);
 void setSlotLevelByKey(int whichSlot, double newLevelByKey);
 void setSlotLevelByVel(int whichSlot, double newLevelByVel);
 void setSlotHpf       (int whichSlot, double newHpfCutoff);
 void setSlotLpf       (int whichSlot, double newLpfCutoff);
 void setSlotRootKey   (int whichSlot, double newRootKey);

 // parameters for the synthesis engine:
 void setXModNumerator   (int    newXModNumerator);
 void setXModDenominator (int    newXModDenominator);
 void setXModStartPhase  (double newXModStartPhase);
 void setXModAttack      (double newXModAttack);
 void setXModRelease     (double newXModRelease);
 void setXModRise        (double newXModRise);
 void setXModAmount      (double newXModAmount);
 void setXModBpm         (double newXModBpm);
 void setXModLoopMode    (bool   newXModLoopModeSwitch);

 void setYModNumerator   (int    newYModNumerator);
 void setYModDenominator (int    newYModDenominator);
 void setYModStartPhase  (double newYModStartPhase);
 void setYModAttack      (double newYModAttack);
 void setYModRelease     (double newYModRelease);
 void setYModRise        (double newYModRise);
 void setYModAmount      (double newYModAmount);
 void setYModBpm         (double newYModBpm);
 void setYModLoopMode    (bool   newYModLoopModeSwitch);

 // the general filter envelope:
 //void switchFilterOn    (bool   shouldBeOn); 
 void setFltMode        (int    newFltMode);
 void setFltFreq        (double newFltFreq);
 void setFltFreqByKey   (double newFltFreqByKey);
 void setFltFreqByVel   (double newFltFreqByVel);
 void setFltReso        (double newFltReso);
 void setFltGain        (double newFltGain);
 void setFltOrder       (int    newFltOrder);

 // the filter envelope-settings:
 void setFltEnvAttack   (double newFltEnvAttack);
 void setFltEnvPeak     (double newFltEnvPeak);
 void setFltEnvPeakByVel(double newFltEnvPeakByVel);
 void setFltEnvDecay    (double newFltEnvDecay);
 void setFltEnvRelease  (double newFltEnvRelease);
 void setFltEnvEnd      (double newFltEnvEnd);
 void setFltEnvDurByKey (double newFltEnvDurByKey);
 void setFltEnvDurByVel (double newFltEnvDurByVel);
 void setFltEnvSlope    (double newFltEnvSlope);

 // the amplitude and vector-mixer settings:
 void setMasterVolume   (double newMasterVolume);  
 void setMidSideMix     (double newMidSideMix); 
 void setX              (double newX);  
 void setY              (double newY);  

 // the amplitude envelope-settings:
 void setAmpEnvAttack   (double newAmpEnvAttack);
 void setAmpEnvPeak     (double newAmpEnvPeak);
 void setAmpEnvPeakByVel(double newAmpEnvPeakByVel);
 void setAmpEnvDecay    (double newAmpEnvDecay);
 void setAmpEnvSustain  (double newAmpEnvSustain);
 void setAmpEnvRelease  (double newAmpEnvRelease);
 void setAmpEnvDurByKey (double newAmpEnvByKey);
 void setAmpEnvDurByVel (double newAmpEnvByVel);
 void setAmpEnvSlope    (double newAmpEnvSlope);

 // the output equalizer-settings:
 void switchOutFilter   (bool   shouldBeOn);
 void setOutHpfFreq     (double newOutHpfFreq);
 void setOutEqFreq     (double newOutEqFreq, int whichBand);
 void setOutEqQ        (double newOutEqQ,    int whichBand);
 void setOutEqGain     (double newOutEqGain, int whichBand);
 void setOutLpfFreq     (double newOutLpfFreq);

 //---------------------------------------------------------------------------
 // event processing:

 void noteOn(long NoteNumber, long Velocity);

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double *outL, double *outR);   
 /**< Calculates the output-samples for both channels and stores them at the
      adresses of *outL and *outR. */

//============================================================================

protected:

 // embedded audio-modules:
 FourierTransformer        fourierTransformer;
 MagicCarpetVoiceSSE2      magicCarpetVoiceArray[MAXVOICES];
 LowpassHighpassStereoSSE2 outputFilter;
 MagicCarpetEqualizerSSE2  outputEqualizer;
 MagicCarpetModulator      xModulator, yModulator;

 // parameters:
 double sampleRate;
 int    numVoices;
 double x,y;        // values of the X,Y-parameters between -1.0...1.0
 double masterAmplitude;
 double midSideMix;
 double x1L, x1R;   // state variables for the differencing filter

 // MIDI-data and voice-management:
 int mostRecentNote, mostRecentNoteVel, mostRecentNoteDetune;
 //....

 bool   outFilterIsOn;

 // allocate the arrays for the actual wavetables (the +6 is needed for
 // the interpolator interpolator):
 double topLeftTable[MAXTABLELENGTH+6];
 double bottomLeftTable[MAXTABLELENGTH+6];
 double topRightTable[MAXTABLELENGTH+6];
 double bottomRightTable[MAXTABLELENGTH+6];

 double xModTable[MAXTABLELENGTH+6];
 double yModTable[MAXTABLELENGTH+6];

 // temporary arrays for magnitude spectrum and wave-tabel to work with:
 double tmpSpectrum[MAXTABLELENGTH];
 double tmpTable[MAXTABLELENGTH];


 int    randomSeed;

 void zeroAllTables();
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE void MagicCarpetSSE2::getSampleFrameStereo(double *outL, double *outR)
{
 double xMapped, yMapped; 
 double leftAmp, rightAmp, topAmp, bottomAmp;
 double topLeftAmp, topRightAmp, bottomLeftAmp, bottomRightAmp;
 double accuL, accuR;
 double left, right, mid, side, sine, cosine;
 int    i;

 // calculate the amplitudes for the 4 slots from the x,y-coordinates
 // (optimize this later!):

 // add the modulators:
 xMapped = x + xModulator.getSample();
 yMapped = y + yModulator.getSample();

 // now map the range -1.0...+1.0 to the range 0.0...1.0:
 xMapped  = 0.5*xMapped + 0.5; 
 yMapped  = 0.5*yMapped + 0.5;

 // now map the range 0.0...1.0 to the range 0.0...PI/2:
 xMapped *= 0.5*PI;
 yMapped *= 0.5*PI;
  // absorb into the mapping above

 // calculate the resulting slot-amplitudes:
 /*
 leftAmp        = fabs(cos(xMapped));
 rightAmp       = fabs(sin(xMapped));
 bottomAmp      = fabs(cos(yMapped));
 topAmp         = fabs(sin(yMapped));
 */
 sinCosApprox(xMapped, &rightAmp, &leftAmp);
 sinCosApprox(yMapped, &topAmp,   &bottomAmp);
 topLeftAmp     = topAmp*leftAmp;
 topRightAmp    = topAmp*rightAmp;
 bottomLeftAmp  = bottomAmp*leftAmp;
 bottomRightAmp = bottomAmp*rightAmp;

 // initialize the output slots with zero in order to accumulate over the 
 // voices:
 accuL = 0.0;
 accuR = 0.0;

 // loop through the (active) voices:
 __m128d outVector = _mm_set1_pd(0.0); // accumulator
 for(i=0; i<numVoices; i++)
 {
  if(!magicCarpetVoiceArray[i].ampEnv.endIsReached())
  {
   // set up the amplitude factors for the slots:
   magicCarpetVoiceArray[i].topLeftAmplitude     = topLeftAmp;
   magicCarpetVoiceArray[i].topRightAmplitude    = topRightAmp;
   magicCarpetVoiceArray[i].bottomLeftAmplitude  = bottomLeftAmp;
   magicCarpetVoiceArray[i].bottomRightAmplitude = bottomRightAmp;

   // let voice i generate its output-sample and accumulate it:
   outVector = _mm_add_pd(outVector, 
                          magicCarpetVoiceArray[i].getSampleVector());
  }
 }

 // apply master equalizers, filters and volume control:
 outVector = _mm_mul_pd(outVector, _mm_set1_pd(masterAmplitude));

 // apply the lowpass/highpass filter:
 outVector = outputFilter.getSampleVector(outVector);

 // apply the 3 band equalizer:
 outVector = outputEqualizer.getSampleVector(outVector);

 // access the two double values in the vector-variable via pointer-trickery: 
 left  = *  (double*)(&outVector);
 right = * ((double*)(&outVector) + 1); 

 // apply mid/side adjustment:
 mid   = SQRT2_INV * (left + right);
 side  = SQRT2_INV * (left - right);
 //mid  *= cos(0.5*PI*midSideMix);
 //side *= sin(0.5*PI*midSideMix);
 //tmp   = 0.5*PI*midSideMix;
 sinCosApprox(0.5*PI*midSideMix, &sine, &cosine);
 mid  *= cosine;
 side *= sine;
 left  = SQRT2_INV * ( mid + side );
 right = SQRT2_INV * ( mid - side);

 //write the final output samples into the memory-slots:
 *outL = left;
 *outR = right; 
}

#endif // MagicCarpetSSE2_h
