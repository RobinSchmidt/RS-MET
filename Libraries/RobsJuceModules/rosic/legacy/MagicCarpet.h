#ifndef MagicCarpet_h
#define MagicCarpet_h

#include "Definitions.h"
#include "FourierTransformer.h"
#include "MoreMath.h"
#include "AudioModule.h"
#include "LowpassHighpassStereo.h"
#include "MagicCarpetEqualizer.h"
#include "MagicCarpetModulator.h"
#include "MagicCarpetVoice.h"
#include "OnePoleFilter.h"

/**

MagicCarpet is a synthesizer based on a mix of special audio-samples. blah...

*/


class MagicCarpet
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 MagicCarpet();   ///< Constructor.
 ~MagicCarpet();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Sets the sample-rate.

 void setNumVoices(int newNumVoices); 
 ///< Sets the maximum number of voices (polyphony).

 int getNumActiveVoices(); 
 ///< Returns the number of currently playnig voices (polyphony).

 //void clearWaveTables(int destinationSlot);

 void invalidateTablePointers(int whichSlot);
 /**< This function is called from the plugIn before a new table is passed to
      set the pointers temporarily to NULL. */

 bool setWaveTables(int     destinationSlot,
                    double* newWaveTableL, 
                    double* newWaveTableR,
                    int     newTableLength);
 /**< Passes a pointers to the wave-tables which should be used for left and 
      right channel. The return-value reporst, if the required memory could be
      allocated.*/

 bool setSpectrum(int     destinationSlot,
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
 void setFltTwoStages   (bool   newFltTwoStagesSwitch);
 void setFltFreq        (double newFltFreq);
 void setFltFreqByKey   (double newFltFreqByKey);
 void setFltFreqByVel   (double newFltFreqByVel);
 void setFltReso        (double newFltReso);
 void setFltGain        (double newFltGain);

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
 void setMasterVolume        (double newMasterVolume);  
 void setMasterVolumeByVoices(double newMasterVolumeByVoices);
 void setMidSideMix          (double newMidSideMix); 
 void setX                   (double newX);  
 void setY                   (double newY);  

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
 void setPitchBend(double transpositionInSemitones);

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double *outL, double *outR);   
 /**< Calculates the output-samples for both channels and stores them at the
      adresses of *outL and *outR. */

//============================================================================

protected:

 // embedded audio-modules:
 FourierTransformer    fourierTransformer;
 MagicCarpetVoice      magicCarpetVoiceArray[MAXVOICES];
 LowpassHighpassStereo outputFilter;
 MagicCarpetEqualizer  outputEqualizer;
 MagicCarpetModulator  xModulator, yModulator;
 OnePoleFilter         xSmoother, ySmoother;

 // parameters:
 double sampleRate;
 int    numVoices;
 double x,y;        // values of the X,Y-parameters between -1.0...1.0
 double masterAmplitude;
 double masterAmplitudeByVoices;
 double midSideMix;
 //double x1L, x1R;   // state variables for the differencing filter

 // MIDI-data and voice-management:
 int mostRecentNote, mostRecentNoteVel, mostRecentNoteDetune;
 //....

 bool   outFilterIsOn;

 // allocate the arrays for the actual wavetables (the +6 is needed for
 // the interpolator):
 double* topLeftTableL;
 double* topLeftTableR;
 double* topRightTableL;
 double* topRightTableR;
 double* bottomLeftTableL;
 double* bottomLeftTableR;
 double* bottomRightTableL;
 double* bottomRightTableR;

 double* xModTable;
 double* yModTable;

 /*
 double topLeftTableL[MAXTABLELENGTH+6];
 double topLeftTableR[MAXTABLELENGTH+6];
 double topRightTableL[MAXTABLELENGTH+6];
 double topRightTableR[MAXTABLELENGTH+6];
 double bottomLeftTableL[MAXTABLELENGTH+6];
 double bottomLeftTableR[MAXTABLELENGTH+6];
 double bottomRightTableL[MAXTABLELENGTH+6];
 double bottomRightTableR[MAXTABLELENGTH+6];

 double xModTable[MAXTABLELENGTH+6];
 double yModTable[MAXTABLELENGTH+6];

 // temporary arrays for magnitude spectrum and wave-tabel to work with:
 double tmpSpectrum[MAXTABLELENGTH];
 double tmpTableL[MAXTABLELENGTH];
 double tmpTableR[MAXTABLELENGTH];
 */

 int    randomSeed;
 int    numActiveVoices; // number of currently playing voices

 //void zeroAllTables();
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE void MagicCarpet::getSampleFrameStereo(double *outL, double *outR)
{
 double xMapped, yMapped; 
 double leftAmp, rightAmp, topAmp, bottomAmp;
 double topLeftAmp, topRightAmp, bottomLeftAmp, bottomRightAmp;
 double accuL, accuR;
 double voiceAmplitudeSum;
 double left, right, mid, side, sine, cosine;
 int    i;

 // calculate the amplitudes for the 4 slots from the x,y-coordinates
 // (optimize this later!):

 // add the modulators:
 xMapped = xSmoother.getSample(x) + xModulator.getSample();
 yMapped = ySmoother.getSample(y) + yModulator.getSample();

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
 // voices and then loop through the (active) voices:
 accuL             = 0.0;
 accuR             = 0.0;
 numActiveVoices   = 0;
 voiceAmplitudeSum = 0.0;
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
   magicCarpetVoiceArray[i].getSampleFrameStereo(&accuL, &accuR,
                                                 &voiceAmplitudeSum);

   // increment the variable which keeps track of the active voices:
   numActiveVoices++;
  }
 }

 left  = accuL;
 right = accuR;

 // apply master equalizers, filters and volume control:
 left  *= masterAmplitude;
 right *= masterAmplitude;

 // apply the amplitude-scaling by the number of voices (weigthed by their 
 // respective amplitudes):
 double voiceScaleFactor = (1.0 + masterAmplitudeByVoices) / 
                           (1.0 + masterAmplitudeByVoices * sqrt(voiceAmplitudeSum) );

 /*
 double voiceScaleFactor = (1.0 + masterAmplitudeByVoices) / 
                           (1.0 + masterAmplitudeByVoices * voiceAmplitudeSum ); */
                         
 left  *= voiceScaleFactor;
 right *= voiceScaleFactor;

 // apply the lowpass/highpass filter:
 outputFilter.getSampleFrameStereo(&left, &right, &left, &right);

 // apply the 3 band equalizer:
 outputEqualizer.getSampleFrameStereo(&left, &right, &left, &right);

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

 // write the final output samples into the memory-slots:
 *outL = left;
 *outR = right; 
}

#endif // MagicCarpet_h
