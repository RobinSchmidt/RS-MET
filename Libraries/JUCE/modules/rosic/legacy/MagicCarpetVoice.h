#ifndef MagicCarpetVoice_h
#define MagicCarpetVoice_h

#include "stdlib.h"
#include "math.h"
#include "AmpEnvRc.h"
//#include "CookbookFilter.h"
#include "MagicCarpetFilter.h"
#include "PitchEnvRc.h"
#include "Definitions.h"
#include "MoreMath.h"
#include "SampleOscillator.h"
#include "EllipticHalfbandFilter.h"

/**

This class realizes a single voice for the MagicCarpet-Synthesizer.

*/

class MagicCarpetVoice  
{

public:

 /** This is an enumeration of the slots for the 6 samples. */
 /*
 enum sampleSlots
 {
  TOP_LEFT_SAMPLE = 0,
  TOP_RIGHT_SAMPLE,
  BOTTOM_LEFT_SAMPLE,
  BOTTOM_RIGHT_SAMPLE,
  X_LFO_SAMPLE,
  Y_LFO_SAMPLE
 };
 */

 //---------------------------------------------------------------------------
 // construction/destruction:

	MagicCarpetVoice();
	~MagicCarpetVoice();

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate      (double newSampleRate); // set the sample rate
 //void setStereo          (bool   Stereo);        // stereo switch
 //void setVol             (double Vol);           // master volume control

 /*
 // key- and velocity sensitivities of the slot-volumes:
 void setTopLeftVolByKey    (double newVolByKey);
 void setTopLeftVolByVel    (double newVolByVel);
 void setTopRightVolByKey   (double newVolByKey);
 void setTopRightVolByVel   (double newVolByVel);
 void setBottomLeftVolByKey (double newVolByKey);
 void setBottomLeftVolByVel (double newVolByVel);
 void setBottomRightVolByKey(double newVolByKey);
 void setBottomRightVolByVel(double newVolByVel);
 */

 // the filter-envelope settings:
 void setFltEnvAttack      (double newFltEnvAttack);
 void setFltEnvPeak        (double newFltEnvPeak);
 void setFltEnvPeakByVel   (double newFltEnvPeakByVel);
 void setFltEnvDecay       (double newFltEnvDecay);
 void setFltEnvSustain     (double newFltEnvSustain);
 void setFltEnvRelease     (double newFltEnvRelease);
 void setFltEnvByKey       (double newFltEnvByKey);
 void setFltEnvByVel       (double newFltEnvByVel);
 void setFltEnvSlope       (double newFltEnvSlope);

 // the amplitude-envelope settings:
 void setAmpEnvAttack      (double newAmpEnvAttack);
 void setAmpEnvPeak        (double newAmpEnvPeak);
 void setAmpEnvPeakByVel   (double newAmpEnvPeakByVel);
 void setAmpEnvDecay       (double newAmpEnvDecay);
 void setAmpEnvSustain     (double newAmpEnvSustain);
 void setAmpEnvRelease     (double newAmpEnvRelease);
 void setAmpEnvByKey       (double newAmpEnvByKey);
 void setAmpEnvByVel       (double newAmpEnvByVel);
 void setAmpEnvSlope       (double newAmpEnvSlope);


 void resetDifferentiator();

 // vector mixing parameter settings:
 /*
 INLINE void setTopLeftVolume     (double newTopLeftVolume);
 INLINE void setTopRightVolume    (double newTopRightVolume);
 INLINE void setBottomLeftVolume  (double newBottomLeftVolume);
 INLINE void setBottomRightVolume (double newBottomRightVolume);
 */

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double *outL, double *outR, 
                                  double *voiceAmplitude);   
 /**< Calculates the output-samples for both channels and adds them at the
      adresses of *outL and *outR. The third output tells the current 
      amplitude of the voice (i.e. the output of the amp-envelope.) */

 //---------------------------------------------------------------------------
 // event processing:

 void noteOn(int NoteNumber, int Velocity);

 //---------------------------------------------------------------------------
 // public member variables (to be accessed by the outlying MagicCarpet 
 // class):

 // MIDI-data:
 int  currentNote, 
      currentNoteVel, 
      currentNoteDetune;
 int  currentNoteAge;  // age of the note in samples (for last note priority
                       // voice assignment)
 bool isPlaying;

 // parameters:
 //double topLeftDetune, topRightDetune, bottomLeftDetune, bottomRightDetune;

 // embedded public audio-modules (they have to be accessed by the outlying
 // MagicCarpet class):
 SampleOscillator  topLeftSource, 
                   topRightSource,
                   bottomLeftSource, 
                   bottomRightSource;
 PitchEnvRc        fltEnv;  
 MagicCarpetFilter filter;
 AmpEnvRc          ampEnv;  

 //===========================================================================

 // nominal source amplitudes (without velocity and key-scaling):
 double topLeftAmplitude;
 double topRightAmplitude;
 double bottomLeftAmplitude;
 double bottomRightAmplitude;

 // key-dependencies:
 double topLeftTuneByKey;
 double topRightTuneByKey;
 double bottomLeftTuneByKey;
 double bottomRightTuneByKey;
 double topLeftLevelByKey;
 double topRightLevelByKey;
 double bottomLeftLevelByKey;
 double bottomRightLevelByKey;
 //double fltFrqByKey;

 double ampEnvPeakByVel;
 double ampEnvDurByKey;
 double ampEnvDurByVel;
 double fltEnvDurByKey;
 double fltEnvDurByVel;

 // velocity-dependencies:
 double topLeftLevelByVel;
 double topRightLevelByVel;
 double bottomLeftLevelByVel;
 double bottomRightLevelByVel;
 double ampByVel;
 //double fltFrqByVel;
 double fltEnvPeakByVel;

 // filter settings:
 bool   filterIsOn;
 double fltFreq;
 double fltFreqByKey;
 double fltFreqByVel;
 double fltFreqScaler;
 //double fltFreqScaled; // includes key and velocity-scaling


 //double fltReso;
 //double fltGain;

 //===========================================================================

protected:

 // embedded protected audio-modules:


 //EllipticHalfbandFilter antiAliasFilterL, antiAliasFilterR;

 // parameters:
 double sampleRate;        // sample-rate
 //int    tableLength;
 bool   stereo;

 // state samples x[n-1] for the differentiators:
 double x1L, x1R;

 // internal functions:
 INLINE void getSourceSampleFrameStereo(double *outL, double *outR);  
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

//----------------------------------------------------------------------------
// audio processing:

INLINE void MagicCarpetVoice::getSampleFrameStereo(double *outL, double *outR, 
                                                   double *voiceAmplitude)
{
 /*
 if(ampEnv.endIsReached()) 
  return; // this voice is dead and adds nothing to the overall output
 */

 //...do stufff....
 double tmpL, tmpR, tmp2L, tmp2R;
 double ampEnvOut;
 double fltEnvOut;

 // for test only:
 tmpL = tmpR = 0.0;

 // generate the source-signal:
 getSourceSampleFrameStereo(&tmpL, &tmpR);

 if( filterIsOn )
 {
  // read the filter-envelope:
  fltEnvOut = fltEnv.getSample();

  // setup the filter:
  filter.setFreq(fltEnvOut * fltFreqScaler * fltFreq);
  filter.calcCoeffs();

  // filter the source-signal:
  filter.getSampleFrameStereo(&tmpL, &tmpR, &tmp2L, &tmp2R);

  // read and apply the amplitude-envelope:
  ampEnvOut = ampEnv.getSample();
  tmp2L *= ampEnvOut;
  tmp2R *= ampEnvOut;

  // write the outputs into the slots (accumulate, because we assume that other
  // voices write into them as well):
  *voiceAmplitude += ampEnvOut;
  *outL           += tmp2L;
  *outR           += tmp2R;
 }
 else
 {
  // read and apply the amplitude-envelope:
  ampEnvOut = ampEnv.getSample();
  tmpL *= ampEnvOut;
  tmpR *= ampEnvOut;

  *voiceAmplitude += ampEnvOut;
  *outL += tmpL;
  *outR += tmpR;
 }

 /*
 // read and apply the amplitude-envelope:
 ampEnvOut = ampEnv.getSample();
 tmp2L *= ampEnvOut;
 tmp2R *= ampEnvOut;

 // write the outputs into the slots (accumulate, because we assume that other
 // voices write into them as well):
 *outL += tmp2L;
 *outR += tmp2R;
 */

 // increment the note age:
 currentNoteAge++;
}

INLINE void MagicCarpetVoice::getSourceSampleFrameStereo(double *outL, 
                                                         double *outR)
{
 double tmpL, tmpR;
 double accuL, accuR;

 // accumulate the outputs of the 4 sources:
 accuL  = 0.0; 
 accuR  = 0.0;
 
 topLeftSource.getSampleFrameStereo(&tmpL, &tmpR);
 accuL += topLeftAmplitude * tmpL;
 accuR += topLeftAmplitude * tmpR;
 topRightSource.getSampleFrameStereo(&tmpL, &tmpR);
 accuL += topRightAmplitude * tmpL;
 accuR += topRightAmplitude * tmpR;
 bottomLeftSource.getSampleFrameStereo(&tmpL, &tmpR);
 accuL += bottomLeftAmplitude * tmpL;
 accuR += bottomLeftAmplitude * tmpR;
 bottomRightSource.getSampleFrameStereo(&tmpL, &tmpR);
 accuL += bottomRightAmplitude * tmpL;
 accuR += bottomRightAmplitude * tmpR;

 *outL = accuL;
 *outR = accuR;

 // apply differencing filter and store the calculates signals in the
 // output-slots:
 //*outL = accuL - x1L;
 //*outR = accuR - x1R;

 // update the state variables of the differencing-filters:
 //x1L = accuL;
 //x1R = accuR; 

}

#endif //  MagicCarpetVoice_h