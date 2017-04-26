#ifndef rosic_SimpleSampler_h
#define rosic_SimpleSampler_h

// rosic-indcludes:
//#include "../infrastructure/rosic_PresetRememberer.h"
#include "../infrastructure/rosic_PolyphonicInstrument.h"
#include "rosic_SimpleSamplerVoice.h"
#include "../filters/rosic_EllipticSubBandFilter.h"

namespace rosic
{

  /**

  SimpleSampler ...

  */

  class SimpleSampler : public PolyphonicInstrument
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SimpleSampler();   

    /** Destructor. */
    virtual ~SimpleSampler();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Overrides setSampleRate (inherited from PolyphonicInstrument) in order to update the 
    sample-rate for the global dc-blokcer as well. */
    //virtual void setSampleRate(double newSampleRate);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates the output-samples for both channels and stores them at the adresses of *outL 
    and *outR. */
    virtual INLINE void getSampleFrameStereo(double *outL, double *outR);   

    //=============================================================================================

    SimpleSamplerVoice*   voiceArray;
    EllipticSubBandFilter antiAliasFilterL, antiAliasFilterR;

  protected:

    SampleBuffer sampleBuffer1;
      // These are the SampleBuffer objects to be used by the samplePlayers.

    SamplePlaybackParameters samplePlaybackParameters1;
      // These are the SamplePlaybackParameters objects to be used by the samplePlayers.

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void SimpleSampler::getSampleFrameStereo(double *outL, double *outR)
  {
    if( voicePointers == NULL )
      return;

    // cast the voicePointer array to the subclass:
    SimpleSamplerVoice* castedVoicePointer;

    double accuL1, accuR1, accuL2, accuR2, left, right;
    double cumulativePower;
    int    i;

    // initialize the output slots with zero in order to accumulate over the 
    // voices and then loop through the (active) voices:
    accuL1          = 0.0;
    accuR1          = 0.0;
    accuL2          = 0.0;
    accuR2          = 0.0;
    numActiveVoices = 0;
    cumulativePower = 0.0;

    // loop through the voices, voicePointers[0] servers only as template for the others and is not 
    // invoked to process audio:
    for(i=1; i<=numPlayableVoices; i++)
    {
      if( voicePointers[i]->isSilent != true ) 
      {
        castedVoicePointer = static_cast<SimpleSamplerVoice*>(voicePointers[i]);

        // let voice i generate its output-sample and accumulate it:
        castedVoicePointer->getSampleFrameStereo(
          &accuL1, &accuR1, &accuL2, &accuR2, &cumulativePower);

        // increment the variable which keeps track of the active voices:
        numActiveVoices++;
      }
    }

    /*
    if( voiceAmplitudeSum < 0.0 )
    {
      DEBUG_BREAK; // this one hickups FL even in release builds
    }
    */


    // anti-alias filtering and decimation:
    left  = antiAliasFilterL.getSampleDirect1(accuL1);
    left  = antiAliasFilterL.getSampleDirect1(accuL2);
    right = antiAliasFilterR.getSampleDirect1(accuR1);
    right = antiAliasFilterR.getSampleDirect1(accuR2);


    // DC blocker:
    //.....


    // apply the amplitude-scaling by the number of voices (weigthed by their 
    // respective amplitudes):
    //voiceAmplitudeSum = max(1.0, voiceAmplitudeSum);
    /*
    double voiceScaleFactor = masterAmplitude * (1.0 + 0.01*masterLevelByVoices) / 
      (1.0 + 0.01*masterLevelByVoices * sqrt(voiceAmplitudeSum) );
    */
    /*
    double voiceScaleFactor; 
    if( voiceAmplitudeSum < 1.0 || masterLevelByVoices == 0.0 )
      voiceScaleFactor = masterAmplitude;
    else
    {
      double w         = 0.01*masterLevelByVoices; // weighting factor
      double factor    = 1.0 / (w *  sqrt(voiceAmplitudeSum));
      voiceScaleFactor = w*factor + (1.0-w);


      voiceScaleFactor *= masterAmplitude;

      //voicScaleFactor = masterAmplitude / (0.01*masterLevelByVoices * sqrt(voiceAmplitudeSum));
    }
    */
    /*
    if( voiceAmplitudeSum < 0.0 || voiceScaleFactor > 2.0 )
    {
      DEBUG_BREAK; // this one hickups FL even in release builds
    }
    */

    double compensator = getCumulativePowerCompensator(cumulativePower);
    left  *= masterAmplitude*compensator;
    right *= masterAmplitude*compensator;


    // write the final output samples into the memory-slots:
    *outL = left;
    *outR = right;
  }

} // end namespace rosic

#endif // SimpleSampler_h
