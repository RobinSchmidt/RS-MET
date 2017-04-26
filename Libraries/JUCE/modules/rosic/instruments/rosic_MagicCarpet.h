#ifndef rosic_MagicCarpet_h
#define rosic_MagicCarpet_h

// rosic-indcludes:
//#include "../infrastructure/rosic_PresetRememberer.h"
#include "../infrastructure/rosic_PolyphonicInstrument.h"
//#include "../modulators/rosic_LowFrequencyOscillator.h"
#include "rosic_MagicCarpetVoice.h"
#include "../filters/rosic_EllipticSubBandFilter.h"
//#include "../filters/rosic_Equalizer.h"
#include "../filters/rosic_EqualizerStereo.h"
//#include "../others/rosic_VectorMixer.h"
#include "../effects/rosic_DelayPhaser.h"

namespace rosic
{

  /**

  MagicCarpet ...

  */

  class MagicCarpet : public PolyphonicInstrument
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    MagicCarpet();   

    /** Destructor. */
    virtual ~MagicCarpet();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    //---------------------------------------------------------------------------------------------
    // audio processing:

    virtual void setSampleRate(double newSampleRate);
    virtual void setBeatsPerMinute(double newBpm);
    virtual INLINE void getSampleFrameStereo(double *outL, double *outR);   

    //=============================================================================================

    MagicCarpetVoice* voiceArray;
    //VectorMixer            vectorMixer;
    //LowFrequencyOscillator xLfo, yLfo;
    EqualizerStereo   equalizer;
    DelayPhaser       delayPhaser;

    // for later: use the vectorMixer above only for th GUI and two additonal vectorMixers for left 
    // and right channel for the actual DSP code (because they are going to be modulatd by the 
    // LFOs).   ???

  protected:

    // SampleBuffer objects to be used by the samplePlayers:
    SampleBuffer sampleBufferTopLeft, sampleBufferTopRight, sampleBufferBottomLeft, 
      sampleBufferBottomRight;

    // SamplePlaybackParameters objects to be used by the samplePlayers:
    SamplePlaybackParameters samplePlaybackParametersTopLeft, samplePlaybackParametersTopRight, 
      samplePlaybackParametersBottomLeft, samplePlaybackParametersBottomRight; 

    // wavetable objects to be used by the LFOs:
    WaveTable xLfoWaveTable, yLfoWaveTable;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void MagicCarpet::getSampleFrameStereo(double *outL, double *outR)
  {
    if( voicePointers == NULL )
      return;

    // cast the voicePointer array to the subclass:
    MagicCarpetVoice* castedVoicePointer;

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
        castedVoicePointer = static_cast<MagicCarpetVoice*>(voicePointers[i]);

        // let voice i generate its output-sample and accumulate it:
        castedVoicePointer->getSampleFrameStereo(
          &accuL1, &accuR1, &accuL2, &accuR2, &cumulativePower);

        // increment the variable which keeps track of the active voices:
        numActiveVoices++;
      }
    }



    // anti-alias filtering and decimation:
    left  = accuL1;
    right = accuR1;


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


    double compensator = getCumulativePowerCompensator(cumulativePower);
    left  *= masterAmplitude*compensator;
    right *= masterAmplitude*compensator;


    // write the final output samples into the memory-slots:
    *outL = left;
    *outR = right;
  }

} // end namespace rosic

#endif // MagicCarpet_h
