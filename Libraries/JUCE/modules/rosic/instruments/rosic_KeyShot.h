#ifndef rosic_KeyShot_h
#define rosic_KeyShot_h

// rosic-indcludes:
//#include "../infrastructure/rosic_PresetRememberer.h"
#include "../infrastructure/rosic_PolyphonicInstrument.h"
#include "rosic_KeyShotVoice.h"
#include "../filters/rosic_EllipticSubBandFilter.h"

namespace rosic
{

  /**

  KeyShot ...

  */

  class KeyShot : public PolyphonicInstrument
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    KeyShot();   

    /** Destructor. */
    virtual ~KeyShot();  

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

    KeyShotVoice*   voiceArray;
    //EllipticSubBandFilter antiAliasFilterL, antiAliasFilterR;

  protected:

    SampleBuffer sampleBuffer;
      // These are the SampleBuffer objects to be used by the samplePlayers.

    SamplePlaybackParameters samplePlaybackParameters;
      // These are the SamplePlaybackParameters objects to be used by the samplePlayers.

    // todo: use arrays of 128 entries here

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void KeyShot::getSampleFrameStereo(double *outL, double *outR)
  {
    if( voicePointers == NULL )
      return;

    // cast the voicePointer array to the subclass:
    KeyShotVoice* castedVoicePointer;

    double accuL, accuR, left, right;
    double cumulativePower;
    int    i;

    // initialize the output slots with zero in order to accumulate over the 
    // voices and then loop through the (active) voices:
    accuL           = 0.0;
    accuR           = 0.0;
    numActiveVoices = 0;
    cumulativePower = 0.0;   // only a dummy

    // loop through the voices, voicePointers[0] servers only as template for the others and is not 
    // invoked to process audio:
    for(i=1; i<=numPlayableVoices; i++)
    {
      if( voicePointers[i]->isSilent != true ) 
      {
        castedVoicePointer = static_cast<KeyShotVoice*>(voicePointers[i]);

        // let voice i generate its output-sample and accumulate it:
        castedVoicePointer->getSampleFrameStereo(&accuL, &accuR, &cumulativePower);

        // increment the variable which keeps track of the active voices:
        numActiveVoices++;
      }
    }

    left  = accuL;
    right = accuR;

    double compensator = getCumulativePowerCompensator(cumulativePower);
    left  *= masterAmplitude*compensator;
    right *= masterAmplitude*compensator;

    // write the final output samples into the memory-slots:
    *outL = left;
    *outR = right;
  }

} // end namespace rosic

#endif // KeyShot_h
