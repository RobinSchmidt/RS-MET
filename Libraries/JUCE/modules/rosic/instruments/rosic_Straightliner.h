#ifndef rosic_Straightliner_h
#define rosic_Straightliner_h

//// rosic-indcludes:
//#include "../infrastructure/rosic_PolyphonicInstrument.h"
//#include "rosic_StraightlinerVoice.h"
//#include "../filters/rosic_EllipticSubBandFilterDirectForm.h"

namespace rosic
{

  /**

  This class respresents the Straightliner synthesizer's core DSP engine

  */

  class Straightliner : public PolyphonicInstrument
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Straightliner();   

    /** Destructor. */
    ~Straightliner();  

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    /** Overrides setSampleRate (inherited from PolyphonicInstrument) in order to update the sample-rate for the global dc-blokcer as 
    well. */
    virtual void setSampleRate(double newSampleRate);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates the output-samples for both channels and stores them at the adresses of *outL and *outR. */
    virtual INLINE void getSampleFrameStereo(double *outL, double *outR);   

    //=====================================================================================================================================

    StraightlinerVoice*             voiceArray;
    rsSubBandFilterMono antiAliasFilterL, antiAliasFilterR;

  protected:
      
    // wavetable objects to be used by the oscillators:
    MipMappedWaveTableStereo waveTableForOsc1, waveTableForOsc2, waveTableForOsc3, 
      waveTableForOsc4;

    OnePoleFilterStereo dcBlocker;   // dc-blocker for the overall output signal

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void Straightliner::getSampleFrameStereo(double *outL, double *outR)
  {
    if( voicePointers == NULL )
      return;

    // cast the voicePointer array to the subclass:
    StraightlinerVoice* castedVoicePointer;

    double accuL1, accuR1, accuL2, accuR2, left, right;
    double cumulativePower;
    int    i;

    // initialize the output slots with zero in order to accumulate over the voices and then loop through the (active) voices:
    accuL1          = 0.0;
    accuR1          = 0.0;
    accuL2          = 0.0;
    accuR2          = 0.0;
    numActiveVoices = 0;
    cumulativePower = 0.0;

    for(i=1; i<=numPlayableVoices; i++)  // voicePointers[0] serves only as template
    {
      if( voicePointers[i]->isSilent != true ) 
      {
        castedVoicePointer = static_cast<StraightlinerVoice*>(voicePointers[i]);
        castedVoicePointer->getSampleFrameStereo(&accuL1, &accuR1, &accuL2, &accuR2, &cumulativePower);
        numActiveVoices++;
      }
    }

    // anti-alias filtering and decimation:
    left  = antiAliasFilterL.getSample(accuL1);
    left  = antiAliasFilterL.getSample(accuL2);
    right = antiAliasFilterR.getSample(accuR1);
    right = antiAliasFilterR.getSample(accuR2);

    // apply DC blocker:
    dcBlocker.getSampleFrameStereo(&left, &right, &left, &right);
    
    // apply cumulative power compensation and mid/side adjustment:
    applyVolumeScalers(&left, &right, cumulativePower);

    // write the final output samples into the memory-slots:
    *outL = left;
    *outR = right; 
  }

} // end namespace rosic

#endif // Straightliner_h
