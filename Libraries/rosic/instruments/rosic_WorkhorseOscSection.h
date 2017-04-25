#ifndef rosic_WorkhorseOscSection_h
#define rosic_WorkhorseOscSection_h

// rosic-indcludes:
//#include "../infrastructure/rosic_PresetRememberer.h"
#include "../generators/rosic_SamplePlayer.h"

namespace rosic
{
  /**

  This class implements the oscillator section of the Workhorse synth. Its purpose is mainly to
  keep track of the osc-section preset.

  */

  class WorkhorseOscSection //: public PresetRememberer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WorkhorseOscSection();

    /** Destructor. */
    virtual ~WorkhorseOscSection();

    //---------------------------------------------------------------------------------------------
    // delegations:

    void setPlaybackFrequencyNominal(double newFrequency);
    void setSampleRate(double newSampleRate);
    void reset();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *outL, double *outR);

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    SamplePlayer samplePlayerTopLeft, samplePlayerTopRight, samplePlayerBottomLeft, 
      samplePlayerBottomRight;

    // FourSourceVectorMixer .... or no, put that into the magicCarpetClass (as this does not 
    // exist per voice

  protected:


  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void WorkhorseOscSection::getSampleFrameStereo(double *outL, double *outR)
  {
    samplePlayerTopLeft.getSampleFrameStereo(outL, outR); // preliminary
  }

} // end namespace rosic

#endif //  rosic_WorkhorseOscSection_h