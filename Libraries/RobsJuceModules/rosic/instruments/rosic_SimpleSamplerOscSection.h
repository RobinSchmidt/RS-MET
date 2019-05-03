#ifndef rosic_SimpleSamplerOscSection_h
#define rosic_SimpleSamplerOscSection_h

//// rosic-indcludes:
//#include "../generators/rosic_SamplePlayer.h"

namespace rosic
{
  /**

  This class implements the oscillator section of the SimpleSampler synth. Its purpose is mainly to
  keep track of the osc-section preset.

  */

  class SimpleSamplerOscSection //: public PresetRememberer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SimpleSamplerOscSection();

    /** Destructor. */
    virtual ~SimpleSamplerOscSection();

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    SamplePlayer samplePlayer1;

  protected:


  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  /*
  INLINE void SimpleSamplerOscSection::getSampleFrameStereo(double *outL, double *outR)
  {
    *outL = 0.0;
    *outR = 0.0;
  }
  */

} // end namespace rosic

#endif //  rosic_SimpleSamplerOscSection_h