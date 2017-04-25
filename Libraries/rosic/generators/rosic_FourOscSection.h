#ifndef rosic_FourOscSection_h
#define rosic_FourOscSection_h

// rosic-indcludes:
//#include "../infrastructure/rosic_PresetRememberer.h"
#include "rosic_OscillatorStereo.h"

namespace rosic
{
  /**

  This class implements the section of four oscillators (of class OscillatorStereo). Its purpose is 
  mainly to keep track of the osc-section preset.

  */

  class FourOscSection //: public PresetRememberer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    FourOscSection();

    /** Destructor. */
    ~FourOscSection();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate for this voice. */
    virtual void setSampleRate(double newSampleRate); 

    //---------------------------------------------------------------------------------------------
    // audio processing:

    //virtual void getSampleFrameStereo(double *outL, double *outR);  
    /**< Calculates a stereo output sample pair at a time. */

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    OscillatorStereo osc1, osc2, osc3, osc4;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  /*
  INLINE void FourOscSection::getSampleFrameStereo(double *outL, double *outR)
  {
    *outL = 0.0;
    *outR = 0.0;
  }
  */

} // end namespace rosic

#endif //  rosic_FourOscSection_h