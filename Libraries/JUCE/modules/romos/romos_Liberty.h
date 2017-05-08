#ifndef romos_Liberty_h
#define romos_Liberty_h

#include "framework/romos_TopLevelModule.h"
#include "framework/romos_NoteEvent.h"

namespace romos
{

  /**

  Liberty ....

  */

  class Liberty   //: public PolyphonicInstrument
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Liberty();   

    /** Destructor. */
    ~Liberty();  

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets the samplerate at which the system receives and produces signals. */
    void setSampleRate(double newSampleRate);

    /** Resets the system into a default state. */
    void reset();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns a pointer to the top-level module. */
    INLINE TopLevelModule* getTopLevelModule() const { return topLevelModule; }

    
    INLINE bool isSilent() { return false; } 

    //-------------------------------------------------------------------------------------------------------------------------------------
    // event handling:

    /** Triggers a note-on for the given key. */
    void noteOn(int key, int velocity);

    /** Triggers a note-off for the given key. */
    void noteOff(int key);

    /** Resets all voices into default state. */
    //void resetAllVoices();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates the output-samples for both channels and stores them at the adresses of *outL and *outR. */
    template <class SampleType>
    INLINE void getSampleFrameStereo(SampleType *inOutL, SampleType *inOutR);  


    /** Produces a block of output samples. */
    template <class SampleType>
    INLINE void getBlockOfSampleFramesStereo(SampleType *inOutL, SampleType *inOutR, int numFrames);


    //=====================================================================================================================================

  protected:

    TopLevelModule *topLevelModule;

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  template <class SampleType>
  INLINE void Liberty::getSampleFrameStereo(SampleType *inOutL, SampleType *inOutR)
  {
    topLevelModule->getSampleFrameStereo(inOutL, inOutR);
  }

  template <class SampleType>
  INLINE void Liberty::getBlockOfSampleFramesStereo(SampleType *inOutL, SampleType *inOutR, int numFrames)
  {
    topLevelModule->getBlockOfSampleFramesStereo(inOutL, inOutR, numFrames);
  }

} // end namespace romos

#endif // Liberty_h
