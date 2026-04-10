#ifndef romos_Liberty_h
#define romos_Liberty_h


/** Liberty is a modular synthesizer */

class Liberty   //: public PolyphonicInstrument
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:

  /** Constructor. */
  Liberty();

  /** Destructor. */
  ~Liberty();

  //-----------------------------------------------------------------------------------------------
  // Setup:

  /** Sets the samplerate at which the system receives and produces signals. */
  void setSampleRate(double newSampleRate);

  /** Not yet implemented. */
  void setOversampling(int newFactor);

  /** Resets the system into a default state. */
  void reset();
  // Should got into Event Handling section

  //-----------------------------------------------------------------------------------------------
  // Inquiry:

  /** Returns a pointer to the top-level module. */
  INLINE TopLevelModule* getTopLevelModule() const { return topLevelModule; }


  INLINE bool isSilent() { return false; }
  // Not yet implemented

  //-----------------------------------------------------------------------------------------------
  // Event handling:

  /** Triggers a note-on for the given key. */
  void noteOn(int key, int velocity);

  /** Triggers a note-off for the given key. */
  void noteOff(int key);

  /** Resets all voices into default state. */
  //void resetAllVoices();

  //-----------------------------------------------------------------------------------------------
  // Audio processing:

  /** Calculates the output-samples for both channels and stores them at the addresses of *outL and 
  *outR. */
  template <class SampleType>
  INLINE void getSampleFrameStereo(SampleType *inOutL, SampleType *inOutR)
  {
    // todo: upsample/downsample if oversampling != 1
    topLevelModule->getSampleFrameStereo(inOutL, inOutR);
  }

  /** Produces a block of output samples. */
  template <class SampleType>
  INLINE void getBlockOfSampleFramesStereo(SampleType *inOutL, SampleType *inOutR, int numFrames)
  {
    // todo: upsample/downsample if oversampling != 1
    topLevelModule->getBlockOfSampleFramesStereo(inOutL, inOutR, numFrames);
  }

protected:

  TopLevelModule* topLevelModule;

  // Stuff for oversampling:
  //using SubBandFilter = RAPT::rsEllipticSubBandFilter<double, double>;
  //SubBandFilter upSamplerL, upSamplerR, downSamplerL, downSamplerR;
  int overSampling = 1;
  // Not yet used. For later. It would be nice to be able to select oversampling on a per-container
  // basis in ToolChain itself. Maybe later. It's complicated.
  // Maybe later use rsOverSampler. And/or maybe try to use SIMD for the up- and downsampling 
  // filters by using rsEllipticSubBandFilter<rsFloat64x2, double>. Maybe later we can also use 
  // other kind of subband filters

  // Under construction:
  //void populateModuleTypeRegistry();  // not yet used - rename to addCustomModulesToFactory
  //ModuleFactory moduleTypeRegistry;   // not yet used 
  // ToDo: Document what this was supposed to be good for. I don't remember. Ah! I think, the idea
  // was that Liberty can add some more modules to the library of available built-in modules on 
  // top of what's already there in romos::TopLevelModule or something like that. Maybe even 
  // pre-built non-atomic modules, i.e. pre-defined containers. Verify and document this!

};

//-------------------------------------------------------------------------------------------------
// inlined functions:
/*
template <class SampleType>
INLINE void Liberty::getSampleFrameStereo(SampleType *inOutL, SampleType *inOutR)
{
  topLevelModule->getSampleFrameStereo(inOutL, inOutR);
}

template <class SampleType>
INLINE void Liberty::getBlockOfSampleFramesStereo(SampleType *inOutL, SampleType *inOutR, 
  int numFrames)
{
  topLevelModule->getBlockOfSampleFramesStereo(inOutL, inOutR, numFrames);
}
*/


#endif // Liberty_h
