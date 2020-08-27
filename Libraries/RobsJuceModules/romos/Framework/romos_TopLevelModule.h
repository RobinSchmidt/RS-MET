#ifndef romos_TopLevelModule_h
#define romos_TopLevelModule_h


// maybe we need a class TopLevelModuleObserver to take care of focusedModuleChanged callbacks
// and/or removal of submodules...but maybe this appaly to any kind of modules

/** This class is used to represent the special top-level, or root module. This usually represents 
a complete instrument or effect unit. It is a CompositeModule with a bit of additional 
functionality that applies only to an instrument/effect as a whole.


It also stores some global and per-voice information such as system sample-rate, the tempo in BPM, 
the notes that are currently active (for each voice), etc. Some modules access these informations.  
...nah not true - we now have a singleton for this

\todo: maybe use block-processing here, also implement global oversampling, etc.  */

class TopLevelModule : public ContainerModule
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Disconnects the input pins of our output modules (if connected) - needed for clean-up when 
  recalling a patch. */
  void disconnectAudioOutputModules();

  /** Overriden to avoid changing the polyphony (a toplevel-module is always monophonic). */
  virtual void setPolyphonic(bool /*shouldBePolyphonic*/)
  {
    // do nothing
  }

  /** Overriden to avoid changing the polyphony (a toplevel-module is always monophonic). */
  virtual void setPolyphonicRecursively(bool /*shouldBePolyphonic*/)
  {
    // do nothing 
  }

  /** Overriden to avoid adding of audio inputs to the toplevel module. */
  virtual Module* addAudioInputModule(rosic::rsString name = rosic::rsString(), int x = 1, 
    int y = 1, bool sortModuleArrayAfterInsertion = true)
  {
    return NULL;
  }

  /** Overriden to avoid adding of audio outputs to the toplevel module. */
  virtual Module* addAudioOutputModule(rosic::rsString /*name*/ = rosic::rsString(), int /*x*/ = 1, 
    int /*y*/ = 1,  bool /*sortModuleArrayAfterInsertion*/ = true)
  {
    return NULL;
  }



  /** Overriden to avoid deletion of I/O modules and additionaly delete connections that go 
  directly from inputs to outputs. */
  virtual void deleteChildModule(Module *moduleToDelete, 
    bool updateHasDelayedConnectionFlag = true);

  /** Overriden to avoid attempts to re-connect our outside connections (which - by definition - 
  don't exist). */
  virtual void sortChildModuleArray();

  //-----------------------------------------------------------------------------------------------
  // Inquiry:

  virtual std::string getTypeName() const override { return "TopLevelModule"; }
   // if you change this, be sure to make a corresponding change in 
   // ModulePropertiesEditorHolder::createPropertiesEditorForSelectedModule

  virtual bool isTopLevelModule() const override { return true; }


  //-----------------------------------------------------------------------------------------------
  // Callbacks:


  //-----------------------------------------------------------------------------------------------
  // Processing:

  /** Calculates the output-samples for both channels and stores them at the adresses of *outL and 
  *outR. */
  template <class SampleType>
  INLINE void getSampleFrameStereo(SampleType *inOutL, SampleType *inOutR);

  /** Produces a block of output samples. */
  template <class SampleType>
  INLINE void getBlockOfSampleFramesStereo(SampleType *inOutL, SampleType *inOutR, int numFrames);


protected:


  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  TopLevelModule();

  /** Destructor. */
  virtual ~TopLevelModule();

  friend class ModuleFactory;    // old
  friend class ModuleFactory;


  /** Does initialization stuff for input buffers and input-modules (memory allocation, setting up
  the pointers, etc.). To be called from the constructor. */
  void initInputBuffersAndModules();

  /** Clean up routine for our input buffers and modules. To be called from the destructor. */
  void cleanUpInputBuffersAndModules();

  /** Allcoates the memory for our input buffers. */
  void allocateInputBuffers();

  /** Frees the memory for our input buffers. */
  void deleteInputBuffers();

  /** Clears the inpu buffers (sets all values to zero). */
  void clearInputBuffers();

  /** Sets up our top-level input modules by assigning their pointers to our input buffers. */
  void setupInputModules();

  double *inL;
  double *inR;


  // \todo maybe use a special processing function that uses a flat list/array of all modules that 
  // represents the whole instrument as one monolith - that is, without any structuring via 
  // containers. this array can contain just pointers to the modules that are inside the containers
  // - no actual duplication is necessary. doing it like this avoids the container-overhead. We 
  // should make sure that this array is always properly sorted - but that's straightforward 
  // -> just traverse the tree and include the I/O modules -  when moving a module inside a 
  // container, the container can request its toplevel module, access the flat array there and 
  // replace a subarray with the re-ordered subarray we can assign the source/target module 
  // pointers in the connections for container modules to the container's respective I/O modules
  // -> as these are first class identity modules, it should work for the DSP side, for the GUI 
  // side, we must take care that the connection still returns the container-module when the GUI 
  // requests the pointer (maybe have another function getApparentSourceModule or something)

};

//-------------------------------------------------------------------------------------------------
// inlined functions:
/*
class TopLevelModuleInfo : public ModuleTypeInfo
{
public:
  BiquadTypeInfo() {
    shortName    = "TopLevelModule";
    fullName     = "TopLevelModule";
    description  = "Biquad (2 pole, 2 zero) filter. Realizes y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]";
    category     = "Filters";
    createModule =  []()->Module* { return new Biquad; };
  }
};
*/

template <class SampleType>
INLINE void TopLevelModule::getSampleFrameStereo(SampleType *inOutL, SampleType *inOutR)
{
  // for debug:
  int offset1 = processingStatus.getBufferSize();
  int offset2 = getRequiredOutputBufferSizePerPin();

  *inL = (double)*inOutL;
  *inR = (double)*inOutR;
  double *doublePointerL = audioOutputs;
    //doublePointerR = audioOutputs + processingStatus.getBufferSize();  // later    
  double *doublePointerR = audioOutputs + getRequiredOutputBufferSizePerPin();


  processSampleFrame();

  *inOutL = (SampleType)*(audioOutputs);
  //*inOutR = (SampleType) *(audioOutputs + processingStatus.getBufferSize());
  *inOutR = (SampleType)*(audioOutputs + getRequiredOutputBufferSizePerPin());

  voiceAllocator.resetTriggerFlags();
}

template <class SampleType>
INLINE void TopLevelModule::getBlockOfSampleFramesStereo(SampleType *inOutL, 
  SampleType *inOutR, int numFrames)
{
  if(numFrames < 1)
    return;


  // for debug:
  int offset1 = processingStatus.getBufferSize();
  int offset2 = getRequiredOutputBufferSizePerPin();

  int maxBlockSize = romos::processingStatus.getBufferSize();
  int blockStart   = 0;

  // when this block-processing function is called immediately after receiving a note-event, we 
  // must process the 1st sample-frame separately using the per-sample processing function to make
  // sure that the voiceAllocator's note-on/off flags are reset to false after the 1st sample:
  if(voiceAllocator.isAnyTriggerFlagSet())
  {
    getSampleFrameStereo(inOutL, inOutR);
    blockStart = 1;
  }

  while(blockStart < numFrames)
  {
    int blockSize = numFrames - blockStart;
    if(blockSize > maxBlockSize)
      blockSize = maxBlockSize;

    // establish input block:
    SampleType *samplePointerL = inOutL + blockStart;
    SampleType *samplePointerR = inOutR + blockStart;
    double     *doublePointerL = inL;
    double     *doublePointerR = inR;
    for(int n = 0; n < blockSize; n++)
    {
      *doublePointerL = (double)*samplePointerL;
      *doublePointerR = (double)*samplePointerR;
      samplePointerL++;
      samplePointerR++;
      doublePointerL++;
      doublePointerR++;
    }

    // process:
    processBlockOfSamples(blockSize);

    // retrieve output block:
    samplePointerL = inOutL + blockStart;
    samplePointerR = inOutR + blockStart;
    doublePointerL = audioOutputs;
    //doublePointerR = audioOutputs + processingStatus.getBufferSize();  // later
    doublePointerR = audioOutputs + getRequiredOutputBufferSizePerPin();
    for(int n = 0; n < blockSize; n++)
    {
      *samplePointerL = (SampleType)*doublePointerL;
      *samplePointerR = (SampleType)*doublePointerR;
      samplePointerL++;
      samplePointerR++;
      doublePointerL++;
      doublePointerR++;
    }

    blockStart += blockSize;
    int dummy = 0;
  }
}


#endif 
