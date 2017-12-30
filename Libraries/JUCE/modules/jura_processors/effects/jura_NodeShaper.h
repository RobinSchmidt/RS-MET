#ifndef jura_NodeShaper_h
#define jura_NodeShaper_h

/** This class  */

class NodeShaperAudioModule : public ModulatableAudioModule
{

public:

  //---------------------------------------------------------------------------------------------
  // \name Construction/Destruction

  NodeShaperAudioModule(CriticalSection *newPlugInLock);

  virtual ~NodeShaperAudioModule();

  AudioModuleEditor* createEditor() override;

  //---------------------------------------------------------------------------------------------
  // \name Setup

  virtual void setSampleRate(double newSampleRate) override
  {
    //wrappedNodeShaper->setSampleRate(newSampleRate);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int i = 0; i < numChannels; i++)
      for(int n = 0; n < numSamples; n++)
        inOutBuffer[i][n] = mapper.getValue(inOutBuffer[i][n]);
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    *left  = mapper.getValue(*left);
    *right = mapper.getValue(*right);
  }

protected:

  void createParameters();

  RAPT::rsInterpolatingFunction<double> mapper;

  friend class NodeShaperModuleEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NodeShaperAudioModule)
};



//=================================================================================================

class NodeShaperModuleEditor : public AudioModuleEditor
{

public:

  //---------------------------------------------------------------------------------------------
  // \name Construction/Destruction

  /** Constructor. */
  NodeShaperModuleEditor(NodeShaperAudioModule* newNodeShaperAudioModule);

  /** Destructor. */
  virtual ~NodeShaperModuleEditor();

  //---------------------------------------------------------------------------------------------
  // callbacks:
  
  virtual void resized() override;

protected:

  virtual void createWidgets();

  virtual void updateWidgetsAccordingToState() override;

  NodeShaperAudioModule *nodeShaperModule;

  rsNodeBasedFunctionEditor *nodeEditor;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NodeShaperModuleEditor)
};

#endif 
