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
    wrappedNodeShaper->setSampleRate(newSampleRate);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {

  }

  virtual void processStereoFrame(double *left, double *right) override
  {

  }

protected:

  void createParameters();

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NodeShaperAudioModule)
};

//=================================================================================================

class NodeShapeModuleEditor : public AudioModuleEditor
{

public:

  //---------------------------------------------------------------------------------------------
  // \name Construction/Destruction

  /** Constructor. */
  NodeShapeModuleEditor(NodeShaperAudioModule* newNodeShaperAudioModule);

  /** Destructor. */
  virtual ~NodeShapeModuleEditor();

  //---------------------------------------------------------------------------------------------
  // callbacks:
  
  virtual void resized() override;

protected:

  virtual void createWidgets();

  virtual void updateWidgetsAccordingToState() override;

  NodeShaperAudioModule *nodeShaperModule;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NodeShapeModuleEditor)
};

#endif 
