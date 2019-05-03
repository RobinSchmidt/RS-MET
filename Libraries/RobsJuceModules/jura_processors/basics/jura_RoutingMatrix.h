#ifndef jura_RoutingMatrix_h
#define jura_RoutingMatrix_h

class RoutingMatrixAudioModule : public AudioModule
{

  friend class RoutingMatrixModuleEditor;

public:

  RoutingMatrixAudioModule(CriticalSection *newPlugInLock, 
    rosic::RoutingMatrix *newRoutingMatrixToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override {}

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::RoutingMatrix *wrappedRoutingMatrix;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class RoutingMatrixModuleEditor : public AudioModuleEditor
{

public:

  RoutingMatrixModuleEditor(CriticalSection *newPlugInLock, 
    RoutingMatrixAudioModule* newRoutingMatrixAudioModule);

  virtual void resized();

protected:

  int numInputs, numOutputs;
  juce::Array<RDraggableNumber*> matrixFields;
  juce::Array<RTextField*>       rowLabels;
  juce::Array<RTextField*>       columnLabels;

  juce_UseDebuggingNewOperator;
};


#endif 
