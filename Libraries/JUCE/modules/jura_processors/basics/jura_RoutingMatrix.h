#ifndef jura_RoutingMatrix_h
#define jura_RoutingMatrix_h

class RoutingMatrixAudioModule : public AudioModule
{

  friend class RoutingMatrixModuleEditor;

public:

  RoutingMatrixAudioModule(CriticalSection *newPlugInLock, 
    rosic::RoutingMatrix *newRoutingMatrixToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::RoutingMatrix *wrappedRoutingMatrix;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

#endif 
