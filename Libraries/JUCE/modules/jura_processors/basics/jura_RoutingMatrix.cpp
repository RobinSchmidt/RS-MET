
//-------------------------------------------------------------------------------------------------
// construction/destruction:

RoutingMatrixAudioModule::RoutingMatrixAudioModule(CriticalSection *newPlugInLock, 
  rosic::RoutingMatrix *newRoutingMatrixToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newRoutingMatrixToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedRoutingMatrix = newRoutingMatrixToWrap;
  moduleName = juce::String("RoutingMatrix");

  // create and initialize the automatable parameters:
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void RoutingMatrixAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedRoutingMatrix == NULL )
    return;

  double value  = parameterThatHasChanged->getValue();
  int    index  = getIndexOfParameter(parameterThatHasChanged);
  int    input  = index / wrappedRoutingMatrix->getNumOutputs();
  int    output = index % wrappedRoutingMatrix->getNumOutputs();

  wrappedRoutingMatrix->setMatrixEntry(input, output, value);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void RoutingMatrixAudioModule::initializeAutomatableParameters()
{
  if( wrappedRoutingMatrix == NULL )
    return;

  AutomatableParameter* p;
  for(int i=0; i<wrappedRoutingMatrix->getNumInputs(); i++)
  {
    for(int o=0; o<wrappedRoutingMatrix->getNumOutputs(); o++)
    {
      juce::String name = juce::String("M_") + juce::String(i+1) + 
        juce::String("_") + juce::String(o+1);
      p = new AutomatableParameter(lock, name, -1.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
      addObservedParameter(p);
    }
  }

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

