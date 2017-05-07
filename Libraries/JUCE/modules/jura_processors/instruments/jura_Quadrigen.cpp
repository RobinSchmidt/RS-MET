
//-------------------------------------------------------------------------------------------------
// construction/destruction:

QuadrigenAudioModule::QuadrigenAudioModule(CriticalSection *newPlugInLock, rosic::Quadrigen *quadrigenToWrap)
: AudioModule(newPlugInLock)
{
  jassert(quadrigenToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedQuadrigen = quadrigenToWrap;
  editor           = NULL;
  moduleName       = juce::String("Quadrigen");
  setActiveDirectory(getApplicationDirectory() + juce::String("/QuadrigenPresets"));

  matrixModule = new RoutingMatrixAudioModule(lock, &wrappedQuadrigen->mixMatrix);
  matrixModule->setModuleName(juce::String("RoutingMatrix"));
  addChildAudioModule(matrixModule);

  acquireLock();

  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
  {
    // create a bypass-module for each slot:
    rosic::BypassModule* bypassCoreModule = 
      static_cast<rosic::BypassModule*> (wrappedQuadrigen->getGeneratorModule(i)); // this dynamic_cast causes bugs in the release version
    jura::BypassAudioModule *audioModule = new jura::BypassAudioModule(lock, bypassCoreModule);
    generatorModules[i] = audioModule;
    addChildAudioModule(generatorModules[i]);

    // allocate memory to store the states internally:
    oscillatorStereoStates[i]        = new XmlElement(juce::String("OscillatorStereo"));
  }

  initializeAutomatableParameters();

  releaseLock();
}

QuadrigenAudioModule::~QuadrigenAudioModule()
{
  //acquireLock();
  mutex.enter();
  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
  {
    delete oscillatorStereoStates[i]; 
  }
  mutex.exit();
  //releaseLock();
}

//-------------------------------------------------------------------------------------------------
// setup:

void QuadrigenAudioModule::setEditor(QuadrigenModuleEditor *newEditor)
{
  acquireLock();
  editor = newEditor;
  releaseLock();
}

void QuadrigenAudioModule::setGeneratorAlgorithm(int slotIndex, int newAlgorithmIndex)
{
  acquireLock();
  if( wrappedQuadrigen == NULL )
  {
    releaseLock();
    return;
  }
  if( slotIndex < 0 || slotIndex >= rosic::Quadrigen::numGeneratorSlots )
  {
    releaseLock();
    return;
  }

  // store the state of the old generator to be replaced:
  int oldAlgorithmIndex = wrappedQuadrigen->getGeneratorAlgorithmIndex(slotIndex);
  switch( oldAlgorithmIndex )
  {
  case rosic::Quadrigen::OSCILLATOR_STEREO: 
    {
      jura::OscillatorStereoAudioModule *audioModule = 
        static_cast<jura::OscillatorStereoAudioModule*> (generatorModules[slotIndex]);
      delete oscillatorStereoStates[slotIndex];
      oscillatorStereoStates[slotIndex] = audioModule->getStateAsXml(juce::String::empty, false);
    } break;

    //....................tbc...............

  }

  // delete the old generator (and its (sub)editor, if present):
  if( editor != NULL )
    editor->removeChildEditorInSlot(slotIndex);
  generatorModules[slotIndex]->removeAllStateWatchers();  
  removeChildAudioModule(generatorModules[slotIndex], true);
  wrappedQuadrigen->setGeneratorAlgorithm(slotIndex, newAlgorithmIndex);

  // create the new generator and restore its state from any previous use:
  switch( newAlgorithmIndex )
  {
  case rosic::Quadrigen::OSCILLATOR_STEREO: 
    {
      rosic::OscillatorStereoModule *core = 
        static_cast<rosic::OscillatorStereoModule*> (wrappedQuadrigen->getGeneratorModule(slotIndex));
      jura::OscillatorStereoAudioModule *audioModule = new jura::OscillatorStereoAudioModule(lock, core);
      audioModule->setModuleName(juce::String("OscillatorStereo") + juce::String(slotIndex+1));
      audioModule->setStateFromXml(*oscillatorStereoStates[slotIndex], juce::String::empty, true);
      generatorModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;

  default: // bypass by default (i.e. value out of range)
    {
      rosic::BypassModule *core = 
        static_cast<rosic::BypassModule*> (wrappedQuadrigen->getGeneratorModule(slotIndex));
      jura::BypassAudioModule *audioModule = new jura::BypassAudioModule(lock, core);
      generatorModules[slotIndex] = audioModule;
      addChildAudioModule(audioModule);
    } break;
  }

  // let the editor (if present) create an appropriate child-editor:
  if( editor != NULL )
    editor->createEditorForSlot(slotIndex, newAlgorithmIndex);

  releaseLock();
}

//-------------------------------------------------------------------------------------------------
// automation:

void QuadrigenAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  acquireLock();
  if( wrappedQuadrigen == NULL )
  {
    releaseLock();
    return;
  }

  /*
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedQuadrigen->setDryWet(  value); break;
  case   1: wrappedQuadrigen->setWetLevel(value); break;
  case   2: triggerInterval = value;
  } // end of switch( parameterIndex )
  */

  releaseLock();
}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* QuadrigenAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  acquireLock();
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  if( wrappedQuadrigen != NULL )
  {
    // store the slot-generator assignments:
    for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
    {
      xmlState->setAttribute(juce::String("Slot")+juce::String(i+1), 
        generatorAlgorithmIndexToString(wrappedQuadrigen->getGeneratorAlgorithmIndex(i)) );
    }
  }
  releaseLock();
  return xmlState;
}

void QuadrigenAudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,    
                                           bool markAsClean)
{
  acquireLock();
  if( wrappedQuadrigen != NULL )
  {
    // recall the slot-generator assignments:
    for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
    {
      setGeneratorAlgorithm(i, stringToGeneratorAlgorithmIndex( 
        xmlState.getStringAttribute( juce::String("Slot")+juce::String(i+1), "Mute")));
    }
  }
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  releaseLock();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void QuadrigenAudioModule::initializeAutomatableParameters()
{
  acquireLock();

  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  juce::Array<double> defaultValues;

  // this pointer will be used to temporarily store the addresses of the created Parameter-objects:
  Parameter* p;

  // #00:
  p = new Parameter(lock, "DryWet", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR); 
  addObservedParameter(p);

  // #01:
  p = new Parameter(lock, "WetLevel", -36.0, 6.0, 0.01, 0.0, Parameter::LINEAR); 
  addObservedParameter(p);

  // #02:
  p = new Parameter(lock, "TriggerInterval", 0.0, 64.0, 1.0, 8.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);

  releaseLock();
}

juce::String QuadrigenAudioModule::generatorAlgorithmIndexToString(int index)
{
  switch( index )
  {
  case rosic::Quadrigen::MUTE:                 return juce::String("Mute");
  case rosic::Quadrigen::OSCILLATOR_STEREO:    return juce::String("OscillatorStereo");

  default:                                     return juce::String("Mute");
  }
}

int QuadrigenAudioModule::stringToGeneratorAlgorithmIndex(const juce::String &algoString)
{
  if( algoString == juce::String("Mute")   )            return rosic::Quadrigen::MUTE;
  if( algoString == juce::String("OscillatorStereo") )  return rosic::Quadrigen::OSCILLATOR_STEREO;

  return rosic::Quadrigen::MUTE;
}