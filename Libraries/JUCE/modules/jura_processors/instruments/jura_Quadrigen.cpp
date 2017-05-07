
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

//=================================================================================================

QuadrigenModuleEditor::QuadrigenModuleEditor(CriticalSection *newPlugInLock, 
  QuadrigenAudioModule* newQuadrigenAudioModule) 
  : AudioModuleEditor(newQuadrigenAudioModule)
{
  setHeadlineStyle(MAIN_HEADLINE);
  jassert(newQuadrigenAudioModule != NULL ); // you must pass a valid module here
  quadrigenModuleToEdit = newQuadrigenAudioModule;

  acquireLock();

  // remember for toggling:
  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++ )
  {
    if( quadrigenModuleToEdit == NULL && quadrigenModuleToEdit->wrappedQuadrigen == NULL )
      oldAlgorithmIndices[i] = quadrigenModuleToEdit->wrappedQuadrigen->getGeneratorAlgorithmIndex(i);
    else
      oldAlgorithmIndices[i] = rosic::Quadrigen::MUTE;
  }

  matrixEditor = new RoutingMatrixModuleEditor(lock, quadrigenModuleToEdit->matrixModule);
  matrixEditor->setHeadlineStyle(AudioModuleEditor::NO_HEADLINE);
  addChildEditor(matrixEditor);

  // create the sub-editors:
  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
    moduleEditors[i] = NULL;

  // attach to the underlying audiomodule to be edited:
  quadrigenModuleToEdit->setEditor(this);

  initializeColourScheme();
  updateWidgetsAccordingToState();

  releaseLock();
}

QuadrigenModuleEditor::~QuadrigenModuleEditor()
{
  // detach from the underlying audiomodule to be edited:
  acquireLock();
  quadrigenModuleToEdit->setEditor(NULL);
  releaseLock();
}

//-------------------------------------------------------------------------------------------------
// setup:

void QuadrigenModuleEditor::initializeColourScheme()
{

}

//-------------------------------------------------------------------------------------------------
// callbacks:

void QuadrigenModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  acquireLock();
  if( quadrigenModuleToEdit == NULL )
  {
    releaseLock();
    return;
  }
  if( quadrigenModuleToEdit->wrappedQuadrigen == NULL )
  {
    releaseLock();
    return;
  }

  //...

  quadrigenModuleToEdit->markStateAsDirty();
  releaseLock();
}

void QuadrigenModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  acquireLock();
  if( quadrigenModuleToEdit == NULL )
  {
    releaseLock();
    return;
  }
  if( quadrigenModuleToEdit->wrappedQuadrigen == NULL )
  {
    releaseLock();
    return;
  }

  rosic::Quadrigen* core = quadrigenModuleToEdit->wrappedQuadrigen;

  //....

  if( quadrigenModuleToEdit != NULL )
    quadrigenModuleToEdit->markStateAsDirty();

  releaseLock();
}

void QuadrigenModuleEditor::rSliderValueChanged(RSlider *rSliderThatHasChanged)
{
  //...

  if( quadrigenModuleToEdit != NULL )
    quadrigenModuleToEdit->markStateAsDirty();
}


void QuadrigenModuleEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void QuadrigenModuleEditor::mouseDown(const MouseEvent &e)
{

}

void QuadrigenModuleEditor::updateWidgetsAccordingToState()
{
  acquireLock();
  if( quadrigenModuleToEdit == NULL )
  {
    releaseLock();
    return;
  }
  if( quadrigenModuleToEdit->wrappedQuadrigen == NULL )
  {
    releaseLock();
    return;
  }

  AudioModuleEditor::updateWidgetsAccordingToState();

  // delete old and create new editors:
  removeChildEditors();
  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
    createEditorForSlot(i, quadrigenModuleToEdit->wrappedQuadrigen->getGeneratorAlgorithmIndex(i));

  releaseLock();
}

void QuadrigenModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);

  fillRectWithBilinearGradient(g, globalRectangle, editorColourScheme.topLeft, 
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);

  g.setColour(editorColourScheme.outline);
  g.drawRect(globalRectangle, 2);

  int x = globalRectangle.getX(); // + middleRectangle.getWidth()/2;
  int y = globalRectangle.getY();
  //drawBitmapFontText(g, x+4, y+4, juce::String(T("Global Settings")), BigFont::getInstance(), 
  //  editorColourScheme.labelTextColour);

  acquireLock();
  if( quadrigenModuleToEdit == NULL )
  {
    releaseLock();
    return;
  }
  if( quadrigenModuleToEdit->wrappedQuadrigen == NULL )
  {
    releaseLock();
    return;
  }

  rosic::Quadrigen* core = quadrigenModuleToEdit->wrappedQuadrigen;
  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
  {
    fillRectWithBilinearGradient(g, slotRectangles[i], editorColourScheme.topLeft, editorColourScheme.topRight, 
      editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
    g.drawRect(slotRectangles[i], 2);

    x = slotRectangles[i].getX(); 
    y = slotRectangles[i].getY();
    juce::String headlineString = juce::String(i+1) + juce::String(" - ") + 
      quadrigenModuleToEdit->generatorAlgorithmIndexToString(core->getGeneratorAlgorithmIndex(i));

    //drawBitmapFontText(g, x+4, y+4, headlineString, &boldFont16px, 
    //  editorColourScheme.headline);
    // old version

    jassertfalse; // something doesn't work with accessing the 16px font instance here
                  // -> check out why
    drawBitmapFontText(g, x+4, y+4, headlineString, &BitmapFontRoundedBoldA10D0::instance, 
      editorColourScheme.headline); 
    // this works (with 10px font) - but we want a big font here ...figure that out...

    //drawBitmapFontText(g, x+4, y+4, headlineString, &BitmapFontRoundedBoldA16D0::instance, 
    //  editorColourScheme.headline);
    // why does this not work? (linker error)

    //drawBitmapFontText(g, x+4, y+4, headlineString, &boldFont16px, editorColourScheme.headline);
    // undeclare identifier

  }

  releaseLock();
}

void QuadrigenModuleEditor::resized()
{
  acquireLock();

  AudioModuleEditor::resized();
  int x = 0;
  int y = getHeadlineBottom();
  int w = getWidth();
  int h = getHeight();

  // setup the sizes of the child editors:
  x = globalRectangle.getRight();
  y = globalRectangle.getY();
  w = (getWidth()-x-4)/2;
  h = globalRectangle.getHeight()/2;
  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
  {
    int x2 = x;
    if( i == 1 || i == 3 )
      x2 += w;

    int y2 = y;
    if( i >=2 )
      y2 += h;

    slotRectangles[i].setBounds(x2, y2,    w, h);
    if( moduleEditors[i] != NULL )
      moduleEditors[i]->setBounds(x2, y2+24, w, h-24);
  }

  releaseLock();
}

//-------------------------------------------------------------------------------------------------
// intenal helper functions:

int QuadrigenModuleEditor::getSlotIndexAtPixelPosition(int x, int y)
{
  for(int i = 0; i < rosic::Quadrigen::numGeneratorSlots; i++ )
  {
    if( slotRectangles[i].contains(x, y) == true )
      return i;
  }
  return -1;
}

void QuadrigenModuleEditor::setGeneratorAlgorithm(int slotIndex, int newAlgorithmIndex)
{
  if( slotIndex < 0 || slotIndex >= rosic::Quadrigen::numGeneratorSlots )
  {
    jassertfalse;
    return;
  }

  acquireLock();
  if( quadrigenModuleToEdit == NULL )
  {
    releaseLock();
    return;
  }

  quadrigenModuleToEdit->setGeneratorAlgorithm(slotIndex, newAlgorithmIndex);
  // the quadrigenModuleToEdit has a pointer to 'this' editor and will take care of updating
  // the GUI (deleting and creating the appropriate child-editor)

  releaseLock();
}

/*
void QuadrigenModuleEditor::createEditorsIfNeeded()
{
acquireLock();
for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
createEditorForSlotIfNeeded(i);
releaseLock();
}

void QuadrigenModuleEditor::createEditorForSlotIfNeeded(int slotIndex)
{

}
*/

void QuadrigenModuleEditor::removeChildEditors()
{
  acquireLock();
  for(int i=0; i<rosic::Quadrigen::numGeneratorSlots; i++)
    removeChildEditorInSlot(i);
  releaseLock();
}

void QuadrigenModuleEditor::removeChildEditorInSlot(int slotIndex)
{
  acquireLock();
  int i = slotIndex;
  if( slotIndex >= 0 && slotIndex < rosic::Quadrigen::numGeneratorSlots 
    && moduleEditors[i] != NULL  )
  {
    moduleEditors[i]->invalidateModulePointer();  // so it doesn't dereference it in the destructor
    removeChildEditor(moduleEditors[i], true);    // deletes the object also
    moduleEditors[i] = NULL;
  }
  releaseLock();
}

void QuadrigenModuleEditor::createEditorForSlot(int slotIndex, int algorithmIndex)
{
  acquireLock();
  if( quadrigenModuleToEdit == NULL )
  {
    releaseLock();
    return;
  }
  if( quadrigenModuleToEdit->wrappedQuadrigen == NULL )
  {
    releaseLock();
    return;
  }

  jassert( moduleEditors[slotIndex] == NULL ); 
  // you should delete the old editor and null the pointer before creating a new one

  switch( algorithmIndex )
  {
  case rosic::Quadrigen::OSCILLATOR_STEREO:
  {
    OscillatorStereoAudioModule *audioModule = static_cast<OscillatorStereoAudioModule*>
      (quadrigenModuleToEdit->getGeneratorAudioModule(slotIndex));
    moduleEditors[slotIndex] = new OscillatorStereoEditor(lock, audioModule); // mmm..should it be OscillatorStereoModuleEditor?
                                                                              // uncomment again later
  } break;


  // ------> INSERT NEW CASE-MARK HERE WHEN ADDING A NEW ALGORITHM <------


  case rosic::Quadrigen::MUTE:
  {
    MuteAudioModule *audioModule = static_cast<MuteAudioModule*>
      (quadrigenModuleToEdit->getGeneratorAudioModule(slotIndex));
    moduleEditors[slotIndex] = new MuteModuleEditor(lock, audioModule); 
  } break;
  default:
  {
    BypassAudioModule *audioModule = static_cast<BypassAudioModule*>
      (quadrigenModuleToEdit->getGeneratorAudioModule(slotIndex));
    moduleEditors[slotIndex] = new BypassModuleEditor(lock, audioModule); 
  }
  }

  moduleEditors[slotIndex]->setHeadlineStyle(Editor::NO_HEADLINE);    
  moduleEditors[slotIndex]->setLinkPosition(AudioModuleEditor::INVISIBLE);
  moduleEditors[slotIndex]->setDescriptionField(infoField, true);
  addChildEditor(moduleEditors[slotIndex]);
  resized();  // to adjust the bounds of the new editor
  repaint();  // to redraw to headline
  moduleEditors[slotIndex]->updateWidgetsAccordingToState();


  releaseLock();
}