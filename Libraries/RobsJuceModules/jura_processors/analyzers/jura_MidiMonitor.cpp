
//-------------------------------------------------------------------------------------------------
// construction/destruction:

MidiMonitorAudioModule::MidiMonitorAudioModule(CriticalSection *newPlugInLock) 
  : AudioModuleWithMidiIn(newPlugInLock)
{
  setModuleTypeName("MidiMonitor");
  initializeAutomatableParameters();
}

AudioModuleEditor* MidiMonitorAudioModule::createEditor(int type)
{
  return new MidiMonitorModuleEditor(lock, this);
}


//-------------------------------------------------------------------------------------------------
// event handling:

void MidiMonitorAudioModule::handleMidiMessage(juce::MidiMessage message)
{
  // add a string representing the the midi-message when the message is one of the obvserved 
  // types:
  messageStringLock.enter();
  if( midiFilter.isMessageSuitable(message) )
  {
    midiMessageString += midiMessageToString(message, true);

    // don't let the string grow indefinitely:
    int maxLength = 5000;
    int length    = midiMessageString.length();
    if( length > maxLength )
      midiMessageString = midiMessageString.substring(length-maxLength+1);

    sendChangeMessage();
  }
  messageStringLock.exit();
}

//-------------------------------------------------------------------------------------------------
// automation:

void MidiMonitorAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case  0: 
    {
      midiFilter.noteOn  = value >= 0.5;
      midiFilter.noteOff = value >= 0.5;
    }
    break;
  case  1: midiFilter.controller      = value >= 0.5;      break;
  case  2: midiFilter.pitchWheel      = value >= 0.5;      break;
  case  3: midiFilter.programChange   = value >= 0.5;      break;
  case  4: midiFilter.aftertouch      = value >= 0.5;      break;
  case  5: midiFilter.channelPressure = value >= 0.5;      break;
  case  6: midiFilter.sysEx           = value >= 0.5;      break;
  case  7: midiFilter.metaEvent       = value >= 0.5;      break;
  case  8: midiFilter.transport       = value >= 0.5;      break;
  case  9: midiFilter.songPosition    = value >= 0.5;      break;
  case 10: midiFilter.machineControl  = value >= 0.5;      break;
  case 11: midiFilter.activeSense     = value >= 0.5;      break;
  case 12: midiFilter.clock           = value >= 0.5;      break;
  case 13: midiFilter.other           = value >= 0.5;      break;
  } // end of switch( parameterIndex )
  markStateAsDirty();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void MidiMonitorAudioModule::initializeAutomatableParameters()
{
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "Notes", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "Controllers", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #002:
  p = new AutomatableParameter(lock, "PitchWheel", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #003:
  p = new AutomatableParameter(lock, "ProgramChange", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #004:
  p = new AutomatableParameter(lock, "Aftertouch", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #005:
  p = new AutomatableParameter(lock, "ChannelPressure", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #006:
  p = new AutomatableParameter(lock, "SystemExclusive", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #007:
  p = new AutomatableParameter(lock, "MetaEvents", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #008:
  p = new AutomatableParameter(lock, "Transport", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #009:
  p = new AutomatableParameter(lock, "SongPosition", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #010:
  p = new AutomatableParameter(lock, "MachineControl", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #011:
  p = new AutomatableParameter(lock, "ActiveSense", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #012:
  p = new AutomatableParameter(lock, "Clock", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // #013:
  p = new AutomatableParameter(lock, "Others", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  // make a call to parameterChanged for each parameter in order to set up the midiFilter to 
  // reflect the initial values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

// construction/destruction:

MidiMonitorModuleEditor::MidiMonitorModuleEditor(CriticalSection *newPlugInLock, MidiMonitorAudioModule* newMidiMonitorAudioModule) 
  : AudioModuleEditor(newMidiMonitorAudioModule)
{
  // assign the pointer to the rosic::MidiMonitor object to be used as aduio engine:
  jassert(newMidiMonitorAudioModule != NULL ); // you must pass a valid module here
  midiMonitorModuleToEdit = newMidiMonitorAudioModule;
  midiMonitorModuleToEdit->addChangeListener(this);

  outputDisplay = new RTextEditor( juce::String(("OutputDisplay")) );
  outputDisplay->setMultiLine(true, true); 
  outputDisplay->setReadOnly(true);
  outputDisplay->setScrollbarsShown(false);
  addAndMakeVisible( outputDisplay );

  addWidget( eventFilterLabel = new RTextField( juce::String(("Event Filter"))) );
  eventFilterLabel->setDescription(juce::String(("Choose, which types of events you want to see")));
  eventFilterLabel->setDescriptionField(infoField);
  eventFilterLabel->setJustification(Justification::centred);

  addWidget( noteButton = new RButton(juce::String(("Notes"))) );
  noteButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("Notes")) );
  noteButton->setDescription(juce::String(("Show note-on/-off events")));
  noteButton->setDescriptionField(infoField);
  noteButton->setClickingTogglesState(true);

  addWidget( controllerButton = new RButton(juce::String(("Controllers"))) );
  controllerButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("Controllers")) );
  controllerButton->setDescription(juce::String(("Show control-change events")));
  controllerButton->setDescriptionField(infoField);
  controllerButton->setClickingTogglesState(true);

  addWidget( pitchWheelButton = new RButton(juce::String(("Pitch Wheel"))) );
  pitchWheelButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("PitchWheel")) );
  pitchWheelButton->setDescription(juce::String(("Show pitch-wheel events")));
  pitchWheelButton->setDescriptionField(infoField);
  pitchWheelButton->setClickingTogglesState(true);

  addWidget( programChangeButton = new RButton(juce::String(("Program Changes"))) );
  programChangeButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("ProgramChange")) );
  programChangeButton->setDescription(juce::String(("Show program-change events")));
  programChangeButton->setDescriptionField(infoField);
  programChangeButton->setClickingTogglesState(true);

  addWidget( aftertouchButton = new RButton(juce::String(("Aftertouch"))) );
  aftertouchButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("Aftertouch")) );
  aftertouchButton->setDescription(juce::String(("Show aftertouch events")));
  aftertouchButton->setDescriptionField(infoField);
  aftertouchButton->setClickingTogglesState(true);

  addWidget( channelPressureButton = new RButton(juce::String(("Channel Pressure"))) );
  channelPressureButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("ChannelPressure")) );
  channelPressureButton->setDescription(juce::String(("Show channel pressure events")));
  channelPressureButton->setDescriptionField(infoField);
  channelPressureButton->setClickingTogglesState(true);

  addWidget( sysExButton = new RButton(juce::String(("System Exclusive"))) );
  sysExButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("SystemExclusive")) );
  sysExButton->setDescription(juce::String(("Show system exclusive events")));
  sysExButton->setDescriptionField(infoField);
  sysExButton->setClickingTogglesState(true);

  addWidget( metaEventButton = new RButton(juce::String(("Meta Events"))) );
  metaEventButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("MetaEvents")) );
  metaEventButton->setDescription(juce::String(("Show meta events")));
  metaEventButton->setDescriptionField(infoField);
  metaEventButton->setClickingTogglesState(true);

  addWidget( transportButton = new RButton(juce::String(("Transport"))) );
  transportButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("Transport")) );
  transportButton->setDescription(juce::String(("Show transport control events")));
  transportButton->setDescriptionField(infoField);
  transportButton->setClickingTogglesState(true);

  addWidget( songPositionButton = new RButton(juce::String(("Song Position"))) );
  songPositionButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("SongPosition")) );
  songPositionButton->setDescription(juce::String(("Show song position events")));
  songPositionButton->setDescriptionField(infoField);
  songPositionButton->setClickingTogglesState(true);

  addWidget( machineControlButton = new RButton(juce::String(("Machine Control"))) );
  machineControlButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("MachineControl")) );
  machineControlButton->setDescription(juce::String(("Show machine control events")));
  machineControlButton->setDescriptionField(infoField);
  machineControlButton->setClickingTogglesState(true);

  addWidget( activeSenseButton = new RButton(juce::String(("Active Sense"))) );
  activeSenseButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("ActiveSense")) );
  activeSenseButton->setDescription(juce::String(("Show active sense events")));
  activeSenseButton->setDescriptionField(infoField);
  activeSenseButton->setClickingTogglesState(true);

  addWidget( clockButton = new RButton(juce::String(("Clock"))) );
  clockButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("Clock")) );
  clockButton->setDescription(juce::String(("Show clock events")));
  clockButton->setDescriptionField(infoField);
  clockButton->setClickingTogglesState(true);

  addWidget( otherButton = new RButton(juce::String(("Others"))) );
  otherButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("Others")) );
  otherButton->setDescription(juce::String(("Show other events")));
  otherButton->setDescriptionField(infoField);
  otherButton->setClickingTogglesState(true);

  addWidget( clearButton = new RButton(juce::String(("Clear"))) );
  //clearButton->assignParameter( midiMonitorModuleToEdit->getParameterByName(("Clear")) );
  clearButton->setDescription(juce::String(("Clear screen")));
  clearButton->setDescriptionField(infoField);
  clearButton->setClickingTogglesState(false);
  clearButton->addRButtonListener(this);

  // set up the widgets:
  updateWidgetsAccordingToState();
  updateScreen();

  setSize(600, 440);
}

MidiMonitorModuleEditor::~MidiMonitorModuleEditor()
{
  midiMonitorModuleToEdit->removeChangeListener(this);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void MidiMonitorModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == clearButton )
    clearScreen(true);
}

void MidiMonitorModuleEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  if( objectThatHasChanged == midiMonitorModuleToEdit )
    updateScreen();
}

void MidiMonitorModuleEditor::resized()
{
  presetSectionPosition = INVISIBLE;
  AudioModuleEditor::resized();
  int x = 0;
  int y = getHeadlineBottom();
  int w = getWidth()-136;
  int h = infoField->getY()-y;

  outputDisplay->setBounds(x, y+4, w-8, h-4);

  x = outputDisplay->getRight();
  w = getWidth()-x;

  eventFilterLabel->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  noteButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  controllerButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  pitchWheelButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  programChangeButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  aftertouchButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  channelPressureButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  sysExButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  metaEventButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  transportButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  songPositionButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  machineControlButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  activeSenseButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  clockButton->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  otherButton->setBounds(x+4, y+4, w-8, 16);

  y = outputDisplay->getBottom();
  w = w/2;
  x = outputDisplay->getRight() + w;
  clearButton->setBounds(x-32, y-32, 64, 24);

  updateScreen(); // wihtout that, re-opening of the GUI will have a wrong caret-position when 
                  // the screen is full (in VSTHost)
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void MidiMonitorModuleEditor::clearScreen(bool clearAlsoMessageString)
{
  outputDisplay->setText(juce::String(), false);
  if( clearAlsoMessageString == true && midiMonitorModuleToEdit != NULL)
  {
    midiMonitorModuleToEdit->messageStringLock.enter();
    midiMonitorModuleToEdit->midiMessageString = juce::String();
    midiMonitorModuleToEdit->messageStringLock.exit();
  }
}

void MidiMonitorModuleEditor::updateScreen()
{
  if( midiMonitorModuleToEdit == NULL )
    return;

  midiMonitorModuleToEdit->messageStringLock.enter();

  // use a local string to pass to the editor because the editor's setText() method will update
  // the text asynchronously (albeit in the same thread as this) - so we cannot be sure that this 
  // thread still has the lock inside setText()...and we need the lock because the AudioModule is   
  // constantly messing with its messageString:
  juce::String localString = midiMonitorModuleToEdit->midiMessageString;

  int length = localString.length();
  outputDisplay->setText(localString, false);
  outputDisplay->setCaretPosition(length);

  midiMonitorModuleToEdit->messageStringLock.exit();
}


