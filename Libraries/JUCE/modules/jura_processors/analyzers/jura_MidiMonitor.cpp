
//-------------------------------------------------------------------------------------------------
// construction/destruction:

MidiMonitorAudioModule::MidiMonitorAudioModule(CriticalSection *newPlugInLock) 
  : AudioModule(newPlugInLock)
{
  moduleName = juce::String(("MidiMonitor"));
  initializeAutomatableParameters();
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