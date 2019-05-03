
MidiMessageFilter::MidiMessageFilter()
{
  noteOn          = true;
  noteOff         = true;
  sysEx           = true;
  programChange   = true;
  pitchWheel      = true;
  aftertouch      = true;
  channelPressure = true;
  controller      = true;
  metaEvent       = true;
  activeSense     = false;
  transport       = true;
  clock           = false;
  songPosition    = true;
  machineControl  = true;
  other           = true;
}

MidiMessageFilter::~MidiMessageFilter()
{

}

bool MidiMessageFilter::isMessageSuitable(const juce::MidiMessage &m)
{  
  int hours, minutes, seconds, frames; // needed for the  MMC goto message

  if( m.isNoteOn() )
    return noteOn;
  else if( m.isNoteOff() )
    return noteOff;
  else if( m.isSysEx() )
    return sysEx;
  else if( m.isProgramChange() )
    return programChange;
  else if( m.isPitchWheel() )
    return pitchWheel;
  else if( m.isAftertouch() )
    return aftertouch;
  else if( m.isChannelPressure() )
    return channelPressure;
  else if( m.isController() )
    return controller;
  else if( m.isMetaEvent() )
    return metaEvent;
  else if( m.isActiveSense() )
    return activeSense;
  else if( m.isMidiStart() || m.isMidiStop() || m.isMidiContinue() )
    return transport;
  else if( m.isMidiClock() )
    return clock;
  else if( m.isSongPositionPointer() )
    return songPosition;
  else if( m.isMidiMachineControlMessage() || m.isMidiMachineControlGoto(hours, minutes, seconds, frames) )
    return machineControl;  
  else
    return other;
}
