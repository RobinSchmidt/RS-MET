
void rsMidiMessageDispatcher::handleMidiMessage(MidiMessage message)
{
  ScopedLock scopedLock(*lock);
  if( message.isNoteOn() )
    noteOn(message.getNoteNumber(), message.getVelocity());
  else if( message.isNoteOff() )
    noteOff(message.getNoteNumber());
  else if( message.isAllNotesOff() )
    allNotesOff();
  else if( message.isController() )
    setMidiController(message.getControllerNumber(), (float) message.getControllerValue());
  else if( message.isPitchWheel() )
    setPitchBend(message.getPitchWheelValue());
  else if (message.isAftertouch())
    setAfterTouch(message.getAfterTouchValue());
  else if (message.isChannelPressure())
    setChannelPressure(message.getChannelPressureValue());
}


//=================================================================================================

void rsVoiceManager::noteOn(int noteNumber, int velocity)
{

}

void rsVoiceManager::noteOff(int noteNumber)
{

}

void rsVoiceManager::setPitchBend(int pitchBendValue)
{

}