
void rsMidiMessageDispatcher::handleMidiMessage(MidiMessage message)
{
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

void rsVoiceManager::setMaxNumVoices(int newNumber)
{
  maxNumVoices    = newNumber;
  numVoices       = RAPT::rsMin(numVoices,       maxNumVoices);
  numActiveVoices = RAPT::rsMin(numActiveVoices, maxNumVoices);
  activeVoices.resize(maxNumVoices);
  idleVoices.resize(maxNumVoices);
  voiceStates.resize(maxNumVoices);
  reset();
}

void rsVoiceManager::noteOn(int key, int vel)
{
  if(vel == 0)
    noteOff(key);
  else if(numActiveVoices < numVoices)
  {
    // maybe factor out into voiceToUse = activateLastIdleVoice()
    int idleIndex  = getNumIdleVoices() - 1;
    int voiceToUse = idleVoices[idleIndex];     // us the last from the array of idle voices
    idleVoices[idleIndex] = -1;                 // -1 is code for an invalid voice index
    activeVoices[numActiveVoices] = voiceToUse; // append it to the array of active voices
    numActiveVoices++;

    // maybe factor out into triggerVoice(voiceIndex, key, vel):
    voiceStates[voiceToUse].frequency = RAPT::rsPitchToFreq(double(key)); // preliminary
    voiceStates[voiceToUse].normalizedVelocity = double(vel) / 127.0;
    voiceStates[voiceToUse].key = key;
    voiceStates[voiceToUse].isHeld = true;
  }
  else
    stealVoice(key, vel);
}
// needs tests

void rsVoiceManager::noteOff(int key)
{
  for(int i = 0; i < numActiveVoices; i++)
  {
    int j = activeVoices[i];
    if(voiceStates[j].key == key)
    {
      // we have found the voice that needs to be put into release state


      voiceStates[j].isHeld = false;
      //releasingVoices.push_back(j);


      int dummy = 0;
    }
  }



}

void rsVoiceManager::setPitchBend(int pitchBendValue)
{

}

void rsVoiceManager::perSampleUpdate()
{

}

void rsVoiceManager::reset()
{
  numActiveVoices = 0;
  for(int i = 0; i < maxNumVoices; i++)
  {
    idleVoices  [i] =  jmax(numVoices-1-i, -1);  // 0th voice is in last slot for first grab
    activeVoices[i] = -1;                        // code for invalid voice index
    voiceStates[i].reset();
  }
}

void rsVoiceManager::stealVoice(int key, int vel)
{
  switch(stealMode)
  {
  case StealMode::oldest:
  {

    // ...

  } break;


  default:
  {
    jassertfalse;  // stealMode is not correctly set up
  }

  }


  int dummy = 0;
}