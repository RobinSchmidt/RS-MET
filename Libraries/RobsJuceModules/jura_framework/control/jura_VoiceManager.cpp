
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
    triggerVoice(activateAndGetLastIdleVoice(), key, vel);
  else
    stealVoice(key, vel);
}
// needs tests

void rsVoiceManager::noteOff(int key)
{
  for(int i = 0; i < numActiveVoices; i++) {
    int j = activeVoices[i];
    if(voiceStates[j].key == key) {
      releaseVoice(j);
      return; }} // should be ok to return because we assume that a given key can only be held
                 // in one voice at a time
}
// needs tests

void rsVoiceManager::setPitchBend(int pitchBendValue)
{

}
// may not be needed - we handle PitchBend as moudlator via the modulation system. this is more
// flexible and laso more convenient to implement

/*
void rsVoiceManager::perSampleUpdatePreRender()
{

}

void rsVoiceManager::perSampleUpdatePostRender()
{

}
*/

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

int rsVoiceManager::activateAndGetLastIdleVoice()
{
  int idleIndex  = getNumIdleVoices() - 1;
  jassert(idleIndex >= 0 && idleIndex < getNumIdleVoices());
  int voiceToUse = idleVoices[idleIndex];     // use the last from the array of idle voices
  jassert(voiceToUse >= 0 && voiceToUse < getNumVoices());
  idleVoices[idleIndex] = -1;                 // -1 is code for an invalid voice index
  activeVoices[numActiveVoices] = voiceToUse; // append it to the array of active voices
  numActiveVoices++;
  return voiceToUse;
}

void rsVoiceManager::triggerVoice(int voiceIndex, int key, int vel)
{
  voiceStates[voiceIndex].pitch  = getPitchForKey(key);
  voiceStates[voiceIndex].vel01  = double(vel) / 127.0;
  voiceStates[voiceIndex].key    = key;
  voiceStates[voiceIndex].isHeld = true;
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

void rsVoiceManager::releaseVoice(int i)
{
  voiceStates[i].isHeld = false;
  //releasingVoices.push_back(i);

  // more to do....

  int dummy = 0;
}

template<class T>
void removeElement(T* x, int length, int index)
{
  for(int i = index; i < length-1; i++)
    x[i] = x[i+1];
}
// move to RAPT::rsArrayTools

void rsVoiceManager::deactivateVoice(int activeIndex)
{
  int voiceIndex = activeVoices[activeIndex];
  removeElement(&activeVoices[0], numActiveVoices, activeIndex); 
  idleVoices[getNumIdleVoices()] = voiceIndex;
  numActiveVoices--;
  activeVoices[numActiveVoices] = -1;
}
// needs test
