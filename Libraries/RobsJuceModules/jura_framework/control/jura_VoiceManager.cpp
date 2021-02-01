
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
  else if(message.isAftertouch())
    setAfterTouch(message.getAfterTouchValue());
  else if(message.isChannelPressure())
    setChannelPressure(message.getChannelPressureValue());
}
// ToDo: handle reset messages - provide a midiReset callback for that. the voice manager should 
// override it to call its reset function... isAllNotesOff, isAllSoundOff
// maybe handle these in an override in rsVoiceManager - oh, and also handle sustain

// isResetAllControllers, 
// isSustainPedalOn/Off, isSostenutoPedalOn/Off, isSoftPedalOn/Off...but they are controllers and
// should perhaps handled in the controller callback - they don't warrant to introduce their own
// callbacks


//=================================================================================================

void rsVoiceManager::setMaxNumVoices(int newNumber)
{
  maxNumVoices    = newNumber;
  numVoices       = RAPT::rsMin(numVoices,       maxNumVoices);
  numActiveVoices = RAPT::rsMin(numActiveVoices, maxNumVoices);
  activeVoices.resize(maxNumVoices);
  idleVoices.resize(maxNumVoices);
  releasingVoices.reserve(maxNumVoices);
  voiceStates.resize(maxNumVoices);
  killCounters.resize(maxNumVoices);
  reset();
}

/*
void rsVoiceManager::handleMidiMessage(const juce::MidiMessage& msg, 
  rsMidiMessageHandler::MidiHandleInfo* info)
{
  if( msg.isNoteOn() )
    noteOn(msg.getNoteNumber(), msg.getVelocity());
  else if( msg.isNoteOff() )
    noteOff(msg.getNoteNumber());
  else if( msg.isPitchWheel() )
    setPitchBend(msg.getPitchWheelValue());
}
*/

int rsVoiceManager::noteOnReturnVoice(int key, int vel)
{
  if(vel == 0)
    return noteOffReturnVoice(key);
  else if(numActiveVoices < numVoices)  {
    int i = activateAndGetLastIdleVoice();
    triggerVoice(i, key, vel);
    return i; }
  else
    return stealVoice(key, vel);
}

int rsVoiceManager::noteOffReturnVoice(int key)
{
  for(int i = 0; i < numActiveVoices; i++) {
    int j = activeVoices[i];
    if(voiceStates[j].key == key) {
      releaseVoice(j);
      return j; }}                    // should be ok to return now, see comment below
  return -1;
  // When we found a voice that's holding the key, we release it and directly return because we 
  // assume that a given key can only be held in one voice at a time. Could there be situations 
  // where this is not true, i.e. we have several voices simultaneously palying the same note? If 
  // so, what should we do on noteOff? Release all voices that hold the key? Or only the oldest or
  // newest?
}
// ToDo: if sustain is active, we should not release it - instead just set the isHeld flag false.
// as soon as we receive a setSustainOff message, we loop through all active voices and check their
// isHeld flag and if it's false, we release the voice

void rsVoiceManager::setPitchBend(int pitchBendValue)
{

  int dummy = 0;
}


/*
void rsVoiceManager::perSampleUpdatePreRender()
{

}

void rsVoiceManager::perSampleUpdatePostRender()
{

}
*/

void rsVoiceManager::findAndKillFinishedVoices()
{
  jassert(voicesBuffer != nullptr);  // must be assigend or else access violations will occur

  double *vb = voicesBuffer;
  for(size_t i = 0; i < releasingVoices.size(); i++)
  {
    int    vi = releasingVoices[i];              // voice index
    double va = RAPT::rsMax(vb[2*i], vb[2*i+1]); // voice output amplitude
    killCounters[vi]--;                          // countdown
    if(va > killThreshold)                       // use strictly greater to allow zero threshold
      killCounters[vi] = killTimeSamples;        // reset counter, if output is above threshold
    if(killCounters[vi] == 0)
      deactivateVoice(vi);
  }

  // what if the amplitude of the voice does not go to zero?

  // implementation of va = .. needs to be changed when numChannels != 2 ...maybe do this later
}

void rsVoiceManager::reset()
{
  numActiveVoices = 0;
  for(int i = 0; i < maxNumVoices; i++)
  {
    activeVoices[i] = -1;                        // code for invalid voice index
    idleVoices  [i] =  jmax(numVoices-1-i, -1);  // 0th voice is in last slot for first grab
    voiceStates[i].reset();
    killCounters[i] =  0;
  }
  newestVoice = 0;
  releasingVoices.clear();
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
  //voiceStates[voiceIndex].isHeld = true;
  newestVoice = voiceIndex;
}

int rsVoiceManager::stealVoice(int key, int vel)
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

  return 0;
}

void rsVoiceManager::releaseVoice(int i)
{
  //voiceStates[i].isHeld = false;
  if(killMode == KillMode::immediately)
  {
    deactivateVoice(i);
  }
  else
  {
    releasingVoices.push_back(i);
    killCounters[i] = killTimeSamples;
  }
}


template<class T>
void removeElement(T* x, int length, int index)
{
  for(int i = index; i < length-1; i++)
    x[i] = x[i+1];
}
// move to RAPT::rsArrayTools - maybe change order of length and index for better consistency with
// findIndexOf

void rsVoiceManager::deactivateVoice(int voiceIndex)
{
  int activeIndex = RAPT::rsArrayTools::findIndexOf(&activeVoices[0], voiceIndex, numActiveVoices);
  removeElement(&activeVoices[0], numActiveVoices, activeIndex); 
  idleVoices[getNumIdleVoices()] = voiceIndex;
  numActiveVoices--;
  activeVoices[numActiveVoices] = -1;

  bool dbg = RAPT::rsRemoveFirstOccurrence(releasingVoices, voiceIndex);
  //jassert(dbg); // voices that are deactivated are supposed to have been in release-mode before
  // ...but maybe that doesn't need to be the case? Are there situations when a voice is directly
  // shut off without going into the release phase before?
}
// needs test

/*

ToDo:

-what about thread safety? would we use a mutex?
-maybe we should also store the most recent note an use that to update the monophonic modValue
 where applicable

see also:
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=558446 Polyphonic voice assignment and retrigger

midi-specs:
https://www.midi.org/specifications-old/item/table-1-summary-of-midi-message
https://www.midi.org/specifications

maybe we should respond to midi reset with calling our reset

*/