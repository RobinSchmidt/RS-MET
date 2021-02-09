
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
  releasingVoices.resize(maxNumVoices);
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

template<class T>
bool areEntriesUnique(T* x, int N)
{
  for(int n = 0; n < N; n++) {
    if(RAPT::rsArrayTools::contains(&x[n+1], N-n-1, x[n]))
      return false; }
  return true;
}

bool rsVoiceManager::isInSaneState()
{
  bool ok = true;
  ok &= numVoices <= maxNumVoices;
  ok &= numActiveVoices <= numVoices;
  ok &= numReleasingVoices <= numActiveVoices;
  ok &= areEntriesUnique(&activeVoices[0], numActiveVoices);
  ok &= areEntriesUnique(&releasingVoices[0], numReleasingVoices);

  // ...what else can we check?

  return ok;
}

int rsVoiceManager::noteOnReturnVoice(int key, int vel)
{
  if(vel == 0)
    return noteOffReturnVoice(key);

  // Reuse an existing voice, if desired and possible:
  if(retriggerMode == RetriggerMode::reuseOldVoice) {
    for(int i = 0; i < numActiveVoices; i++) {  // factor out into getVoiceWithKey and reuse it
      int j = activeVoices[i];                  // in noteOffReturnVoice - we do the same thing 
      if(voiceStates[j].key == key) {
        triggerVoice(j, key, vel);
        return j;  }}}

  // No voice was retriggered. We either trigger an idle voice or steal an active voice:
  int k = -1;
  if(numActiveVoices < numVoices)
    k = activateAndGetLastIdleVoice();
  else
    k = stealVoice(key, vel);
  triggerVoice(k, key, vel);
  return k;
}


int rsVoiceManager::noteOffReturnVoice(int key)
{
  int k = noVoice;
  for(int i = 0; i < numActiveVoices; i++) 
  {
    int j = activeVoices[i];
    if(voiceStates[j].key == key)
    {
      releaseVoice(j);
      k = j;
    }
  }
  return k;

  /*
  // old:
  for(int i = 0; i < numActiveVoices; i++) {
    int j = activeVoices[i];
    if(voiceStates[j].key == key) {
      releaseVoice(j);
      return j; }}                    // should be ok to return now, see comment below
  return noVoice;
  */
  // When we found a voice that's holding the key, we release it and directly return because we 
  // assume that a given key can only be held in one voice at a time. Could there be situations 
  // where this is not true, i.e. we have several voices simultaneously playing (as in holding) the
  // same note? If so, what should we do on noteOff? Release all voices that hold the key? Or only
  // the oldest or newest? See here:
  // http://www.martin-finke.de/blog/articles/audio-plugins-016-polyphony/
  // Martin says that yes, that situation could occur. But i think, that's wrong: a noteOff will 
  // be received between the 1st and 2nd noteOn because the keyboard player has to lift the key
  // before they can press it again and that lifting will trigger the noteOff. ...but maybe in a 
  // piano-roll editor, the user could enter overlapping notes on the same key
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
  jassert(voicesBuffer != nullptr);         // must be assigend or access violations will occur
  double *vb = voicesBuffer;

  for(int i = numReleasingVoices-1; i >= 0; i--)
  {
    int vi = releasingVoices[i];            // voice index
    jassert(vi >= 0 && vi < numVoices);
    double va = getVoiceAmplitude(vi);      // voice output amplitude
    killCounters[vi]--;                     // countdown
    if(va > killThreshold)                  // use strictly greater to allow zero threshold
      killCounters[vi] = killTimeSamples;   // reset counter, if output is above threshold
    if(killCounters[vi] == 0)
      deactivateVoice(vi);  
  }
  // It's important to run the loop downward because within the loop, the releasingVoices array 
  // gets manipulated (in the deactivateVoice call) but it's also accessed in vi = ... and that 
  // works only when traversing the array in reverse order.


  /*
  // this is botched - we manipulate the releasingVoices array in deactivateVoice but also use it 
  // for retrieving vi. 
  int numReleasingTmp = numReleasingVoices; // bcs member gets updated in the loop
  for(int i = 0; i < numReleasingTmp; i++)
  {
    int    vi = releasingVoices[i];         // voice index
    jassert(vi >= 0 && vi < numVoices);
    double va = getVoiceAmplitude(vi);      // voice output amplitude
    killCounters[vi]--;                     // countdown
    if(va > killThreshold)                  // use strictly greater to allow zero threshold
      killCounters[vi] = killTimeSamples;   // reset counter, if output is above threshold
    if(killCounters[vi] == 0)
      deactivateVoice(vi);                  // updates member numReleasingVoices, hence the tmp
  }
  */

  // What if the amplitude of the voice does not go to zero because the user has not connected an 
  // amp-env or that env does not go to zero at the end? The voice will never get killed...but
  // maybe that's ok, i.e. considered user error
}

void rsVoiceManager::reset()
{
  numActiveVoices    = 0;
  numReleasingVoices = 0;
  newestVoice        = 0;
  for(int i = 0; i < maxNumVoices; i++)
  {
    activeVoices[i]    = -1;                        // code for invalid voice index
    releasingVoices[i] = -1;
    idleVoices[i]      =  jmax(numVoices-1-i, -1);  // 0th voice is in last slot for first grab
    killCounters[i]    =  0;
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
  //voiceStates[voiceIndex].isHeld = true;
  newestVoice = voiceIndex;
}

int rsVoiceManager::stealVoice(int key, int vel)
{  
  jassert(numActiveVoices > 0); // steal a voice when none is active? something is fishy!

  int i = 0;  // index in the activeVoices array
  switch(stealMode)
  {
  case StealMode::oldest:  { i = 0;                 } break;
  case StealMode::newest:  { i = numActiveVoices-1; } break;
  case StealMode::noSteal: { i = noVoice;           } break;
  }

  if(i == noVoice) {
    jassertfalse;   // stealMode is not correctly set up
    return noVoice; }

  int k = activeVoices[i];

  // Shift the content of the activeVoices array and put the stolen voice at the end of the array 
  // to keep it sorted by age:
  for(int j = i+1; j < numVoices; j++)
    activeVoices[j-1] = activeVoices[j];
  activeVoices[numVoices-1] = k;

  // Update content of voiceStates[k]:
  triggerVoice(k, key, vel);

  // Remove the voice from the releasingVoicesarray as well, if it's present there:
  i = RAPT::rsArrayTools::findIndexOf(&releasingVoices[0], k, numReleasingVoices);
  if(i > -1)
  {
    for(int j = i+1; j < numReleasingVoices; j++)
      releasingVoices[j-1] = releasingVoices[j];
    numReleasingVoices--;
  }


  //RAPT::rsRemoveFirstOccurrence(releasingVoices, k);

  return k;
}
// ToDo: In "oldest" mode we do not take into account, if the voice is releasing or not. But we 
// probably should: releasing voices should be stolen first. But maybe that should be another mode:
// oldestInRelease. It should first try to find the oldest releasing voice and only when no voice
// is releasing use the other active voices. We can implement this by checking if the 
// releasingVoices array is empty

void rsVoiceManager::releaseVoice(int i)
{
  //voiceStates[i].isHeld = false;
  if(killMode == KillMode::immediately)
  {
    deactivateVoice(i);
  }
  else
  {
    //releasingVoices.push_back(i);

    // What if the voice is already in the releasingVoices array? Could that actually happen?

    if(!RAPT::rsArrayTools::contains(&releasingVoices[0], numReleasingVoices, i))
    {
      releasingVoices[numReleasingVoices] = i;  
      numReleasingVoices++; 
      killCounters[i] = killTimeSamples;
      int dummy = 0;
    }  
    //killCounters[i] = killTimeSamples;
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
  jassert(voiceIndex >= 0 && voiceIndex < numVoices);

  int activeIndex = RAPT::rsArrayTools::findIndexOf(&activeVoices[0], voiceIndex, numActiveVoices);
  jassert(activeIndex >= 0 && activeIndex < numActiveVoices);
  removeElement(&activeVoices[0], numActiveVoices, activeIndex); 
  idleVoices[getNumIdleVoices()] = voiceIndex;
  numActiveVoices--;
  activeVoices[numActiveVoices] = -1;

  int releaseIndex = RAPT::rsArrayTools::findIndexOf(&releasingVoices[0], voiceIndex, 
    numReleasingVoices);
  if(releaseIndex != -1) {
    removeElement(&releasingVoices[0], numReleasingVoices, releaseIndex); 
    numReleasingVoices--;
    releasingVoices[numReleasingVoices] = -1; }


  //bool dbg = RAPT::rsRemoveFirstOccurrence(releasingVoices, voiceIndex);
  //jassert(dbg); // voices that are deactivated are supposed to have been in release-mode before
  // ...but maybe that doesn't need to be the case? Are there situations when a voice is directly
  // shut off without going into the release phase before?
}


/*

Bug:

in ToolChain:
-set maxNumVoices to 4, trigger 5 notes, wait until all are killed, play more notes - the new notes
 are cut off early, as if they are immeditely killed on noteOff - is something wrong about the 
 killSamples array? ah - no - it's not immediate - it seems to be cut off after killSamples. It's 
 as if the killCounter is not reset properly when the voice is still producing audio
-it also seems one voice (3) never gets killed ..oh - another time it was voice 2
-looks like the va (voice-amplitude) doesn't go down to zero and instead setlles to a fixed nonzero
 value - is the voiceBuffer not updated correctly

ToDo:

-what about thread safety? would we use a mutex?
-maybe we should also store the most recent note an use that to update the monophonic modValue
 where applicable

see also:
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=558446 
"Polyphonic voice assignment and retrigger"

http://www.martin-finke.de/blog/articles/audio-plugins-016-polyphony/

midi-specs:
https://www.midi.org/specifications-old/item/table-1-summary-of-midi-message
https://www.midi.org/specifications

maybe we should respond to midi reset with calling our reset


implement sustain, soft and sustenuto pedals: https://en.wikipedia.org/wiki/Piano_pedals

Sustain:
"It raises all the dampers off the strings so that they keep vibrating after the player releases 
the key. "
so: no voice renters release on noteOff until the pedal is turned off again. when it's turned off,
then all notes that have previously received noteOff, should enter release

Soft:
"On the pianos of the late eighteenth to early nineteenth centuries, the pianist could shift from 
the normal three-string (tre corde) position to one in which either two strings (due corde) or 
only one (una corda) would be struck, depending on how far the player depressed the pedal."

Sostenuto: 
"The pedal holds up only dampers that were already raised at the moment that it was depressed"
...so it's a partial sustain affecting only those notes already held, future notes are not 
affected

-instead of killing a voice immediately when it falls below the kill-threshold (for some amount of 
 time), let it enter a short fade-out phase of a few milliseconds to avoid the potential low-level
 clicks that we would otherwise get
-maybe use short or even char instead of int where values are supposed to fit in that range
 (like for midi keys and velocities) - that may save a little memory - we'll see

see romos::VoiceAllocator and rosic::PolyphonicInstrument/Voice for ideas for implementing 
certain behaviors, code and formulas. for example, there's this polyphonic glide feature 
("PolyGlide"?) if i remember correctly, also: glide back, if note is released - stuff like that

*/