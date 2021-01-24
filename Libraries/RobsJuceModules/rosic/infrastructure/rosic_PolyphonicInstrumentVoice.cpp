
//-------------------------------------------------------------------------------------------------
// construction/destruction:

PolyphonicInstrumentVoice::PolyphonicInstrumentVoice()
{
  // init member variables:
  level          = 0.0;
  levelByKey     = 0.0;
  levelByVel     = 0.0;
  pitchBendRange = 12.0;
  ampRampTime    = 10.0;
  glideTime      = 100.0;
  sampleRate     = 44100.0;
  glideIsActive  = true;
  tuningTable    = NULL;

  // this call will init a few more members:
  reset();
}

PolyphonicInstrumentVoice::~PolyphonicInstrumentVoice()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void PolyphonicInstrumentVoice::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
}

void PolyphonicInstrumentVoice::setLevel(double newLevel)
{
  level = newLevel;

  if( !isReleasing )
    prepareForAmplitudeRamp(getCurrentNoteKey(), getCurrentNoteVelocity(), true);
  else
    prepareForAmplitudeRamp(noteBeingReleased, noteBeingReleasedVel, true);

  //prepareForAmplitudeRamp(getCurrentNoteKey(), getCurrentNoteVelocity(), true);
}

void PolyphonicInstrumentVoice::setLevelByKey(double newLevelByKey)
{
  levelByKey = newLevelByKey;

  if( !isReleasing )
    prepareForAmplitudeRamp(getCurrentNoteKey(), getCurrentNoteVelocity(), true);
  else
    prepareForAmplitudeRamp(noteBeingReleased, noteBeingReleasedVel, true);
}

void PolyphonicInstrumentVoice::setLevelByVel(double newLevelByVel)
{
  levelByVel = newLevelByVel;

  if( !isReleasing )
    prepareForAmplitudeRamp(getCurrentNoteKey(), getCurrentNoteVelocity(), true);
  else
    prepareForAmplitudeRamp(noteBeingReleased, noteBeingReleasedVel, true);


  //prepareForAmplitudeRamp(getCurrentNoteKey(), getCurrentNoteVelocity(), true);
}

void PolyphonicInstrumentVoice::setMasterTuneA4(double newTuneA4)
{
  if( tuningTable != NULL )
    tuningTable->setMasterTuneA4(newTuneA4);

  mutex.lock();
  if( noteList.empty() )
    targetFrequency = newTuneA4;
  else
    targetFrequency = getNoteFrequency( noteList.front().getKey() );
  mutex.unlock();

  targetPitch = RAPT::rsFreqToPitch(targetFrequency);
}

void PolyphonicInstrumentVoice::setPitchWheelRange(double newRange)
{
  if( newRange >= 0.0 )
    pitchBendRange = newRange;
  pitchBendFactor = RAPT::rsPitchOffsetToFreqFactor(pitchBendRange*pitchBend);
}

void PolyphonicInstrumentVoice::setGlideTime(double newGlideTime)
{
  glideTime   = newGlideTime;
  ampRampTime = max(glideTime, 10.0);
}

void PolyphonicInstrumentVoice::setGlideMode(bool shouldGlide)
{
  glideIsActive = shouldGlide;
}

void PolyphonicInstrumentVoice::setBeatsPerMinute(double /*newBeatsPerMinute*/)
{
  // this function here is empty and serves only as dummy for the PolyphonicInstrument class which
  // may call this - if you need sync-functionality you will need to override this in your 
  // subclass.
}

void PolyphonicInstrumentVoice::setTuningTable(rosic::TuningTable *newTable)
{
  tuningTable = newTable;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double PolyphonicInstrumentVoice::getLevel()
{
  return level;
}

double PolyphonicInstrumentVoice::getLevelByKey()
{
  return levelByKey;
}

double PolyphonicInstrumentVoice::getLevelByVel()
{
  return levelByVel;
}

double PolyphonicInstrumentVoice::getMasterTuneA4()
{
  if( tuningTable != NULL )
    return tuningTable->getMasterTuneA4();
  else
    return 440.0;
}

double PolyphonicInstrumentVoice::getPitchWheelRange()
{
  return pitchBendRange;
}

double PolyphonicInstrumentVoice::getGlideTime()
{
  return glideTime;
}

bool PolyphonicInstrumentVoice::isInGlideMode()
{
  return glideIsActive;
}

int PolyphonicInstrumentVoice::getCurrentNoteKey()
{
  int result;

  mutex.lock();
  if( noteList.empty() )
    result = -1;
  else
    result = noteList.front().getKey();
  mutex.unlock();

  return result;
}

int PolyphonicInstrumentVoice::getCurrentNoteVelocity()
{
  int result;

  mutex.lock();
  if( noteList.empty() )
    result = -1;
  else
    result = noteList.front().getVelocity();
  mutex.unlock();

  return result;
}

bool PolyphonicInstrumentVoice::hasNoteInList(int noteKeyToCheckFor)
{
  bool result = false;

  mutex.lock();
  if( noteList.empty() )
    result = false;
  else
  {
    // search for the note through the list via a list-iterator:
    list<MidiNoteEvent>::iterator i;
    for(i = noteList.begin(); i != noteList.end(); i++)
    {
      if( i->getKey() == noteKeyToCheckFor )
        result = true; // the note was found in the list
    }
  }
  mutex.unlock();

  return result;
}

//-------------------------------------------------------------------------------------------------
// event handling:

void PolyphonicInstrumentVoice::noteOn(int newKey, int newVelocity, int newDetune)
{
  mutex.lock();

  if( newVelocity == 0 )
    noteOff(newKey);
  else // velocity was not zero
  {
    // check if the note-list is empty - if so, trigger a new note. if not, glide to the new note:
    if( noteList.empty() || glideIsActive == false )
      triggerNote(newKey, newVelocity, newDetune);
    else
      glideToNote(newKey, newVelocity, newDetune);

    // add the new note to our list:
    MidiNoteEvent newNote(newKey, newVelocity, newDetune);
    noteList.push_front(newNote);
  }

  // setup some bookkeeping variables:
  currentNoteAge    = 0;
  isSilent          = false;
  //adjustAmplitude();

  mutex.unlock();
}

void PolyphonicInstrumentVoice::noteOff(int newKey)
{
  mutex.lock();

  // when the list is already empty on note-off, we are already in a relase phase - nothing more
  // to do here:
  if( noteList.empty() )
  {
    mutex.unlock();
    return;
  }

  MidiNoteEvent releasedNote(newKey, 0);
  noteList.remove(releasedNote);
  // this removal relies on the operator '==' defined in MidiNoteEvent which equates two
  // instances when they have the same value for their 'key' member

  // Check if the note-list is empty now. If so, trigger a release, otherwise slide to the note
  // at the beginning of the list (this is the most recent one which is still in the list):
  if( noteList.empty() )
  {
    triggerRelease(releasedNote.getKey(), releasedNote.getVelocity());
  }
  else
    glideToNote(noteList.front().getKey(), noteList.front().getVelocity(),
      (int) noteList.front().getDetune());

  mutex.unlock();
}

void PolyphonicInstrumentVoice::allNotesOff()
{
  mutex.lock();

  while( !noteList.empty() )
  {
    MidiNoteEvent releasedNote = noteList.back();
    noteList.remove(releasedNote);
    triggerRelease(releasedNote.getKey(), releasedNote.getVelocity());
  }

  mutex.unlock();
}

void PolyphonicInstrumentVoice::triggerNote(int newKey, int newVelocity, int newDetune)
{
  /*
  //targetLev
  targetLevel             = 0.0;   // will be set in prepareForAmpRamp
  targetAmplitude         = 1.0;   // will be set in prepareForAmpRamp
  currentLevel            = 0.0; // redundant with currentAmplitude
  currentAmplitude        = 1.0; // resetting this here leads to clicks
  levelIncPerSample       = 0.0;   // redundant with levelIncPerSample
  ampFactorPerSample      = 1.0;  // will be set in prepareForAmpRamp
  remainingAmpRampSamples = 0;   // will be set in prepareForAmpRamp
  */

  targetPitch             = (double) newKey + 0.01 * (double) newDetune;
  currentPitch            = targetPitch;
  targetFrequency         = getNoteFrequency(currentPitch);
  currentFrequency        = targetFrequency;
  pitchIncPerSample       = 0.0;
  freqFactorPerSample     = 1.0;
  remainingGlideSamples   = 0;

  isSilent               = false;
  isReleasing            = false;
  noteBeingReleased      = -1;
  noteBeingReleasedVel   = 0;

  //remainingAmpRampSamples = roundToInt(0.001*ampRampTime*sampleRate);
  //prepareForAmplitudeRamp((double)newKey, (double)newVelocity, false);
  prepareForAmplitudeRamp((double)newKey, (double)newVelocity, true);
}

void PolyphonicInstrumentVoice::glideToNote(int newKey, int newVelocity, int newDetune)
{
  // \todo: wrap this stuff in a member function prepareForPitchGlide similar to
  // prepareForAmplitudeRamp

  targetPitch   = (double) newKey + 0.01 * (double) newDetune;

  // is this still necesarry, now that we directly update currentPitch as member variable?:
  if( tuningTable != NULL )
    currentPitch = RAPT::rsFreqToPitch(currentFrequency, tuningTable->getMasterTuneA4() );
  else
    currentPitch = RAPT::rsFreqToPitch(currentFrequency, 440.0);

  double pitchDelta     = targetPitch - currentPitch;

  targetFrequency         = getNoteFrequency(targetPitch);
  remainingGlideSamples   = roundToInt(0.001*glideTime*sampleRate);

  if( remainingGlideSamples < 1 )
  {
    // immediate glide:
    currentFrequency    = targetFrequency;
    freqFactorPerSample = 1.0;
  }
  else
  {
    // calculate frequency multiplier per sample (multiplication is done in getSample)
    pitchIncPerSample   = pitchDelta / (double) remainingGlideSamples;
    freqFactorPerSample = RAPT::rsPitchOffsetToFreqFactor(pitchIncPerSample);
  }

  isReleasing           = false;
  noteBeingReleased     = -1;
  noteBeingReleasedVel  = 0;


  //remainingAmpRampSamples = roundToInt(0.001*ampRampTime*sampleRate);
  prepareForAmplitudeRamp((double)newKey, (double)newVelocity, true);
}

void PolyphonicInstrumentVoice::triggerRelease(int noteToBeReleased, int noteToBeReleasedVel)
{
  isReleasing          = true;
  noteBeingReleased    = noteToBeReleased;
  noteBeingReleasedVel = noteToBeReleasedVel;
}

void PolyphonicInstrumentVoice::setPitchBend(double newPitchBend)
{
  pitchBend       = newPitchBend;
  pitchBendFactor = RAPT::rsPitchOffsetToFreqFactor(pitchBendRange*pitchBend);
}

void PolyphonicInstrumentVoice::reset()
{
  mutex.lock();
  noteList.clear();
  mutex.unlock();

  targetLevel                   = 0.0;
  targetAmplitude               = 1.0;
  currentLevel                  = 0.0;
  currentAmplitude              = 1.0;
  levelIncPerSample             = 0.0;
  ampFactorPerSample            = 1.0;
  remainingAmpRampSamples       = 0;

  targetPitch                   = 69.0;
  currentPitch                  = 69.0;
  targetFrequency               = 440.0;
  currentFrequency              = 440.0;
  pitchIncPerSample             = 0.0;
  freqFactorPerSample           = 1.0;
  remainingGlideSamples         = 0;

  pitchBend                     = 0.0;
  pitchBendFactor               = 1.0;
  currentPitchWithPitchBend     = 69.0;
  currentFrequencyWithPitchBend = 440.0;

  currentNoteAge                = 0;
  noteBeingReleased             = -1;
  noteBeingReleasedVel          = 0;
  isSilent                      = true;
  isReleasing                   = false;
}

//-------------------------------------------------------------------------------------------------
// others:

double PolyphonicInstrumentVoice::getNoteFrequency(int noteNumber)
{
  if( tuningTable != NULL )
    return tuningTable->getFrequency(noteNumber);
  else
    return RAPT::rsPitchToFreq((double) noteNumber);
}

double PolyphonicInstrumentVoice::getNoteFrequency(double noteNumber)
{
  if( tuningTable != NULL )
    return tuningTable->getFrequency(noteNumber);
  else
    return RAPT::rsPitchToFreq(noteNumber);
}

void PolyphonicInstrumentVoice::prepareForAmplitudeRamp(double newKey, double newVel, 
  bool shouldRamp)
{
  if( newKey < 0 || newVel <= 0 )
    shouldRamp = false; // function was called while no note was being played

  targetLevel     = level + levelByVel*(newVel-64.0)/63.0 + levelByKey*(newKey-64.0)/63.0;
  targetAmplitude = RAPT::rsDbToAmp(targetLevel);



  //if( shouldRamp == false || remainingAmpRampSamples < 1 )
  if( shouldRamp == false )
  {
    remainingAmpRampSamples = 0;
    currentLevel            = targetLevel;
    currentAmplitude        = targetAmplitude;
    levelIncPerSample       = 0.0;
    ampFactorPerSample      = 1.0;
  }
  else
  {
    // calculate level increment and amplitude multiplier per sample - multiplication and update
    // of members currentLevel and currentAmplitude is then done in getSample:
    remainingAmpRampSamples = roundToInt(0.001*ampRampTime*sampleRate);
    //double levelDelta  = targetLevel - currentLevel;
    double levelDelta  = targetLevel - RAPT::rsAmpToDb(currentAmplitude);
    levelIncPerSample  = levelDelta / (double) remainingAmpRampSamples;
    ampFactorPerSample = RAPT::rsDbToAmp(levelIncPerSample);
  }
}
