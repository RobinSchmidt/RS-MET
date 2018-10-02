//#include "rosic_PolyphonicInstrument.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

PolyphonicInstrument::PolyphonicInstrument(int numVoicesToAllocate)
{
  numAllocatedVoices      = numVoicesToAllocate;
  numPlayableVoices       = numVoicesToAllocate-1;
  numActiveVoices         = 0;
  voicePointers           = NULL;
  masterLevel             = 0.0;
  masterAmplitude         = 1.0;
  masterLevelByVoices     = 0.0;
  midSideRatio            = 0.5;
  midScale                = sqrt(0.5);
  sideScale               = sqrt(0.5);
  mostRecentNote          = -1;
  mostRecentNoteVel       = 0;
  mostRecentNoteDetune    = 0;
  fillPendingNoteOffArray(false);
  sustainIsActive         = false;
}

PolyphonicInstrument::~PolyphonicInstrument()
{
  //if( voiceArray != NULL )
    //delete[] voiceArray;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void PolyphonicInstrument::setSampleRate(double newSampleRate)
{
  if( voicePointers == NULL )
    return;
  for(int i=0; i<numAllocatedVoices; i++)
    voicePointers[i]->setSampleRate(newSampleRate);
}

void PolyphonicInstrument::setNumPlayableVoices(int newNumVoices)
{
  if( voicePointers == NULL )
    return;

  if( newNumVoices > 0 && newNumVoices <= numAllocatedVoices-1 )
    numPlayableVoices = newNumVoices;

  // reset all the voices:
  for(int i=0; i<numAllocatedVoices; i++)
    voicePointers[i]->reset();
}

int PolyphonicInstrument::getNumActiveVoices()
{
  return numActiveVoices;
}

bool PolyphonicInstrument::isSilent()
{
  bool result = true;
  for(int i=1; i<numAllocatedVoices; i++) // start at index 1 - ignore template-voice
    result &= voicePointers[i]->isSilent;
  return result;
}

void PolyphonicInstrument::setMasterLevel(double newMasterLevel)
{ 
  masterLevel     = newMasterLevel;
  masterAmplitude = RAPT::rsDbToAmp(masterLevel);
}

void PolyphonicInstrument::setVoiceLevelByKey(double newVoiceLevelByKey)
{ 
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setLevelByKey(newVoiceLevelByKey);
  }
}

void PolyphonicInstrument::setVoiceLevelByVel(double newVoiceLevelByVel)
{ 
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setLevelByVel(newVoiceLevelByVel);
  }
}

void PolyphonicInstrument::setMasterLevelByVoices(double newMasterLevelByVoices)
{ 
  masterLevelByVoices = newMasterLevelByVoices;
}

void PolyphonicInstrument::setMidSideRatio(double newMidSideRatio)
{
  midSideRatio = newMidSideRatio;
  double x = 0.5 * PI * midSideRatio;
  sinCos(x, &sideScale, &midScale);
}

void PolyphonicInstrument::setMasterTuneA4(double newTuneA4)
{
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setMasterTuneA4(newTuneA4);
  }
}

void PolyphonicInstrument::setPitchWheelRange(double newRange)
{
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setPitchWheelRange(newRange);
  }
}

void PolyphonicInstrument::setGlideTime(double newGlideTime)
{
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setGlideTime(newGlideTime);
  }
}

void PolyphonicInstrument::setGlideMode(bool shouldUseGlide)
{
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setGlideMode(shouldUseGlide);
  }
}

void PolyphonicInstrument::setBeatsPerMinute(double newBeatsPerMinute)
{
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setBeatsPerMinute(newBeatsPerMinute);
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int PolyphonicInstrument::getNumPlayableVoices()
{
  return numPlayableVoices;
}

double PolyphonicInstrument::getMasterLevel()
{
  return masterLevel;
}

double PolyphonicInstrument::getVoiceLevelByKey()
{
  if( voicePointers != NULL )
    return voicePointers[0]->getLevelByKey();
  else
    return 0.0;
}

double PolyphonicInstrument::getVoiceLevelByVel()
{
  if( voicePointers != NULL )
    return voicePointers[0]->getLevelByVel();
  else
    return 0.0;
}

double PolyphonicInstrument::getMasterLevelByVoices()
{
  return masterLevelByVoices;
}

double PolyphonicInstrument::getMidSideRatio()
{
  return midSideRatio;
}

double PolyphonicInstrument::getMasterTuneA4()
{
  if( voicePointers != NULL )
    return voicePointers[0]->getMasterTuneA4();
  else
    return 440.0;
}

double PolyphonicInstrument::getPitchWheelRange()
{
  if( voicePointers != NULL )
    return voicePointers[0]->getPitchWheelRange();
  else
    return 12.0;
}

double PolyphonicInstrument::getGlideTime()
{
  if( voicePointers != NULL )
    return voicePointers[0]->getGlideTime();
  else
    return 50.0;
}

bool PolyphonicInstrument::isInGlideMode()
{
  if( voicePointers != NULL )
    return voicePointers[0]->isInGlideMode();
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// event processing:

void PolyphonicInstrument::noteOn(int newNoteNumber, int newVelocity, int newDetune)
{
  if( voicePointers == NULL )
    return;

  mostRecentNote       = newNoteNumber;
  mostRecentNoteVel    = newVelocity;
  mostRecentNoteDetune = newDetune;

  // loop through the voices, voicePointers[0] servers only as template for the others and is not 
  // invoked to process audio:

  int i = 0;

  if( mostRecentNoteVel == 0 )   // a note-off event occured
    noteOff(newNoteNumber);
  else                         // a note-on event occured
  {
    // look at all the voices, to see if ALL of them have reached their ends -
    // only in this case, the modulator-phases are reset:
    bool aVoiceIsPlaying  = false;

    for(i=1; i<=numPlayableVoices; i++)
    {
      if( !voicePointers[i]->isSilent )
        aVoiceIsPlaying = true;
    }
    if( aVoiceIsPlaying == false )
    {
      // reset modulators....
    }

    // loop through the voices to find a free one:
    int  oldestNoteAge = 0;
    int  oldestVoice   = 1;       // voice with the oldest note
    for(i=1; i<=numPlayableVoices; i++)
    {
      // check, if voice i is currently releasing the note, which is coming in,
      // if so, use that voice again:
      // check if the voice i is free:
      //if( voicePointers[i]->getCurrentNoteKey() == mostRecentNote )
      if( voicePointers[i]->isReleasing && voicePointers[i]->noteBeingReleased == mostRecentNote )
      {
        // voice is used for the same note (which is currently releasing) again:
        voicePointers[i]->noteOn(mostRecentNote, mostRecentNoteVel);
        return;  // jump out of the function
      }
      // check if the voice i is free:
      else if(voicePointers[i]->isSilent)
      {
        // voice i is free and can be used for the new note:
        voicePointers[i]->noteOn(mostRecentNote, mostRecentNoteVel);
        return;  // jump out of the function
      }
      // keep track of the oldest note for the case, when no voice is free for the new note 
      // (last note priority assignment):
      else if( voicePointers[i]->currentNoteAge > oldestNoteAge )
      {
        oldestNoteAge = voicePointers[i]->currentNoteAge;
        oldestVoice   = i;
      }
    } // end of for-loop

    // no free voice has been found - set the voice with the oldest note
    // to the new note:
    voicePointers[oldestVoice]->noteOn(mostRecentNote, mostRecentNoteVel);

  } // end of else

}

void PolyphonicInstrument::noteOff(int newNoteNumber)
{
  if( newNoteNumber < 0 || newNoteNumber > 127 )
  {
    DEBUG_BREAK; // invalid note-number
    return;
  }

  if( sustainIsActive == true )
  {
    // if sustain is active we just set a flag - the actual message will then be sent later when 
    // the sustain switch is released again:
    pendingNoteOffs[newNoteNumber] = true;
  }
  else
  {
    // loop through the voices and send the note-off to all those voices which have a note with the 
    // current key in their note-list:
    for(int i=1; i<=numPlayableVoices; i++)
    {
      while( voicePointers[i]->hasNoteInList(newNoteNumber) )
        voicePointers[i]->noteOff(newNoteNumber);
    }
  }
}

void PolyphonicInstrument::allNotesOff()
{
  if( voicePointers == NULL )
    return;
  for(int i=1; i<=numPlayableVoices; i++)
    voicePointers[i]->allNotesOff();
}

void PolyphonicInstrument::setMidiController(int controllerNumber, int controllerValue)
{
  if( controllerNumber == 64 )
  {
    if( controllerValue <= 63 )
      setSustain(false);
    else
      setSustain(true);
  }
  else if( controllerNumber == 66 ) // Sustenuto
  {
    /*
    if( controllerValue <= 63 )
      setSustain(true);
    else
      setSustain(false);
    */
  }
  else if( controllerNumber == 123 ) // all notes off
  {
    allNotesOff();
  }
  /*
  else
    AutomatableModule::setMidiController(controllerNumber, controllerValue);
  */
}

void PolyphonicInstrument::setPitchBend(double newPitchBend)
{
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->setPitchBend(newPitchBend);
  }
}

void PolyphonicInstrument::setSustain(bool shouldBeSustained)
{
  sustainIsActive = shouldBeSustained;

  // when the sustain switch was switched off, we need to discharge all the collected note-off
  // messages now:
  if( sustainIsActive == false )
    dischargePendingNoteOffs();
}

void PolyphonicInstrument::resetAllVoices()
{
  if( voicePointers != NULL )
  {
    for(int i=0; i<numAllocatedVoices; i++)
      voicePointers[i]->reset();
  }
}

/*
//-------------------------------------------------------------------------------------------------
// automation:

void PolyphonicInstrument::parameterChanged(Parameter* parameterThatHasChanged)
{

}
*/

//-------------------------------------------------------------------------------------------------
// internal:

void PolyphonicInstrument::fillPendingNoteOffArray(bool valueToFillWith)
{
  for(int i=0; i<128; i++)
    pendingNoteOffs[i] = valueToFillWith;
}

void PolyphonicInstrument::dischargePendingNoteOffs()
{
  for(int i=0; i<128; i++)
  {
    if( pendingNoteOffs[i] == true )
      noteOff(i);
    pendingNoteOffs[i] = false;
  }
}

