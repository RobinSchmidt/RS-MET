//#include "rosic_AcidSequencer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AcidSequencer::AcidSequencer()
{
  sampleRate    = 44100.0;  
  bpm           = 140.0;
  activePattern = 0;
  running       = false;
  countDown     = 0;
  step          = 0;
  sequencerMode = OFF;
  driftError    = 0.0;
  modeChanged   = false;

  for(int k=0; k<=12; k++)
    keyPermissible[k] = true;

  //bpm = 30;  // very slow, for debugging the animation
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AcidSequencer::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
}

void AcidSequencer::setMode(int newMode)
{
  if( newMode >= 0 && newMode < NUM_SEQUENCER_MODES )
  {
    sequencerMode = newMode;
    modeChanged   = true;
  }
}

void AcidSequencer::setKeyPermissible(int key, bool shouldBePermissible)
{
  if( key >= 0 && key <= 12 )
    keyPermissible[key] = shouldBePermissible;
}

void AcidSequencer::toggleKeyPermissibility(int key)
{
  if( key >= 0 && key <= 12 )
    keyPermissible[key] = !keyPermissible[key];
}

void AcidSequencer::copyActiveToClipboard()
{
  AcidPattern* p = &patterns[activePattern];
  int L = p->getNumSteps();
  clipboardLength = L;
  for(int i = 0; i < L; i++)
    clipboard.notes[i] = p->notes[i];
  for(int i = L; i < p->getMaxNumSteps(); i++)
    clipboard.notes[i] = AcidNote();
}

void AcidSequencer::pasteClipboardToActive()
{
  AcidPattern* p = &patterns[activePattern];
  int L = clipboardLength;
  p->setNumSteps(L);
  for(int i = 0; i < L; i++)
     p->notes[i] = clipboard.notes[i];
  for(int i = L; i < p->getMaxNumSteps(); i++)
    p->notes[i] = AcidNote();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

AcidPattern* AcidSequencer::getPattern(int index)
{
  if( index < 0 || index >= numPatterns )
    return NULL;
  else
    return &patterns[index];
}

bool AcidSequencer::modeWasChanged()
{
  bool result = modeChanged;
  modeChanged = false;
  return result;
  // mmm...wouldn't we need mutexes here? the mode changes from the GUI and modeWasChanged
  // is called from the audio-thread - otherwise note-hangs could happen?
}

bool AcidSequencer::isKeyPermissible(int key)
{
  if( key >= 0 && key <= 12 )
    return keyPermissible[key];
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// event handling:

void AcidSequencer::start()
{
  // set up members such that we will trap in the else-branch in the next call to getNote():
  running    = true;


  //countDown  = -1;    // is this wrong?

  // new:
  double secondsToNextStep = RAPT::rsBeatsToSeconds(0.25, bpm);
  double samplesToNextStep = secondsToNextStep * sampleRate;
  countDown                = roundToInt(samplesToNextStep);
  // code is duplicated in getNote() - try to refactor!


  step       = 0;
  driftError = 0.0;
}

void AcidSequencer::stop()
{
  running = false;
}

//-------------------------------------------------------------------------------------------------
// others:

/*=================================================================================================

Ideas:
-Pattern export: Have an "Export" button on the Sequencer GUI that allows to export the current 
 pattern to midi file. User can re-import it in the DAW, turn off the sequencer and then edit the 
 midi data in finer detail than possible in the internal sequencer. It becomes possible to edit 
 the step-lengths individually and to apply partial accents to notes via the velocity (the internal
 DSP allows continuous accent values between 0..1 rather than just "off" or "on"
-Undo: Any modifcation of the current pattern should be undoable, perhaps with a user adjustbale 
 depth (maybe defaulting to 8 or 10) 
-Circular shift of all aspects at once.
-Mouse-select of to-be-transformed (aspects of) events instead of whole rows. Could also apply a 
 cirular Shift of parts of the pattern, e.g. only steps 4-12 are affected, steps 1-3 and 13-16 stay
 in place
-Adjustable number of steps - rather than the hardcoded 16 steps, the user selects, how many steps
 a pattern has. For example, 24 steps could be useful for triolic patterns.
-Perhaps it could make sense to store a bunch different patterns and let the user switch between 
 them by clicks on a button. Maybe there should be copy/paste functionality such that the current
 pattern can be copied into a temporary "clipboard" and then pasted into another pattern. There is
 some infrastructure for that in place already.


*/