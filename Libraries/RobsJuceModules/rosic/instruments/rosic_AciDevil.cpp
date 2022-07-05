//#include "rosic_AciDevil.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AciDevil::AciDevil()
{
  ampScaler        =     1.0;
  oscFreq          =   440.0;
  sampleRate       = 44100.0;
  level            =   -12.0;
  levelByVel       =    12.0;
  accent           =     0.0;
  slideTime        =    60.0;
  cutoff           =  1000.0;
  envUpFraction    =     2.0/3.0;
  normalAttack     =     3.0;
  accentAttack     =     3.0;
  normalDecay      =  1000.0;
  accentDecay      =   200.0;
  normalAmpRelease =     1.0;
  accentAmpRelease =    50.0;
  accentGain       =     0.0;
  subOscGain       =     0.0;
  clipperGain      =     1.0;
  clipperDC        =     0.0;
  pitchWheelRange  =    12.0;
  pitchWheelFactor =     1.0;
  currentNote      =    -1;
  currentVel       =     0;
  noteOffCountDown =     0;
  slideToNextNote  = false;
  idle             = true;

  setEnvMod(12.0);

  oscillator.setWaveTable1(&waveTable1);
  oscillator.setWaveForm1(MipMappedWaveTable::SAW303);
  oscillator.setWaveTable2(&waveTable2);
  oscillator.setWaveForm2(MipMappedWaveTable::SQUARE303);

  subOscillator.setWaveTable1(&subWaveTable1);
  subOscillator.setWaveForm1(MipMappedWaveTable::SAW303);
  subOscillator.setWaveTable2(&subWaveTable2);
  subOscillator.setWaveForm2(MipMappedWaveTable::SQUARE303);

  //mainEnv.setNormalizeSum(true);
  mainEnv.setNormalizeSum(false);

  ampEnv.setAttack(0.0);
  ampEnv.setDecay(1230.0);
  ampEnv.setSustainLevel(0.0);
  ampEnv.setRelease(0.5);
  ampEnv.setTauScale(1.0);

  pitchSlewLimiter.setTimeConstant(60.0);
  ampDeClicker.setTimeConstant(2.0);

  rc1.setTimeConstant(0.0);
  rc2.setTimeConstant(15.0);

  // tweakables:

  /*
  // this parameter-set matches the ABL:
  oscillator.setPulseWidth(43.5);
  subOscillator.setPulseWidth(43.5);
  filter.inputAllpass.setCutoff(70.0);
  filter.inputHighpass.setCutoff(200.0);
  filter.feedbackHighpass.setCutoff(150.0);
  */

  /*
  // this parameter-set matches autodafe's TB303_square_C100_R0_F2_87Hz.wav:
  oscillator.setPulseWidth(47.5);
  subOscillator.setPulseWidth(47.5);
  filter.inputAllpass.setCutoff(75.0);
  filter.inputHighpass.setCutoff(40.0);
  filter.feedbackHighpass.setCutoff(150.0);  // copied from ABL
  */

  /*
  // this parameter-set roughly matches the sample SquareNoResoC2.wav by rv0:
  oscillator.setPulseWidth(45.0);
  subOscillator.setPulseWidth(45.0);
  filter.inputFilter.setAllpassFrequency(60.0);
  filter.inputFilter.setHighpassCutoff(95.0);
  filter.inputFilter.setLowpassCutoff(20000.0);
  filter.feedbackHighpass.setCutoff(150.0);  // copied from ABL
  */


  // measurement for square:
  // sawdown(symmetric) -> tanh(a*x+b) with a = db2amp(36.9), b = -4.37 -> HP@60.8
  // -> filter/dist -> HP@46.9
  oscillator.setPulseWidth(50.0);
  subOscillator.setPulseWidth(50.0);
  hp1.setMode(rsOnePoleFilterDD::HIGHPASS_MZT);
  hp1.setCutoff(60.8);
  hp2.setMode(rsOnePoleFilterDD::HIGHPASS_MZT);
  hp2.setCutoff(46.9);
  filter.setFeedbackHighpassCutoff(150.0);      // copied from ABL

  setSampleRate(sampleRate);
}

AciDevil::~AciDevil()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AciDevil::setSampleRate(double newSampleRate)
{
  mainEnv.setSampleRate         (       newSampleRate);
  ampEnv.setSampleRate          (       newSampleRate);
  pitchSlewLimiter.setSampleRate((float)newSampleRate);
  ampDeClicker.setSampleRate(    (float)newSampleRate);
  rc1.setSampleRate(             (float)newSampleRate);
  rc2.setSampleRate(             (float)newSampleRate);
  sequencer.setSampleRate(              newSampleRate);

  hp1.setSampleRate           (  oversampling*newSampleRate);
  hp2.setSampleRate           (  oversampling*newSampleRate);
  oscillator.setSampleRate    (  oversampling*newSampleRate);
  subOscillator.setSampleRate (  oversampling*newSampleRate);
  filter.setSampleRate        (  oversampling*newSampleRate);
}

void AciDevil::setPitchBend(double newPitchBend)
{
  pitchWheelFactor = RAPT::rsPitchOffsetToFreqFactor(pitchWheelRange*newPitchBend);
}

void AciDevil::setMasterLevel(double newLevel)
{
  level     = newLevel;
  ampScaler = RAPT::rsDbToAmp(level);
}

void AciDevil::setAccent(double newAccent)
{
  accent = 0.01 * newAccent;
}

void AciDevil::setSlideTime(double newSlideTime)
{
  if( newSlideTime >= 0.0 )
  {
    slideTime = newSlideTime;
    pitchSlewLimiter.setTimeConstant((float)(0.2*slideTime));  // \todo: tweak the scaling constant
  }
}
void AciDevil::setCutoff(double newCutoff)
{
  cutoff = newCutoff;
  calculateEnvModScalerAndOffset();
}

void AciDevil::setEnvMod(double newEnvMod)
{
  envMod           = newEnvMod;

  /*
  double upRatio   = pitchOffsetToFreqFactor(      envUpFraction *envMod);
  double downRatio = pitchOffsetToFreqFactor(-(1.0-envUpFraction)*envMod);
  envScaler        = upRatio - downRatio;
  if( envScaler != 0.0 ) // avoid division by zero
    envOffset = - (downRatio - 1.0) / (upRatio - downRatio);
  else
    envOffset = 0.0;
   */

  calculateEnvModScalerAndOffset();
}

//------------------------------------------------------------------------------------------------------------
// others:

void AciDevil::noteOn(int noteNumber, int velocity, double /*detune*/)
{
  if( sequencer.modeWasChanged() )
    allNotesOff();

  if( sequencer.getSequencerMode() != AcidSequencer::OFF )
  {
    if( velocity == 0 )
    {
      sequencer.stop();
      releaseNote(currentNote);
      currentNote = -1;
      currentVel  = 0;
    }
    else
    {
      sequencer.start();
      noteOffCountDown = INT_MAX;
      slideToNextNote  = false;
      currentNote      = noteNumber;
      currentVel       = velocity;
    }
    idle = false;
    return;
  }

  if( velocity == 0 ) // velocity zero indicates note-off events
  {
    MidiNoteEvent releasedNote(noteNumber, 0);
    noteList.remove(releasedNote);
    if( noteList.empty() )
    {
      currentNote = -1;
      currentVel  = 0;
    }
    else
    {
      currentNote = noteList.front().getKey();
      currentVel  = noteList.front().getVelocity();
    }
    releaseNote(noteNumber);
  }
  else // velocity was not zero, so this is an actual note-on
  {
    // check if the note-list is empty (indicating that currently no note is playing) - if so,
    // trigger a new note, otherwise, slide to the new note:
    if( noteList.empty() )
      triggerNote(noteNumber, velocity >= 100);
    else
      slideToNote(noteNumber, velocity >= 100);

    currentNote = noteNumber;
    currentVel  = 64;

    // and we need to add the new note to our list, of course:
    MidiNoteEvent newNote(noteNumber, velocity);
    noteList.push_front(newNote);
  }
  idle = false;
}

void AciDevil::allNotesOff()
{
  noteList.clear();
  ampEnv.noteOff();
  currentNote = -1;
  currentVel  = 0;
}

void AciDevil::triggerNote(int noteNumber, bool hasAccent)
{
  // retrigger osc and reset filter buffers only if amplitude is near zero (to avoid clicks):
  if( ampEnv.endIsReached() )
  {
    oscillator.resetPhase();
    filter.reset();
  }

  if( hasAccent )
  {
    accentGain = accent;
    setMainEnvDecay(accentDecay);
    ampEnv.setRelease(accentAmpRelease);
  }
  else
  {
    accentGain = 0.0;
    setMainEnvDecay(normalDecay);
    ampEnv.setRelease(normalAmpRelease);
  }

  oscFreq = RAPT::rsPitchToFreq(noteNumber);
  pitchSlewLimiter.setState(oscFreq);
  mainEnv.trigger();
  ampEnv.noteOn(true, noteNumber, 64);
  idle = false;
}

void AciDevil::slideToNote(int noteNumber, bool hasAccent)
{
  oscFreq = RAPT::rsPitchToFreq(noteNumber);

  if( hasAccent )
  {
    accentGain = accent;
    setMainEnvDecay(accentDecay);
    ampEnv.setRelease(accentAmpRelease);
  }
  else
  {
    accentGain = 0.0;
    setMainEnvDecay(normalDecay);
    ampEnv.setRelease(normalAmpRelease);
  }
  idle = false;
}

void AciDevil::releaseNote(int /*noteNumber*/)
{
  // check if the note-list is empty now. if so, trigger a release, otherwise slide to the note
  // at the beginning of the list (this is the most recent one which is still in the list). this
  // initiates a slide back to the most recent note that is still being held:
  if( noteList.empty() )
  {
    //filterEnvelope.noteOff();
    ampEnv.noteOff();
  }
  else
  {
    // initiate slide back:
    oscFreq     = RAPT::rsPitchToFreq(currentNote);
  }
}

void AciDevil::setMainEnvDecay(double newDecay)
{
  mainEnv.setDecayTimeConstant(newDecay);
  updateNormalizer1();
  updateNormalizer2();
}

void AciDevil::calculateEnvModScalerAndOffset()
{
  bool useMeasuredMapping = true;
  if( useMeasuredMapping == true )
  {
    // define some constants that arise from the measurements:
    const double c0   = 3.138152786059267e+002;  // lowest nominal cutoff
    const double c1   = 2.394411986817546e+003;  // highest nominal cutoff
    const double oF   = 0.048292930943553;       // factor in line equation for offset
    const double oC   = 0.294391201442418;       // constant in line equation for offset
    const double sLoF = 3.773996325111173;       // factor in line eq. for scaler at low cutoff
    const double sLoC = 0.736965594166206;       // constant in line eq. for scaler at low cutoff
    const double sHiF = 4.194548788411135;       // factor in line eq. for scaler at high cutoff
    const double sHiC = 0.864344900642434;       // constant in line eq. for scaler at high cutoff

    // do the calculation of the scaler and offset:
    double e   = RAPT::rsLinToLin(envMod, 0.0, 80.0, 0.0, 1.0);
    double c   = RAPT::rsExpToLin(cutoff, c0,  c1,   0.0, 1.0);
    double sLo = sLoF*e + sLoC;
    double sHi = sHiF*e + sHiC;
    envScaler  = (1-c)*sLo + c*sHi;
    envOffset  =  oF*c + oC;
  }
  else
  {
    double upRatio   = RAPT::rsPitchOffsetToFreqFactor(      envUpFraction *envMod);
    double downRatio = RAPT::rsPitchOffsetToFreqFactor(-(1.0-envUpFraction)*envMod);
    envScaler        = upRatio - downRatio;
    if( envScaler != 0.0 ) // avoid division by zero
      envOffset = - (downRatio - 1.0) / (upRatio - downRatio);
    else
      envOffset = 0.0;
  }
}

void AciDevil::updateNormalizer1()
{
  n1 = LeakyIntegrator::getNormalizer(mainEnv.getDecayTimeConstant(), rc1.getTimeConstant(),
    sampleRate);
  n1 = 1.0; // test
}

void AciDevil::updateNormalizer2()
{
  n2 = LeakyIntegrator::getNormalizer(mainEnv.getDecayTimeConstant(), rc2.getTimeConstant(),
    sampleRate);
  n2 = 1.0; // test
}



/*


Ideas:

-Let each note have its own length (or length-scaler)
-Let each note have an adjustable amount of accent rather than just a binary on/off switch
-Let each note define its own scaler for the slide time
-Provide more waveforms. maybe morphable ones, pulse-width control, etc.
-Let the filter have a ResoByCutoff parameter. At the moment, the resonance at high cutoffs is a
 bit too low - we need to turn it up as function of cutoff a little bit. Maybe get a TD3 and make
 measurements to mimic it behavior.
-Let the filter have nonlinearities that make it growl, scream, etc.


*/