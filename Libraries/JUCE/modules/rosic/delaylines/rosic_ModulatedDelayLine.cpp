//#include "rosic_ModulatedDelayLine.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulatedDelayLine::ModulatedDelayLine(int maximumDelayInSamples)
: FractionalDelayLine(maximumDelayInSamples)
{
  delayModDepth = 0.0;
  ampModDepth   = 0.0;
  feedforward   = 1.0;
  feedback      = 0.0;
  blend         = 0.0;
  pan           = 0.0;
  gL            = sqrt(0.5);
  gR            = sqrt(0.5);
  g             = 1.0;
  delayModCycle = 0.5;
  ampModCycle   = 0.5;
  delayModSync  = false;  
  ampModSync    = false;  

  interpolatorL.setInterpolationMethod(Interpolator::WARPED_ALLPASS);
  interpolatorR.setInterpolationMethod(Interpolator::WARPED_ALLPASS);

  interpolatorL.setInterpolationMethod(Interpolator::LINEAR);  // for debug
  interpolatorR.setInterpolationMethod(Interpolator::LINEAR);

  delayOscL.setFrequency(7.0);
  delayOscR.setFrequency(7.0);
  delayOscR.setStartPhase(-0.5*PI);

  ampOscL.setFrequency(4.0);
  ampOscR.setFrequency(4.0);
  ampOscR.setStartPhase(-0.5*PI);

  triggerOscillators();

  setDelayTime(0.02);     // 20 ms


  /*
  // for test:
  md  = 0.0;
  ma  = 1.0;
  ff  = 1.0;
  fb  = 0.9;
  bl  = 1.0;
  pan = 0.0;
  gL  = 1.0;
  gR  = 1.0;
  g   = 1.0;
  */
}

ModulatedDelayLine::~ModulatedDelayLine()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void ModulatedDelayLine::setTempoInBPM(double newTempoInBPM)
{
  FractionalDelayLine::setTempoInBPM(newTempoInBPM);
  setupDelayModulationFrequency();
  setupAmplitudeModulationFrequency();
}

void ModulatedDelayLine::setSampleRate(double newSampleRate)
{
  FractionalDelayLine::setSampleRate(newSampleRate);
  equalizer.setSampleRate(sampleRate);
  delayOscL.setSampleRate(sampleRate);
  delayOscR.setSampleRate(sampleRate);
  ampOscL.setSampleRate(sampleRate);
  ampOscR.setSampleRate(sampleRate);
}

void ModulatedDelayLine::setDelayModulationCycleLength(double newCycleLength)
{
  if( newCycleLength > 0.0 )
  {
    delayModCycle = newCycleLength;
    setupDelayModulationFrequency();
  }
}

void ModulatedDelayLine::setDelayModulationDepth(double newDepth)
{
  delayModDepth = clip(newDepth, 0.0, 1.0);
}

void ModulatedDelayLine::setDelayModulationPhaseLeft(double newPhase)
{
  delayOscL.setStartPhase(RAPT::rsDegreeToRadiant(newPhase));
}

void ModulatedDelayLine::setDelayModulationPhaseRight(double newPhase)
{
  delayOscR.setStartPhase(RAPT::rsDegreeToRadiant(newPhase));
}

void ModulatedDelayLine::setDelayModulationSyncMode(bool shouldBeSynced)
{
  delayModSync = shouldBeSynced;
  setupDelayModulationFrequency();
}

void ModulatedDelayLine::setAmplitudeModulationCycleLength(double newCycleLength)
{
  if( newCycleLength > 0.0 )
  {
    ampModCycle = newCycleLength;
    setupAmplitudeModulationFrequency();
  }
}

void ModulatedDelayLine::setAmplitudeModulationDepth(double newDepth)
{
  ampModDepth = clip(newDepth, 0.0, 1.0);
}

void ModulatedDelayLine::setAmplitudeModulationPhaseLeft(double newPhase)
{
  ampOscL.setStartPhase(RAPT::rsDegreeToRadiant(newPhase));
}

void ModulatedDelayLine::setAmplitudeModulationPhaseRight(double newPhase)
{
  ampOscR.setStartPhase(RAPT::rsDegreeToRadiant(newPhase));
}

void ModulatedDelayLine::setAmplitudeModulationSyncMode(bool shouldBeSynced)
{
  ampModSync = shouldBeSynced;
  setupAmplitudeModulationFrequency();
}

//-------------------------------------------------------------------------------------------------
// others:

void ModulatedDelayLine::resetBuffers()
{
  FractionalDelayLine::clearDelayBuffer();
  equalizer.reset();
}

void ModulatedDelayLine::triggerDelayModOscillators()
{
  delayOscL.trigger();
  delayOscR.trigger();
}

void ModulatedDelayLine::triggerAmpModOscillators()
{
  ampOscL.trigger();
  ampOscR.trigger();
}

void ModulatedDelayLine::triggerOscillators()
{
  delayOscL.trigger();
  delayOscR.trigger();
  ampOscL.trigger();
  ampOscR.trigger();
}

void ModulatedDelayLine::setupDelayModulationFrequency()
{
  double cycleLengthInSeconds = delayModCycle;
  if( delayModSync == true )
    cycleLengthInSeconds = RAPT::rsBeatsToSeconds(delayModCycle, bpm);
  delayOscL.setFrequency(1.0/cycleLengthInSeconds);
  delayOscR.setFrequency(1.0/cycleLengthInSeconds);
}

void ModulatedDelayLine::setupAmplitudeModulationFrequency()
{
  double cycleLengthInSeconds = ampModCycle;
  if( ampModSync == true )
    cycleLengthInSeconds = RAPT::rsBeatsToSeconds(ampModCycle, bpm);
  ampOscL.setFrequency(1.0/cycleLengthInSeconds);
  ampOscR.setFrequency(1.0/cycleLengthInSeconds);
}