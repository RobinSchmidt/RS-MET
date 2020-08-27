//#include "romos_ModulationModules.h"
//#include "../Algorithms/romos_Interpolation.h"

//-------------------------------------------------------------------------------------------------

void EnvelopeADSR::initialize()
{
  //initInputPins(7, "Att", "AtSh", "Dec", "DcSh", "Sus", "Rel", "RlSh");
  //initOutputPins(1, "Env");

  initInputPins({ "Att", "AtSh", "Dec", "DcSh", "Sus", "Rel", "RlSh" });
  //initOutputPins({ "Env" });
  initOutputPins({ "" });

  // todo: maybe have a trigger input "Trig" - how does it actually currently get triggered?
  // then, we probably also need to have a "Gate" input in order

  // or: maybe have both MidiADSR and GateADSR maybe with a common baseclass
  // or: have only a gate input and remember here the previous value of the gate and trigger when
  // it was 0 before and is 1 now


  inputPins[4].setDefaultValue(1.0); 
  // wrap this into a function setPinDefaultValue(int pinIndex, double newDefaultValue) or
  // better setPinDefaultValue(std::string pinName, double value)
}

void EnvelopeADSR::processWithoutTriggerFlagCheck(Module *module, const double *Att,
  const double *AtSh, const double *Dec, const double *DcSh, const double *Sus, const double *Rel,
  const double *RlSh, double *out, const int voiceIndex)
{
  EnvelopeADSR  *env     = static_cast<EnvelopeADSR*> (module);
  unsigned long *counter = &(env->counters[voiceIndex]);
  double         yStart  = env->startValues[voiceIndex];

  double sustain   = *Sus;
  double timeScale = 1.0;   // todo: make input
  double scl = timeScale * processingStatus.getSystemSampleRate();
  unsigned long attackSamples  = (unsigned long) ::round(*Att * scl);
  unsigned long decaySamples   = (unsigned long) ::round(*Dec * scl);
  unsigned long releaseSamples = (unsigned long) ::round(*Rel * scl);

  typedef RAPT::rsNodeBasedFunction<double> NBF;
  double a;
  double p;  // proportion of passedLength/fullLength of the phase ( = I(x), page 184, Eq.3-30)
  if(*counter < attackSamples)
  {
    // attack phase
    //a = *AtSh;  // apply formula
    a = NBF::linVsExpFormulaScaler(*AtSh);
    p = (double)*counter / (double)attackSamples;
    if(a == 0.0)
      *out = yStart + (1.0 - yStart) * p;
    else
      *out = yStart + (1.0 - yStart) * (1.0 - exp(p*a)) / (1.0 - exp(a)); // encapsulate that function
    *counter += 1;
  }
  else if(*counter < attackSamples + decaySamples)
  {
    // decay phase
    //a = *DcSh;
    a = NBF::linVsExpFormulaScaler(*DcSh);
    p = (double)(*counter-attackSamples)  / (double)(decaySamples);
    if(a == 0.0)
      *out = 1.0 + (sustain-1.0) * p;
    else
      *out = 1.0 + (sustain-1.0) * (1.0 - exp(p*a)) / (1.0 - exp(a));
    *counter += 1;
  }
  else if(voiceAllocator.isNoteOn(voiceIndex) && *counter == attackSamples + decaySamples)
  {
    // the second part of the conditional is needed in order to not jump up to the sustain level when a new attack is started while
    // the envelope was in its release phase

    // sustain phase
    *out = sustain;
  }
  else if(*counter < attackSamples + decaySamples + releaseSamples)
  {
    // release phase
    //a = *RlSh;
    a = NBF::linVsExpFormulaScaler(*RlSh);
    p = (double)(*counter-attackSamples-decaySamples)  / (double)(releaseSamples);
    if(a == 0.0)
      *out = yStart - yStart * p;
    else
      *out = yStart - yStart * (1.0 - exp(p*a)) / (1.0 - exp(a));
    *counter += 1;
  }
  else
  {
    // end is reached
    *out = 0.0;
  }
}
INLINE void EnvelopeADSR::process(Module *module, double *in1, double *in2, double *in3,
  double *in4, double *in5, double *in6, double *in7, double *out, int voiceIndex)
{
  EnvelopeADSR  *env     = static_cast<EnvelopeADSR*> (module);

  // in case a note-on/off was triggered at this sample, we compute a preliminary value and use this value as value to start from - this
  // makes sure to correctly handle double-attacks and early releases
  if(voiceAllocator.getNoteOffTriggerFlag(voiceIndex) == true)
  {
    // we check note-off flag first to ensure that it attacks again when note-on and -off occur simultanously

    double timeScale = 1.0;
    unsigned long attackSamples  = (unsigned long) ::round(*in1 * timeScale * processingStatus.getSystemSampleRate());
    unsigned long decaySamples   = (unsigned long) ::round(*in3 * timeScale * processingStatus.getSystemSampleRate());
    if(env->counters[voiceIndex] >= attackSamples)
      env->startValues[voiceIndex] = *in5; // set it temporarily to the sustain value

    processWithoutTriggerFlagCheck(module, in1, in2, in3, in4, in5, in6, in7, out, voiceIndex);
    env->startValues[voiceIndex] = *out;
    env->counters[voiceIndex]    = attackSamples + decaySamples;
  }
  if(voiceAllocator.getNoteOnTriggerFlag(voiceIndex) == true)
  {
    processWithoutTriggerFlagCheck(module, in1, in2, in3, in4, in5, in6, in7, out, voiceIndex);
    env->startValues[voiceIndex] = *out;
    env->counters[voiceIndex]    = 0;
  }

  processWithoutTriggerFlagCheck(module, in1, in2, in3, in4, in5, in6, in7, out, voiceIndex);
}
void EnvelopeADSR::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  counters[voiceIndex]    = 0;
  startValues[voiceIndex] = 0.0;
}
void EnvelopeADSR::allocateMemory()
{
  AtomicModule::allocateMemory();
  counters    = new unsigned long[getNumVoices()];
  startValues = new double[getNumVoices()];
}
void EnvelopeADSR::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] counters;     counters    = nullptr; // use deleteAndSetNull function
  delete[] startValues;  startValues = nullptr;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_7(EnvelopeADSR);

//-------------------------------------------------------------------------------------------------
