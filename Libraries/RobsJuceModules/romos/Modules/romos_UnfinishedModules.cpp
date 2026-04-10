
//-------------------------------------------------------------------------------------------------

void PhasorPitchDithered::initialize()
{
  initInputPins({ "Freq", "Min", "Max" });
  initOutputPins({ "" });
  inputPins[2].setDefaultValue(1); // Max is 1 by default
}
INLINE void PhasorPitchDithered::process(Module* module, double* in1, double* in2, double* in3, 
  double* out, int voiceIndex)
{
  //PhasorPitchDithered *phasor = static_cast<PhasorPitchDithered*> (module);
  //*out = *in2 + (*in3 - *in2) * phasor->phases[voiceIndex];  // generate output signal
  //updatePhase(phasor, *in1, voiceIndex);
}

INLINE void PhasorPitchDithered::updatePhase(PhasorPitchDithered* phasor, double freq, 
  int voiceIndex)
{
  //// Increment and wraparound:
  //phasor->phases[voiceIndex] += freq * processingStatus.getSystemSamplePeriod();
  //while(phasor->phases[voiceIndex] >= 1.0)
  //  phasor->phases[voiceIndex] -= 1.0;
  //while(phasor->phases[voiceIndex] <  0.0)
  //  phasor->phases[voiceIndex] += 1.0;
}
void PhasorPitchDithered::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  phases[voiceIndex] = 0.0;  // Introduce a startphase later (as GUI parameter)
}
void PhasorPitchDithered::allocateMemory()
{
  AtomicModule::allocateMemory();
  phases = new double[getNumVoices()];
}
void PhasorPitchDithered::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] phases;
  phases = nullptr;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(PhasorPitchDithered);

