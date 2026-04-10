
//-------------------------------------------------------------------------------------------------

void PhasorPitchDithered::initialize()
{
  initInputPins({ "Freq" });
  initOutputPins({ "" });                     // Maybe call it Phase or 0..1?
  //inputPins[0].setDefaultValue(1000);         // 1 kHz. Not sure. Maybe 440 is better? Or none?
}
INLINE void PhasorPitchDithered::process(Module* module, double* in1, double* out, int voiceIndex)
{
  PhasorPitchDithered *phasor = static_cast<PhasorPitchDithered*> (module);
  *out = phasor->oscs[voiceIndex].getSamplePhasor(phasor->rangeClosed);
}
void PhasorPitchDithered::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  oscs[voiceIndex].reset(rangeClosed);
}
void PhasorPitchDithered::allocateMemory()
{
  AtomicModule::allocateMemory();
  oscs = new RAPT::rsPitchDitherOsc<double>[getNumVoices()];
}
void PhasorPitchDithered::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] oscs;
  oscs = nullptr;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(PhasorPitchDithered);

