
//-------------------------------------------------------------------------------------------------

void PhasorPitchDithered::initialize()
{
  initInputPins({ "Freq" });
  initOutputPins({ "" });                     // Maybe call it Phase or 0..1?
  //inputPins[0].setDefaultValue(1000);         // 1 kHz. Not sure. Maybe 440 is better? Or none?

  // Maybe we should follow the convention that when a module has just a single output, this output
  // remains unnamed such that the rendering of the box on the GUI can be smaller.
}
INLINE void PhasorPitchDithered::process(Module* module, double* in1, double* out, int voiceIndex)
{
  PhasorPitchDithered *phasor = static_cast<PhasorPitchDithered*> (module);
  double newFreq = *in1;
  if(newFreq != phasor->freqs[voiceIndex])
  {
    double sampleRate  = processingStatus.getSystemSampleRate();
    double cycleLength = sampleRate / newFreq;
    phasor->oscs[voiceIndex].setMeanCycleLength(cycleLength, phasor->rangeClosed);
    phasor->freqs[voiceIndex] = newFreq;
  }
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
  oscs  = new RAPT::rsPitchDitherOsc<double>[getNumVoices()];
  freqs = new double[getNumVoices()];
}
void PhasorPitchDithered::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] oscs;  oscs  = nullptr;
  delete[] freqs; freqs = nullptr;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(PhasorPitchDithered);

// ToDo:
//
// - It currently doesn't make use of the "Freq" input. I think, we should maintain a "freq" member
//   and in process, compare the freq member against the input and when they differ, we need to
//   update the freq member and call setMeanCycleLength() on the embedded DSP object. The purpose 
//   of maintaining the freq member is to avoid calling setMeanCycleLength() for every sample even
//   when the freq didn't change because a call to it is costly. ...done. Maybe we could have two 
//   different modes controlled by a GUI switch where we swicth between calling 
//   setMeanCycleLength() and setMeanCycleLengthNoUpdate(). The former should be used if we want to
//   do FM of fast modulations whereas otherwise we could get away with the latter which is 
//   potentially cheaper.
