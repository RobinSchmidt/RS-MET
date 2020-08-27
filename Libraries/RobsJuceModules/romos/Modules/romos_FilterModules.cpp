
//-------------------------------------------------------------------------------------------------

void FirstOrderLowpass::initialize()
{
  initInputPins({ "In", "Cutoff" });
  initOutputPins({ "Out" });
}
INLINE void FirstOrderLowpass::process(Module *module, double *in1, double *in2, double *out,
  int voiceIndex)
{
  FirstOrderLowpass *filter = static_cast<FirstOrderLowpass*> (module);
  double *buf               = filter->buffers + 2*voiceIndex;

  // coefficient computation (via impulse invariant transform):
  double frequency = *in2;
  double omega     = frequency * processingStatus.getFreqToOmegaFactor();
  double x         = exp(-omega);
  double a1 = x;
  double b0 = 1.0-x;
  double b1 = 0.0;

  // filtering:
  *out   = b0 * *in1 + b1 * buf[0] + a1 * buf[1]; // b0*x[n] + b1*x[n-1] + a1*y[n-1]
  buf[0] = *in1;
  buf[1] = *out;
  // optimize away b1, ..also b0, a1
}
void FirstOrderLowpass::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  double *buf = buffers + 2*voiceIndex;
  buf[0] = buf[1] = 0.0;
}
void FirstOrderLowpass::allocateMemory()
{
  AtomicModule::allocateMemory();
  buffers = new double[2*getNumVoices()]; // x[n-1], y[n-1]
}
void FirstOrderLowpass::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] buffers;
  buffers = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(FirstOrderLowpass);

//-------------------------------------------------------------------------------------------------

void FirstOrderFilter::initialize()
{
  initInputPins({ "In", "b0", "b1", "a1" });
  initOutputPins({ "Out" });
}
INLINE void FirstOrderFilter::process(Module *module, double *in1, double *in2, double *in3, 
  double *in4, double *out, int voiceIndex)
{
  FirstOrderFilter *filter = static_cast<FirstOrderFilter*> (module);
  double *buf              = filter->buffers + 2*voiceIndex;
  *out   = *in2 * *in1 + *in3 * buf[0] - *in4 * buf[1]; // b0*x[n] + b1*x[n-1] - a1*y[n-1]
  buf[0] = *in1;
  buf[1] = *out;
}
void FirstOrderFilter::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  double *buf = buffers + 2*voiceIndex;
  buf[0] = buf[1] = 0.0;
}
void FirstOrderFilter::allocateMemory()
{
  AtomicModule::allocateMemory();
  buffers = new double[2*getNumVoices()]; // x[n-1], y[n-1]
}
void FirstOrderFilter::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] buffers;
  buffers = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_4(FirstOrderFilter);

//-------------------------------------------------------------------------------------------------

void Biquad::initialize()
{
  initInputPins({ "In", "b0", "b1", "b2", "a1", "a2" });
  initOutputPins({ "Out" });
}
INLINE void Biquad::process(Module *module, double *in1, double *in2, double *in3, double *in4, 
  double *in5, double *in6, double *out, int voiceIndex) // rename inputs
{
  Biquad *biquad = static_cast<Biquad*> (module);
  double *buf    = biquad->buffers + 4*voiceIndex;
  *out    = (*in2 * *in1 + TINY) + (*in3 * buf[0] + *in4 * buf[1]) 
    - (*in5 * buf[2] + *in6 * buf[3]); // get rid of +TINY 
  buf[3]  = buf[2];
  buf[2]  = *out;
  buf[1]  = buf[0];
  buf[0]  = *in1;
}
void Biquad::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  double *buf = buffers + 4*voiceIndex;
  buf[0] = buf[1] = buf[2] = buf[3] = 0.0;
}
void Biquad::allocateMemory()
{
  AtomicModule::allocateMemory();
  buffers = new double[4*getNumVoices()]; // x[n-1], x[n-2], y[n-1], y[n-2]
}
void Biquad::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] buffers;
  buffers = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_6(Biquad);

//-------------------------------------------------------------------------------------------------

void BiquadDesigner::initialize()
{
  initInputPins({ "Freq", "Q", "Gain" });
  initOutputPins({ "b0", "b1", "b2", "a1", "a2" });
  addParameter(rosic::rsString("Mode"), "Bypass");
  //parameterChanged(0); // nah - it calls resetStateForAllVoices() - which is not allowed before 
  // memory is allocated by the factory
  mode = BYPASS;
}
INLINE void BiquadDesigner::process(Module *module, double *in1, double *in2, double *in3, 
  double *outs, int voiceIndex)
{
  BiquadDesigner *designer = static_cast<BiquadDesigner*> (module);
  double *oldParams = designer->oldParameters + 3*voiceIndex;
  double *oldOuts   = designer->oldOutputs    + 5*voiceIndex;

  // maybe use struct biquadParamsFQG (freq/Q/gain), struct biquadCoeffs, struct biquadState

  // recalculate coefficients if necessary (factor out):
  if(*in1 != oldParams[0] || *in2 != oldParams[1] || *in3 != oldParams[2])
  {
    switch(designer->mode)
    {
    case LOWPASS_6_BILINEAR:            romos::biquadLowpassCoeffsBilinear1(outs, *in1); break;
    case HIGHPASS_6_BILINEAR:           romos::biquadHighpassCoeffsBilinear1(outs, *in1); break;
    case LOW_SHELF_1_BILINEAR:          romos::biquadLowShelfCoeffsBilinear1(outs, *in1, *in3); break;
    case HIGH_SHELF_1_BILINEAR:         romos::biquadHighShelfCoeffsBilinear1(outs, *in1, *in3); break;
    case ALLPASS_1_BILINEAR:            romos::biquadAllpassCoeffsBilinear1(outs, *in1); break;
    case LOWPASS_12_BILINEAR:           romos::biquadLowpassCoeffsBilinear2(outs, *in1, *in2); break;
    case HIGHPASS_12_BILINEAR:          romos::biquadHighpassCoeffsBilinear2(outs, *in1, *in2); break;
    case BANDPASS_CONST_SKIRT_BILINEAR: romos::biquadBandpassConstSkirtCoeffs(outs, *in1, *in2); break;
    case BANDPASS_CONST_PEAK_BILINEAR:  romos::biquadBandpassConstPeakCoeffs(outs, *in1, *in2); break;
    case BANDREJECT_BILINEAR:           romos::biquadBandrejectCoeffs(outs, *in1, *in2); break;
    case PEAK_BILINEAR:                 romos::biquadPeakCoeffs(outs, *in1, *in2, *in3); break;
    case LOW_SHELF_2_BILINEAR:          romos::biquadLowShelfCoeffsBilinear2(outs, *in1, *in2, *in3); break;
    case HIGH_SHELF_2_BILINEAR:         romos::biquadHighShelfCoeffsBilinear2(outs, *in1, *in2, *in3); break;
    case ALLPASS_2_BILINEAR:            romos::biquadAllpassCoeffsBilinear2(outs, *in1, *in2); break;
    default:                            romos::biquadBypassCoeffs(outs);
    }
    oldParams[0] = *in1;
    oldParams[1] = *in2;
    oldParams[2] = *in3;
  }
  else
  {
    // otherwise use same outputs as before:
    memcpy(outs, oldOuts, 5*sizeof(double));  // try direct asssignment - check which is faster
  }

  // remember outputs, in case we hit the else-branch next time:
  memcpy(oldOuts, outs, 5*sizeof(double));
}
void BiquadDesigner::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  double *oldParams = oldParameters + 3*voiceIndex;
  double *oldOuts   = oldOutputs    + 5*voiceIndex;
  oldParams[0] = oldParams[1] = oldParams[2] = 0.0; // div by zero? -> we should tolerate this
  oldOuts[0] = oldOuts[1] = oldOuts[2] = oldOuts[3] = oldOuts[4] = 0.0;
}
void BiquadDesigner::parameterChanged(int index)
{
  rosic::rsString m = parameters[0].value;

  if(     m == "Lowpass, 6 dB/oct, BLT")       mode = LOWPASS_6_BILINEAR;
  else if(m == "Highpass, 6 dB/oct, BLT")      mode = HIGHPASS_6_BILINEAR;
  else if(m == "Low Shelf, 1st order, BLT")    mode = LOW_SHELF_1_BILINEAR;
  else if(m == "High Shelf, 1st order, BLT")   mode = HIGH_SHELF_1_BILINEAR;
  else if(m == "Allpass, 1st order, BLT")      mode = ALLPASS_1_BILINEAR;
  else if(m == "Lowpass, 12 dB/oct, BLT")      mode = LOWPASS_12_BILINEAR;
  else if(m == "Highpass, 12 dB/oct, BLT")     mode = HIGHPASS_12_BILINEAR;
  else if(m == "Bandpass, const. skirt, BLT")  mode = BANDPASS_CONST_SKIRT_BILINEAR;
  else if(m == "Bandpass, const. peak, BLT")   mode = BANDPASS_CONST_PEAK_BILINEAR;
  else if(m == "Bandreject, BLT")              mode = BANDREJECT_BILINEAR;
  else if(m == "Peak, BLT")                    mode = PEAK_BILINEAR;
  else if(m == "Low Shelf, 2nd order, BLT")    mode = LOW_SHELF_2_BILINEAR;
  else if(m == "High Shelf, 2nd order, BLT")   mode = HIGH_SHELF_2_BILINEAR;
  else if(m == "Allpass, 2nd order, BLT")      mode = ALLPASS_2_BILINEAR;
  else                                         mode = BYPASS;
  // maybe use std::map or rosic::KeyValueMap to do the back-and-forth

  // other types to come:  impulse-invariant designs, prescribed nyquist-gain, resonator/FOF, 
  // 2nd order oscillator, Lowpass/Highpass chain, Low/High Shelf chain

  resetStateForAllVoices(); // enforces coefficient re-calculation in next call to process
}
void BiquadDesigner::allocateMemory()
{
  AtomicModule::allocateMemory();
  oldParameters = new double[3*getNumVoices()]; // frequency, gain, Q
  oldOutputs    = new double[5*getNumVoices()]; // b0, b1, b2, a1, a2
}
void BiquadDesigner::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] oldParameters;
  oldParameters = NULL;
  delete[] oldOutputs;
  oldOutputs = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(BiquadDesigner);

//-------------------------------------------------------------------------------------------------

void LadderFilter::initialize()
{
  initInputPins({ "In", "Freq", "Reso", "AutoGain" });
  initOutputPins({ "Out" });
  addParameter("Mode", "Lowpass, 24 dB/oct");
  addParameter("SaturationMode", "No Saturation");
  filterMode     = LP_24;
  saturationMode = NO_SATURATION;

  // maybe include a GUI parameter to select the saturation function (tanh, asinh, atan, hardclip, 
  // softclip, etc.)
}
INLINE void LadderFilter::process(Module *module, double *in1, double *in2, double *in3, 
  double *in4, double *out, int voiceIndex)
{
  LadderFilter *filter = static_cast<LadderFilter*> (module);
  double *c         = filter->coeffs        + 7*voiceIndex;
  double *y         = filter->outputs       + 5*voiceIndex;
  double *oldParams = filter->oldParameters + 3*voiceIndex;

  // recalculate coefficients if necessary:
  if(*in2 != oldParams[0] || *in3 != oldParams[1] || *in4 != oldParams[2])
  {
    ladderCoeffs(c, filter->filterMode, *in2, *in3, *in4);
    oldParams[0] = *in2;
    oldParams[1] = *in3;
    oldParams[2] = *in4;
  }

  // maybe use a function pointer later that is assigned in parameterChanged - check which 
  // version is more efficient
  switch(filter->saturationMode)
  {
  case NO_SATURATION:
  {
    y[0] = *in1 - c[1]*y[4];
    y[1] = y[0] + c[0]*(y[0]-y[1]);
    y[2] = y[1] + c[0]*(y[1]-y[2]);
    y[3] = y[2] + c[0]*(y[2]-y[3]);
    y[4] = y[3] + c[0]*(y[3]-y[4]);
  };
  break;
  case LAST_STAGE:
  {
    y[0] =            *in1 - c[1]*y[4];
    y[1] =            y[0] + c[0]*(y[0]-y[1]);
    y[2] =            y[1] + c[0]*(y[1]-y[2]);
    y[3] =            y[2] + c[0]*(y[2]-y[3]);
    y[4] = RAPT::rsTanhApprox(y[3] + c[0]*(y[3]-y[4]));
  };
  break;
  case FEEDBACK:
  {
    y[0] = *in1 - RAPT::rsTanhApprox(c[1]*y[4]);
    y[1] = y[0] + c[0]*(y[0]-y[1]);
    y[2] = y[1] + c[0]*(y[1]-y[2]);
    y[3] = y[2] + c[0]*(y[2]-y[3]);
    y[4] = y[3] + c[0]*(y[3]-y[4]);
  };
  break;
  case EACH_STAGE:
  {
    y[0] =            *in1 - c[1]*y[4];
    y[1] = RAPT::rsTanhApprox(y[0] + c[0]*(y[0]-y[1]));
    y[2] = RAPT::rsTanhApprox(y[1] + c[0]*(y[1]-y[2]));
    y[3] = RAPT::rsTanhApprox(y[2] + c[0]*(y[2]-y[3]));
    y[4] = RAPT::rsTanhApprox(y[3] + c[0]*(y[3]-y[4]));
  };
  break;
  }

  *out = c[2]*y[0] + c[3]*y[1] + c[4]*y[2] + c[5]*y[3] + c[6]*y[4];
}
void LadderFilter::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  double *ptr = outputs + 5*voiceIndex;
  ptr[0] = ptr[1] = ptr[2] = ptr[3] = ptr[4] = 0.0;
  ptr = coeffs + 7*voiceIndex;
  ptr[0] = ptr[1] = ptr[2] = ptr[3] = ptr[4] = ptr[5] = ptr[6] = 0.0;
  ptr = oldParameters + 3*voiceIndex;
  ptr[0] = ptr[1] = ptr[2] = 0.0;
}
void LadderFilter::parameterChanged(int index)
{
  rosic::rsString m = parameters[0].value;
  if(     m == "Lowpass, 6 dB/oct")      filterMode = LP_6;
  else if(m == "Lowpass, 12 dB/oct")     filterMode = LP_12;
  else if(m == "Lowpass, 18 dB/oct")     filterMode = LP_18;
  else if(m == "Lowpass, 24 dB/oct")     filterMode = LP_24;
  else if(m == "Highpass, 6 dB/oct")     filterMode = HP_6;
  else if(m == "Highpass, 12 dB/oct")    filterMode = HP_12;
  else if(m == "Highpass, 18 dB/oct")    filterMode = HP_18;
  else if(m == "Highpass, 24 dB/oct")    filterMode = HP_24;
  else if(m == "Bandpass, 6/6 dB/oct")   filterMode = BP_6_6;
  else if(m == "Bandpass, 6/12 dB/oct")  filterMode = BP_6_12;
  else if(m == "Bandpass, 12/6 dB/oct")  filterMode = BP_12_6;
  else if(m == "Bandpass, 12/12 dB/oct") filterMode = BP_12_12;
  else if(m == "Bandpass, 6/18 dB/oct")  filterMode = BP_6_18;
  else if(m == "Bandpass, 18/6 dB/oct")  filterMode = BP_18_6;
  else                                   filterMode = FLAT;

  m = parameters[1].value;
  if(m == "No Saturation")  saturationMode = NO_SATURATION;
  else if(m == "Last Stage")  saturationMode = LAST_STAGE;
  else if(m == "Feedback")  saturationMode = FEEDBACK;
  else if(m == "Each Stage")  saturationMode = EACH_STAGE;

  resetStateForAllVoices(); // enforces coefficient re-calculation in next call to process
}
void LadderFilter::allocateMemory()
{
  AtomicModule::allocateMemory();
  outputs       = new double[5*getNumVoices()];  // y0[n], y1[n], y2[n], y3[n], y4[n]
  coeffs        = new double[7*getNumVoices()];  // a1, k, c0, c1, c2, c3, c4
  oldParameters = new double[3*getNumVoices()];  // frequency, resonance, autogain
}
void LadderFilter::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] outputs;
  outputs = NULL;
  delete[] coeffs;
  coeffs = NULL;
  delete[] oldParameters;
  oldParameters = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_4(LadderFilter);
