//#include "romos_SoundGeneratorModules.h"

//-------------------------------------------------------------------------------------------------

void WhiteNoise::initialize()
{
  initOutputPins({ "" });
  //addParameter(rosic::rsString("Seed"), 0.0);
  addParameter("Seed", 0.0);
  parameterChanged(0);
  // add shape, min, max
}
INLINE void WhiteNoise::process(Module *module, double *out, int voiceIndex)
{
  WhiteNoise *noiseGen = static_cast<WhiteNoise*> (module);
  unsigned long *state = noiseGen->state + voiceIndex;
  *state  = 1664525 * (*state) + 1013904223;      // mod implicitely by integer overflow
  *out = -1.0 + ((2.0/4294967296.0) * (*state));  // min + (max-min)/429.... * (*state)
}
void WhiteNoise::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  state[voiceIndex] = (unsigned long)parameters[0].value.asInt();
  //state[voiceIndex] = (unsigned long) round( *(inputPins[0].outputPointer + voiceIndex * inputPins[0].outputVoiceStride) );
}
void WhiteNoise::parameterChanged(int index)
{
  // nothing to do - we initialize the PNRG state directly from parameters[0] in resetVoiceState
}
void WhiteNoise::allocateMemory()
{
  AtomicModule::allocateMemory();
  state = new unsigned long[getNumVoices()]; // PRNG state per voice, use self defined uint32 later
}
void WhiteNoise::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] state;
  state = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(WhiteNoise);

//-------------------------------------------------------------------------------------------------

void Phasor::initialize()
{
  initInputPins({ "Freq", "Min", "Max" });
  initOutputPins({ "" });
  inputPins[2].setDefaultValue(1); // Max is 1 by default
}
INLINE void Phasor::process(Module* module, double* in1, double* in2, double* in3, double* out, 
  int voiceIndex)
{
  Phasor *phasor = static_cast<Phasor*> (module);
  *out = *in2 + (*in3 - *in2) * phasor->phases[voiceIndex];  // generate output signal
  updatePhase(phasor, *in1, voiceIndex);
}

INLINE void Phasor::updatePhase(Phasor* phasor, double freq, int voiceIndex)
{
  // increment and wraparound:
  phasor->phases[voiceIndex] += freq * processingStatus.getSystemSamplePeriod();
  while(phasor->phases[voiceIndex] >= 1.0)
    phasor->phases[voiceIndex] -= 1.0;
  while(phasor->phases[voiceIndex] <  0.0)
    phasor->phases[voiceIndex] += 1.0;
}

void Phasor::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  phases[voiceIndex] = 0.0;  // introduce a startphase later (as GUI parameter)
}
void Phasor::allocateMemory()
{
  AtomicModule::allocateMemory();
  phases = new double[getNumVoices()];
}
void Phasor::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] phases;
  phases = nullptr;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(Phasor);

//-------------------------------------------------------------------------------------------------

void SineOscillator::initialize()
{
  initInputPins({ "Freq", "Amp", "Phase" });
  initOutputPins({ "" });
  inputPins[1].setDefaultValue(1); // Amp is 1 by default
}
INLINE void SineOscillator::process(Module* module, double* Freq, double* Amp, double* Phase, 
  double* out, int voiceIndex)
{
  //SineOscillator* sinOsc = static_cast<SineOscillator*> (module);
  Phasor *phasor = static_cast<Phasor*> (module);
  *out = *Amp * sin(2*PI* (phasor->phases[voiceIndex] + *Phase));
  updatePhase(phasor, *Freq, voiceIndex);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(SineOscillator);

//-------------------------------------------------------------------------------------------------

void BandlimitedImpulseTrain::initialize()
{
  initInputPins({ "Freq", "Phase" });
  initOutputPins({ "Out" });
  fixedPhaseOffset = 0.0;
}

// define this chunk of the prcoess function as macro because it's reused in BiSawOscillator:
#define COMPUTE_LOCALS_FOR_BLIT                                                       \
  static const double maxNumHarmonics = 100000.0;                                     \
  double numHarmonics;                                                                \
  double absFreq = fabs(*in1);                                                        \
  if( absFreq < NANO )                                                                \
    numHarmonics = maxNumHarmonics;                                                   \
  else                                                                                \
    numHarmonics = floor(0.5 * processingStatus.getSystemSampleRate() / absFreq);     \
  double ampScaler;                                                                   \
  if( numHarmonics > 0 )                                                              \
    ampScaler = 0.5 * RAPT::rsSign(*in1) / numHarmonics;                              \
  else                                                                                \
    ampScaler = 0.0;                                                                  \
  double omega    = *in1 * processingStatus.getFreqToOmegaFactor();                   \
  double absOmega = absFreq * processingStatus.getFreqToOmegaFactor();                \
  double theta = blit->phases[voiceIndex] + TWO_PI * (*in2 + blit->fixedPhaseOffset); \

INLINE void BandlimitedImpulseTrain::process(Module *module, double *in1, double *in2, double *out, 
  int voiceIndex)
{
  // for performance tests:
  //*in1 = 440.0;

  BandlimitedImpulseTrain *blit = static_cast<BandlimitedImpulseTrain*> (module);
  COMPUTE_LOCALS_FOR_BLIT; // computes numHarmonics, ampScaler, omega and theta
  *out = ampScaler * computeUnscaledBlitValue(theta, numHarmonics);
  incrementPhases(blit, voiceIndex, omega);
}
INLINE void BandlimitedImpulseTrain::incrementPhases(BandlimitedImpulseTrain *blit, 
  const int voiceIndex, const double omega)
{
  blit->phases[voiceIndex] += omega;
  while(blit->phases[voiceIndex] >= TWO_PI)
    blit->phases[voiceIndex] -= TWO_PI;
  while(blit->phases[voiceIndex] <  0.0)
    blit->phases[voiceIndex] += TWO_PI;
}
INLINE double BandlimitedImpulseTrain::computeUnscaledBlitValue(const double theta, 
  const double numHarmonics)
{
  double halfTheta = 0.5 * theta;
  double numerator, denominator;

  static const double margin = 0.000001;  // must be really small, otherwise the cosine branch produces garbage
  denominator = sin(halfTheta);
  if(fabs(denominator) > margin)
    numerator = sin(halfTheta * (2.0*numHarmonics+1.0));
  else
  {
    numerator   = cos(halfTheta * (2.0*numHarmonics+1.0)) * (2.0*numHarmonics+1.0);
    denominator = cos(halfTheta);
  }

  return (numerator/denominator) - 1.0;
}
void BandlimitedImpulseTrain::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  phases[voiceIndex] = fixedPhaseOffset;  // introduce a startphase later (as GUI parameter)
}
void BandlimitedImpulseTrain::allocateMemory()
{
  AtomicModule::allocateMemory();
  phases = new double[getNumVoices()];
}
void BandlimitedImpulseTrain::freeMemory()
{
  AtomicModule::freeMemory();
  delete[] phases;
  phases = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(BandlimitedImpulseTrain);

//-------------------------------------------------------------------------------------------------

//int      BlitIntegratorInitialStates::maxNumHarmonics;
double** BlitIntegratorInitialStates::stateValues;
BlitIntegratorInitialStates::BlitIntegratorInitialStates()
{
  stateValues = nullptr;    
  //createStateValueTables(); // commented to fix memory leak
}
BlitIntegratorInitialStates::~BlitIntegratorInitialStates()
{
  //deleteStateValueTables(); // commented to fix memory leak
}
void BlitIntegratorInitialStates::createStateValueTables()
{
  if(stateValues == nullptr)
  {
    int numTables = maxNumHarmonics;

    stateValues = new double*[maxNumHarmonics];
    for(int i = 0; i < maxNumHarmonics; i++)
    {
      int numStates = i+1;
      stateValues[i] = new double[numStates];
    }


    for(int tableIndex = 0; tableIndex < numTables; tableIndex++)
    {
      int numPhases    = tableIndex+1;
      int numHarmonics = tableIndex;

      for(int phaseIndex = 0; phaseIndex < numPhases; phaseIndex++)
      {

        stateValues[tableIndex][phaseIndex] = 0.0;

        for(int k = 1; k < numHarmonics; k++)
          stateValues[tableIndex][phaseIndex] += (1.0 / (k*PI)) * sin(2*k*PI * phaseIndex/(numPhases-1.0));


        int dummy = 0;
      }


      /*
      if( numHarmonics == 4 )
      {
        double *indices = new double[maxNumHarmonics+1];
        rosic::fillWithIndex(indices, maxNumHarmonics+1);
        Plotter::plotData(numPhases, indices, stateValues[tableIndex]);
        delete[] indices;
      }
      */
    }

  }
}
void BlitIntegratorInitialStates::deleteStateValueTables()
{
  if(stateValues != NULL)
  {
    for(int i = 0; i < maxNumHarmonics; i++)
      delete[] stateValues[i];
    delete[] stateValues;
    stateValues = NULL;
  }
}

//-------------------------------------------------------------------------------------------------

BlitIntegratorInitialStates BlitSaw::initialStates;
void BlitSaw::initialize()
{
  BandlimitedImpulseTrain::initialize();
  fixedPhaseOffset = 0.0;
}
INLINE void BlitSaw::applyFilters(Module *module, double *out, int voiceIndex, double oscOmega)
{
  BlitSaw *saw = static_cast<BlitSaw*> (module);

  // apply leaky integrator (coefficient computation inlined):
  double yLPF  = saw->oldIntegratorOutputs[voiceIndex];
  double coeff = getLeakyIntegratorCoefficient(oscOmega);
  yLPF  = *out + coeff * yLPF;
  //yLPF *= coeff; // test
  *out  = yLPF;    // other test
  //*out  = -2.0*yLPF; 

  saw->oldIntegratorOutputs[voiceIndex] = yLPF; // state update
}
INLINE double BlitSaw::getLeakyIntegratorCoefficient(double oscOmega)
{
  //return 1.0;
  //return exp(- 10.0 * processingStatus.getFreqToOmegaFactor());  // at constant frequency
  return exp(-0.01*fabs(oscOmega));   // seems to be a good compromise between waveshape and 
                                      // dc-fadeout time
  //return exp(-0.1*fabs(oscOmega));
}

INLINE void BlitSaw::process(Module *module, double *in1, double *in2, double *out, int voiceIndex)
{
  // test:
  //*in1 = 440.0;

  BandlimitedImpulseTrain::process(module, in1, in2, out, voiceIndex);

  if(voiceAllocator.getNoteOnTriggerFlag(voiceIndex) == true)
    ((BlitSaw*)module)->resetIntegratorState(module, voiceIndex, *in2, *out, *in1);

  applyFilters(module, out, voiceIndex, *in1 * processingStatus.getFreqToOmegaFactor());
}
void BlitSaw::resetVoiceState(int voiceIndex)
{
  BandlimitedImpulseTrain::resetVoiceState(voiceIndex);

  resetIntegratorState(this, voiceIndex, 0.0, 1.0, 440.0);
    // will be reset with proper oscOmega value on note-on
}
void BlitSaw::resetIntegratorState(Module *module, int voiceIndex, double startPhase, double blitOut, double oscFreq)
{
  // To initialize the integrator's state, we first prescribe, which value the first output sample 
  // should have. The desired first output sample is computed by considering the continuous time 
  // sawtooth, defined as: saw(t) = 0.5 - mod(t, 1) which is a downward sawtooth between +0.5 
  // and -0.5 with unit period (unit period is appropriate because the startPhase is considered to 
  // be normalized to the interval 0...1 too). Having found this value, we set up the integrator's 
  // state such that desiredFirstSample = blitOutputSample + c * integratorState. Solving for 
  // integratorState gives the intial state value: 
  // double desiredFirstSample = 0.5 - fmod(startPhase + fixedPhaseOffset, 1.0);
  //...that all doens't seem to work - maybe we need indeed table for desiredFirstSample vs 
  // startPhase - for each value of numHarmonics, one table is needed - so we perhaps need an upper 
  // bound for numHarmonics ...perhaps 2000 or so (allows for full spectra (up to 20 kHz) for 
  // fundamentals down to 10 Hz. 


  double desiredFirstSample = getDesiredFirstSample(oscFreq, startPhase);
  double integratorCoeff    = getLeakyIntegratorCoefficient(oscFreq * processingStatus.getFreqToOmegaFactor());
  double integratorState    = (desiredFirstSample - blitOut) / integratorCoeff;

  oldIntegratorOutputs[voiceIndex] = integratorState;
}

double BlitSaw::getDesiredFirstSample(double frequency, double startPhase)
{
  //double desiredFirstSample = 0.5 - fmod(startPhase + fixedPhaseOffset, 1.0); 
    // old version - doesn't work well


  // algorithm to obtain the desired 1st sample:
  // -compute the cycle-length in samples
  // -round to next integer cycle length
  // -create a cycle of a bandlimited impulse train with the corresponding frequency and desired start phase
  // -perform integration over that cycle
  // -compute mean of integrated cycle
  // -use negative mean as desired start sample


  int    cycleLength = (int)floor(processingStatus.getSystemSampleRate() / frequency);
  double roundedFreq = processingStatus.getSystemSampleRate() / cycleLength;


  static const double maxNumHarmonics = 100000.0;
  double numHarmonics;
  double absFreq = fabs(frequency);
  if(absFreq < NANO)
    numHarmonics = maxNumHarmonics;
  else
    numHarmonics = floor(0.5 * processingStatus.getSystemSampleRate() / absFreq);
  double ampScaler;
  if(numHarmonics > 0)
    ampScaler = 0.5 * RAPT::rsSign(frequency) / numHarmonics;
  else
    ampScaler = 0.0;
  double omega    = frequency * processingStatus.getFreqToOmegaFactor();
  double absOmega = absFreq * processingStatus.getFreqToOmegaFactor();
  double theta    = TWO_PI * (startPhase);


  double firstBlitOut = ampScaler * BandlimitedImpulseTrain::computeUnscaledBlitValue(theta, numHarmonics);


  double blitOut;       // output of the BLIT
  double sum1 = 0.0;    // BLIT running sum
  double sum2 = 0.0;    // running sum over the running sum (to get the mean over the integrated cycle later on)
  double desiredFirstSample = 0.0;
  double integratorCoeff    = getLeakyIntegratorCoefficient(frequency * processingStatus.getFreqToOmegaFactor());
  integratorCoeff           = 1.0;
  double integratorState    = (desiredFirstSample - firstBlitOut) / integratorCoeff;

  // todo: find a closed form expression for this sum:
  for(int n = 0; n < cycleLength; n++)
  {
    blitOut          = ampScaler * BandlimitedImpulseTrain::computeUnscaledBlitValue(theta, numHarmonics);
    integratorState += integratorCoeff * blitOut;
    sum1             = integratorState;
    sum2            += sum1;
    theta           += omega;
    theta            = RAPT::rsWrapToInterval(theta, 0.0, TWO_PI);
  }


  double sawMean     = sum2 / cycleLength;
  desiredFirstSample = -sawMean;

  return desiredFirstSample;
}

void BlitSaw::allocateMemory()
{
  BandlimitedImpulseTrain::allocateMemory();
  oldIntegratorOutputs = new double[getNumVoices()];
}
void BlitSaw::freeMemory()
{
  BandlimitedImpulseTrain::freeMemory();
  delete[] oldIntegratorOutputs; oldIntegratorOutputs = NULL;
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(BlitSaw);

//-------------------------------------------------------------------------------------------------

void DualBlitSaw::initialize()
{
  initInputPins({ "Freq", "Phase", "Offset", "Mix" });
  initOutputPins({ "Out" });
  fixedPhaseOffset = 0.5;
}
INLINE void DualBlitSaw::process(Module *module, double *in1, double *in2, double *in3, 
  double *in4, double *out, int voiceIndex)
{
  DualBlitSaw *blit = static_cast<DualBlitSaw*> (module);

  // this macro computes numHarmonics, ampScaler, omega and theta and assigs them to local variables:
  COMPUTE_LOCALS_FOR_BLIT;

  // cretae the first blit:
  *out = ampScaler * computeUnscaledBlitValue(theta, numHarmonics);

  // add the second blit:
  theta  = blit->phases[voiceIndex] + 2.0 * PI * (*in2 + *in3 + blit->fixedPhaseOffset);
  *out  += *in4 * ampScaler * computeUnscaledBlitValue(theta, numHarmonics);


  incrementPhases(blit, voiceIndex, omega);


  if(voiceAllocator.getNoteOnTriggerFlag(voiceIndex) == true)
    ((DualBlitSaw*)module)->resetIntegratorState(module, voiceIndex, *in2, *out, *in1 * processingStatus.getFreqToOmegaFactor(), *in3, *in4);

  applyFilters(module, out, voiceIndex, *in1 * processingStatus.getFreqToOmegaFactor());
}
void DualBlitSaw::resetVoiceState(int voiceIndex)
{
  BlitSaw::resetVoiceState(voiceIndex);
}
void DualBlitSaw::resetIntegratorState(Module *module, int voiceIndex, double startPhase, double blitOut, double oscOmega,
  double phaseOffset, double secondBlitAmplitude)
{
  double desiredSaw1 =  0.5 - fmod(startPhase + fixedPhaseOffset, 1.0);
  double desiredSaw2 = (0.5 - fmod(startPhase + fixedPhaseOffset + phaseOffset, 1.0)) * secondBlitAmplitude;

  double desiredFirstSample = desiredSaw1 + desiredSaw2;


  // add the contribution of the 2nd saw to the desired 1st sample...


  double integratorCoeff    = getLeakyIntegratorCoefficient(oscOmega);
  double integratorState    = (desiredFirstSample - blitOut) / integratorCoeff;

  oldIntegratorOutputs[voiceIndex] = integratorState;

  //oldIntegratorOutputs[voiceIndex] = 0; // test

}

void DualBlitSaw::allocateMemory()
{
  BlitSaw::allocateMemory();
}
void DualBlitSaw::freeMemory()
{
  BlitSaw::freeMemory();
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_4(DualBlitSaw);
