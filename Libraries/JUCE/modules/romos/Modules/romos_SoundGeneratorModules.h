#ifndef romos_SoundGeneratorModules_h
#define romos_SoundGeneratorModules_h

#include "../Framework/romos_ModuleAtomic.h"
#include "romos_ModuleDefinitionMacros.h"

namespace romos
{

  /** Generates noise that is spectrally white and has uniform distribution between -1 and +1. 
  Parameters:
   -Seed (initial state for the pseudo random number generator)
  Inputs:
   none
  Outputs:
   0: the noise signal

  \todo
  -include a shape parameter (number of added noise-gens - 2: triangular, etc.)
  -include min/max parameters
  -include (boolean) variance-normalization parameter (if true, divide output by sqrt(shape) before scaling and offsetting)
  */
  class WhiteNoise : public ModuleWithParameters
  {
    CREATE_COMMON_DECLARATIONS_0(WhiteNoise);
  public:
    virtual void resetVoiceState(int voiceIndex);
    virtual void parameterChanged(int index);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    unsigned long *state; // actually, we should use a self-defined uint32 to make sure it's machine independent
  };


  /** Generates a linear ramp that ascends from 0 to 1 in a given time the resets back to zero, ascends again, etc. 
  Parameters:
   none
  Inputs:
   0: Frq -> Frequency in Hz (1/period)
   later: include a startphase input and a retrigger input
  Outputs:
   0: the ramp signal
  */
  class PeriodicLinearRamp : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_1(PeriodicLinearRamp);
  public:
    virtual void resetVoiceState(int voiceIndex);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    double *phases;
  };


  
  /** Generates an impulse train that is bandlimited by the Nyquist frequency (i.e. an anti-aliased impulse-train) using a closed form 
  expression for a summation of a sine series.
  Parameters:
   none
  Inputs:
   0: Frequency of the fundamental
   1: Phase (an offset for the instantaneous phase)
  Outputs:
   0: the bandlimited impulse train
  References: 
   -http://www.music.mcgill.ca/~gary/307/week5/bandlimited.html
  */
  class BandlimitedImpulseTrain : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_2(BandlimitedImpulseTrain);
  public:
    virtual void resetVoiceState(int voiceIndex);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();

    INLINE static void incrementPhases(BandlimitedImpulseTrain *blit, const int voiceIndex, const double omega);

    INLINE static double computeUnscaledBlitValue(const double theta, const double numHarmonics);



    double *phases;

    double fixedPhaseOffset;

    // intermediate variables, to be used locally in process:
    //double numHarmonics, ampScaler, omega, theta;

  };



  class BlitIntegratorInitialStates
  {

  public:

    BlitIntegratorInitialStates();
    virtual ~BlitIntegratorInitialStates();

    static void createStateValueTables();
    static void deleteStateValueTables();


    static double **stateValues; // 1st index: numHarmonics, 2nd index: startPhase

    static int maxNumHarmonics;

  };

  /** This is similar to the BandlimitedImpulseTrain module, but it incorporates a leaky integrator after the impulse train so as to produce
  as sawtooth wave.
  Parameters:
   none
  Inputs:
   0: Frequency of the fundamental
   1: Phase (an offset for the instantaneous phase)
  Outputs:
   0: the sawtooth wave
  \todo:
   -there are still problems with the initial value for the leaky integrator for certain startPhase/frequency combinations
    ->perhaps related to initlaizing such that the desired first output is that of a non-bandlimited saw?
    ->further investigations necessary
  */
  class BlitSaw : public BandlimitedImpulseTrain
  {
    CREATE_COMMON_DECLARATIONS_2(BlitSaw);
  public:
    virtual void resetVoiceState(int voiceIndex);
    virtual void resetIntegratorState(Module *module, int voiceIndex, double startPhase, double blitOut, double oscOmega);
    INLINE static void applyFilters(Module *module, double *out, int voiceIndex, double oscOmega);
    INLINE static double getLeakyIntegratorCoefficient(double oscOmega);

    double getDesiredFirstSample(double frequency, double startPhase);

  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    double *oldIntegratorOutputs; // y[n-1] for the integrator


    static BlitIntegratorInitialStates initialStates;
  };



  /** This is similar to the BlitSaw module, but it generates two sawtooths with some phase-offset and weighting factor for the 
  second sawtooth. When the offset is +- 0.5 and the weight is -1.0, a square wave is generated. Keeping the weight at -1.0 and changing 
  the offset, a pulse-wave with variable pulse-width is obtained.
  Parameters:
   none
  Inputs:
   0: Freq -> frequency of the fundamental
   1: Phase -> offset for the instantaneous phase of the 1st (master) sawtooth
   2: Offset -> phase offset of the 2nd sawtooth with respect to them 1st
   3: Mix -> weight for the second sawtooth
  Outputs:
   0: the output wave
  */
  class DualBlitSaw : public BlitSaw // rename to BlitOscillator, get rid of the BlitSaw (it's a special case of this one)
  {
    CREATE_COMMON_DECLARATIONS_4(DualBlitSaw);
  public:
    virtual void resetVoiceState(int voiceIndex);
    virtual void resetIntegratorState(Module *module, int voiceIndex, double startPhase, double blitOut, double oscOmega, 
                                      double phaseOffset, double secondBlitAmplitude);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
  };





  // DigitalNoise (randomly switch between 1 and -1) Ins: MeanSwitchFrequency
  // SyncModOsc: oscillator is synced to a master-osc and has it's start-phase modulated by another osc 

  // SawPulseOsc: add two BLITs (one positive, one negative) and integrate the sum. the negtive BLIT's amplitude can be scaled
  // from 0...1: at 1: we have a pulse-wave, at 0: a saw-wave, so the scaling blends between saw and pulse
  // pulse-width parameter adjusts offset between positive and negative impulse
  // always subtract the DC before integrating, maybe on reset initialize the integrator so as to avoid DC transients on note-on
  // scale the resulting waveform so as to preserve the power for all settings of the Width and Blend parameters
  // -> this can serve as general purpose oscillator in VA synths
  // -later: maybe provide a "Sync" input (implement sync via MinBLEP) and a "Phase" input (for phase modulation)
  //...maybe call it BiSawOsc instead - inputs: offset (0.5: half-period), w2: amplitude of second saw, -1: for pulse)

  // SineOscBank, QuadratureOscBank -> let user enter relative frequencies on the GUI


} 

#endif 
