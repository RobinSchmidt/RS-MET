#ifndef rosic_EffectModules_h
#define rosic_EffectModules_h

namespace rosic
{

  /**

  This file defines wrapper classes that wrap some core signal-processing objects into
  Module objects to facilitate their use in a (semi) modular framework such as Quadrifex.

  */

  class BypassModule : public Module
  {
  public:
    virtual void setSampleRate(double /*newSampleRate*/) { }
    virtual void processSampleFrame(double* /*inOutL*/, double* /*inOutR*/) { }
    virtual void reset() { }
  };

  class MuteModule : public Module
  {
  public:
    virtual void setSampleRate(double /*newSampleRate*/) { }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { *inOutL = *inOutR = 0.0; }
    virtual void reset() { }
  };

  class BitCrusherModule : public Module, public BitCrusher
  {
  public:
    virtual void setSampleRate(double /*newSampleRate*/) { }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { BitCrusher::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { BitCrusher::reset();       }
  };

  class ChorusModule : public Module, public Chorus
  {
  public:
    virtual void setSampleRate(double newSampleRate) { Chorus::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(double newTempo)      { Chorus::setTempoInBPM(newTempo);      }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Chorus::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { Chorus::reset();                 }
    virtual void trigger()     { Chorus::resetOscillatorPhases(); }
  };

  class CombBankModule : public Module, public CombBank
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { CombBank::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { CombBank::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { CombBank::reset(); }
  };

  class CombResonatorStereoModule : public Module, public CombResonatorStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { CombResonatorStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { CombResonatorStereo::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { CombResonatorStereo::reset(); }
  };

  class CombStereoizerModule : public Module, public CombStereoizer
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { CombStereoizer::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { CombStereoizer::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { CombStereoizer::reset();                 }
  };

  class DualTwoPoleFilterModule : public Module, public DualTwoPoleFilter
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { DualTwoPoleFilter::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { DualTwoPoleFilter::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { DualTwoPoleFilter::reset();                 }
  };

  class EqualizerModule : public Module, public EqualizerStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate) { EqualizerStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { EqualizerStereo::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()  { EqualizerStereo::reset();                }
  };

  class FlangerModule : public Module, public Flanger
  {
  public:
    virtual void setSampleRate(double newSampleRate) { Flanger::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(double newTempo)      { Flanger::setTempoInBPM(newTempo);      }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Flanger::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { Flanger::reset();                 }
    virtual void trigger()     { Flanger::resetOscillatorPhases(); }
  };

  class FormantShifterModule : public Module, public FormantShifterStereo
  {
  public:
    FormantShifterModule(int maxBlockSize) : FormantShifterStereo(maxBlockSize) {}
    virtual void setSampleRate(double newSampleRate)
    { FormantShifterStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { FormantShifterStereo::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { FormantShifterStereo::reset(); }
  };

  class FourPoleFilterModule : public Module, public FourPoleFilter
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { FourPoleFilter::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { FourPoleFilter::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { FourPoleFilter::reset();                 }
  };

  class FrequencyShifterStereoModule : public Module, public FrequencyShifterStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { FrequencyShifterStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { FrequencyShifterStereo::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { FrequencyShifterStereo::reset();                 }
  };

  class HarmonicsModule : public Module, public Harmonics
  {
  public:
    virtual void setSampleRate(double newSampleRate) { Harmonics::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Harmonics::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { Harmonics::reset();                 }
  };

  class LadderFilterModule : public Module, public LadderFilterOld
  {
  public:
    virtual void setSampleRate(double newSampleRate) { LadderFilterOld::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { LadderFilterOld::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { LadderFilterOld::reset();                 }
  };

  class LimiterModule : public Module, public Limiter
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { Limiter::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Limiter::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { Limiter::reset(); }
  };

  class ModulatedAllpassModule : public Module, public ModulatedAllpass
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { ModulatedAllpass::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { ModulatedAllpass::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { ModulatedAllpass::reset();                 }
  };

  class NoiseGateModule : public Module, public NoiseGate
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { NoiseGate::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { NoiseGate::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { NoiseGate::reset(); }
  };

  class NoisifierModule : public Module, public Noisifier
  {
  public:
    virtual void setSampleRate(double newSampleRate) { Noisifier::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Noisifier::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { Noisifier::reset();                 }
  };

  class PhaserModule : public Module, public Phaser
  {
  public:
    virtual void setSampleRate(double newSampleRate) { Phaser::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(double newTempo)      { Phaser::setTempoInBPM(newTempo);      }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Phaser::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { Phaser::reset();                 }
    virtual void trigger()     { Phaser::resetOscillatorPhases(); }
  };

  class PhaseStereoizerModule : public Module, public PhaseStereoizer
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { PhaseStereoizer::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { PhaseStereoizer::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { PhaseStereoizer::reset();                 }
  };

  class PingPongEchoModule : public Module, public PingPongEcho
  {
  public:
    virtual void setSampleRate(double newSampleRate) { PingPongEcho::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(double newTempo)      { PingPongEcho::setTempoInBPM(newTempo);      }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { PingPongEcho::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { PingPongEcho::reset();                 }

  };

  class PitchShifterModule : public Module, public PitchShifterGrainAdaptive
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { PitchShifterGrainAdaptive::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { PitchShifterGrainAdaptive::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { PitchShifterGrainAdaptive::reset(); }
  };

  class rsReverbModule : public Module, public rsReverb
  {
  public:
    virtual void setSampleRate(double newSampleRate) { rsReverb::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { rsReverb::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { rsReverb::reset();                 }
  };

  class RingModulatorModule : public Module, public RingModulatorStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { RingModulatorStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { RingModulatorStereo::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { RingModulatorStereo::reset(); }
  };

  class SingleSidebandModulatorModule : public Module, public SingleSidebandModulatorStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { SingleSidebandModulatorStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { SingleSidebandModulatorStereo::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { SingleSidebandModulatorStereo::reset(); }
  };

  class SimpleDelayModule : public Module, public FractionalDelayLineStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate) { FractionalDelayLineStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { FractionalDelayLineStereo::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { FractionalDelayLineStereo::clearDelayBuffers(); }
  };

  class SineOscillatorModule : public Module, public SineOscillator
  {
  public:
    virtual void setSampleRate(double newSampleRate) { SineOscillator::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { *inOutL = *inOutR = SineOscillator::getSample(); }
    virtual void reset() { SineOscillator::trigger(); }
  };

  class SlewRateLimiterStereoModule : public Module, public SlewRateLimiterStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate) { SlewRateLimiterStereo::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { SlewRateLimiterStereo::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { SlewRateLimiterStereo::reset();                 }
  };

  class SlopeFilterModule : public Module, public SlopeFilter
  {
  public:
    virtual void setSampleRate(double newSampleRate) { SlopeFilter::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { SlopeFilter::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { SlopeFilter::reset();                 }
  };

  class SoftKneeExpanderModule : public Module, public SoftKneeExpander
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { SoftKneeExpander::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { SoftKneeExpander::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { SoftKneeExpander::reset(); }
  };

  class SoftKneeCompressorModule : public Module, public SoftKneeCompressor
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { SoftKneeCompressor::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { SoftKneeCompressor::getSampleFrameStereo(inOutL, inOutR); }
    virtual void reset() { SoftKneeCompressor::reset(); }
  };

  class StereoPanModule : public Module, public StereoPan
  {
  public:
    virtual void setSampleRate(double /*newSampleRate*/) { }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { StereoPan::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset() {  }
  };

  class StereoWidthModule : public Module, public StereoWidth
  {
  public:
    virtual void setSampleRate(double /*newSampleRate*/) { }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { StereoWidth::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset() {  }
  };

  class TremoloModule : public Module, public Tremolo
  {
  public:
    virtual void setSampleRate(double newSampleRate) { Tremolo::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(double newTempo)      { Tremolo::setTempoInBPM(newTempo);      }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Tremolo::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()   { Tremolo::resetOscillatorPhases(); }
    virtual void trigger() { Tremolo::resetOscillatorPhases(); }
  };

  class TwoPoleFilterModule : public Module, public TwoPoleFilter
  {
  public:
    virtual void setSampleRate(double newSampleRate)
    { TwoPoleFilter::setSampleRate(newSampleRate); }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { TwoPoleFilter::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { TwoPoleFilter::reset();                 }
  };

  class VibratoModule : public Module, public Vibrato
  {
  public:
    virtual void setSampleRate(double newSampleRate) { Vibrato::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(double newTempo)      { Vibrato::setTempoInBPM(newTempo);      }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { Vibrato::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { Vibrato::reset();                 }
    virtual void trigger()     { Vibrato::resetOscillatorPhases(); }
  };

  class WahWahModule : public Module, public WahWah
  {
  public:
    virtual void setSampleRate(double newSampleRate) { WahWah::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(double newTempo)      { WahWah::setTempoInBPM(newTempo);      }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { WahWah::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { WahWah::reset();                 }
    virtual void trigger()     { WahWah::resetOscillatorPhases(); }
  };

  class WaveShaperModule : public Module, public WaveShaper
  {
  public:
    virtual void setSampleRate(double /*newSampleRate*/) {  }
    virtual void processSampleFrame(double *inOutL, double *inOutR)
    { WaveShaper::getSampleFrameStereo(inOutL, inOutR);         }
    virtual void reset()       { WaveShaper::reset();                 }
  };

} // end namespace rosic

#endif
