#ifndef rosic_SamplerEffects_h
#define rosic_SamplerEffects_h

// ToDo: rename to SamplerEffects.h/cpp, factor out a file SamplerEffectCores.h/cpp

namespace rosic {
namespace Sampler {

//=================================================================================================

/** Class to represent parameters of effect processors and modulators. */

class Parameter
{

public:

  // Lifetime:
  //Parameter() {}  // maybe needed later to construct arrays
  Parameter(Opcode opcode)
  {
    this->opcode = opcode;
    //this->value  = SfzCodeBook::getInstance()->getDefaultValue(opcode); // maybe do this later
  }

  // Setup:
  void setValue(float newValue) { value = newValue; }

  // Inquiry:
  Opcode getOpcode() const { return opcode; }
  float getValue() const { return value; }

  // provide operator to convert to float to help avoidng boilerplate later

protected:

  Opcode opcode = Opcode::Unknown;
  float value = 0.f;        // maybe rename to nominalValue, maybe init to 666
  // float modulatedValue;  // ...later, for modulation
};

//=================================================================================================

/** Baseclass for effect processors that can be applied to layers while they are the played back.
Subclasses can be various kinds of filters, equalizers, waveshapers, effects, etc. */

class Effect
{
public:

  using uchar = unsigned char;

  // Lifetime:

  // Setup:
  void addParameter(Opcode opcode); // maybe should return an integer parameter index?
  void setParameter(Opcode opcode, float value);

  // Inquiry:
  DspType getType() { return type; }
  int getNumParameters() const { return (int) params.size(); }

  //virtual const Parameter* getParameter() const { return nullptr; } // maybe relax constness later
  // check, if that's the right form of constness - we want to disable modfifying the parameter
  // itself but maybe reassign the pointer - yep, looks right:
  // https://stackoverflow.com/questions/21476869/constant-pointer-vs-pointer-to-constant

  // Processing:
  virtual void prepareToPlay(uchar key, uchar vel, double sampleRate) = 0;  // maybe it needs to receive the key, vel

  virtual void processFrame(float* L, float* R) = 0;
  virtual void processBlock(float* L, float* R, int N) = 0;

  /** Resets all parameters to their default values. The caller needs to pass an index because 
  those default values may depend on that index. For example eq1_freq has 50, eq2_freq has 500 
  and eq3_freq has 5000 as default. So, the index is the index of a dsp of the particular given 
  type within the dsp chain. i.e. in the case of the eqs 1,2,3. It's not the index at which point 
  it sits in the dspChain. Only dsps of the same type count. */
  virtual void resetSettings(int index);

protected:

  std::vector<Parameter> params;
  DspType type = DspType::Unknown;
  uchar key = 0, vel = 0;

  // For response to midi control:
  // MidiStatus midiStatus;
  // bool dirty = false;
  // clean, transitional
  // uchar voice; // maybe needed for cutoff_polyaft etc. - should be passed into prepareToPlay
  // uchar chan;  // maybe needed for cutoff_chanaft

  
//private:  // doesn't compile because subclasses complain...hmm...
  
  /*
  // make uncopyable:
  SignalProcessor() = default;
  SignalProcessor(const SignalProcessor&) = delete;
  SignalProcessor & operator=(const SignalProcessor&) = delete;
  */
  // ToDo: define a macro for that, see:
  // https://stackoverflow.com/questions/2173746/how-do-i-make-this-c-object-non-copyable
};

/** Baseclass for modulators that can be applied to parameters of signal processors. Subclasses
can be envelopes, LFOs, etc. */
class Modulator
{
public:
  virtual float getSample() = 0;
  // todo: processBlock, prepareToPlay
};


//=================================================================================================
// these should go into ModulatorCores:

/** Envelope generator for the sampler. */

class rsSamplerEnvGen
{

public:

  /** Time parameters are in seconds. */
  void setup(float delay, float start, float attack, float hold, float decay, float sustain,
    float release, float sampleRate);

  float getSample();

  //void processBlock(float* samples, int numSamples); 
  void reset() { sampleCount = 0; out = 0.f; }

protected:

  RAPT::rsUint32 sampleCount;                          // sample counter
  RAPT::rsUint32 delay, attack, hold, decay, release;  // time parameters in samples
  float start, sustain;                                // level parameters 
  float out;                                           // output signal (stored between calls)

};
// -What value is the env supposed to produce during the delay stage? the start value or zero?

//=================================================================================================

/** Low frequency oscillator for the sampler. */

class rsSamplerLowFreqOsc
{

public:

  void setup(float freq, float delay, float fade, float sampleRate);


protected:

  float pos;                   // normalized position in the wave in 0..1
  float inc;                   // per sample increment for pos
  RAPT::rsUint32 delay, fade;  // delay and fade-in time in samples

};
// -Maybe we can get rid of delay by initializing pos to -delay and the implementation of "fade"
//  returns zero for pos < 0
// -Introduce a start-phase variable
// -Maybe have a function pointer to a function that produces the actual waveform


//=================================================================================================

// ToDo: maybe let the preprocessor generate this boilerplate using a macro:

class Amplifier : public Effect
{
public:
  Amplifier();
  void prepareToPlay(uchar key, uchar vel, double fs) override;
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N) override;
protected:
  AmplifierCore core;
};

class Filter : public Effect
{
public:
  Filter();
  void prepareToPlay(uchar key, uchar vel, double fs) override;
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N) override;
  static FilterCore::Type convertTypeEnum(FilterType sfzType);
protected:
  FilterCore core;
};

class Equalizer : public Effect
{
public:
  Equalizer();
  void prepareToPlay(uchar key, uchar vel, double fs) override;
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N) override;
protected:

  FilterCore core;
  // ToDo: use a more efficient implementation - it needs to support only a biquad mode. Maybe 
  // use TDF1. A patch can use a lot of eq bands, so we may need many eqs, so we should be more
  // frugal with memory than for the filter opcode where there is typically only one or maybe two
  // per note. ...but we probably should support shelving modes, too
};

class WaveShaper : public Effect
{

public:

  WaveShaper()
  {
    type = DspType::WaveShaper;
    params.reserve(3);
    addParameter(Opcode::distortN_shape);
    addParameter(Opcode::distortN_drive);
    addParameter(Opcode::distortN_dc);
  }
  void prepareToPlay(uchar key, uchar vel, double fs) override
  {
    core.setup((DistortShape)(int)params[0].getValue(), params[1].getValue(),
      params[2].getValue(), 1.f, 0.f, 0.f);
    int dummy = 0;
  }
  void processFrame(float* L, float* R) override { core.processFrame(L, R); }
  void processBlock(float* L, float* R, int N) override
  {
    for(int n = 0; n < N; n++)
      processFrame(&L[n], &R[n]);
  }

protected:

  WaveshaperCore core;

};

// before writing too much boilerplate, fix the API:
// -processFrame should work either with 2 floats or rsFloat32x4...or maybe make a class
//  rsFloat32x2...could be just a synonym for rsFloat32x4 (because that can be simdified), but 
//  the name makes clear that we use only 2 of the 4 available slots
// -can we avoid the need for the boilerplate? ...or at least reduce the amount? maybe with 
//  similar strategies as in romos, using macros?


//=================================================================================================

/** A class for storing a pool of Effect objects...tbc... */

class EffectPool
{

public:

  EffectPool();

  ~EffectPool();

  /** Allocates the effects by resizing our vectors (which contain direct objects). */
  void allocateEffects();
  // todo: Let it have an argument that somehow specifies, how many of each type should be 
  // allocated. Maybe that could be a reference to the sfz-data itself or something derived from it


  /** A client can request an effect of the given type. If an effect of the desired type is 
  available, a pointer to it will returned. If not, a nullptr will be returned. The "client" will 
  typically be a RegionPlayer and call grabEffect on noteOn when assembling the effect chain. When
  no effect of the desired type is available anymore, the calling code should probably forego the 
  whole RegionPlayer. If the region can't be played correctly due to lack of resources, it should 
  not play at all. What it certainly should not do is to just replace the non-available effect by a
  bypass dummy effect because that could have really bad consequences: imagine a missing 
  attenuation effect. Regions are always played back either correctly or not at all but never 
  wrongly. */
  Effect* grabEffect(DspType type);

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositEffect(Effect* p);

  /** This is currently only meant to facilitate unit testing overload conditions. In such tests,
  we want a well defined and small number of filters to be available so we can simulate conditions
  where the engine is running out of filters in a controlled way. */
  void setMaxNumFilters(int newMax) { filters.init(newMax); }


  int getNumUsedFilters() const { return filters.getNumUsedItems(); }
  // also to facilitate testing


protected:

  rsObjectPool<Amplifier>  amplifiers;
  rsObjectPool<Filter>     filters;
  rsObjectPool<Equalizer>  equalizers;
  rsObjectPool<WaveShaper> waveShapers;

};
// Maybe templatize, rename funcs to generic grabItem, repositItem. But maybe not - maybe we don't
// want to store a flat array of SignalProcessors of any kind but maintain different arrays for 
// different kinds of processors? or maybe keep the processors somehow sorted by type?

//=================================================================================================

/** Structure to consolidate the different kinds of DSP resources */ 

struct DspResourcePool
{
  EffectPool effectPool;
  //ModulatorPool modulatorPool;
  //ConnectionPool connectionPool;
};










}} // namespaces

#endif