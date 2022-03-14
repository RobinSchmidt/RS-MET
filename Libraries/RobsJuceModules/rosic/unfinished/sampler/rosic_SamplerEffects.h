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

  // Misc:

  /** Initializes the modulated value to the nominal value. */
  void initModulatedValue() { modulatedValue = value; }



protected:

  Opcode opcode = Opcode::Unknown;
  float value = 0.f;               // maybe rename to nominalValue, maybe init to 666

  float modulatedValue = 666.f;
  // Uncomment that later for modulation, maybe use floatStereo because we want to support 
  // different modulator values for left and right. Think of an LFO mapped to cutoff with a phase
  // offset between left and right channel. That will easily give nice stereo movement. However,
  // not for all parameters does it make sense to have different left/right values. Think, for 
  // example, of a pan or width parameter - there is no meaningful concept of applying different
  // pan values to left and right channels. ...not yet sure, how to handle that... Maybe have two
  // generic outputs, i.e. channel1/channel2 instead of left/right and how they are used is opcode
  // defined - some (like cutoff) may interpret them as left/right, others (like pan) may just 
  // ignore the second channel. dunno...
};

//=================================================================================================

/** Baseclass for Effect and Modulator factoring out some common stuff...tbc... */

class Processor // maybe find a better name
{

public:

  using uchar = unsigned char;

  // Setup:
  void addParameter(Opcode opcode); // maybe should return an integer parameter index?
  void setParameter(Opcode opcode, float value);

  /** Resets all parameters to their default values. The caller needs to pass an index because 
  those default values may depend on that index. For example eq1_freq has 50, eq2_freq has 500 
  and eq3_freq has 5000 as default. So, the index is the index of a dsp of the particular given 
  type within the dsp chain. i.e. in the case of the eqs 1,2,3. It's not the index at which point 
  it sits in the dspChain. Only dsps of the same type count. */
  virtual void resetSettings(int index);
  // rename to resetParameters or setParametersToDefaults

  // Inquiry:
  OpcodeType getType() { return type; }
  int getNumParameters() const { return (int) params.size(); }

  /** Returns a pointer to our paramter object that "listens" to the given opcode if we have such a
  Parameter here, nullptr otherwise. The method is not const because the caller may modify the
  returned Parameter, i.e. call non-const methods on it. This is used for modulation...tbc... */
  Parameter* getParameter(Opcode op);

  //virtual const Parameter* getParameter() const { return nullptr; } // maybe relax constness later
  // check, if that's the right form of constness - we want to disable modfifying the parameter
  // itself but maybe reassign the pointer - yep, looks right:
  // https://stackoverflow.com/questions/21476869/constant-pointer-vs-pointer-to-constant



  // Processing:
  virtual void prepareToPlay(uchar key, uchar vel, double sampleRate) = 0;  

protected:

  std::vector<Parameter> params;
  OpcodeType type = OpcodeType::Unknown;

  uchar key = 0, vel = 0; 
  // needed for key/vel tracking? maybe should be in baseclass? envelopes will need it, too

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

//-------------------------------------------------------------------------------------------------

/** Baseclass for effect processors that can be applied to layers while they are the played back.
Subclasses can be various kinds of filters, equalizers, waveshapers, effects, etc. */

class Effect : public Processor  // maybe get rid -> unify Effect and Modulator
{
public:

  virtual void processFrame(float* L, float* R) = 0;
  virtual void processBlock(float* L, float* R, int N) = 0;

};

//-------------------------------------------------------------------------------------------------

/** Baseclass for modulators that can be applied to parameters of signal processors. Subclasses
can be envelopes, LFOs, etc. */
class Modulator : public Processor
{
  // maybe don't derive from Effect but instead derive both Effect and Modulator from a common
  // baseclass (maybe Processor) which contains only params and type

// Maybe treat modulators uniformly with effects...this may actually simplify the implementation
// *and* make it more flexible. We want stereo modulation anyway. Maybe let the LFO's load 
// samples as well using a syntax like lfoN_wave="SingleCycle/TwoSines.wav"
// https://sfzformat.com/opcodes/lfoN_wave

public:

  virtual float getSample() = 0;


  void updateModValue()
  {
    modValue = getSample();
  }
  float modValue = 0.f;


  // todo: processBlock, prepareToPlay
};
// To unify Effect and Modulator, we should use processFrame in modulators too to compute the 
// output. Then we may get rid of Modulator::getSample() and updateModValue(). The modValue does
// not need to be stored as member of the modulator - it can be held in a local variable in the 
// per-sample processing. The purely virtual Effect::processFrame and Effect::processBlock would 
// be moved into Processor and the Effect subclass would become obsolete. When doing it that way,
// modulators would become inherently capable of producing stereo output because processFrame is 
// layed out for producing stereo signals anyway. Not having to deal with these two subclasses 
// will simplify the code in other ways as well.


//=================================================================================================

class ModulationConnector  // maybe rename to ModulationWire or just Modulation or ModConnector
{

public:



  /** Computes and returns a contribution to be added to a modulated value that comes from this
  connection ...tbc.... */
  float getContribution(float modulatorOutput, float unmodulatedValue)
  {
    switch(mode)
    {
    case ModMode::absolute: return depth * modulatorOutput;
    case ModMode::relative: return depth * modulatorOutput * unmodulatedValue;
    default: RAPT::rsError("Unknown ModType");
    }
  }
  // ToDo:
  // -Try to find better name
  // -Factor out the computations - we use similar computations in the mod-system in jura -> 
  //  consolidate the code
  // -Maybe move to .cpp to avoid inlining

  inline void initTarget()
  {
    RAPT::rsAssert(targetParam != nullptr);
    targetParam->initModulatedValue();
  }

  void setSource(Modulator* newSource) { source = newSource; }


  void setTarget(Processor* newTargetProcessor, Parameter* newTargetParameter) 
  { 
    targetProc  = newTargetProcessor; 
    targetParam = newTargetParameter;
  }

  void setDepth(float newDepth) { depth = newDepth; }

  void setMode(ModMode newMode) { mode = newMode; }

  void reset()
  {
    source      = nullptr;
    targetProc  = nullptr; 
    targetParam = nullptr;
    depth = 0.f; 
    mode  = ModMode::absolute; // maybe use ModMode::unknown as default
  }
  // maybe move to cpp




private:

  Modulator* source      = nullptr;  // e.g. EnvGen, LowFreqOsc, etc.
  Processor* targetProc  = nullptr;  // e.g. Filter, Amplifier, etc.
  Parameter* targetParam = nullptr;  // Parameter objects are members of a Processor
  float depth  = 0.f;                // Strength of modulation
  ModMode mode = ModMode::absolute;  // maybe the default should depend on the target?

  // Maybe we also should keep a pointer to the target Processor not only the Parameter?
  // ...could be useful...we'll see
};

//=================================================================================================
// these should go into ModulatorCores:

/** Envelope generator for the sampler. */

class EnvGen : public Modulator
{

public:

  /** Time parameters are in seconds. */
  void setup(float delay, float start, float attack, float hold, float decay, float sustain,
    float release, float sampleRate);

  //float getSample();

  //void processBlock(float* samples, int numSamples); 
  void reset() { sampleCount = 0; out = 0.f; }

  void prepareToPlay(uchar key, uchar vel, double sampleRate) override
  {
    RAPT::rsError("Not yet implemented");
  }

  float getSample() override
  {
    RAPT::rsError("Not yet implemented");
    return 0.f;
  }


protected:

  RAPT::rsUint32 sampleCount;                          // sample counter
  RAPT::rsUint32 delay, attack, hold, decay, release;  // time parameters in samples
  float start, sustain;                                // level parameters 
  float out;                                           // output signal (stored between calls)

};
// -What value is the env supposed to produce during the delay stage? the start value or zero?

//=================================================================================================

/** Low frequency oscillator for the sampler. */

class LowFreqOsc : public Modulator
{

public:

  LowFreqOsc();
  void prepareToPlay(uchar key, uchar vel, double sampleRate) override;
  float getSample() override;

  //void setup(float freq, float delay, float fade, float sampleRate); // move into a core class



protected:

  // move into a core class:
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
  // ToDo: use a more efficient implementation and call it EqualizerCore - it needs to support 
  // only a biquad mode. Maybe use TDF1. A patch can use a lot of eq bands, so we may need many
  // eqs, so we should be more frugal with memory than for the filter opcode where there is 
  // typically only one or maybe two per note. ...but we probably should support shelving modes,
  // too. Maybe the EQ could have a type, too and also provide some of the types that filter has.
  // Then the sfz author could choose to pick eqN instead of filN, knowing that the implementation
  // is cheaper. ...but maybe it should provide keytrack and veltrack, too - but maybe not 
  // controller response ...what about MS-processing? maybe have another opcode big_eqN ...how
  // about EngineersFilter? maybe engfilN_type, engfilN_freq, etc.
};

class WaveShaper : public Effect
{
public:
  WaveShaper();
  void prepareToPlay(uchar key, uchar vel, double fs) override;
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N);
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

class EffectPool  // get rid
{

public:

  //EffectPool();
  //~EffectPool();

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
  Effect* grabEffect(OpcodeType type);
  // rename to grabProcessor

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositEffect(Effect* p);

  /** This is currently only meant to facilitate unit testing overload conditions. In such tests,
  we want a well defined and small number of filters to be available so we can simulate conditions
  where the engine is running out of filters in a controlled way. */
  //void setMaxNumFilters(int newMax) { filters.init(newMax); }


  //int getNumUsedFilters() const { return filters.getNumUsedItems(); }
  // also to facilitate testing


//protected:

  rsObjectPool<Amplifier>  amplifiers;
  rsObjectPool<Filter>     filters;
  rsObjectPool<Equalizer>  equalizers;
  rsObjectPool<WaveShaper> waveShapers;

};
// Maybe templatize, rename funcs to generic grabItem, repositItem. But maybe not - maybe we don't
// want to store a flat array of SignalProcessors of any kind but maintain different arrays for 
// different kinds of processors? or maybe keep the processors somehow sorted by type?

//=================================================================================================

class ModulatorPool // get rid
{

public:

  //ModulatorPool();
  //~ModulatorPool();

  void allocateModulators();
  Modulator* grabModulator(OpcodeType type);
  void repositModulator(Modulator* p);


protected:

  rsObjectPool<EnvGen>     envGens;
  rsObjectPool<LowFreqOsc> lowFreqOscs;

};
// Maybe consolidate the ModulatorPool and EffectPool into a single ProcessorPool class that 
// contains both types of processors. The distinction seems to just lead to code duplication.





//=================================================================================================

/** Structure to consolidate the different kinds of DSP resources */ 

class DspResourcePool
{

public:

  DspResourcePool();

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
  Effect* grabEffect(OpcodeType type);

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositEffect(Effect* p);


  void allocateModulators();
  Modulator* grabModulator(OpcodeType type);
  void repositModulator(Modulator* p);

  void allocateConnectors()
  {
    connectorPool.init(512);  // preliminary
  }
  ModulationConnector* grabConnector() { return connectorPool.grabItem(); }
  void repositConnector(ModulationConnector* c) { connectorPool.repositItem(c); }


  /** This is currently only meant to facilitate unit testing overload conditions. In such tests,
  we want a well defined and small number of filters to be available so we can simulate conditions
  where the engine is running out of filters in a controlled way. */
  void setMaxNumFilters(int newMax) { effectPool.filters.init(newMax); }


  int getNumUsedFilters() const { return effectPool.filters.getNumUsedItems(); }
  // also to facilitate testing

  int getNumIdleConnectors() const { return connectorPool.getNumIdleItems(); }



protected:

  EffectPool    effectPool;
  ModulatorPool modulatorPool;
  rsObjectPool<ModulationConnector> connectorPool;
};
// maybe turn this class into a sort of facade such that client code doesn't need to rech through
// into the member - have functions like grab/reposit Effect/Modulator/Connector. Maybe get rid of
// classes effectPool, modulatorPool - move implementations of its members into this class. data 
// members should be rsObjectPool<Effect> effectPool; rsObjectPool<Modulator> modulatorPool;
// -> more consistency/uniformity










}} // namespaces

#endif