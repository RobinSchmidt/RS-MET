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

  void applyModulation(float m)
  {
    modulatedValue += m;
  }



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

class Processor // maybe find a better name maybe SignalProcessor
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
  void setParametersToDefaults(int index);
  // rename to resetParameters or setParametersToDefaults, maybe make non-virtual

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
  virtual void prepareToPlay(uchar key, uchar vel, double sampleRate);
  // Maybe make non-virtual. ...but maybe some subclasses want to override it? So far, none but 
  // maybe later - but then we can make it virtual only if needed


  virtual void processFrame(float* L, float* R) = 0;


  virtual void processBlock(float* L, float* R, int N)
  {
    for(int n = 0; n < N; n++)
      processFrame(&L[n], &R[n]);
  }



  /** Subclasses should override this in order to compute the algorithm's internal coefficients 
  from the user parameters. This will get called on in prepareToPlay and possibly also during 
  parameter modulation */
  virtual void updateCoeffs(double sampleRate) = 0;
  // Maybe it doesn't need to be purely virtual. Some very simple Processors may not have to do
  // anything ...but that will be the exception, i guess so the vast majority will need to override
  // it anyway

  // maybe it should take the sampleRate as parameter - then we don't nee a member for it - we 
  // could store it in the PlayStatus. actually, we could perhaps also directly access the 
  // PlayStatus object ourselves...but maybe it's more efficient to obtain it once in processFrame
  // and then pass it as parameter to all Processors
  // Maybe sampleRate should be float -> make benchmarks and numeric accuracy tests with both 
  // variants. Generally double -> float conversions tend to be somewhat expensive, i think.

  /** Subclasses should override this function to reset their internal state variables (if any), 
  such as a filter's past in- and output samples, an oscillator's phase etc. */
  virtual void resetState() {}

protected:

  std::vector<Parameter> params;
  OpcodeType type = OpcodeType::Unknown;
  uchar key = 0, vel = 0;                 // Used for key/vel tracking
  bool dirty = true; 
  // Maybe some more advanced processors need a finer dirtification control - use a char of 8 flags
  // to set certain parts of the coeff-set clean/dirty

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

/*
class Effect : public Processor  // maybe get rid -> unify Effect and Modulator
{
public:

  //virtual void processFrame(float* L, float* R) = 0;
  //virtual void processBlock(float* L, float* R, int N) = 0;

  // Maybe processFrame should not be purely virtual. It's just more boilerplate. In most cases
  // we will want to use processBlock exclusively anyway. We could have a default processFrame 
  // function that just calls processBlock with N=1
  // Maybe the block-processing should use also buffers for the modulated parameters maintained in
  // PlayStatus which can be re-used for all the different processors? Otherwise, I think, each 
  // Parameter object would need its own buffer which would make Parameter quite heavyweight.

};
*/

//-------------------------------------------------------------------------------------------------

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
    default:                return 0.f;
    }
  }
  // ToDo:
  // -Try to find better name
  // -Factor out the computations - we use similar computations in the mod-system in jura -> 
  //  consolidate the code
  // -Maybe move to .cpp to avoid inlining
  // -Maybe the multiplicative mode should literally multiply and not "translate" it into addition.
  //  The main reason for this trickery in jura's mod system was the fact that the ordering of the 
  //  modulators was hidden under the hood there - but here in sfz files, it's visible, so it may 
  //  be viable to make modulation dependent on the order of the connections...however, showing it 
  //  as mod-matrix would then again hide the order unless the matrix entries also show it somehow.
  //  It could be matrix of mod-widgets with a slider for amount, box for mode, draggable nuber for 
  //  position in the "chain". When we do this, we perhaps need to replace the getContribution
  //  function with an applyContribution function that directly applies the modulation to a 
  //  pre-existing accumulator value passed by reference instead of returning a value that is 
  //  supposed to be (additively) accumulated by the caller.

  /* obsolete
  inline void initTarget()
  {
    RAPT::rsAssert(targetParam != nullptr);
    targetParam->initModulatedValue();
  }
  */

  void setSource(Processor* newSource) { source = newSource; }

  void setSourceIndex(int newIndex) { sourceIndex = newIndex; }


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
    sourceIndex = -1;
    depth = 0.f; 
    mode  = ModMode::absolute; // maybe use ModMode::unknown as default
  }
  // maybe move to cpp

  int getSourceIndex() const { return sourceIndex; }


  Parameter* getTargetParam() const { return targetParam; }

private:

  Processor* source      = nullptr;  // e.g. EnvGen, LowFreqOsc, etc.
  Processor* targetProc  = nullptr;  // e.g. Filter, Amplifier, etc.
  // ToDo: replace source/targetProc with integers sourceIndex/targetProcIndex


  Parameter* targetParam = nullptr;  // Parameter objects are members of a Processor
  // ToDo: maybe replace targetParam with targetParamIndex...but maybe not because doing so would 
  // require an additional indirection in the per-sample code - we would have to go to the params
  // array of the target Processor...on the other hand, it may make the connector memory footprint 
  // smaller (we store a 32 bit int instead of a 64 bit pointer - maybe even 16 bit int would be 
  // enough)

  int sourceIndex = -1;
  // maybe get rid - if this is needed in the SamplePlayer, it should just keep an array of 
  // integers modSourceIndices for that purpose...hmm..or maybe it's ok this way..it's actually
  // more economic to avoid having an additional array there



  float depth  = 0.f;                // Strength of modulation
  ModMode mode = ModMode::absolute;  // maybe the default should depend on the target?

  // Maybe we also should keep a pointer to the target Processor not only the Parameter?
  // ...could be useful...we'll see
};

//=================================================================================================
// these should go into ModulatorCores:

/** Envelope generator for the sampler. */

class EnvGen : public Processor
{

public:

  /** Time parameters are in seconds. */
  void setup(float delay, float start, float attack, float hold, float decay, float sustain,
    float release, float sampleRate);


  void resetState() override 
  { 
    sampleCount = 0; out = 0.f;
    //core.resetState(); 
  }

  void processFrame(float* L, float* R) override
  {
    *L = *R = 0.f;  // preliminary
  }


  //void reset() { sampleCount = 0; out = 0.f; }

  /*
  void prepareToPlay(uchar key, uchar vel, double sampleRate) override
  {
    RAPT::rsError("Not yet implemented");
  }

  float getSample() override
  {
    RAPT::rsError("Not yet implemented");
    return 0.f;
  }
  */
  // get rid

  void updateCoeffs(double sampleRate) override { /* core.updateCoeffs(); */ }


protected:

  // move into a core class:
  RAPT::rsUint32 sampleCount;                          // sample counter
  RAPT::rsUint32 delay, attack, hold, decay, release;  // time parameters in samples
  float start, sustain;                                // level parameters 
  float out;                                           // output signal (stored between calls)

};
// -What value is the env supposed to produce during the delay stage? the start value or zero?

//=================================================================================================

/** Low frequency oscillator for the sampler. */

class LowFreqOsc : public Processor
{

public:

  LowFreqOsc();


  //void prepareToPlay(uchar key, uchar vel, double sampleRate) override;
  //float getSample() override;

  //void setup(float freq, float delay, float fade, float sampleRate); // move into a core class


  void processFrame(float* L, float* R) override
  {
    *L = *R = 0.1f;  // preliminary
  }


  void updateCoeffs(double sampleRate) override { /* core.updateCoeffs(); */ }

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

class Amplifier : public Processor
{
public:
  Amplifier();
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N) override;
  void updateCoeffs(double sampleRate) override;
protected:
  AmplifierCore core;
};

class Filter : public Processor
{
public:
  Filter();
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N) override;
  void updateCoeffs(double sampleRate) override;
  void resetState() override { core.resetState(); }
  static FilterCore::Type convertTypeEnum(FilterType sfzType);
protected:
  FilterCore core;
};

class Equalizer : public Processor
{
public:
  Equalizer();
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N) override;
  void updateCoeffs(double sampleRate) override;
  void resetState() override { core.resetState(); }
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

class WaveShaper : public Processor
{
public:
  WaveShaper();
  void processFrame(float* L, float* R) override;
  void processBlock(float* L, float* R, int N);
  void updateCoeffs(double sampleRate) override;
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
  Processor* grabEffect(OpcodeType type);
  // rename to grabProcessor

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositEffect(Processor* p);

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


  void allocateModulators();
  Processor* grabModulator(OpcodeType type);
  void repositModulator(Processor* p);


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
  Processor* grabEffect(OpcodeType type);

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositEffect(Processor* p);


  void allocateModulators();
  Processor* grabModulator(OpcodeType type);
  void repositModulator(Processor* p);

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