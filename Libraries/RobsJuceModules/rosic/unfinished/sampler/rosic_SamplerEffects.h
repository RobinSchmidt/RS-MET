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
    //this->value  = SfzCodeBook::getInstance()->opcodeDefaultValue(opcode, 1);
    //initModulatedValue();
    // the getInstance call triggers an assert
  }
  // move to cpp

  // Setup:

  void setOpcode(Opcode newOpcode) {  opcode = newOpcode; }
  
  void setValue(float newValue) { value = newValue; }  // maybe rename to setNominalValue

  // Inquiry:
  Opcode getOpcode() const { return opcode; }
  float getValue() const { return value; }  // rename to getNominalValue

  // provide operator to convert to float to help avoidng boilerplate later

  // Misc:

  /** Initializes the modulated value to the nominal value. */
  void initModulatedValue() { modulatedValue = value; }

  void applyModulation(float m) { modulatedValue += m; }
  // use stereo input later

  /** Returns the modulated value. The function name is kept short to reduce the verbosity of the
  boilerplate in the Processors. */
  inline float mv() const { return modulatedValue; }
  // Maybe use the conversion operator instead to make the code even shorter...but this may be a
  // bit ambiguous because the caller may think that this refers to the unmodulated value


protected:

  Opcode opcode = Opcode::Unknown;
  float value = 0.f;               // maybe rename to nominalValue, maybe init to 666
  float modulatedValue = 666.f;
  // Maybe use rsFloatStereo later because we want to support different modulator values for left
  // and right. Think of an LFO mapped to cutoff with a phase offset between left and right 
  // channel. That will easily give nice stereo movement. However, not for all parameters does it
  // make sense to have different left/right values. Think, for example, of a pan or width 
  // parameter - there is no meaningful concept of applying different pan values to left and right
  // channels. ...not yet sure, how to handle that... Maybe have two generic outputs, i.e. 
  // channel1/channel2 instead of left/right and how they are used is opcode defined - some (like
  // cutoff) may interpret them as left/right, others (like pan) may just ignore the second 
  // channel. dunno...
  // Maybe the Parameter should also have a dirty flag to allow for a more fine-grained control
  // of the modulation optimization.
};

//=================================================================================================

/** Baseclass for Effect and Modulator factoring out some common stuff...tbc... */

class Processor // maybe find a better name maybe SignalProcessor
{

public:

  using uchar = unsigned char;

  // Setup:

  /** Adds a parameter that is controlled by the given Opcode. */
  void addParameter(Opcode opcode); // maybe should return an integer parameter index?


  void removeParameter(Opcode opcode);


  /** Replaces the oldOpcode with the given newOpcode. This is used for re-assigning existing 
  parameter objects to new opcodes. This is needed when we have different modules that are 
  functionally the same but have different opcode names (like ampeg_, fileg_, pitcheg_, egN_).
  See the constructor of EnvGenAmp for how this works and why we may want to do this. It returns 
  the index of the parameter in our params array which was replaced array or -1, if no param 
  listening to the given oldOpcode was found. */
  int replaceOpcode(Opcode oldOpcode, Opcode newOpcode);




  void setParameter(Opcode opcode, float value);

  /** Resets all parameters to their default values. The caller needs to pass an index because 
  those default values may depend on that index. For example eq1_freq has 50, eq2_freq has 500 
  and eq3_freq has 5000 as default. So, the index is the index of a dsp of the particular given 
  type within the dsp chain. i.e. in the case of the eqs 1,2,3. It's not the index at which point 
  it sits in the dspChain. Only dsps of the same type count. */
  void setParametersToDefaults(int index);
  // rename to resetParameters or setParametersToDefaults, maybe make non-virtual

  /** Sets our dirty flag. Used for optimizing modulation: updateCoeffs only needs to be called,
  if the dirty flag is true. */
  void setDirty() { dirty = true; }

  // Inquiry:
  OpcodeType getType() { return type; }

  int getNumParameters() const { return (int) params.size(); }

  bool isDirty() const { return dirty; }

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

  /** Returns the index of the Parameter object inside this Processor that listens to the given 
  opcode  or -1 if no such Parameter is found */
  int findParameter(Opcode opcode);



  std::vector<Parameter> params;
  OpcodeType type = OpcodeType::Unknown;
  uchar key = 0, vel = 0;                 // Used for key/vel tracking
  bool dirty = true; 
  // Maybe some more advanced processors need a finer dirtification control - use a char of 8 flags
  // to set certain parts of the coeff-set clean/dirty

  // float dryGain = 0.f;
  // float wetGain = 1.f;
  // should by default be 0,1 for effects and 1,1 for generators and x,1 for modulators (dryGain
  // doesn't matter for modulators because they receive no input signal anayway, i.e. all zeros).
  // Maybe have a Parameter for dry/wet mix in all Processors.

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

//=================================================================================================

class ModulationConnector  // maybe rename to ModulationWire or just Modulation or ModConnector
{

public:



  /** Computes and returns a contribution to be added to a modulated value that comes from this
  connection ...tbc.... */
  float getContribution(float modulatorOutput, float unmodulatedValue)
  {
    float m = modulatorOutput;
    float u = unmodulatedValue;
    float d = depth;

    switch(mode)
    {
    case ModMode::absolute: return d * m;
    case ModMode::relative: return d * m * u;
    case ModMode::cents:
    {
      float k = RAPT::rsPitchOffsetToFreqFactor(d * m / 100.f); // factor to be applied
      float c = k*u - u;  // Contribution. Verify this!
      return c;           // What about applying the depth here? Maybe in another mode
    }

    default:
    {
      RAPT::rsError("Unknown modulation mode");
      return 0.f;
    }

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
  // -see https://sfzformat.com/opcodes/_mod , https://sfzformat.com/opcodes/varNN_mod for how the 
  //  ARIA engine handles additive vs multiplicative modulation. It specifies it per modulation 
  //  target, not per connection

  /* obsolete
  inline void initTarget()
  {
    RAPT::rsAssert(targetParam != nullptr);
    targetParam->initModulatedValue();
  }
  */

  void setSource(Processor* newSource) { source = newSource; }
  // maybe rename to setSourceProcessor

  void setSourceArrayIndex(int newIndex) { sourceIndex = newIndex; }
  // maybe get rid and instead give a 2nd int parameter to setSource


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


  Processor* getSourceProcessor() { return source; }

  const Processor* getSourceProcessor() const { return source; }

  /** Returns the index of the source in the modSources array in SamplePlayer. This is set up in
  SamplePlayer::assembleRoutableModulations and retrieved in SamplePlayer::handleModulations. It's
  used to figure out, where the moudlation signal of the source was buffered. This should not be 
  confused with the sourceIndex stored in ModulationRouting - it's a different indeand means 
  something different. */
  int getSourceArrayIndex() const { return sourceIndex; }



  Parameter* getTargetParam() { return targetParam; }
  const Parameter* getTargetParam() const { return targetParam; }
  // maybe rename to getTargetParameter


  Processor* getTargetProcessor() { return targetProc; }

  const Processor* getTargetProcessor() const { return targetProc; }


  /*
  bool hasMatchingSource(ModulationConnector& c) const
  {
    return c.source == source;
  }
  bool hasMatchingTarget(const ModulationConnector& c) const
  {
    return c.targetParam == targetParam;
    // If targetParam pointers match, it implies that targetProc pointers must also match (because
    // the Processor objects own their Parameter objects), so we don't need to check both. 
  }
  bool hasMatchingEndpoints(const ModulationConnector& c) const
  {
    return hasMatchingSource(c) && hasMatchingTarget(c);
  }
  */


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
  // rename to sourceArrayIndex
  // maybe get rid - if this is needed in the SamplePlayer, it should just keep an array of 
  // integers modSourceIndices for that purpose...hmm..or maybe it's ok this way..it's actually
  // more economic to avoid having an additional array there



  float depth  = 0.f;                // Strength of modulation
  //ModMode mode = ModMode::absolute;  // maybe the default should depend on the target?
  ModMode mode = ModMode::unknown;   // Default mode depends on target, so we leave it blank here.

  // Maybe we also should keep a pointer to the target Processor not only the Parameter?
  // ...could be useful...we'll see
};

//=================================================================================================
// these should go into ModulatorCores:

/** Envelope generator for the sampler. */

/*
class EnvGen : public Processor
{

public:

  // Time parameters are in seconds. 
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

  //void updateCoeffs(double sampleRate) override { // core.updateCoeffs();  }


protected:

  // move into a core class:
  RAPT::rsUint32 sampleCount;                          // sample counter
  RAPT::rsUint32 delay, attack, hold, decay, release;  // time parameters in samples
  float start, sustain;                                // level parameters 
  float out;                                           // output signal (stored between calls)

};
*/
// -What value is the env supposed to produce during the delay stage? the start value or zero?

//=================================================================================================

// ToDo: maybe let the preprocessor generate this boilerplate using a macro:

class LowFreqOsc : public Processor
{
public:
  LowFreqOsc();
  void processFrame(float* L, float* R) override { core.processFrame(L, R); }
  void updateCoeffs(double sampleRate) override;
  void resetState() override { core.resetState(); }
protected:
  LowFreqOscCore core;
};

class LowFreqOscAmp : public LowFreqOsc
{
public:
  LowFreqOscAmp();
};

class LowFreqOscFil : public LowFreqOsc
{
public:
  LowFreqOscFil();
};





/** The general extended ADSR envelope generator controlled by the SFZ2 opocodes egN_attack, 
egN_decay, egN_sustain, egN_release, etc.. */
class EnvGen : public Processor
{
public:
  EnvGen();

  float getRelease()  const { return core.getRelease();  }
  float getEnd()      const { return core.getEnd();      }
  bool  hasFinished() const { return core.hasFinished(); }

  void processFrame(float* L, float* R) override { core.processFrame(L, R); }
  void updateCoeffs(double sampleRate) override;
  void resetState() override { core.resetState(); }

  void noteOff() { core.noteOff(); }

protected:
  EnvGenCore core;
};

/** A specialization of the general extended ADSR EG for the amplitude controlled by the SFZ1
opcodes ampeg_attack, ampeg_decay, etc.. */
class EnvGenAmp : public EnvGen
{
public:
  EnvGenAmp();
};

/** A specialization of the general extended ADSR EG for the amplitude controlled by the SFZ1
opcodes fileg_attack, fileg_decay, etc.. */
class EnvGenFil : public EnvGen
{
public:
  EnvGenFil();
};

/*
class EnvGenPitch : public EnvGen
{
public:
  EnvGenPitch();
};
*/


class PlayStatus; 
// Why do we need this forward declaration here? Try to re-organize the code to make it 
// unnecessary! maybe it's enough to change the include order of the .h files

/** A class to generate modulation signals from MIDI controllers */ 
class MidiController : public Processor
{
public:

  MidiController();
  void processFrame(float* L, float* R) override;
  void updateCoeffs(double sampleRate) override;

  /** Sets up the PlayStatus object to be used to grab the raw midi controller data from. Should be
  called soon after construction. */
  void setPlayStatusToUse(PlayStatus* ps) { playStatus = ps; }


protected:

  PlayStatus* playStatus = nullptr; 
  // Needs to be assigned shortly after creation. Is used to read out the current value of the midi
  // controller with the index that is set up by our index-parameter.

  int ctrlIndex = -1;           // midi controller number to listen to, e.g. 74 for cutoff
  RAPT::rsUint8 neutralVal = 0; // not yet used

  // ToDo: Maybe write a class MidiControllerCore. Maybe that class should contain the pointer to
  // the PlayStatus - or maybe not - but if we include smoothing later, the core would be the right
  // place to implement that.
};




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

/** Structure to consolidate the different kinds of DSP resources */ 

class DspResourcePool
{

public:

  DspResourcePool(PlayStatus* playStatusToUse);

  /**  */
  //void setPlayStatusToUse(PlayStatus* ps);

  //PlayStatus* playStatus = nullptr;


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

  /** Dispatches between grabEffect and grabModulator. */
  Processor* grabProcessor(OpcodeType type);

  // Maybe we should just have grabProcessor/repositProcessor - we don't really need to distinguish
  // between the two kinds here anymore. But first, get rid of EffectPool, ModulatorPool

  void allocateConnectors()
  {
    connectorPool.init(512);  // preliminary
  }
  ModulationConnector* grabConnector() { return connectorPool.grabItem(); }
  void repositConnector(ModulationConnector* c) 
  { 
    c->reset();  // not really neceassry
    connectorPool.repositItem(c); 
  }


  /** This is currently only meant to facilitate unit testing overload conditions. In such tests,
  we want a well defined and small number of filters to be available so we can simulate conditions
  where the engine is running out of filters in a controlled way. */
  void setMaxNumFilters(int newMax) { filters.init(newMax); }


  int getNumUsedFilters() const { return filters.getNumUsedItems(); }
  // also to facilitate testing

  int getNumIdleConnectors() const { return connectorPool.getNumIdleItems(); }



protected:

  // Effects:
  rsObjectPool<Amplifier>  amplifiers;
  rsObjectPool<Filter>     filters;
  rsObjectPool<Equalizer>  equalizers;
  rsObjectPool<WaveShaper> waveShapers;

  // Modulators:
  rsObjectPool<EnvGen>    freeEnvGens;
  rsObjectPool<EnvGenAmp> ampEnvGens;
  rsObjectPool<EnvGenFil> filEnvGens;

  rsObjectPool<LowFreqOsc>    freeLowFreqOscs;
  rsObjectPool<LowFreqOscAmp> ampLowFreqOscs;
  rsObjectPool<LowFreqOscFil> filLowFreqOscs;

  rsObjectPool<MidiController> midiControllers;


  // Connectors:
  rsObjectPool<ModulationConnector> connectorPool;
  // maybe this isn't needed - connections may be held as direct objects in the players - no
  // pointers needed


  PlayStatus* playStatus = nullptr;

};
// maybe turn this class into a sort of facade such that client code doesn't need to rech through
// into the member - have functions like grab/reposit Effect/Modulator/Connector. Maybe get rid of
// classes effectPool, modulatorPool - move implementations of its members into this class. data 
// members should be rsObjectPool<Effect> effectPool; rsObjectPool<Modulator> modulatorPool;
// -> more consistency/uniformity










}} // namespaces

#endif