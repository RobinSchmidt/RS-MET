#ifndef rosic_SamplerProcessors_h
#define rosic_SamplerProcessors_h

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

/** Class where all the boilerplate for making DSP processors available in the sampler goes. */

class rsSamplerEffects  // maybe make this a namespace and/or maybe get rid
{

public:

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
    void processFrame(float* L, float* R) override; //{ core.processFrame(L, R); }
    void processBlock(float* L, float* R, int N) override;
 

    static FilterCore::Type convertTypeEnum(FilterType sfzType)
    {
      // Conversion of filter type enum values used in the sfz data and those used in the dsp core.
      // Maybe we should try to avoid the translation step between the core-enum and sfz-enum by 
      // using a single enum for both. I'm not yet sure, if that's practical - we'll see. It could
      // turn out to be problematic when we want to use bit-twiddling of the enum-values to switch 
      // between different filter topologies in the core while the sfz-type number must allow for
      // lossless roundtrip with a float.
      using TC = FilterCore::Type;      // enum used in the DSP core
      using TO = FilterType;            // enum used in the sfz opcode
      switch(sfzType)
      {
      case TO::lp_6:   return TC::FO_Lowpass;
      case TO::hp_6:   return TC::FO_Highpass;

      case TO::lp_12:  return TC::BQ_Lowpass;        // ToDo: Use SVF as default implementation
      case TO::hp_12:  return TC::BQ_Highpass;       // for 2nd order filters...maybe...
      case TO::bp_6_6: return TC::BQ_Bandpass_Skirt; 
      case TO::br_6_6: return TC::BQ_Bandstop; 


      //case TO::lp_12: return TC::SVF_Lowpass_12;
      //case TO::hp_12: return TC::SVF_Highpass_12;
      }
      //RAPT::rsError("Unknown filter type in convertTypeEnum.");
      //return TC::Unknown;
      // This may actually happen when the user only defines cutoff but not fil_type. In this 
      // case, sfz prescribes that the default mode is lpf_2p.

      return TC::BQ_Lowpass;
    }
    // todo: avoid this conversion - use the same enum in both, the sfz codebook and the 
    // FilterCore just like we do with the waveshaper's distortion shapes

  protected:

    FilterCore core;
  };

  class Equalizer : public Effect
  {

  public:

    Equalizer()
    {
      type = DspType::Equalizer;
      params.reserve(3);
      addParameter(Opcode::eqN_gain);
      addParameter(Opcode::eqN_freq);
      addParameter(Opcode::eqN_bw); 
    }
    void prepareToPlay(uchar key, uchar vel, double fs) override 
    { 
      core.setupGainFreqBw(
        FilterCore::Type::BQ_Bell,
        params[0].getValue(),
        params[1].getValue() * float(2*PI/fs),
        params[2].getValue()
      );
      core.resetState();
    }
    void processFrame(float* L, float* R) override { core.processFrame(L, R); }
    void processBlock(float* L, float* R, int N) override 
    {
      for(int n = 0; n < N; n++)
        processFrame(&L[n], &R[n]);
    }

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

};

//=================================================================================================
// Should go to SamplerTools and maybe later to RAPT or somewhere in the rosic library.

/** Implements a pool of objects from which items can be grabbed and later returned. To "grab" an
item means to retrieve a pointer to it using grabItem(). The caller can then use the item fro as 
long as it needs it and when it is finished, it returns the item into the pool by calling
repositItem(T*) with the pointer. This will make the item available again such that it may be 
handed out again to other clients (or to the same client again). */

template<class T>
class rsObjectPool
{

public:

  /** Initializes the pool by allocating the given number of items. */
  void init(int numItems);

  /** Returns a pointer to an item or a nullptr if no items are available anymore. */
  T* grabItem();

  /** Reposits an item into the pool such that it becomes available again. Returns the index in 
  our array where the item is located or -1, if the item could not be found. The latter condition
  indicates an error in the client code: you may have tried to return an item that you did 
  previously grabbed from the pool or you may have re-alloacted the pool while there were still
  items  in use from the old one. Typically, it is irrelevant for the client, at which index the
  object lives, but the error detection may be useful. */
  int repositItem(T* item);

  /** Like repositItem(T*) but accepts also pointers that were not declared as pointing to T. 
  This version is needed when the client uses a pointer to a baseclass of T. */
  int repositItem(void* item) { return repositItem((T*)item); }
  // todo: maybe templalize it on some type B and make an assertion that B is baseclass of T, if 
  // that is possible. That would be a bit more restrictive than accepting void and perhaps
  // catch more errors at compile time


  int getNumItems() const { return numUsed + numIdle; }

  int getNumUsedItems() const { return numUsed; }

  int getNumIdleItems() const { return numIdle; }
  


  /** Checks, if the class invariants are satisfied. */
  bool isInConsistentState();


protected:

  std::vector<T> items; 
  /**< Array of the actual items. We store them here as direct objects and hand out pointers
  to them to the client. 
  \todo: Maybe it should later use rsRetainedArray (a class currently in development in the 
  Prototypes) to enable dynamic growth of the pool without invalidation of the pointers at the 
  cost of discontiguous object storage (a cost which is only incurred when growth is actually 
  happening - the first pre-allocated chunk, the size of which can be chosen at initialization, 
  will be stored contiguously). Currently, dynamic growth is not supported. */

  std::vector<char> used;  
  /**< Parallel array of flags indicating, if item[i] is currently used by anyone. If 0 (false),
  the item is avaible and can be handed out. If true, it can't **/

  int numUsed = 0;  /**< Number of used items */
  int numIdle = 0;  /**< Number of items available for handing out. */

};
// todo: 
// -move to rapt
// -write unit tests
// -maybe switch std::vector to rsRetainedArray
// -in this implementation, the items are assumed to be indsitigu
// -maybe don't use a std::vector for the items but rather a raw pointer
// -maybe introduce another level on indirection to enable grab/reposit in constant or at least
//  log(N) time: always grab objects from the end of a list of available objects, when objects
//  are reposited, they are appended to the end...hmm...i have to think about it some more - we 
//  want to avoid having to iterate through the whole used[] array to find an idle item and we also
//  want to avoid iterating through the list on reposit to find the location where we need to set
//  the used flag back to false

template<class T>
void rsObjectPool<T>::init(int numItems)
{
  RAPT::rsAssert(numUsed == 0 && isInConsistentState());
  items.resize(numItems);
  used.resize(numItems);  // will be initialized to zero (that's important!)
  numUsed = 0;
  numIdle = (int) used.size();
}

template<class T>
T* rsObjectPool<T>::grabItem()
{
  for(size_t i = 0; i < used.size(); i++) {
    if(!used[i]) {
      used[i] = 1;
      numUsed++;
      numIdle--;
      return &(items[i]); }}
  return nullptr;
}

template<class T>
int rsObjectPool<T>::repositItem(T* item)
{
  for(size_t i = 0; i < items.size(); i++) {
    if(item == &(items[i]) ) {
      used[i] = 0;
      numUsed--;
      numIdle++;
      return (int)i; }}
  return -1;
}
// -grab/reposit have currently linear complexity in the number of stored items. ToDo: try to 
//  reduce it to log or even const - keep this class then as prototype
// -maybe use range-based loops to later make the switch to rsRetainedArray easier

template<class T>
bool rsObjectPool<T>::isInConsistentState()
{
  bool ok = items.size() == used.size();
  int count = 0;
  for(const auto& it : used) {
    if( it != 0 )
      count++; }
  ok &= count == numUsed;
  ok &= numIdle == (int) used.size() - numUsed;
  return ok;
}


//=================================================================================================

/** A class for storing a pool of SignalProcessor objects...tbc... */

class SignalProcessorPool  // rename to EffectPool
{

public:

  SignalProcessorPool();

  ~SignalProcessorPool();

  /** Allocates the processors. */
  void allocateProcessors();
  // todo: Let it have an argument that somehow specifies, how many of each type should be 
  // allocated. Maybe that could be a reference to the sfz-data itself or something derived from it




  /** A client can request a effect of the given type. If a effect of the desired type is 
  available, a pointer to it will returned. If not, a nullptr will be returned. The "client" will 
  typically be a RegionPlayer and call grabEffect on noteOn. When no effect of the desired 
  type is available anymore, the calling code should probably forego the whole RegionPlayer. If
  the region can't be played correctly due to lack of resources, it should not play at all. What
  it certainly should not do is to just replace the non-available effect by a bypass dummy 
  effect because that could have really bad consequences: imagine a missing attenuation 
  effect. Regions are always played back either correctly or not at all but never wrongly. */
  Effect* grabProcessor(DspType type);  // rename to grabEffect

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositProcessor(Effect* p);
  // maybe rename it to repositProcessor


  /** This is currently only meant to facilitate unit testing overload conditions. In such tests,
  we want a well defined and small number of filters to be available so we can simulate conditions
  where the engine is running out of filters in a controlled way. */
  void setMaxNumFilters(int newMax) { filters.init(newMax); }


  int getNumUsedFilters() const { return filters.getNumUsedItems(); }
  // also to facilitate testing


protected:

  using SP = rsSamplerEffects;  // rename or remove
  rsObjectPool<SP::Amplifier>  amplifiers;
  rsObjectPool<SP::Filter>     filters;
  rsObjectPool<SP::Equalizer>  equalizers;
  rsObjectPool<SP::WaveShaper> waveShapers;

};
// Maybe templatize, rename funcs to generic grabItem, repositItem. But maybe not - maybe we don't
// want to store a flat array of SignalProcessors of any kind but maintain different arrays for 
// different kinds of processors? or maybe keep the processors somehow sorted by type?

//=================================================================================================

/** Structure to consolidate the different kinds of DSP resources */ 

struct DspResourcePool
{
  SignalProcessorPool processorPool;  // rename to EffectPool
  //ModulatorPool modulatorPool;
  //ConnectionPool connectionPool;
};










}} // namespaces

#endif