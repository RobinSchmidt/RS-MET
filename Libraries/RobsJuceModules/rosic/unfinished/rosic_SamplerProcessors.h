#ifndef rosic_SamplerProcessors_h
#define rosic_SamplerProcessors_h

namespace rosic {
namespace Sampler {

//=================================================================================================

/** Class to represent parameters of DSP processors and modulators. */
class Parameter
{

public:

  // Lifetime:
  //Parameter() {}  // maybe needed later to construct arrays
  Parameter(Opcode opcode)
  {
    this->opcode = opcode;
    //this->value  = OpcodeDefaultValues[opcode];
  }

  // Setup:
  void setValue(float newValue) { value = newValue; }

  // Inquiry:
  //const char* getName() const { return name; }
  Opcode getOpcode() const { return opcode; }
  float getValue() const { return value; }
  //float getDefaultValue() const { return defaultValue; }

protected:

  Opcode opcode = Opcode::Unknown;
  float value = 0.f; 
};

/** Baseclass for signal processors that can be applied to layers while they are the played back.
Subclasses can be various kinds of filters, equalizers, waveshapers, effects, etc. */
class SignalProcessor
{
public:

  // Lifetime:

  // Setup:
  void addParameter(Opcode opcode); // maybe should return an integer parameter index?
  void setParameter(Opcode opcode, float value);

  // Inquiry:
  SignalProcessorType getType() { return type; }
  int getNumParameters() const { return (int) params.size(); }

  //virtual const Parameter* getParameter() const { return nullptr; } // maybe relax constness later
  // check, if that's the right form of constness - we want to disable modfifying the parameter
  // itself but maybe reassign the pointer - yep, looks right:
  // https://stackoverflow.com/questions/21476869/constant-pointer-vs-pointer-to-constant

  // Processing:
  virtual void prepareToPlay(double sampleRate) = 0;
  virtual void processFrame(rsFloat64x2& inOut) = 0;
  virtual void processBlock(rsFloat64x2* inOut, int N) = 0;

  //virtual void resetState() = 0;    // maybe remove
  //virtual void resetSettings() = 0; // ditto

protected:

  SignalProcessorType type = SignalProcessorType::Unknown;
  std::vector<Parameter> params;

  
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
  virtual double getSample() = 0;
  virtual void resetState() = 0;
  virtual void resetSettings() = 0;
  // todo: processBlock
};

//=================================================================================================
// The DSP classes below are meant to be used in the sampler, but they are nevertheless implemented
// as pure DSP classes without the infrastructural aspect, i.e. without being a subclass of 
// Sampler::SignalProcessor. This has been done in order to facilitate dragging them out of the 
// Sampler sub-namespace to make them available in other contexts as well. The infrastructure is 
// provided by boilerplate classes that derive from SignalProcessor and the actual core DSP class 
// via multiple inheritance. Maybe move this pure DSP code elsewhere...

/** A multimode filter that implements not only different filter frequency response types (like
lowpass, highbpass, bandpass, etc.) but even completely differently structured filters. Depending
on what mode has been chosen, the internal state and coefficient data may be interpreted in
different ways.... */

class rsSamplerFilter   // rename to FilterCore - the rs is not needed within this namespace
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsSamplerFilter()
  {

  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  enum class Type // maybe don't use all-caps for the types
  {
    // Biquad filter modes:
    Bypass,
    Lowpass_6,   // reanme to BQD_Lowpass_6 etc
    Highpass_6,

    // State variable filter modes:
    SVF_Lowpass_12,
    SVF_Lowpass_24,
    SVF_Highpass_12,
    SVF_Highpass_24,

    // Ladder filter modes:
    LDR_Lowpass_6,
    LDR_Lowpass_12,
    LDR_Lowpass_18,
    LDR_Lowpass_24
  };
  // hmm...
  // -maybe split the type into two parts: topology (svf, ladder, etc.), mode (lpf, hpf, etc.)
  // -maybe we can use a bitfield: 2 bits for the topology, 6 bits for the type within the selected
  //  topology, makes a total of 8 bits to specify the mode in a structured way. This is a 
  //  potential space optimization that may be done later. If done, it should not affect the API.


  void setup(Type type, float cutoff, float resonance);
  //void initCoeffs();
  //void updateCoeffs();

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  void processFrame(float& L, float& R);
  void resetState();


protected:


  //-----------------------------------------------------------------------------------------------
  /** \name Data Structures. We use structs to define the coefficient sets and internal states
  of the different filter topologies and then make a union from the structs to represent either of
  these. Then, we declare a member of the union type to store our data. */

  using TCoef = float;
  using TSig  = RAPT::rsVector2D<float>;  // for stereo
  // todo: maybe templatize this class and use float for TPar, and rsfloat32x2 for TSig in the 
  // sampler

  struct OnePoleVars    // ToDo: remove the "Vars" from the name
  {
    TSig  x1, y1;       // state
    TCoef b0, b1, a1;   // coeffs
    void resetState() { x1 = y1 = TSig(0); }
    void initCoeffs() { b0 = TCoef(1); b1 = a1 = TCoef(0); }
    TSig getSample(const TSig& in)
    {
      TSig y = b0*in + b1*x1 + a1*y1; // compute output
      x1 = in; y1 = y;                // update state
      return y;
    }
  };
  struct BiquadVars             // biquad filter, using DF2 (todo: try TDF1 -> smaller state)
  {
    TSig  x1, x2, y1, y2;       // state
    TCoef b0, b1, b2, a1, a2;   // coeffs
    void resetState() { x1 = x2 = y1 = y2 = TSig(0); }
    void initCoeffs() { b0 = TCoef(1); b1 = b2 = a1 = a2 = TCoef(0); }
    TSig getSample(const TSig& in)
    {
      TSig y = b0*in + b1*x1 + b2*x2 - a1*y1 - a2*y2;  // compute output
      x2 = x1; x1 = in; y2 = y1; y1 = y;               // update state
      return y;
    }
  };
  struct StateVars              // state variable filter (using ZDF)
  {
    TSig  s1, s2;               // state
    TCoef R;                    // damping(?)
  };
  struct LadderVars             // ladder filter (using UDF)
  {
    TSig  y1, y2, y3, y4;       // state
    TCoef a, k;                 // filter and feedback coeffs
    TCoef c0, c1, c2, c3, c4;   // output gains
  };
  union FilterVars              // Try to find a better name
  {
    FilterVars() {}             // without it, msc complains

    OnePoleVars p1z1;           // 1-pole 1-zero
    BiquadVars  bqd;
    StateVars   svf;
    LadderVars  ldr;
  };
  // ToDo: implement reset/getSample etc also in StateVars, etc. all thes structs should provide 
  // the sample API, but implement it in a way that is suitable to the given filter topology.


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  Type type = Type::Bypass;
  FilterVars vars;           // rename to something better

  // ToDo: maybe include also a state-vector filter (maybe rename to state phasor filter to avoid
  // name clash in abbreviation)
};

//=================================================================================================

/**  */

class rsSamplerWaveShaper
{

public:

  enum class Shape
  {
    None,
    Tanh,
    HardClip,
    SoftClipCubic,
    SoftClipHexic
    // ...etc.
  };

  void setup(Shape shape, float preGain, float dcOffset, float postGain, float shapePar1,
    float shapePar2)
  {
    this->shape     = shape;
    this->preGain   = preGain;
    this->dcOffset  = dcOffset;
    this->postGain  = postGain;
    this->shapePar1 = shapePar1;
    this->shapePar2 = shapePar2;
  }

  void processFrame(rsFloat64x2& inOut) 
  {
    float L = (float) inOut[0];
    float R = (float) inOut[1];
    L = postGain * tanh(preGain * L + dcOffset);
    R = postGain * tanh(preGain * R + dcOffset);
    inOut[0] = L;
    inOut[1] = R;
  }
  // under construction - optimize this!
  // -include a switch based on shape
  // -use a branchless tanh approximation directly operating on rsFloat64x2
  // -maybe switch to rsFloat32x4 for audio-samples...but maybe that's not advantageous
  // -implement a buffer-based variant

protected:

  Shape shape = Shape::None;
  float preGain, dcOffset, postGain;
  float shapePar1, shapePar2;            // meaning depends on chosen shape

  // what about oversampling..but maybe that should be reserved to an advanced variant? or maybe
  // that should apply to a region as a whole?
};

//=================================================================================================

/** Class where all the boilerplate for making DSP processors available in the sampler goes. */

class rsSamplerProcessors
{

public:

  class Filter : public SignalProcessor
  {

  public:

    Filter() 
    { 
      type = SignalProcessorType::Filter;
      params.reserve(3);
      addParameter(Opcode::FilterType);
      addParameter(Opcode::FilterCutoff);
      addParameter(Opcode::FilterResonance);
      // Having to pass a magic number to reserve() is bad and error-prone -> try to find a better
      // way. The number of parameters is actually known at compile time. Maybe use std::array 
      // instead of std::vector...hmm...but the number varies between the subclasses, so the array
      // could not be a baseclass member then...hmm...
    }
    void prepareToPlay(double fs) override 
    { 
      using TypeCore   = rsSamplerFilter::Type; // enum used in the DSP core
      using TypeOpcode = FilterType;            // enum used in the sfz opcode
      TypeCore type = TypeCore::Lowpass_6;      // preliminary
      core.setup(
        type, 
        params[1].getValue() * float(2*PI/fs),
        params[2].getValue());
      core.resetState();
    }
    // ToDo:
    // -try to avoid the translation step between the core-enum and sfz-enum - use a single
    //  enum for both
    // -not ideal that this code depends on the order, how we add the params in the constructor

    void processFrame(rsFloat64x2& inOut) override 
    {
      // ToDo: avoid these conversions - settle for a uniform sample format throughout the sampler,
      // perhaps rsFloat32x4 for compatibility with the float coeffs that are used in the DSP
      // core:
      float L = (float) inOut[0];
      float R = (float) inOut[1];
      core.processFrame(L, R);
      inOut[0] = L;
      inOut[1] = R;
    }
    void processBlock(rsFloat64x2* inOut, int N) override {}


    //void resetState() override { core.resetState(); }
    //void resetSettings() override { core.initCoeffs(); }

  protected:

    rsSamplerFilter core;

  };


  class WaveShaper : public SignalProcessor
  {

  public:

    WaveShaper() 
    { 
      type = SignalProcessorType::WaveShaper; 
      params.reserve(2);
      addParameter(Opcode::DistShape);
      addParameter(Opcode::DistDrive);
    }
    void prepareToPlay(double fs) override
    {
      using ShapeCore = rsSamplerWaveShaper::Shape;
      // using ShapeOpcode ...
      ShapeCore shape = ShapeCore::Tanh;      // preliminary
      core.setup(shape, params[1].getValue(), 0.f, 1.f, 0.f, 0.f);
      int dummy = 0;
    }
    void processFrame(rsFloat64x2& inOut) override { core.processFrame(inOut); }
    void processBlock(rsFloat64x2* inOut, int N) override {}
    //void resetState() override {}
    //void resetSettings() override {}

  protected:

    rsSamplerWaveShaper core;

  };

  // before writing too much boilerplate, fix the API:
  // -processFrame should work either with 2 floats or rsFloat32x4...or maybe make a class
  //  rsFloat32x2...could be just a synonym for rsFloat32x4 (because that can be simdified), but 
  //  the name makes clear that we use only 2 of the 4 available slots
  // -can we avoid the need for the boilerplate? ...or at least reduce the amount? maybe with 
  //  similar strategies as in romos, using macros?

};

//=================================================================================================

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
  /**< Array of the actual items. */

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

class SignalProcessorPool
{

public:

  SignalProcessorPool();

  ~SignalProcessorPool();

  /** Allocates the processors. */
  void allocateProcessors();
  // todo: Let it have an argument that somehwo specifies, how many of each type should be 
  // allocated. Maybe that could be a reference to the sfz-data itself or something derived from it




  /** A client can request a processor of the given type. If a processor of the desired type is 
  available, a pointer to it will returned. If not, a nullptr will be returned. The "client" will 
  typically be a RegionPlayer and call grabProcessor on noteOn. When no processor of the desired 
  type is available anymore, the calling code should probably forego the whole RegionPlayer. If
  the region can't be played correctly due to lack of resources, it should not play at all. What
  it certainly should not do is to just replace the non-available processor by a bypass dummy 
  processor because that could have really bad consequences: imagine a missing attenuation 
  processor. Regions are always played back either correctly or not at all but never wrongly. */
  SignalProcessor* grabProcessor(SignalProcessorType type);

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositProcessor(SignalProcessor* p);
  // maybe rename it to repositProcessor


protected:

  using SP = rsSamplerProcessors;

  rsObjectPool<SP::Filter>     filters;
  rsObjectPool<SP::WaveShaper> waveShapers;

  //std::vector<SP::Filter>     filters;
  //std::vector<SP::WaveShaper> waveShapers;
  // todo: use rsObjectPool instead of std::vector

};
// Maybe templatize, rename funcs to generic grabItem, repositItem. But maybe not - maybe we don't
// want to store a flat array of SignalProcessors of any kind but maintain different arrays for 
// different kinds of processors? or maybe keep the processors somehow sorted by type?

//=================================================================================================

/** Structure to consolidate the different kinds of DSP resources */ 

struct DspResourcePool
{
  SignalProcessorPool processorPool;
  //ModulatorPool modulatorPool;
  //ConnectionPool connectionPool;
};










}} // namespaces

#endif