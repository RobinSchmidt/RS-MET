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
  Opcode getOpcode() const { return opcode; }
  float getValue() const { return value; }

  // provide operator to convert to float to help avoidng boilerplate later

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
  DspType getType() { return type; }
  int getNumParameters() const { return (int) params.size(); }

  //virtual const Parameter* getParameter() const { return nullptr; } // maybe relax constness later
  // check, if that's the right form of constness - we want to disable modfifying the parameter
  // itself but maybe reassign the pointer - yep, looks right:
  // https://stackoverflow.com/questions/21476869/constant-pointer-vs-pointer-to-constant

  // Processing:
  virtual void prepareToPlay(double sampleRate) = 0; // maybe it needs a flag, if input is stereo
  virtual void processFrame(float* L, float* R) = 0;
  virtual void processBlock(float* L, float* R, int N) = 0;

  /** Resets all parameters to their default values. The caller needs to pass an index because 
  those default values may depend on that index. For example eq1_freq has 50, eq2_freq has 500 
  and eq3_freq has 5000 as default. So, the index is the index of a dsp of the particular given 
  type within the dsp chain. i.e. in the case of the eqs 1,2,3. It's not the index at which point 
  it sits in the dspChain. Only dsps of the same type count. */
  virtual void resetSettings(int index);

protected:

  DspType type = DspType::Unknown;
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
  virtual float getSample() = 0;
  // todo: processBlock, prepareToPlay
};

//=================================================================================================
// The DSP classes below are meant to be used in the sampler, but they are nevertheless implemented
// as pure DSP classes without the infrastructural aspect, i.e. without being a subclass of 
// Sampler::SignalProcessor. This has been done in order to facilitate dragging them out of the 
// Sampler sub-namespace to make them available in other contexts as well. The infrastructure is 
// provided by boilerplate classes that derive from SignalProcessor and the actual core DSP class 
// via multiple inheritance. Maybe move this pure DSP code elsewhere...


//=================================================================================================

class AmplifierCore
{

public:

  /** Sets up the internal algo parameters in terms of user parameters defined in the sfz spec.
  The following descriptions of the parameters are taken from: https://sfzformat.com/legacy/ with
  some additions of my own signified by [square brackets].
  (ToDo: verify, if this is implemented correctly - at the moment, my formulas are a guess)

  volume: -144 to 6 dB
  The volume for the region, in decibels.

  pan: -100..+100
  The panoramic position for the region. If a mono sample is used, pan value defines the position 
  in the stereo image where the sample will be placed. When a stereo sample is used, the pan value 
  [controls] the relative amplitude of one channel respect the other. A value of zero means 
  centered, negative values move the panoramic to the left, positive to the right.

  width: -100% to 100%
  Only operational for stereo samples, width defines the amount of channel mixing applied to play 
  the sample. A width value of 0 makes a stereo sample play as if it were mono (adding both 
  channels and compensating for the resulting volume change). A value of 100 will make the stereo 
  sample play as original. Any value in between will mix left and right channels with a part of the
  other, resulting in a narrower stereo field image. Negative width values will reverse 
  [i.e. swap? not negate?] left and right channels.

  position: -100% to 100%
  Only operational for stereo samples, position defines the position in the stereo field of a 
  stereo signal, after channel mixing as defined in the width opcode. A value of zero means 
  centered, negative values move the panoramic to the left, positive to the right. Examples:
  width=0 position=-100 will mix both channels and play the result at left, width=50 position=30 
  will make the stereo image narrower and play it slightly right. */
  void setup(float volume, float pan, float width, float position);

  // ToDo: 
  // -Maybe provide a different parametrization where the volume is expressed as linear gain, 
  //  so we may also model polarity inversions. Or maybe realize this here via an additional 
  //  boolean flag "invertPolarity"

  /** Returns the channel mixing matrix coefficients that determine the scale factors by which both
  input channels (L,R) go to the output channels. */
  void getChannelMixCoeffs(float* L2L, float* R2L, float* L2R, float* R2R)
  { *L2L = gLL; *R2L = gLR; *L2R = gRL; *R2R = gRR; }


  void processFrame(float* L, float* R)
  {
    float t = *L;            // temporary
    *L = gLL * *L + gLR * *R;
    *R = gRL *  t + gRR * *R;
  }


protected:

  float gLL = 1.f, gLR = 0.f, gRL = 0.f, gRR = 1.f;

};

//=================================================================================================

/** A multimode filter that implements not only different filter frequency response types (like
lowpass, highbpass, bandpass, etc.) but even completely differently structured filters. Depending
on what mode has been chosen, the internal state and coefficient data may be interpreted in
different ways.... */

class FilterCore
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  //FilterCore() { }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  enum class Type // maybe don't use all-caps for the types
  {
    Bypass,

    // 1st order filter modes:
    // FO_Start = 1 << 24,
    FO_Lowpass,
    FO_Highpass,
    FO_LowShelf,
    FO_HighShelf,
    FO_Allpass,

    // Biquad filter modes:
    // BQ_Start = 2 << 24,
    BQ_Lowpass,
    BQ_Highpass,
    BQ_Bandpass_Skirt,
    BQ_Bandpass_Peak,
    BQ_Bandstop,
    BQ_Bell,
    BQ_LowShelf,
    BQ_HighShelf,
    BQ_Allpass,

    // State variable filter modes:
    // SVF_Start = 3 << 24,
    SVF_Lowpass_12,
    SVF_Lowpass_24,
    SVF_Highpass_12,
    SVF_Highpass_24,
    //...

    // Ladder filter modes:
    // LDR_Start = 4 << 24,
    LDR_Lowpass_6,
    LDR_Lowpass_12,
    LDR_Lowpass_18,
    LDR_Lowpass_24,
    //...


    Unknown
  };
  // hmm...
  // -maybe split the type into two parts: topology (svf, ladder, etc.), mode (lpf, hpf, etc.)
  // -maybe we can use a bitfield: 2 bits for the topology, 6 bits for the type within the selected
  //  topology, makes a total of 8 bits to specify the mode in a structured way. This is a 
  //  potential space optimization that may be done later. If done, it should not affect the API.
  //  ...maybe use the 1st 8 bits for topology - retrieve via shift+mask when needed

  /** Sets the filter up in terms of cutoff (or center) frequency as normalized radian frequency 
  and resonance in decibels. This parametrization is suitable when used to implement the filter
  opcodes in sfz. */
  void setupCutRes(Type type, float cutoff, float resonance);

  /** Sets the filter up in terms of gain (in decibels), frequency (normalized radian) and 
  bandwidth (in octaves). This parametrization is suitable when used for the equalizer opcodes in 
  sfz. */
  void setupGainFreqBw(Type type, float gain, float omega, float bw);

  //void initCoeffs();
  //void updateCoeffs();

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  void processFrame(float* L, float* R);
  void resetState();


protected:


  //-----------------------------------------------------------------------------------------------
  /** \name Data Structures. We use structs to define the coefficient sets and internal states
  of the different filter topologies and then make a union from the structs to represent either of
  these. Then, we declare a member of the union type to store our data. */

  using TCoef = float;
  using TSig  = RAPT::rsVector2D<float>;  
  // for stereo, preliminary. maybe use rsSimdVector<float, 2> if possible, else 
  // rsSimdVector<float, 4>...hmm...but that may increase the size of the struct. Maybe keep using
  // rsVector2D but use simd within the computation, if possible - but let's not do premature
  // optimizations...


  // todo: maybe templatize this class and use float for TPar, and rsfloat32x2 for TSig in the 
  // sampler

  struct OnePoleImpl
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
  struct BiquadImpl             // biquad filter, using DF2 (todo: try TDF1 -> smaller state)
  {
    TSig  x1, x2, y1, y2;       // state
    TCoef b0, b1, b2, a1, a2;   // coeffs
    void resetState() { x1 = x2 = y1 = y2 = TSig(0); }
    void initCoeffs() { b0 = TCoef(1); b1 = b2 = a1 = a2 = TCoef(0); }
    TSig getSample(const TSig& in)
    {
      TSig y = b0*in + b1*x1 + b2*x2 + a1*y1 + a2*y2;  // compute output
      x2 = x1; x1 = in; y2 = y1; y1 = y;               // update state
      return y;
    }
  };
  struct SvfImpl                // state variable filter (using ZDF)
  {
    TSig  s1, s2;               // state
    TCoef R;                    // damping(?) - what obout the freq-scaling and weights?
    // ...stuff to do...
  };
  struct LadderImpl             // ladder filter (using UDF)
  {
    TSig  y1, y2, y3, y4;       // state
    TCoef a, k;                 // filter and feedback coeffs
    TCoef c0, c1, c2, c3, c4;   // output gains
    // ...stuff to do...
  };
  union FilterImpl
  {
    FilterImpl() {}             // without it, msc complains - why?
    OnePoleImpl fo;             // first order
    BiquadImpl  bqd;            // biquad
    SvfImpl     svf;            // state variable filter
    LadderImpl  ldr;            // ladder
    //PhasorImpl psr;           // phasor filter
  };
  // ToDo: implement reset/getSample etc. also in StateVars, etc. all these structs should provide 
  // the same API, but implement it in a way that is suitable to the given filter topology.


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  Type type = Type::Bypass;
  FilterImpl impl;

  // ToDo: maybe include also a state-vector filter (maybe rename to phasor filter to avoid name 
  // clash in abbreviation with svf)
};

//=================================================================================================

/**  */

class WaveshaperCore
{

public:

  /*
  enum class Shape
  {
    None,
    Tanh,
    HardClip,
    SoftClipCubic,
    SoftClipHexic
    // ...etc.
  };
  */

  using Shape = DistortShape;

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

  void processFrame(float* L, float* R) 
  {
    *L = preGain * *L + dcOffset;
    *R = preGain * *R + dcOffset;
    using namespace RAPT;
    switch(shape)
    {
    case Shape::tanh: { *L = tanh(*L);   *R = tanh(*R);   } break; // todo: use rsTanh based on exp
    case Shape::clip: { *L = rsClip(*L); *R = rsClip(*R); } break;
    default: 
      break;
    }
    *L *= postGain;
    *R *= postGain;
  }
  // under construction
  // -when we remove the RAPT from rsClip we get a compile error: "no matching overload found" ->
  //  check if there's another non-templated rsClip somehwere else and remove it, if so
  // -use a branchless tanh approximation directly operating on rsFloat64x2
  // -maybe switch to rsFloat32x4 for audio-samples...but maybe that's not advantageous
  // -implement a block-based variant

protected:

  Shape shape = Shape::linear;
  float preGain, dcOffset, postGain;
  float shapePar1, shapePar2;            // meaning depends on chosen shape

  // what about oversampling..but maybe that should be reserved to an advanced variant? or maybe
  // that should apply to a region as a whole?
};

//=================================================================================================

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

class rsSamplerProcessors  // maybe make this a namespace
{

public:



  class Amplifier : public SignalProcessor
  {

  public:

    Amplifier() 
    { 
      type = DspType::Amplifier;
      params.reserve(4);
      addParameter(Opcode::volumeN);
      addParameter(Opcode::panN);
      addParameter(Opcode::widthN);
      addParameter(Opcode::positionN);
      // Having to pass a magic number to reserve() is bad and error-prone -> try to find a better
      // way. The number of parameters is actually known at compile time. Maybe use std::array 
      // instead of std::vector...hmm...but the number varies between the subclasses, so the array
      // could not be a baseclass member then...hmm...Maybe SignalProcessor could have an int as
      // template parameter - but no: then i can't treat them uniformly in the DSP chain
    }
    void prepareToPlay(double fs) override 
    { 
      core.setup(
        params[0].getValue(),
        params[1].getValue(),
        params[2].getValue(),
        params[3].getValue());
      // ToDo:
      // -get rid of the getValue() calls by allowing the params to convert to float
      // -it's not ideal that this code depends on the order, how we add the params in the 
      //  constructor - try to avoid that
    }
    void processFrame(float* L, float* R) override 
    {
      core.processFrame(L, R);
    }

    void processBlock(float* L, float* R, int N) override 
    {
      for(int n = 0; n < N; n++)
        processFrame(&L[n], &R[n]);
    }


  protected:

    AmplifierCore core;

  };



  class Filter : public SignalProcessor
  {

  public:

    Filter() 
    { 
      type = DspType::Filter;
      params.reserve(3);
      addParameter(Opcode::filN_type);
      addParameter(Opcode::cutoffN);
      addParameter(Opcode::resonanceN);   // in sfz, this is a gain in dB
      //addParameter(Opcode::FilterBandwidth);
    }
    void prepareToPlay(double fs) override 
    { 
      FilterType sfzType = (FilterType)(int)params[0].getValue();
      FilterCore::Type coreType = convertTypeEnum(sfzType);
      core.setupCutRes(
        coreType,
        params[1].getValue() * float(2*PI/fs),
        params[2].getValue());
      core.resetState();
    }
    void processFrame(float* L, float* R) override { core.processFrame(L, R); }
    void processBlock(float* L, float* R, int N) override 
    {
      for(int n = 0; n < N; n++)
        processFrame(&L[n], &R[n]);
    }

  protected:

    FilterCore::Type convertTypeEnum(FilterType sfzType)
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

    FilterCore core;
  };

  class Equalizer : public SignalProcessor
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
    void prepareToPlay(double fs) override 
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



  class WaveShaper : public SignalProcessor
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
    void prepareToPlay(double fs) override
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
  // todo: Let it have an argument that somehow specifies, how many of each type should be 
  // allocated. Maybe that could be a reference to the sfz-data itself or something derived from it




  /** A client can request a processor of the given type. If a processor of the desired type is 
  available, a pointer to it will returned. If not, a nullptr will be returned. The "client" will 
  typically be a RegionPlayer and call grabProcessor on noteOn. When no processor of the desired 
  type is available anymore, the calling code should probably forego the whole RegionPlayer. If
  the region can't be played correctly due to lack of resources, it should not play at all. What
  it certainly should not do is to just replace the non-available processor by a bypass dummy 
  processor because that could have really bad consequences: imagine a missing attenuation 
  processor. Regions are always played back either correctly or not at all but never wrongly. */
  SignalProcessor* grabProcessor(DspType type);

  /** This function should be called by the client when it doesn't need the processor anymore, For
  example, because the region for which it was used has stopped playing. The client returns the 
  processor to the pool so it becomes available again for playing other notes. */
  void repositProcessor(SignalProcessor* p);
  // maybe rename it to repositProcessor


  /** This is currently only meant to facilitate unit testing overload conditions. In such tests,
  we want a well defined and small number of filters to be available so we can simulate conditions
  where the engine is running out of filters in a controlled way. */
  void setMaxNumFilters(int newMax) { filters.init(newMax); }


  int getNumUsedFilters() const { return filters.getNumUsedItems(); }
  // also to facilitate testing


protected:

  using SP = rsSamplerProcessors;
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
  SignalProcessorPool processorPool;
  //ModulatorPool modulatorPool;
  //ConnectionPool connectionPool;
};










}} // namespaces

#endif