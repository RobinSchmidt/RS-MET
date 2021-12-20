#ifndef rosic_SamplerProcessors_h
#define rosic_SamplerProcessors_h

namespace rosic {
namespace Sampler {

//=================================================================================================

/** Enumeration of the different signal processor types that may be used in the definition of
instruments. What kinds of processors are used within a region is implicitly determined by the sfz 
opcodes, e.g. the presence of a FilterCutoff opcode dictates the presence of a filter within the 
respective region. In order to facilitating to build the DSP chain for a region player,
we also need an explicit representation of the DSP processor types. */
/*
enum SignalProcessorType
{
  SamplePlayer,

  // The modulators:
  // AmpEnv, FilterEnv, PitchEnv, AmpLFO, ...

  // The actual DSP processors:
  Filter,
  WaveShaper,

  Unknown
};
// moved into rosic_SamplerData.h
*/

/** Class to represent parameters of DSP processors and modulators. */
class Parameter
{

public:

  // Lifetime:
  //Parameter() {}  // maybe needed later to construct arrays
  Parameter(const char* name, float defaultValue)
  {
    this->name         = name;
    this->defaultValue = defaultValue;
    this->value        = defaultValue;
  }

  // Setup:
  void setValue(float newValue) { value = newValue; }
  // maybe later: setName, setDefaultValue

  // Inquiry:
  const char* getName() const { return name; }
  float getValue() const { return value; }
  float getDefaultValue() const { return defaultValue; }


protected:

  float value        = 0.f; 
  float defaultValue = 0.f;
  const char* name   = "\0";  // shall be assigned to some fixed global list of strings


  Opcode opcode = Opcode::Unknown;
  // todo: move PlaybackSetting::Type out of rsSamplerData, maybe rename to Opcode



  // what about min/max and mapping?

};

/** Baseclass for signal processors that can be applied to layers while they are the played back.
Subclasses can be various kinds of filters, equalizers, waveshapers, effects, etc. */
class SignalProcessor
{
public:

  // Inquiry:
  SignalProcessorType getType() { return type; }
  virtual int getNumParameters() const { return 0; }
  virtual const Parameter* getParameter() const { return nullptr; } // maybe relax constness later
  // check, if that's the right form of constness - we want to disable modfifying the parameter
  // itself but maybe reassign the pointer - yep, looks right:
  // https://stackoverflow.com/questions/21476869/constant-pointer-vs-pointer-to-constant





  // Processing:
  virtual void processFrame(rsFloat64x2& inOut) = 0;
  virtual void processBlock(rsFloat64x2* inOut, int N) = 0;
  virtual void resetState() = 0;
  virtual void resetSettings() = 0;

protected:

  // Setup:
  virtual void addParameter(const char* name, float defaultValue); // maybe should return an integer parameter index?



  SignalProcessorType type = SignalProcessorType::Unknown;
  //bool busy = false;
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
  void initCoeffs();


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

    BiquadVars bqd;
    StateVars  svf;
    LadderVars ldr;
  };
  // ToDo: implement reset/getSample etc also in StateVars, etc. all thes structs should provide 
  // the sample API, but implement it in a way that is suitable to the given filter topology.


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  Type type = Type::Bypass;
  FilterVars vars;

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
      addParameter("Cutoff",    1000.f);
      addParameter("Resonance", 0.f);
      // Maybe we should use sfz opcode names here and also use the ranges and default values
      // defined there. Maybe we should somewhere have a global enum with opcode identifiers, then
      // a table with their names, etc.
    }

    void processFrame(rsFloat64x2& inOut) override {}
    void processBlock(rsFloat64x2* inOut, int N) override {}
    void resetState() override { core.resetState(); }
    void resetSettings() override { core.initCoeffs(); }

  protected:

    rsSamplerFilter core;

  };


  class WaveShaper : public SignalProcessor
  {

  public:

    WaveShaper() { type = SignalProcessorType::WaveShaper; }
    void processFrame(rsFloat64x2& inOut) override {}
    void processBlock(rsFloat64x2* inOut, int N) override {}
    void resetState() override {}
    void resetSettings() override {}

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
  std::vector<SP::Filter>     filters;
  std::vector<SP::WaveShaper> waveShapers;

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