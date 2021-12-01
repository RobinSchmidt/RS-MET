#ifndef rosic_SamplerProcessors_h
#define rosic_SamplerProcessors_h

namespace rosic {
namespace Sampler {


//=================================================================================================

/** Baseclass for signal processors that can be applied to layers while they are the played back.
Subclasses can be various kinds of filters, equalizers, waveshapers, effects, etc. */

class SignalProcessor
{
public:
  virtual void processFrame(rsFloat64x2& inOut) = 0;
  virtual void processBlock(rsFloat64x2* inOut, int N) = 0;
  virtual void resetState() = 0;
  virtual void resetSettings() = 0;
};


//=================================================================================================
// The DSP classes below are meant to be used in the sampler, but they are nevertheless implemented
// as pure DSP classes without the infrastructural aspect, i.e. without being a subclass of 
// Sampler::SignalProcessor. This has been done in order to facilitate dragging them out of the 
// Sampler sub-namespace to make them available in other contexts as well. The infrastructure is 
// provided by boilerplate classes that derive from SignalProcessor and the actual core DSP class 
// via multiple inheritance.

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
    Lowpass_6,
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
  //  topology, makes a total of 8 bits to specify the mode in a structured way


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


}} // namespaces

#endif