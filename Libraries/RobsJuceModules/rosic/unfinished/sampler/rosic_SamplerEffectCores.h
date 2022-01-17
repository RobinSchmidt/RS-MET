#ifndef rosic_SamplerEffectCores_h
#define rosic_SamplerEffectCores_h

namespace rosic {
namespace Sampler {

//=================================================================================================
// The DSP classes below are meant to be used in the sampler, but they are nevertheless implemented
// as pure DSP classes without the infrastructural aspect, i.e. without being a subclass of 
// Sampler::Effect. This has been done in order to facilitate dragging them out of the Sampler 
// sub-namespace to make them available in other contexts as well. The infrastructure is 
// provided by boilerplate classes that derive from Effect and the actual core DSP class 
// having a member called "core".


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
  // -Maybe get rid and use the same enum as the sfzCodeBook. The conversion may just be verbose 
  //  cruft.

  /** Sets the filter up in terms of cutoff (or center) frequency as normalized radian frequency 
  and resonance in decibels. This parametrization is suitable when used to implement the filter
  opcodes in sfz. */
  void setupCutRes(Type type, float cutoffOmega, float resonance);

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







}
}


#endif
