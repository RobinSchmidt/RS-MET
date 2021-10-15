#ifndef rosic_SamplerProcessors_h
#define rosic_SamplerProcessors_h

namespace rosic
{


//=================================================================================================

/** A multimode filter that implements not only different filter frequency response types (like 
lowpass, highbpass, bandpass, etc.) but even completely differently structured filters. Depending 
on what mode has been chosen, the internal state and coefficient data may be interpreted in 
different ways.... */

class rsSamplerFilter
{

public:

  rsSamplerFilter()
  {

  }

  enum class Mode
  {
    // Biquad filter modes:
    BYPASS,
    LPF_6,
    HPF_6,

    // State variable filter modes:
    SVF_LPF_12,
    SVF_LPF_24,
    SVF_HPF_12,
    SVF_HPF_24,

    // Ladder filter modes:
    LDR_LPF_6,
    LDR_LPF_12,
    LDR_LPF_18,
    LDR_LPF_24
  };


protected:


  //-----------------------------------------------------------------------------------------------
  /** \name Data Structures */

  using TCoef = float;
  using TSig  = RAPT::rsVector2D<float>;  // for stereo
  // todo: maybe templatize this class and use float for TPar, and rsfloat32x2 for TSig in the 
  // sampler

  struct BiquadCoeffs
  {
    TCoef b0, b1, b2, a1, a2;
  };
  struct BiquadState
  {
    TSig x1, x2, y1, y2;
  };

  struct StateVarCoeffs
  {

  };
  struct StateVarState
  {

  };

  struct LadderCoeffs
  {
    TCoef a, k;                 // filter and feedback coeffs
    TCoef c0, c1, c2, c3, c4;   // output gains
  };
  struct LadderState
  {
    TSig y1, y2, y3, y4;
  };



  // make unions from these different structs

  union Coeffs
  {
    BiquadCoeffs   bqd;
    LadderCoeffs   ldr;
    StateVarCoeffs svf;
  };
  union State
  {
    BiquadState   bqd;
    LadderState   ldr;
    StateVarState svf;
  };


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  Mode   mode = Mode::BYPASS;
  //Coeffs coeffs;
  //State  state;


  // ToDo: maybe include also a state-vector filter (maybe rename to state phasor filter to avoid
  // name clash in abbreviation)

};

//=================================================================================================


}
#endif