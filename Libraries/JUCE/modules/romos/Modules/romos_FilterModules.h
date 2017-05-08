#ifndef romos_FilterModules_h
#define romos_FilterModules_h

#include "../Framework/romos_ModuleAtomic.h"
#include "../Algorithms/romos_FilterDesign.h"
#include "../Framework/romos_WorkArea.h"
#include "romos_ModuleDefinitionMacros.h"

namespace romos
{


  /** A first order lowpass filter in direct form 1 with 1 pole and 1 zero realizing the difference eqaution: 
  y[n] = b0[n]*x[n] + b1[n]*x[n-1] - a1[n]*y[n-1]. 
  As opposed to the general FirstOrderFilter with input pins for the coefficients, this filter computes its coefficients itself from the
  second input which defines the cutoff frequency. It uses a bilinear transform design equation. */
  class FirstOrderLowpass : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_2(FirstOrderLowpass);
  public:
    virtual void resetVoiceState(int voiceIndex);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    double *buffers;
  };



  /** A first order filter in direct form 1 with 1 pole and 1 zero realizing the difference eqaution: 
  y[n] = b0[n]*x[n] + b1[n]*x[n-1] - a1[n]*y[n-1]. */
  class FirstOrderFilter : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_4(FirstOrderFilter);
  public:
    virtual void resetVoiceState(int voiceIndex);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    double *buffers;
  };

  /** A biquad filter module in direct form 1 realizing the difference equation: 
  y[n] = b0[n]*x[n] + b1[n]*x[n-1] + b2[n]*x[n-2] - a1[n]*y[n-1] - a2[n]*y[n-2]. */
  class Biquad : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_6(Biquad);
  public:
    virtual void resetVoiceState(int voiceIndex);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    double *buffers;
  };

  /** Module to calculate biquad coefficients from a specification. The type is selected on the GUI and there are input pins for 
  characteristic frequency, gain, etc. 

  \todo implement modes:
  // LP6, LP12, HP6, HP12, LP6->HP6, LP6+HP6, BP(const skirt), Bell, 
  // HiShelf1 (1st order), HiShelf2 (2nd order), LoShelf1, LoShelf2
  // Allpass90 (1st order allpas), Allpass180 (2nd order allpass, with Q), BP(const peak), 

  */
  class BiquadDesigner : public ModuleWithParameters
  {
    CREATE_COMMON_DECLARATIONS_3(BiquadDesigner);
  public:
    enum modes
    {
      BYPASS = 0,

      // 1st order bilinear transform designs:
      LOWPASS_6_BILINEAR,
      HIGHPASS_6_BILINEAR,
      LOW_SHELF_1_BILINEAR, 
      HIGH_SHELF_1_BILINEAR,
      ALLPASS_1_BILINEAR,

      // 2nd order bilinear transform designs (RBJ cookbook):
      LOWPASS_12_BILINEAR,  
      HIGHPASS_12_BILINEAR,
      BANDPASS_CONST_SKIRT_BILINEAR,
      BANDPASS_CONST_PEAK_BILINEAR,
      BANDREJECT_BILINEAR,
      PEAK_BILINEAR,
      LOW_SHELF_2_BILINEAR,
      HIGH_SHELF_2_BILINEAR,
      ALLPASS_2_BILINEAR,
   

      //LOWPASS_6_IMPULSE_INVARIANT,
      //LOWPASS_6_PRESCRIBED_NYQUIST_GAIN,
      // etc...
    };

    virtual void resetVoiceState(int voiceIndex);
    virtual void parameterChanged(int index);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    double *oldParameters; // frequency, gain, Q
    double *oldOutputs;    // b0, b1, b2, a1, a2
    int    mode;
  };



  /** A Moog style ladder filter module. */
  class LadderFilter : public ModuleWithParameters
  {
    CREATE_COMMON_DECLARATIONS_4(LadderFilter); // x, Freq, Reso, GainCompensation
  public:

    enum filterModes
    {
      FLAT = 0,
      LP_6,
      LP_12,
      LP_18,
      LP_24,
      HP_6,
      HP_12,
      HP_18,
      HP_24,
      BP_12_12,
      BP_6_18,
      BP_18_6,
      BP_6_12,
      BP_12_6,
      BP_6_6,
      NUM_FILTER_MODES
    };

    enum saturationModes
    {
      NO_SATURATION = 0,
      LAST_STAGE,
      FEEDBACK,
      EACH_STAGE,
      NUM_SATURATION_MODES
    };

    virtual void resetVoiceState(int voiceIndex);
    virtual void parameterChanged(int index);

  protected:
    virtual void allocateMemory();
    virtual void freeMemory();

    double *outputs, *coeffs;
    double *oldParameters; // frequency, resonance, autogain
    int    filterMode, saturationMode;
  };




  // DirectFormFilter: y[n] = sum_j( bj[n] * x[n-j] ) - sum_k( ak[n] * y[n-k] ) - variable number of a, b coefficients
  // BiquadCascade, FirstOrderCascade - versions with identical and different coeffs vor each stage
  // DualBiquad/DualBiquadDesigner
  // LadderFilter, LatticeFilter, DualBiquad (with serial/parallel blend)
  // BiquadCascadeToDirectForm/DirectFormToBiquadCascade/BiquadToLattice/LatticeToBiquad coefficient converters

  // ModalFilterBank: per mode parameters: relative frequency, relative amplitude, relative decay time, phase-offset
  //                  audio inputs: fundamental/reference frequency, decay, maybe freq-offset/shift -> shifts all partials by a constant
  //                  

} 

#endif 
