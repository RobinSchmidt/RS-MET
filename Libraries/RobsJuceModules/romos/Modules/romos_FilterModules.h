#ifndef romos_FilterModules_h
#define romos_FilterModules_h

//-------------------------------------------------------------------------------------------------

/** A first order lowpass filter in direct form 1 with 1 pole and 1 zero realizing the difference 
eqaution:
y[n] = b0[n]*x[n] + b1[n]*x[n-1] - a1[n]*y[n-1].
As opposed to the general FirstOrderFilter with input pins for the coefficients, this filter 
computes its coefficients itself from the second input which defines the cutoff frequency. It uses
a bilinear transform design equation. 
// todo: make it switchable to impulse variant via a gui parameter */

class FirstOrderLowpass : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_2(FirstOrderLowpass);
public:
  virtual void resetVoiceState(int voiceIndex);
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  double *buffers;
};
class FirstOrderLowpassTypeInfo : public ModuleTypeInfo
{
public:
  FirstOrderLowpassTypeInfo() {
    shortName    = "LPF6";
    fullName     = "FirstOrderLowpass";
    description  = "Applies a first order lowpass filter (with 6 dB/oct slope) to the signal";
    category     = "Filters";
    createModule =  []()->Module* { return new FirstOrderLowpass; };
  }
};

//-------------------------------------------------------------------------------------------------

/** A first order filter in direct form 1 with 1 pole and 1 zero realizing the difference 
equation:
y[n] = b0[n]*x[n] + b1[n]*x[n-1] - a1[n]*y[n-1]. */

class FirstOrderFilter : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_4(FirstOrderFilter);
public:
  virtual void resetVoiceState(int voiceIndex);
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  double *buffers;
};
class FirstOrderFilterTypeInfo : public ModuleTypeInfo
{
public:
  FirstOrderFilterTypeInfo() {
    shortName    = "1p1z";
    fullName     = "FirstOrderFilter";
    description  = "One pole, one zero filter. Realizes y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1]";
    category     = "Filters";
    createModule =  []()->Module* { return new FirstOrderFilter; };
  }
};

//-------------------------------------------------------------------------------------------------

/** A biquad filter module in direct form 1 realizing the difference equation:
y[n] = b0[n]*x[n] + b1[n]*x[n-1] + b2[n]*x[n-2] - a1[n]*y[n-1] - a2[n]*y[n-2]. */

class Biquad : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_6(Biquad);
public:
  virtual void resetVoiceState(int voiceIndex);
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  double *buffers;
};
class BiquadTypeInfo : public ModuleTypeInfo
{
public:
  BiquadTypeInfo() {
    shortName    = "Biquad";
    fullName     = "Biquad";
    description  = "Biquad (2 pole, 2 zero) filter. Realizes y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]";
    category     = "Filters";
    createModule =  []()->Module* { return new Biquad; };
  }
};

//-------------------------------------------------------------------------------------------------

/** Module to calculate biquad coefficients from a specification. The type is selected on the GUI 
and there are input pins for characteristic frequency, gain, etc.

\todo implement modes:
// LP6, LP12, HP6, HP12, LP6->HP6, LP6+HP6, BP(const skirt), Bell,
// HiShelf1 (1st order), HiShelf2 (2nd order), LoShelf1, LoShelf2
// Allpass90 (1st order allpas), Allpass180 (2nd order allpass, with Q), BP(const peak), */

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
class BiquadDesignerTypeInfo : public ModuleTypeInfo
{
public:
  BiquadDesignerTypeInfo() {
    shortName    = "BqdDsgn";
    fullName     = "BiquadDesigner";
    description  = "Computes biquad coefficients from a specification";
    category     = "Filters";
    createModule =  []()->Module* { return new BiquadDesigner; };
  }
};

//-------------------------------------------------------------------------------------------------

/** Under construction.
A Moog style ladder filter module. */

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

class LadderFilterTypeInfo : public ModuleTypeInfo
{
public:
  LadderFilterTypeInfo() {
    shortName    = "Ladder";
    fullName     = "LadderFilter";
    description  = "A Moog style ladder filter with multiple modes";
    category     = "Filters";
    createModule =  []()->Module* { return new LadderFilter; };
  }
};

//-------------------------------------------------------------------------------------------------

// DirectFormFilter: y[n] = sum_j( bj[n] * x[n-j] ) - sum_k( ak[n] * y[n-k] ) - variable number of 
// a, b coefficients BiquadCascade, FirstOrderCascade - versions with identical and different 
// coeffs vor each stage DualBiquad/DualBiquadDesigner LadderFilter, LatticeFilter, DualBiquad 
// (with serial/parallel blend) BiquadCascadeToDirectForm/DirectFormToBiquadCascade/
// BiquadToLattice/LatticeToBiquad coefficient converters

// ModalFilterBank: per mode parameters: relative frequency, relative amplitude, 
// relative decay time, phase-offset
// audio inputs: fundamental/reference frequency, decay, 
// maybe freq-offset/shift -> shifts all partials by a constant (->FM)
// maybe use the StateVectorFilter for the implementation

#endif 
