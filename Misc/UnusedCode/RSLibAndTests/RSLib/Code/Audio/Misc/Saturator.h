#ifndef RS_SATURATOR_H
#define RS_SATURATOR_H

namespace RSLib
{

/**

This is a class which applies saturation to a value asssumed to be nonnegative. The user can define
a saturation threshold t <= 1. Below the threshold, the output is just equal to the input. Above 
the threshold, a user-defined saturation function is applied to a shifted and scaled input value
and the output of the saturation function is shifted and scaled again - so, eventually:
y = t + s2 * saturate( (x-t)*s1 ) where s1 chosen such that the input domain t..1 maps to 0..1 and 
s2 is chosen such that the saturation-function's output range 0..1 maps back to t..1. For a typical
saturating behavior, you'll probably want to choose a saturation function that satisfies:
f(0)=0, f'(0)=1, f(inf)=1. f(0)=0 will ensure a matched function value at the threshold, f'(0)=0 
will ensure a matched derivative at the threshold and f(inf)=1 will ensure that the overall 
saturation value is unity. 

The main purpose of this class is to be used inside the rsSaturator class for saturating lower and 
upper halfwave separately, using one halfwave-saturator object for each halfwave (factoring out a 
halfwave-saturator class avoids some code duplication in the rsSaturator class).

*/

class RSLib_API rsHalfWaveSaturator
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsHalfWaveSaturator();


  /** \name Setup */

  /** Sets up the threshold t <= 1 below which the input/output mapping is the identity function.
  Above the threshold, the user-defined saturation function will be used. */
  void setThreshold(double newThreshold);

  /** Sets up the saturation function. You should pass a function-pointer to a function that takes
  a "double" parameter and has a "double" as return value. */
  void setSaturationFunction(double (*newFunction)(double));


  /** \name Audio Processing */

  /** Returns one output sample at a time. */
  RS_INLINE double getSample(double in);


protected:

  /** \name Data */

  double thresh;                  // threshold
  double (*saturate)(double);     // saturation function
  double inScale, outScale;       // scale factors for in- and output

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

RS_INLINE double rsHalfWaveSaturator::getSample(double in)
{
  if(in <= thresh)
    return in;
  else
  {
    in = (in-thresh)*inScale;        // shift/scale saturator input
    in = saturate(in);               // apply saturation
    return outScale*in + thresh;     // scale/shift saturator output
  }
}

//=================================================================================================

/**

This is a class for saturation waveshaping. It allows to define separate thresholds and saturation
functions for the lower and upper halfwave by using two rsHalfWaveSaturator objects for positive 
and negative halfwave separately. See documentation of rsHalfWaveSaturator.

ToDo: 
maybe instead of using function pointers, use the rsParametricSigmoid class

*/

class RSLib_API rsSaturator
{

public:

  /** \name Setup */

  /** Sets the threshold for the nonlinear range for the lower halfwave. */
  void setLowerThreshold(double newThreshold);

  /** Sets the threshold for the nonlinear range for the upper halfwave. */
  void setUpperThreshold(double newThreshold);

  /** Sets lower and upper threshold at once to the same value. */
  void setThresholds(double newThreshold);

  /** Sets the saturation function for the lower halfwave. */
  void setLowerSaturationFunction(double (*newFunction)(double));

  /** Sets the saturation function for the upper halfwave. */
  void setUpperSaturationFunction(double (*newFunction)(double));

  /** Sets lower and upper saturation function at once to the function. */
  void setSaturationFunctions(double (*newFunction)(double));


  /** \name Audio Processing */

  /** Returns one output sample at a time. */
  RS_INLINE double getSample(double in);


protected:

  /** \name Embedded objects */

  rsHalfWaveSaturator loSaturator, upSaturator;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

RS_INLINE double rsSaturator::getSample(double in)
{
  if(in >= 0.0)
    return upSaturator.getSample(in);
  else
    return -loSaturator.getSample(-in);
}

//=================================================================================================

/**

This class extends rsParametricSigmoid by re-expressing the saturation function y = f(x) in various 
ways, like, for example as a signal-dependent factor: y = f(x) = x * a(x) where a(x) = f(x)/x and 
then allowing a sidechain signal s to go into a(x) instead of x itself, such that effectively, we 
compute y = x * a(s) instead of y = x * a(x). There are different modes which amount to different
ways of rewriting y = f(x). The example above corresponds to the multiplicative mode, there's also
an additive mode: y = f(x) = x + a(x), etc.

*/

class RSLib_API rsSidechainSaturator : public rsParametricSigmoid
{

public:

  // application modes:
  enum modes
  {
    BYPASS = 0,      // y = x
    DIRECT1,         // y = f(s)
    DIRECT2,         // y = f(s) - f(s-x)
    ADDITIVE,        // y = x + a(s), a(s) = f(s) - s
    MULTIPLICATIVE,  // y = x * a(s), a(s) = f(s) / s

    NUM_MODES
  };
  // DIRECT2: idea: s is given by s = x + t, where t is an offset that was added to x...
  // todo: invent more modes, write proper comment, documenting the available modes


  /** \name Setup */

  /** Sets the saturation mode */
  void setMode(int newMode);

  /** \name Audio Processing */

  /** Returns one output sample at a time. You have to pass the actual input signal and the 
  sidechain signal. */
  RS_INLINE double getSample(double input, double sidechain);


protected:

  /** \name Data */

  int mode = DIRECT1;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

RS_INLINE double rsSidechainSaturator::getSample(double x, double s)
{
  double as;                                     // a(s)
  double fs = rsParametricSigmoid::getValue(s);  // f(s) 

  switch(mode)
  {
  case DIRECT1:  return fs;
  case DIRECT2:  return fs - rsParametricSigmoid::getValue(s-x);
  case ADDITIVE: 
  {
    as = fs - s;
    return x + as;   // can be collapsed into one line
  }
  case MULTIPLICATIVE: 
  {
    if(s == 0.0)
      as = 1.0;      // we assume the derivative of f(x) to be unity at x=0
    else
      as = fs / s;
    return x * as;
  }
  default: return x;
  };
}

}
#endif
