#ifndef jura_ParameterMappers_h
#define jura_ParameterMappers_h

//=================================================================================================

/** Objects of (subclasses of) this class are used in Parameter to do the mapping between the 
normalized range 0..1 and the actual parameter range. Subclasses can implement different mapping 
functions (shapes) to go from some minimum to some maximum parameter value. */

class JUCE_API rsParameterMapper
{

public:

  /** Constructor. */
  rsParameterMapper() {}

  /** Destructor. */
  virtual ~rsParameterMapper() {}

  /** Override this function in your subclass to map a normalized input value (in the range 0..1) 
  to the corresponding actual parameter value (in the range min..max).  */
  virtual double map(double normalizedValue) const = 0;

  /** Override this function in your subclass to map the actual parameter value (in 
  the range min..max) to the corresponding normalized value (in the range 0..1). It should be the
  inverse function of map. For example, if map is exp then unmap should be log. */
  virtual double unmap(double value) const = 0;

  /** Sets up the range for the (mapped, actual) parameter value. */
  virtual void setRange(double newMin, double newMax)
  {
    jassert(newMin < newMax);  // Must be strictly less for division by (max-min) in unmap
    min = newMin;
    max = newMax;
  }

  /** Returns the minimum value. */
  inline double getMin() const { return min; }

  /** Returns the maximum value. */
  inline double getMax() const { return max; }

protected:

  double min = 0, max = 1;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapper)
};

//=================================================================================================

/** Subclass of rsParameterMapper for linear mapping. This is appropriate for parameters that are
either intrinsically perceived on a linear scale (such as a phase between -180..+180) or a 
perceptually linearized measure of some quantity (such as decibels or semitones). */

class JUCE_API rsParameterMapperLinear : public rsParameterMapper
{
public:
  rsParameterMapperLinear() {}
  double   map(double x) const override { return min + (max-min) * x; }
  double unmap(double y) const override { return (y-min) / (max-min); }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperLinear)
};

//=================================================================================================

/** Subclass of rsParameterMapper for identity mapping. This mapper should be used when the range
must be unrestricted, i.e. go from -inf to +inf. In this case, the linear mapper would produce
NaN. */

class JUCE_API rsParameterMapperIdentity : public rsParameterMapper
{
public:
  rsParameterMapperIdentity() { min = -INF; max = INF; }
  double   map(double x) const override { return x; }
  double unmap(double y) const override { return y; }
  void setRange(double newMin, double newMax) override {} // do nothing
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperIdentity)
};

//=================================================================================================

/** Subclass of rsParameterMapper for exponential mapping. This is appropriate for parameters that 
are perceived linearly on a logarithmic scale, but nevertheless are entered as raw values by the 
user. An example would be a frequency parameter with a range of 20..20000. With exponential mapping, 
equal differences in slider value translate to equal multiplication factors for the parameter 
value. */

class JUCE_API rsParameterMapperExponential : public rsParameterMapper
{
public:
  rsParameterMapperExponential() {}
  double   map(double x) const override { return min * exp(x*(log(max/min))); }
  double unmap(double y) const override 
  {
    return jlimit(0.0, 1.0, log(y/min) / (log(max/min)));
    //return log(y/min) / (log(max/min)); 
  }

  // \todo: optimize by precomputing log(max/min), 1/log(max/min), 1/min

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperExponential)
};

//=================================================================================================

/** A mapper based on the hyperbolic sine function, using y = a * sinh(b*x) where x is a value
between -1 and +1 (derived from the normalized 0..1 parameter p as x = 2*p-1). This function is 
suitable for parameters that should be mapped exponentially but should nevertheless be bipolar.
An example would be a frequency between -20000 and +20000 Hz. You can set up a shape parameter 
which controls the trade-off between precision around zero and high frequency precision (in the 
case of the freq-example). This shape parameter is actually the "b" in the formula. The "a" will 
then be determined by "b" and the max value (i.e. 20000). It currently supports only ranges that 
are symmetric around zero, i.e. min should be -max. */

class JUCE_API rsParameterMapperSinh : public rsParameterMapper
{
public:

  rsParameterMapperSinh(double minValue, double maxValue, double shape)
  {
    b = shape;
    setRange(minValue, maxValue); // updates the a-coeff
  }

  double map(double x) const override
  { 
    x = 2*x - 1;           // 0..1 to -1..+1
    return a * sinh(b*x);  // -max..max
      // maybe generalize to y = a * sinh(b*(x+c)) + d for unsymmetric ranges, etc. 
  } 

  double unmap(double y) const override
  { 
    y = asinh(y/a) / b;    // -1..+1
    return 0.5 * (y+1);    //  0..1
  }

  /** The range must be symmetrical around 0: newMin == -newMax. */
  void setRange(double newMin, double newMax) override
  {
    jassert(newMin == -newMax); // supports currently only 0-centered symmetric mapping
    rsParameterMapper::setRange(newMin, newMax);
    updateCoeffs();
  }

  void setShape(double newShape)
  {
    jassert(newShape > 0);
    b = newShape;
    updateCoeffs();
  }

  // todo: provide a helper function that computes the coefficients a,b, such that:
  //   a * sinh(b*x1) = y1
  //   a * sinh(b*x2) = y2
  // for some user-given input points (x1,y1),(x2,y2). This leads to:
  //   a = y1/sinh(b*x1) = y2/sinh(b*x2) -> y1/y2 = sinh(b*x1) / sinh(b*x2) -> solve for b...

protected:

  void updateCoeffs()
  {
    a = max / sinh(b);
  }

  double a = 1, b = 1;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperSinh)
};

//=================================================================================================

/** Parameter mapper based on the hyperbolic tangent. ...the code actually exactly parallels the 
sinh mapper - maybe we can avoid the duplication by refactoring? ...maybe it needs to use function 
pointers to sinh/asinh and tanh/atanh respectively ...and maybe that can then be generalized.
maybe factor out a class rsParameterMapperBipolar that is defined via the function:

y = a * f(b*x)

for some function f wich can be sinh, tanh, etc. and that is defined by a function pointer. 
subclasses then just set the function-pointer in the constructor (or actually 2 function pointers, 
for forward (sinh/tanh) and backward (asinh/atanh) mapping */

class JUCE_API rsParameterMapperTanh : public rsParameterMapper
{
public:

  rsParameterMapperTanh(double minValue, double maxValue, double shape)
  {
    b = shape;
    setRange(minValue, maxValue); // updates the a-coeff
  }

  double map(double x) const override
  { 
    x = 2*x - 1;           // 0..1 to -1..+1
    return a * tanh(b*x);  // -max..max
      // maybe generalize to y = a * tanh(b*(x+c)) + d for unsymmetric ranges, etc. 
  } 

  double unmap(double y) const override
  { 
    y = atanh(y/a) / b;    // -1..+1
    return 0.5 * (y+1);    //  0..1
  }

  /** The range must be symmetrical around 0: newMin == -newMax. */
  void setRange(double newMin, double newMax) override
  {
    jassert(newMin == -newMax); // supports currently only 0-centered symmetric mapping
    rsParameterMapper::setRange(newMin, newMax);
    updateCoeffs();
  }

  void setShape(double newShape)
  {
    jassert(newShape > 0);
    b = newShape;
    updateCoeffs();
  }

protected:

  void updateCoeffs()
  {
    a = max / tanh(b);
  }

  double a = 1, b = 1;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperTanh)
};

//=================================================================================================

/* Elan's rational mapping function, see here: https://www.desmos.com/calculator/xaklfkriac */

class JUCE_API rsParameterMapperRational : public rsParameterMapper
{
public:

  rsParameterMapperRational(double minValue, double maxValue, double shape)
  {
    tension = shape;
    rsParameterMapper::setRange(minValue, maxValue);
  }

  double curve(double t, double v) const
  {
    double tv = t*v;
    return (tv-v)/(2*tv-t-1);
  }

  /** Shape must be between -1 and +1, negative for log, positive for exp */
  void setShape(double newShape) { tension = newShape;}

  double map(double x) const override 
  {    
    return jmap(curve(tension, x), min, max);
  }
  double unmap(double y) const override 
  { 
    return curve(-tension, jmap(y, min, max, 0.0, 1.0));
  }

protected:

  double tension;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperRational)
};

//=================================================================================================

/* Bipolar version of the rational mapping function. */

class JUCE_API rsParameterMapperRationalBipolar : public rsParameterMapper
{
public:

  rsParameterMapperRationalBipolar(double minValue, double maxValue, double shape)
  {
    tension = shape;
    rsParameterMapper::setRange(minValue, maxValue);
  }

  double s_curve(double t, double v) const
  {
    double tv = t*v;

    if (v < 0.5)
      return (tv-v) / (4*tv-t-1);

    v += -0.5;
    t *= -0.5;
    tv = t*v;
    return (tv-v*0.5) / (4*tv-t-0.5) + 0.5;
  }

  /** Shape must be between -1 and +1, negative for log, positive for exp */
  void setShape(double newShape) { tension = newShape; }

  double map(double x) const override
  {
    return jmap(s_curve(tension, x), min, max);
  }
  double unmap(double y) const override
  {
    return s_curve(-tension, jmap(y, min, max, 0.0, 1.0));
  }

protected:

  double tension;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperRationalBipolar)
};

#endif