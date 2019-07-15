#pragma once

/** Implements a one pole filter for non-uniformly sampled data. This filter class is real-valued, 
so it supports only lowpass mode (todo: add highpass). One-pole bandpasses are necessarrily complex 
valued.....

References:
(1) High-Order Recursive Filtering of Non-Uniformly Sampled Signals for Image and Video Processing
http://inf.ufrgs.br/~eslgastal/NonUniformFiltering/Gastal_Oliveira_EG2015_Non-Uniform_Filtering.pdf

todo: 
-add highpass mode ...and maybe others, too? but how? maybe just highpass = in - lowpass? ...but 
 in this case, this should be better left to client code 
-maybe we can implement this as subclass of the uniform filter? hmm...probably no good idea */




template<class T>
class rsNonUniformOnePole
{

public:


  rsNonUniformOnePole() { reset(); }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setOmega(T newOmega);


  enum class NormalizeMode
  {
    noNormalization,
    spatiallyVariantScaling,
    piecewiseResampling
  };
  // maybe drag outside the class, so it can be used by the complex and high-order version, too
  // rsNonUniformFilterNormalization

  /** Sets the normalization mode for the getSample dispatcher method, which calls one the three
  versions  */
  void setNormalizationMode(NormalizeMode newMode)
  {
    normMode = newMode;
  }
  // todo: figure out, which normalization approach is most suitable in various use cases and add
  // that info to the comments



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes one output sample y[n] at a time from an incoming input sample x[n] and a 
  time-difference dt = t[n] - t[n-1] between the sample instant for x and the previous sample 
  instant for x[n-1]. It dispatches between various variations of the algorithm that use different
  approaches to normalize the output. The variations are getSampleNonNormalized, 
  getSampleSpatiallyInvariantScaled, getSamplePiecewiseResampled. Client code my also call any of 
  these directly, thereby bypassing the disaptcher (for optimization purposes). */
  T getSample(T x, T dt);

  /**  Implements none of the two normalization methods suggested in (1) - it's just
  the raw filter. But this raw filter is the one that actually agrees with the continuous impulse 
  response at the sample instants - so perhaps it's not such a bad choice after all. */
  T getSampleNonNormalized(T x, T dt);

  /** Computes a sample using the "spatially-variant scaling" normalization method from (1), 
  section 3.2.4. In our context, we may interpret the term "spatially-variant" actually as 
  "time-variant". The idea is to treat the input signal as if it were uniformly sampled, but 
  treating the filter as a time-variant system. What varies in this case is the output gain of the
  filter. */
  T getSampleSpatiallyVariantScaled(T x, T dt);
  // rename to TimeVariantScaled

  /** Computes a sample using the "piecewise resampling" normalization method from (1), 
  section 3.2.3. The idea is to compute the output from an uniformly resampled signal that is 
  obtained by linearly interpolating the incoming non-uniformly sampled signal. This is probably 
  more suitable for image processing than for audio processing - although, linear interpolation 
  may be a reasonably choice for control signals. */
  T getSamplePiecewiseResampled(T x, T dt);

  /** Resets the filter state. */
  void reset();


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  T y = 0;         // output (g in the paper)
  T a = 1, b = 0;  // coefficients, a: feedforward, b: feedback (notation as in the paper)
  // in more common dsp conventions a = b0, b = a1, b1 does not exist - maybe rename later
  // maybe use r and p for residue an pole - make its consistent with the complex version

  T x1 = 0; // previous input (f in the paper) - only needed for "piecewise resampling" method
  T s  = 1; // scaler (gamma in the paper) - only needed for "spatially-variant scaling" method
  // todo: maybe split into two classes, one with the x-member, the other with s-member
  //T p = 1; 
  //T q = 1; 

  T w = 0;   // (normalized radian) cutoff frequency ..maybe get rid...

  NormalizeMode normMode = NormalizeMode::noNormalization;
};

//=================================================================================================

/** Implements a complex-valued one-pole filter for non-uniformly sampled data. Higher order 
non-uniform filters (Butterworth, Gaussian, etc.) are created as a parallel connection of these
complex one-pole units via a partial fraction expansion of the transfer function. */

template<class T>
class rsNonUniformComplexOnePole
{

public:

  rsNonUniformComplexOnePole() { reset(); }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the feedforward (i.e. feed-in or input-gain) and feedback coefficients which correspond 
  to the resiude and the pole of this complex one-pole section. */
  void setCoeffs(const std::complex<T>& newFeedin, const std::complex<T>& newFeedback)
  {
    a = newFeedin;
    b = newFeedback;
  }


  void setReferenceOmega(T newOmega)
  {
    wr = newOmega;
  }
  // hmm....this is not yet really used properly - look up in the paper, how to deal with that


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  std::complex<T> getSampleNonNormalized(std::complex<T> x, T dt);

  std::complex<T> getSampleSpatiallyVariantScaled(std::complex<T> x, T dt);

  std::complex<T> getSamplePiecewiseResampled(std::complex<T> x, T dt);



 
  std::complex<T> getSample(std::complex<T> x, T dt)
  {
    //return getSampleNonNormalized(x, dt);
    //return getSampleSpatiallyVariantScaled(x, dt);
    return getSamplePiecewiseResampled(x, dt);

    //return getSamplePiecewiseResampled(x, dt) / rsAbs(s);  // test
  }
 
  // preliminary - todo: implement it as dispatcher - have a normalizationMode variable

  /** Accepts a real input sample and returns the real part of our complex output. The idea here is
  that in a setting with real inputs and outputs, each complex stage has a partner that produces 
  the complex conjugate of its output. In the final output, when all the outputs of the individual 
  one-poles are added together, the imaginary parts will cancel each other out and the real parts 
  will add up to twice the value of the output of the single stage. That means, we do not actually 
  have to compute the output of the partner stage by implementing it as an actual complex one-pole. 
  Instead, we can obtain the output of the two stages as twice the real part of one of those 
  stages. The multiplication by two is not done here - it's left to outside code for optimization 
  purposes (it needs to be done only once for the whole sum over all stages). */
  T getSampleReal(T x, T dt)
  {
    std::complex<T> z = getSample(std::complex<T>(x, T(0)), dt);
    return z.real();
  }

  /** Resets the filter state. */
  void reset();



  std::complex<T> getScaler() const { return s; }

protected:

  std::complex<T> y = T(0);            // maybe use g (as in the paper)
  std::complex<T> a = T(1), b = T(0);  
  // maybe rename to b0, a1 for consistency with other filters - or to r,p for residue,pole

  // maybe factor out into subclass(es):
  std::complex<T> x1 = 0; // for production, keep this in the outlying rsNonUniformFilterIIR

  // remove these for production - we don't use time-variant scaling:
  std::complex<T> s  = 1;
  T wr = 0;  // reference frequency for unit gain in time-variant scaling normalization

};

//=================================================================================================

/** Implements a high order infinite impulse response filter for non-uniformly sampled signals. */

template<class T>  // todo: have template parameters for signal, time and parameter
class rsNonUniformFilterIIR
{

public:


  rsNonUniformFilterIIR();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  ///** Sets the normalized radian frequency omega = 2*PI*frequency/sampleRate. */
  //void setOmega(T newOmega);

  /** Sets the cutoff frequency of the filter in Hz. */
  void setFrequency(T newFreq);

  /** Enumeration of the available filter modes. This is currently only a subset of the types
  available in rsPrototypeDesigner because we are currently restricted to allpole filters. For
  filters with zeros, i need to figure out, how to transform the analog prototype to digital
  because the impulse invariant transform (which is used here) is applicable only to allpole
  filters.
  ...actually, we'll get issues already when creating highpass and bandpass filters
  ..hmm..we'll see  */
  enum class ApproximationMethod
  {
    gaussian,
    bessel,
    butterworth,
    //chebychev,   // needs a ripple parameter
    papoulis,
    halpern,

    elliptic      // experimental - does not yet work
  };

  /** Sets the approximation method for the prototype filter (Butterworth, Bessel, etc.). */
  void setApproximationMethod(ApproximationMethod newMethod);

  /** Sets the order of the (prototype) filter - when we do bandpasses later, the order of the
  actual filter will be twice the prototype order, but that's not yet implemented. */
  void setOrder(int newOrder);


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  // getTransferFunction(z), getMagnitudeAt(w)



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes one output sample y[n] at a time from an incoming input sample x[n] and a 
  time-difference dt = t[n] - t[n-1] between the sample instant for x and the previous sample 
  instant for x[n-1]. dt is in seconds, i.e. 1/sampleRate for uniformly sampled signals. */
  T getSample(T x, T dt)
  {
    dt *= dtScaler;
    std::complex<T> y = fir[0]*x;  // accumulator for parallel (complex) filter outputs
    for(int i = 0; i < order; i++)
      y += onePoles[i].getSamplePiecewiseResampled(x, dt); // this seems to work best
      //y += onePoles[i].getSampleNonNormalized(x, dt);
      //y += onePoles[i].getSample(x, dt);
      //y += onePoles[i].getSampleSpatiallyVariantScaled(x, dt); 
         // nope! this doesn't work! normalizing should not be done per stage!
    return outScaler * y.real();
  }
  // could we also use getSampleReal and use a real accumulator? i think so - try it!
  // as an optimization, we may run only over half of the samples - the other half would just 
  // produce the complex conjugate - and at the end multiply everything by 2 - but the output of
  // the real stage, if present, should not be multiplied by 2

  /** Resets the internal state of the filter. */
  void reset()
  {
    for(int i = 0; i < order; i++)
      onePoles[i].reset();
  }


protected:

  /** Updates all filter coefficients according to the settings. */
  void updateCoeffs();


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  // array of complex one-pole filters:
  static const int maxOrder = 10;
  rsNonUniformComplexOnePole<T> onePoles[maxOrder]; 
  // maybe get rid - has redundancy with r,p - but we then need a y-array - but maybe not - we
  // to use the piecewise resampling method per stage and implement it within the stage

  // settings:
  int order = 1;
  T freq    = 0.25;  // 0.25 for halfband -> w = 2*PI*f/fs = PI/2 for fs = 1
  ApproximationMethod approxMethod = ApproximationMethod::butterworth;

  // experimental:
  T dtScaler = 1.0;
  //T fs = 1.0;


  T outScaler = 1.0;

  // prototype designer and buffers for partial fraction expansion routine:
  RAPT::rsPrototypeDesigner<T> protoDesigner;
  std::complex<T> p[maxOrder];     // prototype poles
  std::complex<T> z[maxOrder];     // prototype zeros
  std::complex<T> r[maxOrder];     // residues
  //std::complex<T> num[1] = { 1 };  // numerator of transfer function
  std::complex<T> num[maxOrder+1]; // numerator of transfer function
  std::complex<T> den[maxOrder+1]; // denominator of transfer function
  std::complex<T> fir[maxOrder+1]; // FIR part (polynomial)
  int muls[maxOrder];              // pole multiplicities (all 1 at the moment)
};