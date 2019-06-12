#pragma once

/** Implements a one pole filter for non-uniformly sampled data  

...very incomplete...implements the lowpass case w=0 in which case all coeffs are real valued

References:
(1) High-Order Recursive Filtering of Non-Uniformly Sampled Signals for Image and Video Processing
http://inf.ufrgs.br/~eslgastal/NonUniformFiltering/Gastal_Oliveira_EG2015_Non-Uniform_Filtering.pdf

todo: 
-add highpass mode ...and maybe others, too?
-maybe we can implement this as subclass of the uniform filter?

*/

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


  T y = 0;         // output (g in the paper)
  T a = 1, b = 0;  // coefficients, a: feedforward, b: feedback (notation as in the paper)
                   // in more common dsp conventions a = b0, b = a1, b1 does not exist


  T x1 = 0; // previous input (f in the paper) - only needed for "piecewise resampling" method
  T s  = 1; // scaler (gamma in the paper) - only needed for "spatially-invariant scaling" method
  // todo: maybe split into two classes, one with the x-member, the other with s-member
  //T p = 1; 
  //T q = 1; 

  T w = 0;   // (normalized radian) cutoff frequency

  NormalizeMode normMode = NormalizeMode::noNormalization;
};