#pragma once

/** Implements a one pole filter for non-uniformly sampled data  

...very incomplete...implements the lowpass case w=0 in which case all coeffs are real valued

References:
(1) High-Order Recursive Filtering of Non-Uniformly Sampled Signals for Image and Video Processing
http://inf.ufrgs.br/~eslgastal/NonUniformFiltering/Gastal_Oliveira_EG2015_Non-Uniform_Filtering.pdf
*/

template<class T>
class rsNonUniformOnePole
{

public:

  rsNonUniformOnePole() { reset(); }

  void setOmega(T newOmega);

  /** Computes one output sample at a time. Takes the input sample x and a time-difference between
  the sample instant for x and the previous sample instant for x[n-1]. Implements the 
  "spatially-invariant scaling" normalization method from (1), section 3.2.4. */
  T getSample(T x, T dt);

  /** Like getSample, but implements the "piecewise resampling" normalization method from 
  (1), section 3.2.3. */
  T getSample2(T x, T dt);

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
};