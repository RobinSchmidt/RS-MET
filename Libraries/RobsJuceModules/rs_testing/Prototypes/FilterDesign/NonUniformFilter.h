#pragma once

/** Implements a one pole filter for non-uniformly sampled data  

...very incomplete...implements the lowpass case w=0 in which case all coeffs are real valued

References:
High-Order Recursive Filtering of Non-Uniformly Sampled Signals for Image and Video Processing
http://inf.ufrgs.br/~eslgastal/NonUniformFiltering/Gastal_Oliveira_EG2015_Non-Uniform_Filtering.pdf
*/

template<class T>
class rsNonUniformOnePole
{

public:

  rsNonUniformOnePole() { reset(); }

  void setOmega(T newOmega);

  T getSample(T x, T dt);

  void reset();


protected:


  T y = 0;         // state/output
  T a = 1, b = 0;  // coefficients

  T s = 1; // scaler
  //T p = 1; 
  //T q = 1; 

  T w = 0;   // (normalized radian) cutoff frequency

  // we 
};