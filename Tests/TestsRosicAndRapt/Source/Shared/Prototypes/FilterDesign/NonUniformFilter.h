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


  T getSample(T x, T dt)
  {
    T bdt = pow(b, dt);

    // update scaler:
    s = a + bdt*s         //  Eq. 16, with w = 0

    // update state:
    y = (a*x + bdt*y) / abs(s);   //  Eq. 13
    // verify, if the scaler really has to be applied when updating the state or if it should be 
    // applied only to the returned output

    // 
    return y;
  }

protected:


  T y;     // state/output
  T a, b;  // coefficients

  T s = 1; // scaler
  //T p = 1; 
  //T q = 1; 

  // we 
};