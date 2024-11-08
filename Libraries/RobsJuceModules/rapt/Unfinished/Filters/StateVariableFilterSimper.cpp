

template<class T>
void rsStateVariableFilterSimper<T>::setupFromBiquad(T b0, T b1, T b2, T a1, T a2)
{
  //rsError("This is under construction. It doesn't work yet!");

  // When done, this function should implement the conversion formulas from biquad to SVF coeffs
  // given here on page 8:
  //
  //   https://cytomic.com/files/dsp/SvfLinearTrapezoidalSin.pdf
  //
  // Note that the paper uses the convention of writing DF biquad transfer functions as:
  //
  //           a0 + a1 z^-1 + a2 z^-2
  //   H(z) = ------------------------
  //           1  - b1 z^-1 - b2 z^-2
  //
  // whereas I usually use the convention:
  //
  //           b0 + b1 z^-1 + b2 z^-2
  //   H(z) = ------------------------
  //           1  + a1 z^-1 + a2 z^-2
  //
  // i.e. with a- and b-coeffs swapped and a sign inversion on the numerator coeffs. The function
  // parameters are supposed to be given in my convention. ...TBC...


  // ToDo:
  // Verify everything in numerical tests. Could the argument of the sqrt become negative? What 
  // then? Maybe use absolute value of argument and then flip the sign of all mixing coeffs? Or 
  // maybe we need to use complex arithmetic? ...but nah! Or well maybe. The formulas use only the
  // quotient and the product of s1 and s2 - if they are complex conjugates, we may end up with
  // real coeffs.


  // Intermediate variables:
  using Complex = std::complex<T>;
  Complex t1 = -1 + a1 + a2;          // Argument of first square root
  Complex t2 = -1 - a1 + a2;          // Argument of second square root
  Complex s1 = sqrt(t1);              // Square root in the denominator of formula for g
  Complex s2 = sqrt(t2);              // Square root in the numerator of formula for g
  Complex qc = s1 / s2;               // Quotient
  Complex pc = s1 * s2;               // Product
  T q = real(qc);                     // imag(qc) should be zero anyway, I guess?
  T p = real(pc);                     // ..same for pc

  // Compute intermediate variables g,k:
  T g = -q;
  T k = (1 + a2) / p;                 //  Solution 1 seems to be the right one

  // Test:
  //g = -g; k = -k; // Test - this is the 2nd solution


  // This assigns our a-coeffcient member variables:
  calcFilterCoeffs(g, k);


  // These formulas use the a1,a2 function parameters, not our member variables:
  m0 = (b0 - b1 + b2) / (1 + a1 - a2);
  m1 = 2*(b0 - b2)    / p;
  m2 = (b0 + b1 + b2) / (1 - a1 - a2);
}


//=================================================================================================
/*

Notes:

- On page 6, there's also the formula "k = 1/Q = 2 - 2*res". Does that mean there could be a 
  "resonance" parameter instead of Q? Figure out! Yes - this seems to be the right formula.
  So we have Q = 1 / (2 - 2*R), R = (1/Q - 2) / (-2) = (2 - 1/Q) / 2


ToDo:

- Figure out and document what the coefficients and intermediate variables mean. Looking at 
  the scribble on the front page on the paper, it seems like k is the feedback factor after the
  1st integrator? And the a1, a2 are the gains of the integrators? And g is affecting them?
  v0 is the input voltage, v1, v2 the voltages after the 1st and 2nd integrator stage 
  representing bandpass and lowpass output? This can be inferred from the mixing coeffs 
  m0,m1,m2. They are 0,1,0 for bandpass and 0,0,1 for lowpass. I guess, the ic1eq, ic2eq in the
  paper (that I have just called i1, i2 here) stand for something like "current into capacitor 1" 
  but somehow "equalized" as in "normalized"? Ah - see here, page 5:
  https://cytomic.com/files/dsp/OnePoleLinearLowPass.pdf
  he talks about "equivalent current"

- Maybe have two template parameters TSig, TPar as in the other filters. I think,
  v0,v1,v2,v3,ic1eq,ic2eq must all be TSig, a1,a2,a3,m1,m2,m3 must be TPar

- Maybe move implementation into .cpp file ...but maybe not.

- Figure out the z-domain transfer function and implement a function 
  getTransferFunctionAt(rsComplex<TPar> z)

- Adjust the order of the modes in the enum to be the same as in the RBJ filter [DONE]...but RBJ
  has two bandpass variants. I think, this here is a const skirt gain bandpass. Maybe to obtain 
  const peak gain behavior, we just need to scale by k = 1/Q? [YES - DONE - seems OK]. The RBJ 
  filters are also missing a "peak" filter in the sense meant here. I think, it's just a 
  resonator? If so, try to introduce it in the RBJ filters as well. Maybe rename the mode to
  "Reson" or "Resonator". Maybe rename "Notch" to "Bandreject" for consistency. Or maybe call them
  Bandstop everywhere.  https://en.wikipedia.org/wiki/Band-stop_filter

- Maybe implement also this stuff:
    https://cytomic.com/files/dsp/SvfLinearTrapAllOutputs.pdf
  and look into the other papers here:
    https://cytomic.com/technical-papers/
  This one:
    https://cytomic.com/files/dsp/SvfLinearTrapezoidalSin.pdf
  has also a part where the g,k,m0,m1,m2 coeffs are computed from direct from biquad coeffs. There
  are also hints for how to compute the transfer function. It's all buried in the Mathematica 
  code, though.




*/