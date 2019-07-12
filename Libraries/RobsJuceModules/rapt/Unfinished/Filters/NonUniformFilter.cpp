template<class T>
void rsNonUniformOnePole<T>::setOmega(T newOmega)
{
  w = newOmega;
  T dummy;       // will be set to zero
  RAPT::rsFirstOrderFilterBase<T, T>::coeffsLowpassIIT(w, &a, &dummy, &b);
}

template<class T>
T rsNonUniformOnePole<T>::getSample(T x, T dt)
{
  typedef NormalizeMode NM;
  switch(normMode)
  {
  case NM::noNormalization:         return getSampleNonNormalized(x, dt);
  case NM::spatiallyVariantScaling: return getSampleSpatiallyVariantScaled(x, dt);
  case NM::piecewiseResampling:     return getSamplePiecewiseResampled(x, dt);
  default: { rsError("unknown normalization mode"); return T(0); }
  }
}

template<class T>
T rsNonUniformOnePole<T>::getSampleNonNormalized(T x, T dt)
{
  T bdt = pow(b, dt);  // b^dt
  y = a*x + bdt*y;
  return y;
}
// maybe inline this function


template<class T>
T rsNonUniformOnePole<T>::getSampleSpatiallyVariantScaled(T x, T dt)
{
  T bdt = pow(b, dt);    // b^dt - optimize: bdt = exp(log(b) * dt) where log(b) is precomputed
  s = a + bdt*s;         // update scaler via Eq. 16 with w = 0 (normalize gain at DC)
  y = a*x + bdt*y;       // update state
  return y / rsAbs(s);   // apply re-normalization
}
/*
// ...could it be that the update of s must be done *after* the output is computed? see Eq 13: 
// gamma_0 = a should used for the 0-th sample - let's try it:
template<class T>
T rsNonUniformOnePole<T>::getSampleSpatiallyVariantScaled(T x, T dt)
{
  T bdt = pow(b, dt);    // b^dt - optimize: bdt = exp(log(b) * dt) where log(b) is precomputed
  y     = a*x + bdt*y;   // update state
  T out = y / rsAbs(s);  // apply re-normalization
  s     = a + bdt*s;     // update scaler via Eq. 16 with w = 0 (normalize gain at DC)
  return out;
}
// ...hmmm - well - the version above with pre-update has the "better looking" step response
*/

template<class T>
T rsNonUniformOnePole<T>::getSamplePiecewiseResampled(T x, T dt)
{
  T bdt = pow(b, dt); // b^dt - express exp(log(b) * dt), precompute log(b)

  // some coeffs r0,r1 that can be precomputed:
  T bm1 = b-1;
  T r0  = bm1*bm1 / (a*b);
  T r1  = a / bm1;

  // compute additional compensation term:
  T R   = (bdt-1)/(r0*dt);
  T Phi = (R-r1*b)*x - (R-r1*bdt)*x1;
  //T test = (R-r1*bdt);  // multiplier for x1 - goes to zero when dt goes to 1

  // update state and return output:
  x1 = x;
  y  = a*x + bdt*y + Phi;
  return y;
}
// maybe allow the user to set a scaler (between 0...1) for Phi - so we can fade between 
// no-normalization (good for impulse response) and full normalization (good for step response)

template<class T>
void rsNonUniformOnePole<T>::reset()
{
  y  = T(0);

  x1 = T(0);

  //x1 = T(1000);  
  // test: doesn't seem to make a difference for imp-resp or step-resp. It seems like the
  // factor (R-r1*bdt) that multiplies x1 is always zero for the first sample -> it goes to zero
  // when dt goes to 1.

  //s = a;              // step response looks wrong...
  s = a / (T(1) - b);   // ...better
  // the general formula is Eq18: s = a / (1 - b * exp(j*wr)) where wr is the reference frequency
  // where we want unit gain, which is zero in this (lowpass) case
  // there's an alternative formula: s = a ...figure out, which one should be used
  // ...comparing results with a uniform 1-pole filter, it seems like s = a/(1-b) is the correct 
  // one - in this case, s comes out as 1 - but probably only in the special case of a lowpass

  // Looking at Eq.15, it seems, this gamma term (our s) is just the sum over all the 
  // impulse-response samples encountered so far when wr = 0. When wr != 0, it's the sum over the 
  // modulated impulse response
}

//=================================================================================================

template<class T>
std::complex<T> rsNonUniformComplexOnePole<T>::getSampleNonNormalized(std::complex<T> x, T dt)
{
  std::complex<T> bdt = pow(b, dt);  // b^dt
  y = a*x + bdt*y;
  return y;
}

template<class T>
std::complex<T> rsNonUniformComplexOnePole<T>::getSampleSpatiallyVariantScaled(
  std::complex<T> x, T dt)
{
  // code copied from rsNonUniformOnePole<T>::getSampleSpatiallyVariantScaled
  // still needs adaptions:

  std::complex<T> bdt = pow(b, dt);  // b^dt - optimize: bdt = exp(log(b) * dt) where log(b) can be precomputed
  s = a + bdt*s;                     // update scaler, Eq. 16, with w = 0 -> todo: generalize for w != 0
  y = a*x + bdt*y;
  return y / rsAbs(s);

  // we may need a function that computes/updates s but doesn't apply it (i.e. doesn't divide by it
  // at the output) - we may want to apply one single scaler for a whole bank of filters....
}

template<class T>
std::complex<T> rsNonUniformComplexOnePole<T>::getSamplePiecewiseResampled(std::complex<T> x, T dt)
{
  // code copied from rsNonUniformOnePole<T>::getSamplePiecewiseResampled but here, all variables
  // are complex - todo: verify in paper, if this is correct


  std::complex<T> bdt = pow(b, dt); // b^dt - express exp(log(b) * dt), precompute log(b)

  // some coeffs r0,r1 that can be precomputed:
  std::complex<T> bm1 = b - T(1);
  std::complex<T> r0  = bm1*bm1 / (a*b);
  std::complex<T> r1  = a / bm1;

  // compute additional compensation term:
  std::complex<T> R   = (bdt-T(1))/(r0*dt);
  std::complex<T> Phi = (R-r1*b)*x - (R-r1*bdt)*x1;

  // update state and return output:
  x1 = x;
  y  = a*x + bdt*y + Phi; // maybe scale Phi by value in 0..1 for "partial renormalization"
  return y;
}

template<class T>
void rsNonUniformComplexOnePole<T>::reset()
{
  x1 = T(0);
  y  = T(0);
  //s = a;
  std::complex<T> j(T(0), T(1));   // imaginary unit
  s = a / (T(1) - b * exp(j*wr));  // ...verify...
}

//=================================================================================================

template<class T>
rsNonUniformFilterIIR<T>::rsNonUniformFilterIIR()
{
  for(int i = 0; i < maxOrder; i++)
    muls[i] = 1;  // pole multiplicities are all 1
  setOrder(8);    // just for test
  updateCoeffs();
}

template<class T>
void rsNonUniformFilterIIR<T>::setFrequency(T newFreq)
{
  freq = newFreq;
  updateCoeffs();
}

template<class T>
void rsNonUniformFilterIIR<T>::setOrder(int newOrder)
{
  rsAssert(order >= 1 && order <= maxOrder, "order out of range");
  order = newOrder;
  protoDesigner.setOrder(order);
  updateCoeffs();
}

template<class T>
void rsNonUniformFilterIIR<T>::setApproximationMethod(ApproximationMethod newMethod)
{
  approxMethod = newMethod;

  // todo: let the prototype designer also use an enum class - maybe even the same - maybe list the
  // allpole types first, so we can simply do a < check to see, if the selected type is allpole
  // ...oh - but enum classes can't be converted to int - or can they explictly? we really want to 
  // avoid the "translation" - it's inelegant
  typedef RAPT::rsPrototypeDesigner<T>::approximationMethods PAM;
  typedef ApproximationMethod AM;
  switch(approxMethod)
  {
  case AM::gaussian:    protoDesigner.setApproximationMethod(PAM::GAUSSIAN);    break;
  case AM::bessel:      protoDesigner.setApproximationMethod(PAM::BESSEL);      break;
  case AM::butterworth: protoDesigner.setApproximationMethod(PAM::BUTTERWORTH); break;
  case AM::papoulis:    protoDesigner.setApproximationMethod(PAM::PAPOULIS);    break;
  case AM::halpern:     protoDesigner.setApproximationMethod(PAM::HALPERN);     break;

  // experimental:
  case AM::elliptic:    protoDesigner.setApproximationMethod(PAM::ELLIPTIC);    break;

  default: rsError("unknown approximation method");  
  }

  updateCoeffs();
}




template<class T>
void rsNonUniformFilterIIR<T>::updateCoeffs()
{
  protoDesigner.getPolesAndZeros(p, z);
  // this gets us only the non-redundant poles and zeros - complex poles always have a conjugate 
  // partner which is not in the array. if a real pole exists, it's in the last position of the
  // array
  // the array of zeros should be all inf - all zeros are at infinity because we have an allpole
  // filter - it's not actually used here (yet?)

  // make the pole-array digestible for the partial fraction expansion routine - the zeros
  // array is not needed in this case - the numerator of our rational function is just unity

  // create the complex conjugate partners (needed for partial fraction expansion):
  int i;
  for(i = (order-1)/2; i >= 0; i--) 
  {
    p[2*i]   = p[i];
    p[2*i+1] = conj(p[i]); 
    z[2*i]   = z[i];
    z[2*i+1] = conj(z[i]); 
  }

  // z is only half-full with inf - fill it up completely:
  //rsArray::fillWithValue(z, order, std::complex<T>(RS_INF(T), T(0)));
  // we should probably handle the zeros in exactly the same way as the poles - then we may later
  // be able to use all filter-types - also those with zeros - but when the number of finite zeros
  // is the same as the number of poles, the partial fraction expansion will produce a polynomial
  // part ...i think, a 0th order polynomial, i.e. just an added constant? may this be implemented
  // by just feeding through some of the filter's input? that would seem plausible...

  // do s-domain lowpass-to-lowpass transform to set up cutoff frequency:
  T k  = T(1);       // maybe make this a member...
  // maybe 1 is not correct in all cases? for Gaussian, Papoulis and Halpern filters, we get wrong 
  // DC gains

  k = T(1) / protoDesigner.getMagnitudeAt(T(0));

  T wc = 2*PI*freq;  // ...and this too?
  rsPoleZeroMapper<T>::sLowpassToLowpass(z, p, &k, z, p, &k, order, wc);
  // ...produces inf - j*nan for the zeros -> fix this!

  // create the sum-form of the denominator:
  rsPolynomial<T>::rootsToCoeffs(p, den, order);

  // todo: do the same for the numerator - we need to figure out the number of finite zeros...
  int nz = protoDesigner.getNumFiniteZeros();
  rsPolynomial<T>::rootsToCoeffs(z, num, nz);


  // do the partial fraction expansion:
  //rsRationalFunction<T>::partialFractionExpansion(num, 0, den, order, p, muls, order, r);
  rsRationalFunction<T>::partialFractionExpansion(num, nz, den, order, p, muls, order, r);
  // does not yet work when nz >= oder (numerator's degree >= denominator's) because the partial
  // fraction expansion will then contain a polynomial ("FIR") part and that's not yet implemented



  // transform analog poles to digital domain by means of impulse-invariant transform
  // see: https://ccrma.stanford.edu/~jos/pasp/Impulse_Invariant_Method.html
  for(i = 0; i < order; i++) {
    //z[i] = exp(z[i]);  // mapped zeros are not needed when dealing with a PFE
    p[i] = exp(p[i]);  // Eq. 9.2 with T=1 (T: sampling interval)
  }
  // todo: factor out into function impulseInvariantAnalogToDigital in class rsPoleZeroMapper

  // set up the one-pole filters:
  for(i = 0; i < order; i++)
    onePoles[i].setCoeffs(k*r[i], p[i]);

  // this seems to fix the problem of 1st order filter step responses shooting at a value higher 
  // than 1 (when piecewise resampling is used) - i don't know, why that works - the paper doesn't
  // say anything about doing such a thing:
  std::complex<T> tmp = T(0);
  for(i = 0; i < order; i++) {
    onePoles[i].reset(); // triggers computation of scaler s_i
    tmp += onePoles[i].getScaler();  
  }
  outScaler = T(1) / tmp.real();
  // but it makes sense: when the transfer function is a sum of terms like r_i / (s - p_i) and we 
  // plug in s = 0, we have a sum -r_i/p_i ...but s_i = r_i / (1-p_i) ...hmmm
}

/*

It seems like the only useful mode is "piecewise resampling" (unless i'm doing something wrong with 
the other modes) - make a production code version that strips off the stuff for the other modes,
includes optimizations and maybe continuous fade between piecewise resampling and no normalization.
This code here may then be moved back into the prototypes section (for further experiments with the
other modes)

Resources:

http://inf.ufrgs.br/~eslgastal/NonUniformFiltering/
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.533.3670&rep=rep1&type=pdf
http://ljk.imag.fr/membres/Brigitte.Bidegaray/Sources/FB10.pdf

*/