// construction/destruction:

template<class T>
rsInfiniteImpulseResponseDesigner<T>::rsInfiniteImpulseResponseDesigner()
{
  sampleRate     = 44100.0;
  frequency      = 1000.0;
  setBandwidth(1.0); // sets up lowerFrequency and upperFrequency
  gain           = rsAmpToDb(T(0.25));
  mode           = LOWPASS;
  prototypeOrder = 2;
}

template<class T>
rsInfiniteImpulseResponseDesigner<T>::~rsInfiniteImpulseResponseDesigner()
{

}

// parameter settings:

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setSampleRate(T newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setMode(int newMode)
{
  if( newMode >= BYPASS && newMode <= PEAK )
    mode = newMode;
  else
    RS_DEBUG_BREAK;  // newMode must be one of the numbers defined in enum 'modes'
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setPrototypeOrder(int newOrder)
{
  if( newOrder >= 1 )
    prototypeOrder = newOrder;
  else
    RS_DEBUG_BREAK;  // newOrder must be at least 1
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setApproximationMethod(int newApproximationMethod)
{
  prototypeDesigner.setApproximationMethod(newApproximationMethod);
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setFrequency(T newFrequency)
{
  frequency = newFrequency;
  calculateLowerAndUpperFrequency();
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setBandwidth(T newBandwidth)
{
  bandwidth = newBandwidth;
  calculateLowerAndUpperFrequency();
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setLowerFrequency(T newLowerFrequency)
{
  if( newLowerFrequency > 0.0 )
    lowerFrequency = newLowerFrequency;
  else
    RS_DEBUG_BREAK;  // negative frequencies are not allowed
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setUpperFrequency(T newUpperFrequency)
{
  if( newUpperFrequency > 0.0 )
    upperFrequency = newUpperFrequency;
  else
    RS_DEBUG_BREAK;  // negative frequencies are not allowed
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setGain(T newGain)
{
  gain = newGain;
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setRipple(T newRipple)
{
  prototypeDesigner.setPassbandRipple(newRipple);
  prototypeDesigner.setPassbandGainRatio(T(1)-T(0.01)*newRipple);
  prototypeDesigner.setStopbandGainRatio(T(0.01)*newRipple);
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::setStopbandRejection(T newStopbandRejection)
{
  prototypeDesigner.setStopbandRejection(newStopbandRejection);
}

// inquiry:

template<class T>
bool rsInfiniteImpulseResponseDesigner<T>::doesModeDoubleTheOrder()
{
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    return true;
  else
    return false;
}

template<class T>
int rsInfiniteImpulseResponseDesigner<T>::getFinalFilterOrder()
{
  if( doesModeDoubleTheOrder() == true )
    return 2*prototypeOrder;
  else
    return prototypeOrder;
}

template<class T>
int rsInfiniteImpulseResponseDesigner<T>::getNumBiquadStages()
{
  int order = getFinalFilterOrder();
  if( rsIsEven(order) )
    return order/2;
  else
    return (order+1)/2;
}

template<class T>
bool rsInfiniteImpulseResponseDesigner<T>::hasCurrentModeBandwidthParameter()
{
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    return true;
  else
    return false;
}

template<class T>
bool rsInfiniteImpulseResponseDesigner<T>::hasCurrentModeGainParameter()
{
  if( mode == LOW_SHELV || mode == HIGH_SHELV || mode == PEAK )
    return true;
  else
    return false;
}

template<class T>
bool rsInfiniteImpulseResponseDesigner<T>::hasCurrentModeRippleParameter()
{
  if( prototypeDesigner.hasCurrentMethodRippleParameter() )
  {
    if( mode != LOW_SHELV || mode != HIGH_SHELV || mode != PEAK )
      return true;
    else
      return false;
  }
  else
    return false;
}

template<class T>
bool rsInfiniteImpulseResponseDesigner<T>::hasCurrentModeRejectionParameter()
{
  if( prototypeDesigner.hasCurrentMethodRejectionParameter() )
  {
    if( mode != LOW_SHELV || mode != HIGH_SHELV || mode != PEAK )
      return true;
    else
      return false;
  }
  else
    return false;
}

// coefficient retrieval:

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::getPolesAndZeros(Complex* poles, Complex* zeros)
{
  // calculate the required order and number of biquads for the filter:
  int finalOrder;
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    finalOrder = 2*prototypeOrder;
  else
    finalOrder = prototypeOrder;

//  int numBiquads;
//  if( isEven(finalOrder) )
//    numBiquads = finalOrder/2;
//  else
//    numBiquads = (finalOrder+1)/2;

  // set up the type for the prototype (lowpass/low-shelv):
  if( mode == LOW_SHELV || mode == HIGH_SHELV || mode == PEAK || mode == BYPASS )
    prototypeDesigner.setPrototypeMode(rsPrototypeDesigner<T>::LOWSHELV_PROTOTYPE);
  else
    prototypeDesigner.setPrototypeMode(rsPrototypeDesigner<T>::LOWPASS_PROTOTYPE);

  T  f1, f2, wd1, wd2, wa1, wa2;
  T  fs = sampleRate;
  int     k;

  // use sanity-checked local frequency variables here:
  f2 = upperFrequency;
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    f1 = lowerFrequency;
    if( f2 > T(0.999*0.5)*fs )
      f2 = T(0.999*0.5)*fs;  // ensure upper frequency < sampleRate/2
    if( f1 > T(0.999)*f2  )
      f1 = T(0.999)*f2;      // ensure lower frequency < upper frequency
  }
  else
  {
    f1 = frequency;
    if( f1 > T(0.999*0.5)*fs )
      f1 = T(0.999*0.5)*fs;  // ensure frequency < sampleRate/2
  }
  // maybe factor this out!

  // prewarp the frequencies to the desired frequencies required for the design of the
  // (unnormalized) analog prototype filter:
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    wd1 = T(2*PI)*f1/fs;           // normalized digital radian frequency 1
    wa1 = T(2)*fs*tan(T(0.5)*wd1); // pre-warped analog radian frequency 1
    wd2 = T(2*PI)*f2/fs;           // normalized digital radian frequency 2
    wa2 = T(2)*fs*tan(T(0.5)*wd2); // pre-warped analog radian frequency 2
  }
  else
  {
    wd1 = T(2*PI)*f1/fs;           // normalized digital radian frequency 1
    wa1 = T(2)*fs*tan(T(0.5)*wd1); // pre-warped analog radian frequency 1
    wd2 = 0.0;                     // unused
    wa2 = 0.0;                     // unused
  }

  // set up the prototype-designer according to the specifications:
  prototypeDesigner.setOrder(prototypeOrder);
  prototypeDesigner.setGain(gain);
  if( mode == LOWPASS || mode == HIGHPASS || mode == BANDPASS || mode == BANDREJECT )
    prototypeDesigner.setReferenceGain(-RS_INF(T));
  else
    prototypeDesigner.setReferenceGain(0.0);

  // Allocate temporary memory (Get rid of this! Preallocate):
  //Complex* protoPoles = new Complex[prototypeOrder];
  //Complex* protoZeros = new Complex[prototypeOrder];
  std::vector<Complex> pp(prototypeOrder), pz(prototypeOrder); // !!!HEAP-ALLOCATION!!!
  Complex* protoPoles = &pp[0];
  Complex* protoZeros = &pz[0];

  // design the analog prototype filter:
  if( mode == HIGH_SHELV )
  {
    // we need to cope with exchange of roles of poles and zeros for high-shelving (inverse)
    // chebychevs because the low-shelv -> high-shelv frequency transform exchanges these roles
    // once again:
    if( prototypeDesigner.getApproximationMethod() == rsPrototypeDesigner<T>::CHEBYCHEV )
    {
      prototypeDesigner.setApproximationMethod(rsPrototypeDesigner<T>::INVERSE_CHEBYCHEV);
      prototypeDesigner.getPolesAndZeros(poles, zeros);
      prototypeDesigner.setApproximationMethod(rsPrototypeDesigner<T>::CHEBYCHEV);
    }
    else if( prototypeDesigner.getApproximationMethod() == rsPrototypeDesigner<T>::INVERSE_CHEBYCHEV )
    {
      prototypeDesigner.setApproximationMethod(rsPrototypeDesigner<T>::CHEBYCHEV);
      prototypeDesigner.getPolesAndZeros(poles, zeros);
      prototypeDesigner.setApproximationMethod(rsPrototypeDesigner<T>::INVERSE_CHEBYCHEV);
    }
    else
      prototypeDesigner.getPolesAndZeros(poles, zeros);
  }
  else
    prototypeDesigner.getPolesAndZeros(poles, zeros);

  // because the PrototypeDesigner returns only one representant for each pair of complex 
  // conjugate poles/zeros, we now create the full set here:
  if( rsIsOdd(prototypeOrder) )
  {
    // copy the real pole/zero to the end:
    protoPoles[prototypeOrder-1] = poles[prototypeOrder/2];
    protoZeros[prototypeOrder-1] = zeros[prototypeOrder/2];
  }
  // for each complex pole/zero create a pair of complex conjugates:
  for(k=0; k<prototypeOrder/2; k++)
  {
    protoPoles[2*k]   = poles[k];
    protoPoles[2*k+1] = conj(poles[k]);
    protoZeros[2*k]   = zeros[k];
    protoZeros[2*k+1] = conj(zeros[k]);
  }

  // s-plane frequency transformation:
  switch( mode )
  {
  case LOWPASS:    rsPoleZeroMapper<T>::sPlanePrototypeToLowpass(   protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  case HIGHPASS:   rsPoleZeroMapper<T>::sPlanePrototypeToHighpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  case BANDPASS:   rsPoleZeroMapper<T>::sPlanePrototypeToBandpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1, wa2); break;
  case BANDREJECT: rsPoleZeroMapper<T>::sPlanePrototypeToBandreject(protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1, wa2); break;
  case LOW_SHELV:  rsPoleZeroMapper<T>::sPlanePrototypeToLowpass(   protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  case HIGH_SHELV:
    {
      // this is ugly - rewrite the prototype design code such that all approximation methods can be treated uniformly:
      if( prototypeDesigner.needsSpecialHighShelvTransform() )
        rsPoleZeroMapper<T>::sPlanePrototypeToHighShelv( protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);
      else
        rsPoleZeroMapper<T>::sPlanePrototypeToHighpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);
    }
    break;
  case PEAK:       rsPoleZeroMapper<T>::sPlanePrototypeToBandpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1, wa2); break;
  default:         rsPoleZeroMapper<T>::sPlanePrototypeToLowpass(   protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  };

  //std::vector<Complex> pDbg, zDbg; // for debugging
  //pDbg = toVector(poles, finalOrder);
  //zDbg = toVector(zeros, finalOrder);

  // transform from the s-domain to the z-domain via the bilinear transform:
  T g; // not used actually
  rsPoleZeroMapper<T>::bilinearAnalogToDigital(poles, finalOrder, zeros, finalOrder, fs, &g);

  //pDbg = toVector(poles, finalOrder);
  //zDbg = toVector(zeros, finalOrder);
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::getBiquadCascadeCoefficients(T *b0, T *b1, 
  T *b2, T *a1, T *a2)
{
  // calculate the required order and number of biquads for the filter:
  int finalOrder, numBiquads;

  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    finalOrder = 2*prototypeOrder;
  else
    finalOrder = prototypeOrder;

  if( rsIsEven(finalOrder) )
    numBiquads = finalOrder/2;
  else
    numBiquads = (finalOrder+1)/2;

  // set up the type for the prototype (lowpass/low-shelv):
  if( mode == LOW_SHELV || mode == HIGH_SHELV || mode == PEAK || mode == BYPASS )
  {
    prototypeDesigner.setPrototypeMode(rsPrototypeDesigner<T>::LOWSHELV_PROTOTYPE);
    if( rsIsCloseTo(gain, 0.0, 0.001) || mode == BYPASS ) // gains of zero yield a 'bypass' filter
    {
      for(int b = 0; b < numBiquads; b++) // msvc gives an "unreachable code" warning here - but
      {                                   // the code is definitely reachable...hmmm
        b0[b] = 1.0;
        b1[b] = 0.0;
        b2[b] = 0.0;
        a1[b] = 0.0;
        a2[b] = 0.0;
        return;
      }
    }
  }
  else
    prototypeDesigner.setPrototypeMode(rsPrototypeDesigner<T>::LOWPASS_PROTOTYPE);

  T  f1, f2, wd1, wd2; // wa1, wa2;
  T  fs = sampleRate;

  // use sanity-checked local frequency variables here:
  f2 = upperFrequency;
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    f1 = lowerFrequency;
    if( f2 > T(0.99*0.5)*fs )
      f2 = T(0.99*0.5)*fs;  // ensure upper frequency < sampleRate/2
    if( f1 > T(0.99)*f2  )
      f1 = T(0.99)*f2;      // ensure lower frequency < upper frequency
  }
  else
  {
    f1 = frequency;
    if( f1 > T(0.99*0.5)*fs )
      f1 = T(0.99*0.5)*fs;  // ensure frequency < sampleRate/2
  }

  // prewarp the frequencies to the desired frequencies required for the design of the 
  // (unnormalized) analog prototype filter:
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    wd1 = T(2*PI)*f1/fs;         // normalized digital radian frequency 1
    //wa1 = 2.0*fs*tan(0.5*wd1);  // pre-warped analog radian frequency 1
    wd2 = T(2*PI)*f2/fs;         // normalized digital radian frequency 2
    //wa2 = 2.0*fs*tan(0.5*wd2);  // pre-warped analog radian frequency 2
  }
  else
  {
    wd1 = T(2*PI)*f1/fs;         // normalized digital radian frequency 1
    //wa1 = 2.0*fs*tan(0.5*wd1);  // pre-warped analog radian frequency 1
    wd2 = 0.0;                  // unused
    //wa2 = 0.0;                  // unused
  }

  std::vector<Complex> poles(finalOrder), zeros(finalOrder); // !!!HEAP-ALLOCATION!!!

  getPolesAndZeros(&poles[0], &zeros[0]); // seems to return zeros in wrong order for elliptic bandpass
  rsFilterCoefficientConverter<T>::polesAndZerosToBiquadCascade(&poles[0], &zeros[0], finalOrder, 
    b0, b1, b2, a1, a2);
  T deWarpedPassbandCenter = T(2) * atan(sqrt(tan(T(0.5)*wd1)*tan(T(0.5)*wd2)));
  normalizeGain(b0, b1, b2, a1, a2, deWarpedPassbandCenter, numBiquads);
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::getDirectFormCoefficients(T *b, T *a)
{
  int numBiquads = getNumBiquadStages();
  T *b0 = new T[numBiquads];
  T *b1 = new T[numBiquads];
  T *b2 = new T[numBiquads];
  T *a1 = new T[numBiquads];
  T *a2 = new T[numBiquads];
  getBiquadCascadeCoefficients(b0, b1, b2, a1, a2);
  rsFilterCoefficientConverter<T>::biquadCascadeToDirectForm(numBiquads, b0, b1, b2, a1, a2, b, a);
  delete[] b0;
  delete[] b1;
  delete[] b2;
  delete[] a1;
  delete[] a2;
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::calculateLowerAndUpperFrequency()
{
  lowerFrequency = frequency / pow(T(2), T(0.5)*bandwidth);
  upperFrequency = lowerFrequency * pow(T(2), bandwidth);
}

template<class T>
void rsInfiniteImpulseResponseDesigner<T>::normalizeGain(T *b0, T *b1, T *b2, 
  T *a1, T *a2, T wc, int numBiquads)
{
  T w = 0.0;
  if( mode == LOWPASS || mode == BANDREJECT || mode == HIGH_SHELV || mode == PEAK ) // redundant
    w = 0.0;           // normalize at DC
  else if( mode == HIGHPASS || mode == LOW_SHELV )
    w = T(PI);            // normalize at Nyquist frequency
  else if( mode == BANDPASS )
    w = wc;            // normalize at passband center frequency

  T normalizeFactor = T(1);
  if( prototypeDesigner.getPrototypeMode() == rsPrototypeDesigner<T>::LOWSHELV_PROTOTYPE )
  {
    if(  prototypeDesigner.getApproximationMethod() == rsPrototypeDesigner<T>::INVERSE_CHEBYCHEV
      || prototypeDesigner.getApproximationMethod() == rsPrototypeDesigner<T>::ELLIPTIC )
    {
      if( rsIsEven(prototypeDesigner.getOrder()) )
      {
        T factor   = T(1)-prototypeDesigner.getPassbandGainRatio();
        T excessDb = -factor * gain;
        normalizeFactor = rsDbToAmp(-excessDb/numBiquads);
      }
    }
  }
  else // prototype is lowpass
  {
    if(  prototypeDesigner.getApproximationMethod() == rsPrototypeDesigner<T>::CHEBYCHEV
      || prototypeDesigner.getApproximationMethod() == rsPrototypeDesigner<T>::ELLIPTIC )
    {
      if( rsIsEven(prototypeDesigner.getOrder()) )
      {
        normalizeFactor = T(1) / rsDbToAmp(prototypeDesigner.getPassbandRipple());
        if( doesModeDoubleTheOrder() )
          normalizeFactor = pow(normalizeFactor, T(1) / prototypeDesigner.getOrder());
        else
          normalizeFactor = pow(normalizeFactor, T(2) / prototypeDesigner.getOrder());
      }
    }
  }

  rsFilterCoefficientConverter<T>::normalizeBiquadStages(b0, b1, b2, a1, a2, w, numBiquads, 
    normalizeFactor);
}


/*

ToDo:

- Get rid of the heap allocations for temporary poles and zeros (i.e. the transformed poles and
  zeros with s-domain and s-to-z transformations applied). Use pre-allocated buffers for that, i.e.
  member variables. This has the additional advantage that we have the poles and zeros available
  for inspection and possibly even for running the filter in complex 1-pole form.

*/