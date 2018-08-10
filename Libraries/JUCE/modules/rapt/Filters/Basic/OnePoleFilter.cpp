// Construction/Destruction:

template<class TSig, class TPar>
rsOnePoleFilter<TSig, TPar>::rsOnePoleFilter()
{
  shelvingGain = 1.0;
  setSampleRate(44100.0);  // sampleRate = 44100 Hz by default
  setMode      (0);        // bypass by default
  setCutoff    (20000.0);  // cutoff = 20000 Hz by default
  reset();                 // reset memorized samples to zero
}

// Setup:

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  sampleRateRec = 1.0 / sampleRate;

  calcCoeffs();
  return;
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setMode(int newMode)
{
  mode = newMode; // 0:bypass, 1:Low Pass, 2:High Pass
  calcCoeffs();
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setCutoff(TPar newCutoff)
{
  if( (newCutoff > 0.0) && (newCutoff <= 20000.0) )
    cutoff = newCutoff;
  else
    cutoff = 20000.0;

  calcCoeffs();
  return;
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setShelvingGain(TPar newGain)
{
  if( newGain > 0.0 )
  {
    shelvingGain = newGain;
    calcCoeffs();
  }
  else
    RS_DEBUG_BREAK; // this is a linear gain factor and must be >= 0.0
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setShelvingGainInDecibels(TPar newGain)
{
  setShelvingGain(rsDB2amp(newGain));
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setCoefficients(TPar newB0, TPar newB1, TPar newA1)
{
  b0 = newB0;
  b1 = newB1;
  a1 = newA1;
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setInternalState(TSig newX1, TSig newY1)
{
  x1 = newX1;
  y1 = newY1;
}

// Inquiry:

template<class TSig, class TPar>
TPar rsOnePoleFilter<TSig, TPar>::getMagnitudeAt(TPar f)
{
  return onePoleMagnitudeAt(b0, b1, -a1, 2*PI*f*sampleRateRec);
    // we use a different sign-convention for the a1 coefficient here ...maybe fix this someday
}

// Misc:

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::calcCoeffs()
{
  // maybe move these to FilterDesignFormulas - factor out, maybe together with biquad formulas
  // BUT: make sure to be clear about the sign-convention for feedback coeffs. here we assume
  // that the feedback coeffs are used with positive sign in the filter update
  // these functions should then take only an omega as input
  switch(mode)
  {
  case LOWPASS_IIT: 
    {
      // formula from dspguide (impulse invariant):
      TPar x = exp(-2.0 * PI * cutoff * sampleRateRec); 
      b0 = 1-x;
      b1 = 0.0;
      a1 = x;
    }
    break;
  case HIGHPASS_MZT:  
    {
      // formula from dspguide (impulse invariant):
      TPar x = exp(-2.0 * PI * cutoff * sampleRateRec);
      b0 =  0.5*(1+x);
      b1 = -0.5*(1+x);  // = -b0 -> optimize
      a1 = x;
    }
    break;
  case ALLPASS_BLT:  
    {
      // formula from DAFX (bilinear):
      TPar t = tan(PI*cutoff*sampleRateRec);
      TPar x = (t-1.0) / (t+1.0);

      b0 = x;
      b1 = 1.0;
      a1 = -x;
    }
    break;
  case LOWSHELV_NMM:
    {
      // formula derived as special case of the Orfanidis equalizers:
      TPar g    = rsDB2amp(shelvingGain);
      TPar wc   = 2*PI*cutoff*sampleRateRec;
      TPar wa   = 2*sampleRate*tan(wc/2);
      TPar gb   = rsSqrt(g);
      TPar beta = rsSqrt( (gb*gb-1)/(g*g-gb*gb) );
      TPar pa   = -beta*wa;
      TPar za   = -g*beta*wa;
      TPar s    = 0.5*sampleRateRec;
      TPar p    = (1+s*pa)/(1-s*pa);
      TPar z    = (1+s*za)/(1-s*za);
      b1        = -z;
      a1        = -p;
      TPar n    = rsSqrt((1+a1*a1-2*a1) / (1+b1*b1-2*b1));
      b0        = n;
      b1       *= n;
      // this seems overly complicated - can't we just derive the coeffs directly from 3 magnitude
      // constraints ...check this out....
    }
    break;
  case HIGHSHELV_NMM:
    {
      // formula derived as special case of the Orfanidis equalizers:
      TPar g    = rsDB2amp(shelvingGain);
      TPar wc   = 2*PI*cutoff*sampleRateRec;
      TPar wa   = 2*sampleRate*tan(wc/2);
      TPar gb   = rsSqrt(g);
      TPar beta = rsSqrt( (gb*gb-1)/(g*g-gb*gb) );
      TPar pa   = -beta*wa;
      TPar za   = -g*beta*wa;
      TPar s    = 0.5*sampleRateRec;
      TPar p    = (1+s*pa)/(1-s*pa);
      TPar z    = (1+s*za)/(1-s*za);
      b1        = -p;
      a1        = -z;
      TPar n    = rsSqrt((1+a1*a1+2*a1) / (1+b1*b1+2*b1));
      b0        = n;
      b1       *= n;
      // \todo get rid of the code duplication
    }
    break;
  case LOWSHELV_BLT:
    {
      // formula from DAFX:
      TPar c = 0.5*(shelvingGain-1.0);
      TPar t = tan(PI*cutoff*sampleRateRec);
      TPar a;
      if( shelvingGain >= 1.0 )
        a = (t-1.0)/(t+1.0);
      else
        a = (t-shelvingGain)/(t+shelvingGain);

      b0 = 1.0 + c + c*a;
      b1 = c + c*a + a;
      a1 = -a;
    }
    break;
  case HIGHSHELV_BLT:
    {
      // formula from DAFX:
      TPar c = 0.5*(shelvingGain-1.0);
      TPar t = tan(PI*cutoff*sampleRateRec);
      TPar a;
      if( shelvingGain >= 1.0 )
        a = (t-1.0)/(t+1.0);
      else
        a = (shelvingGain*t-1.0)/(shelvingGain*t+1.0);

      b0 = 1.0 + c - c*a;
      b1 = a + c*a - c;
      a1 = -a;
    }
    break;
  case LOWPASS_BLT:  
  {
    TPar t = tan(PI*cutoff*sampleRateRec);
    TPar x = (t-1.0) / (t+1.0);
    b0 = 0.5*(1+x);
    b1 = b0;
    a1 = -x;
  }
  break;
  case HIGHPASS_BLT:  
  {
    TPar t = tan(PI*cutoff*sampleRateRec);
    TPar x = (t-1.0) / (t+1.0);
    b0 = 0.5*(1-x);
    b1 = -b0;
    a1 = -x;
  }
  break;
  default: // bypass
    {
      b0 = 1.0;
      b1 = 0.0;
      a1 = 0.0;
    }
    break;
  }
  // \todo provide IIT, BLT, and magnitude-matched (MMT) transform versions
  // get rid of these strange names DAFX, etc.
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::reset()
{
  x1 = 0.0;
  y1 = 0.0;
}
