// Construction/Destruction:

template<class TSig, class TPar>
rsOnePoleFilter<TSig, TPar>::rsOnePoleFilter()
{
  shelvingGain = 1.0;
  setSampleRate(44100.0);  // sampleRate = 44100 Hz by default
  setMode      (0);        // bypass by default
  setCutoff    (20000.0);  // cutoff = 20000 Hz by default
}

// Setup:

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  if(newSampleRate > 0.0)
    freqToOmega = TPar(2*PI)/newSampleRate;
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
  //if( (newCutoff > 0.0) && (newCutoff <= 20000.0) )
  //  cutoff = newCutoff;
  //else
  //  cutoff = 20000.0;
  // ...check, if any code relies on this safety check - if so, restrict the cutoff range on the 
  // higher level


  cutoff = newCutoff;
  calcCoeffs();
  return;
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setShelvingGain(TPar newGain)
{
  rsAssert(newGain >= 0.0); // this is a linear gain factor and must be >= 0.0
  shelvingGain = newGain;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::setShelvingGainInDecibels(TPar newGain)
{
  setShelvingGain(rsDB2amp(newGain));
}

// Inquiry:

template<class TSig, class TPar>
TPar rsOnePoleFilter<TSig, TPar>::getMagnitudeAt(TPar f)
{
  return onePoleMagnitudeAt(this->b0, this->b1, -this->a1, freqToOmega*f);
    // we use a different sign-convention for the a1 coefficient here ...maybe fix this someday
}

// Misc:

template<class TSig, class TPar>
void rsOnePoleFilter<TSig, TPar>::calcCoeffs()
{
  typedef rsFirstOrderFilterBase<TSig, TPar> B; // B for "baseclass"
  TPar w = freqToOmega*cutoff;
  TPar g = shelvingGain;
  switch(mode)
  {
    // handle these cases in baseclass:
  case LOWPASS_IIT:   { B::coeffsLowpassIIT(  w,    &this->b0, &this->b1, &this->a1); } break;
  case HIGHPASS_MZT:  { B::coeffsHighpassMZT( w,    &this->b0, &this->b1, &this->a1); } break;
  case ALLPASS_BLT:   { B::coeffsAllpassBLT(  w,    &this->b0, &this->b1, &this->a1); } break;
  case LOWPASS_BLT:   { B::coeffsLowpassBLT(  w,    &this->b0, &this->b1, &this->a1); } break;
  case HIGHPASS_BLT:  { B::coeffsHighpassBLT( w,    &this->b0, &this->b1, &this->a1); } break;
  case LOWSHELV_BLT:  { B::coeffsLowShelfBLT( w, g, &this->b0, &this->b1, &this->a1); } break;
  case HIGHSHELV_BLT: { B::coeffsHighShelfBLT(w, g, &this->b0, &this->b1, &this->a1); } break;

  // these two need clean-up (and tests):
  case LOWSHELV_NMM:
    {
      // formula derived as special case of the Orfanidis equalizers:
      // TPar g    = rsDB2amp(shelvingGain);
      //TPar g    = shelvingGain;
      //TPar wc   = 2*PI*cutoff*sampleRateRec;
      TPar sampleRate = TPar(2*PI)/freqToOmega;
      TPar wa   = TPar(2)*sampleRate*tan(w/2);
      TPar gb   = rsSqrt(g);
      TPar beta = rsSqrt( (gb*gb-1)/(g*g-gb*gb) );
      TPar pa   = -beta*wa;
      TPar za   = -g*beta*wa;
      TPar s    = TPar(0.5)/sampleRate;
      TPar p    = (1+s*pa)/(1-s*pa);
      TPar z    = (1+s*za)/(1-s*za);
      this->b1  = -z;
      this->a1  = -p;
      TPar n    = rsSqrt((1+this->a1*this->a1-2*this->a1) / (1+this->b1*this->b1-2*this->b1));
      this->b0  = n;
      this->b1 *= n;
      // this seems overly complicated - can't we just derive the coeffs directly from 3 magnitude
      // constraints - should be a simple 3x3 linear system ...check this out....
    }
    break;
  case HIGHSHELV_NMM:
    {
      // formula derived as special case of the Orfanidis equalizers:
      //TPar g    = rsDB2amp(shelvingGain);  // wrong? shelvingGain is already linear
      //TPar g    = shelvingGain;
      //TPar wc   = 2*PI*cutoff*sampleRateRec; // redundant with w
      TPar sampleRate = TPar(2*PI)/freqToOmega;
      TPar wa   = 2*sampleRate*tan(w/2);
      TPar gb   = rsSqrt(g);
      TPar beta = rsSqrt( (gb*gb-1)/(g*g-gb*gb) );
      TPar pa   = -beta*wa;
      TPar za   = -g*beta*wa;
      TPar s    = TPar(0.5)/sampleRate;
      TPar p    = (1+s*pa)/(1-s*pa);
      TPar z    = (1+s*za)/(1-s*za);
      this->b1  = -p;
      this->a1  = -z;
      TPar n    = rsSqrt((1+this->a1*this->a1+2*this->a1) / (1+this->b1*this->b1+2*this->b1));
      this->b0  = n;
      this->b1 *= n;
      // \todo get rid of the code duplication
    }
    break;


  default: { B::coeffsBypass(&this->b0, &this->b1, &this->a1); } break;
  }
}