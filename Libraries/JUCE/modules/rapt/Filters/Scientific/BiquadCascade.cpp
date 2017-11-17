// construction/destruction:

template<class TSig, class TCoef>
rsBiquadCascade<TSig, TCoef>::rsBiquadCascade(int newMaxNumStages)
{
  if( newMaxNumStages >= 1 )
    maxNumStages = newMaxNumStages;
  else
    maxNumStages = 12;

  // allcocate memory for coefficients and buffers:
  a1 = new TCoef[maxNumStages];
  a2 = new TCoef[maxNumStages];
  b0 = new TCoef[maxNumStages];
  b1 = new TCoef[maxNumStages];
  b2 = new TCoef[maxNumStages];
  x1 = new TCoef[maxNumStages];
  x2 = new TCoef[maxNumStages];
  y1 = new TCoef[maxNumStages];
  y2 = new TCoef[maxNumStages];

  initBiquadCoeffs();        
  reset();                   
  numStages = maxNumStages; 
}

template<class TSig, class TCoef>
rsBiquadCascade<TSig, TCoef>::~rsBiquadCascade()
{
  delete[] a1;
  delete[] a2;
  delete[] b0;
  delete[] b1;
  delete[] b2;
  delete[] x1;
  delete[] x2;
  delete[] y1;
  delete[] y2;
}

// parameter settings:

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::setNumStages(int newNumStages)
{
  if( (newNumStages >= 0 ) && (newNumStages <= maxNumStages) )
    numStages = newNumStages;
  else
    DEBUG_BREAK;
  reset();
}

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::setOrder(int newOrder)
{
  if( isEven(newOrder) )
    setNumStages(newOrder/2);
  else
    setNumStages( (newOrder+1)/2 );
}

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::setGlobalGainFactor(TCoef newGain)
{
  b0[0] *= newGain;
  b1[0] *= newGain;
  b2[0] *= newGain;
}

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::copySettingsFrom(rsBiquadCascade *other)
{
  setNumStages(other->getNumStages());

  TCoef* pB0 = other->getAddressB0();
  TCoef* pB1 = other->getAddressB1();
  TCoef* pB2 = other->getAddressB2();
  TCoef* pA1 = other->getAddressA1();
  TCoef* pA2 = other->getAddressA2();
  for(int s = 0; s < numStages; s++)
  {
    b0[s] = pB0[s]; b1[s] = pB1[s]; b2[s] = pB2[s]; 
    a1[s] = pA1[s]; a2[s] = pA2[s]; 
  }
}

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::turnIntoAllpass()
{
  for(int i = 0; i < numStages; i++)
  {
    if( a2[i] == 0.0 )  // biquad stage is actually 1st order
    {
      b0[i] = a1[i];
      b1[i] = 1.0;
      b2[i] = 0.0;
    }
    else
    {

      b0[i]  = a2[i];
      b1[i]  = a1[i];
      b2[i]  = 1.0;
    }
  }

  // normalize gain of each stage to unity:
  for(int i = 0; i < numStages; i++)
  {
    TCoef num = b0[i]*b0[i] + b1[i]*b1[i] + b2[i]*b2[i] 
      + 2.0*(b0[i]*b1[i] + b1[i]*b2[i]) + 2.0*b0[i]*b2[i];
    TCoef den = 1.0 + a1[i]*a1[i] + a2[i]*a2[i] + 2.0*(a1[i] + a1[i]*a2[i]) + 2.0*a2[i];
    TCoef g   = sqrt(den/num);
    b0[i] *= g;
    b1[i] *= g;
    b2[i] *= g;
  }
}

// inquiry:

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::getFrequencyResponse(TCoef* w, std::complex<TCoef>* H, 
  int numBins, int accumulationMode)
{
  FilterAnalyzer::getBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numStages, w, H, numBins, 
    accumulationMode);
}

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::getMagnitudeResponse(TCoef *w, TCoef *magnitudes, int numBins, 
  bool inDecibels, bool accumulate)
{
  FilterAnalyzer::getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numStages, w, magnitudes, 
    numBins, inDecibels, accumulate);
}

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::getMagnitudeResponse(TCoef* frequencies, TCoef sampleRate, 
  TCoef* magnitudes, int numBins, bool inDecibels, bool accumulate)
{
  FilterAnalyzer::getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numStages, frequencies, 
    sampleRate, magnitudes, numBins, inDecibels, accumulate);
}

//-------------------------------------------------------------------------------------------------
// others:

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::initBiquadCoeffs()
{
  for(int i = 0; i < maxNumStages; i++)
  {
    b0[i] = 1.0;
    b1[i] = 0.0;
    b2[i] = 0.0;
    a1[i] = 0.0;
    a2[i] = 0.0;
  }
}

template<class TSig, class TCoef>
void rsBiquadCascade<TSig, TCoef>::reset()
{
  for(int i = 0; i < maxNumStages; i++)
  {
    x2[i] = 0.0;
    x1[i] = 0.0;
    y2[i] = 0.0;
    y1[i] = 0.0;
  }
}
