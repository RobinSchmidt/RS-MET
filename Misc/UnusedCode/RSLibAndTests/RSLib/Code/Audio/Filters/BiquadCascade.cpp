#ifndef RSLib_BiquadCascade_cpp
#define RSLib_BiquadCascade_cpp


// construction/destruction:

rsBiquadCascade::rsBiquadCascade(int newMaxNumStages)
{
  if( newMaxNumStages >= 1 )
    maxNumStages = newMaxNumStages;
  else
    maxNumStages = 12;

  // allcocate memory for coefficients and buffers:
  a1 = new double[maxNumStages];
  a2 = new double[maxNumStages];
  b0 = new double[maxNumStages];
  b1 = new double[maxNumStages];
  b2 = new double[maxNumStages];
  x1 = new double[maxNumStages];
  x2 = new double[maxNumStages];
  y1 = new double[maxNumStages];
  y2 = new double[maxNumStages];

  initBiquadCoeffs();        
  reset();                   
  numStages = maxNumStages; 
}

rsBiquadCascade::~rsBiquadCascade()
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

void rsBiquadCascade::setNumStages(int newNumStages)
{
  if( (newNumStages >= 0 ) && (newNumStages <= maxNumStages) )
    numStages = newNumStages;
  else
    rsError("numStages > maxNumStages");
  reset();
}

void rsBiquadCascade::setOrder(int newOrder)
{
  if( rsIsEven(newOrder) )
    setNumStages(newOrder/2);
  else
    setNumStages( (newOrder+1)/2 );
}

void rsBiquadCascade::setGlobalGainFactor(double newGain)
{
  b0[0] *= newGain;
  b1[0] *= newGain;
  b2[0] *= newGain;
}

void rsBiquadCascade::copySettingsFrom(rsBiquadCascade *other)
{
  setNumStages(other->getNumStages());

  double *pB0 = other->getAddressB0();
  double *pB1 = other->getAddressB1();
  double *pB2 = other->getAddressB2();
  double *pA1 = other->getAddressA1();
  double *pA2 = other->getAddressA2();
  for(int s=0; s<numStages; s++)
  {
    b0[s] = pB0[s]; b1[s] = pB1[s]; b2[s] = pB2[s]; 
    a1[s] = pA1[s]; a2[s] = pA2[s]; 
  }
}

void rsBiquadCascade::turnIntoAllpass()
{
  for(int i=0; i<numStages; i++)
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
  for(int i=0; i<numStages; i++)
  {
    double num = b0[i]*b0[i]+b1[i]*b1[i]+b2[i]*b2[i]+2.0*(b0[i]*b1[i]+b1[i]*b2[i])+2.0*b0[i]*b2[i];
    double den = 1.0        +a1[i]*a1[i]+a2[i]*a2[i]+2.0*(      a1[i]+a1[i]*a2[i])+2.0*      a2[i];
    double g   = sqrt(den/num);
    b0[i] *= g;
    b1[i] *= g;
    b2[i] *= g;
  }
}

// inquiry:

void rsBiquadCascade::getFrequencyResponse(double *w, rsComplexDbl *H, int numBins, 
  int accumulationMode)
{
  rsFilterAnalyzer::getBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numStages, w, H, numBins,
    accumulationMode);
}

void rsBiquadCascade::getMagnitudeResponse(double *w, double *magnitudes, int numBins, 
  bool inDecibels, bool accumulate)
{
  rsFilterAnalyzer::getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numStages, w, magnitudes, 
    numBins, inDecibels, accumulate);
}

void rsBiquadCascade::getMagnitudeResponse(double *frequencies, double sampleRate, 
  double *magnitudes, int numBins, bool inDecibels, bool accumulate)
{
  rsFilterAnalyzer::getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numStages, frequencies, 
    sampleRate, magnitudes, numBins, inDecibels, accumulate);
}


double biquadRingingTime(double a1, double a2, double threshold)
{
  rsComplexDbl p = -a1/2 + rsSqrtC(rsComplexDbl(a1*a1/4.0-a2)); // 1st pole
  return rsLogB(threshold, p.getRadius());
}
// move to somewhere else

double rsBiquadCascade::getRingingTimeEstimate(double threshold)
{
  double rt = 0.0;
  for(int i = 0; i < numStages; i++)
    rt = rsMax(rt, biquadRingingTime(a1[i], a2[i], threshold));
  return rt;
}


/*
double biquadRingingTime(double a1, double a2, double threshold)
{
  // compute poles from coeffs (factor out):
  rsComplexDbl p1, p2, sq;
  sq = rsSqrtC(rsComplexDbl(a1*a1/4.0-a2));
  p1 = -a1/2 + sq;
  p2 = -a1/2 - sq;

  // compute ringing times of the indiviual poles and add them up:
  double rt;
  rt = rsLogB(threshold, p1.getRadius());
  if( a2 != 0.0 )
    rt += rsLogB(threshold, p2.getRadius());

  return rt;
}
// move to somewhere else

double rsBiquadCascade::getRingingTime(double threshold)
{
  double rt = 0.0;
  for(int i = 0; i < numStages; i++)
    rt += biquadRingingTime(a1[i], a2[i], threshold);
  return rt;
  // this doesn't work well for higher order filters, maybe, we shouldn't add the times but choose
  // the maximum
}
*/

// misc:

void rsBiquadCascade::initBiquadCoeffs()
{
  for(int i=0; i<maxNumStages; i++)
  {
    b0[i] = 1.0;
    b1[i] = 0.0;
    b2[i] = 0.0;
    a1[i] = 0.0;
    a2[i] = 0.0;
  }
}

void rsBiquadCascade::reset()
{
  for(int i=0; i<maxNumStages; i++)
  {
    x2[i] = 0.0;
    x1[i] = 0.0;
    y2[i] = 0.0;
    y1[i] = 0.0;
  }
}

//-------------------------------------------------------------------------------------------------
// class rsBiquadCascadeStereo:

rsBiquadCascadeStereo::rsBiquadCascadeStereo(int newMaxNumStages) 
  : rsBiquadCascade(newMaxNumStages)
{
  x1R = new double[maxNumStages];
  x2R = new double[maxNumStages];
  y1R = new double[maxNumStages];
  y2R = new double[maxNumStages];
       
  reset();                   
}

rsBiquadCascadeStereo::~rsBiquadCascadeStereo()
{
  delete[] x1R;
  delete[] x2R;
  delete[] y1R;
  delete[] y2R;
}

void rsBiquadCascadeStereo::reset()
{
  rsBiquadCascade::reset();
  for(int i=0; i<maxNumStages; i++)
  {
    x2R[i] = 0.0;
    x1R[i] = 0.0;
    y2R[i] = 0.0;
    y1R[i] = 0.0;
  }
}

#endif