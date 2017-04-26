#include "rosic_DirectFormFilter.h"
using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

DirectFormFilter::DirectFormFilter(int maximumOrder)
{
  maxOrder = maximumOrder;
  order    = maxOrder;
  w        = new double[order+1];
  a        = new double[order+1];
  b        = new double[order+1];
  initializeCoefficients();
  reset();  
}

DirectFormFilter::~DirectFormFilter()
{
  delete[] w;
  delete[] a;
  delete[] b;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void DirectFormFilter::setCoefficients(double *newCoeffsA, double *newCoeffsB, int newOrder)
{
  initializeCoefficients();  // why that?
  order = newOrder;
  for(int k=0; k<=newOrder; k++)
  {
    a[k] = newCoeffsA[k];
    b[k] = newCoeffsB[k];
  }
}

void DirectFormFilter::setGlobalGainFactor(double newFactor)
{
  for(int k=0; k<=order; k++)
    b[k] *= newFactor;
}

double DirectFormFilter::getMagnitudeResponseAt(double omega)
{
  double ca = 0.0;
  double sa = 0.0;
  double cb = 0.0; 
  double sb = 0.0;
  double sk, ck;
  for(int k=0; k<=order; k++)
  {
    sk  = sin(k*omega);  // \todo optimize by trigonometric recursion
    ck  = cos(k*omega);  
    ca += ck * a[k];
    sa += sk * a[k];
    cb += ck * b[k];
    sb += sk * b[k];
  }
  return sqrt( (cb*cb + sb*sb) / (ca*ca + sa*sa) );
}

void DirectFormFilter::getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, double sampleRate, 
                                            bool inDecibels, bool accumulate)
{
  double H;
  for(int k=0; k<numBins; k++)
  {
    H = getMagnitudeResponseAt(2*PI*frequencies[k] / sampleRate);

    // \todo move this decision function somewher where it can be shared - other filter use that, too
    if( !inDecibels && !accumulate )
      magnitudes[k]  = H;
    else if( !inDecibels && accumulate )
      magnitudes[k] *= H;
    else if( inDecibels && !accumulate )
      magnitudes[k]  = amp2dBWithCheck(H, 0.0000001); // clip at -140 dB
    else if( inDecibels && accumulate )
      magnitudes[k] += amp2dBWithCheck(H, 0.0000001);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void DirectFormFilter::reset()
{
  for(int i=0; i<=maxOrder; i++)
    w[i] = 0.0;
}

void DirectFormFilter::initializeCoefficients()
{
  b[0] = 1.0;
  a[0] = 1.0;
  for(int i=1; i<=maxOrder; i++)
  {
    b[i] = 0.0;
    a[i] = 0.0;
  }
}
