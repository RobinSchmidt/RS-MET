
//-------------------------------------------------------------------------------------------------
// construction/destruction:

FiniteImpulseResponseFilter::FiniteImpulseResponseFilter()
{
  kernelLength = 0;
  h            = NULL;
  //setImpulseResponseLength(401);
  setImpulseResponseLength(801);
}

FiniteImpulseResponseFilter::~FiniteImpulseResponseFilter()
{
  delete[] h;
}

//-------------------------------------------------------------------------------------------------
// setup:

void FiniteImpulseResponseFilter::setSampleRate(double newSampleRate)
{
  designer.setSampleRate(newSampleRate);
  updateCoefficients();
}

void FiniteImpulseResponseFilter::setMode(int newMode)
{
  designer.setMode(newMode);
  updateCoefficients();
}

void FiniteImpulseResponseFilter::setFrequency(double newFrequency)
{
  designer.setFrequency(newFrequency);
  updateCoefficients();
}

void FiniteImpulseResponseFilter::setBandwidth(double newBandwidth)
{
  designer.setBandwidth(newBandwidth);
  updateCoefficients();
}

void FiniteImpulseResponseFilter::setImpulseResponseLength(int newLength)
{
  //rassert( isOdd(newLength) ); // currently, only odd lengths are supported 
                                 // ...nah: a differentiator may have even length

  if( newLength != kernelLength )
  {
    kernelLength = newLength;
    delete[] h;
    h = new double[kernelLength];
    updateCoefficients();
  }
}

void FiniteImpulseResponseFilter::setWindowType(int newWindow)
{
  designer.setWindowType(newWindow);
  updateCoefficients();
}

/*
void FiniteImpulseResponseFilter::setImpulseResponse(double *newImpulseResponse, int newLength)
{
  convolver.setImpulseResponse(newImpulseResponse, newLength);
  kernelLength = newLength;
}
*/

//-------------------------------------------------------------------------------------------------
// inquiry:

Complex FiniteImpulseResponseFilter::getTransferFunctionAt(Complex z)
{
  z          = 1.0/z;  // we need powers of z^(-1)
  Complex zn = 1.0;    // accumulator for z^(-n)
  Complex H  = 0.0;    // accumulator for H(z)
  for(int n=0; n<kernelLength; n++)
  {
    H += h[n] * zn;
    zn *= z;
  }
  return H; // == sum over h[n] * z^(-n)
}

void FiniteImpulseResponseFilter::getMagnitudeResponse(
  double *frequencies, double *magnitudes, int numBins, bool inDecibels, bool accumulate)
{
  Complex j(0.0, 1.0);
  for(int k=0; k<numBins; k++)
  {
    double  w = 2*PI*frequencies[k] / designer.getSampleRate();
    Complex z = expC(j*w);   // z = e^(j*w)
    double  m = getTransferFunctionAt(z).getRadius();

    // The following occurs very similarly in other filter's getMagnitudeResponse functions 
    // -> we should factor it out into a fucntion:
    if(      inDecibels == false && accumulate == false) magnitudes[k]  = m;
    else if( inDecibels == false && accumulate == true ) magnitudes[k] *= m;
    else if( inDecibels == true  && accumulate == false) magnitudes[k]  = RAPT::rsAmpToDb(m);
    else if( inDecibels == true  && accumulate == true ) magnitudes[k] += RAPT::rsAmpToDb(m);
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void FiniteImpulseResponseFilter::updateCoefficients()
{
  designer.getImpulseResponse( h, kernelLength);
  convolver.setImpulseResponse(h, kernelLength);
}
