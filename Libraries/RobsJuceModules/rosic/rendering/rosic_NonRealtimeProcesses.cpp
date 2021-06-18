//#include "rosic_NonRealtimeProcesses.h"
//#include "../math/rosic_Interpolation.h"
//using namespace rosic;

// Synthesis:

void rosic::synthesizeWaveform(double *x, int N, int shape, double frequency, double sampleRate,
                               double phase, bool antiAlias)
{
  double w = 2*PI*frequency/sampleRate;
  RAPT::rsArrayTools::fillWithZeros(x, N);
  switch( shape )
  {
  case SINE:
    {
      if( !(frequency >= sampleRate/2 && antiAlias == true) )
      {
        for(int n=0; n<N; n++)
          x[n] = sin(w*n + phase);
      }
    }
    break;
  case SAW:
    {
      if( antiAlias == false )
      {
        for(int n=0; n<N; n++)
          x[n] = RAPT::rsSawWave(w*n + phase);
      }
      else
      {
        int k = 1;
        while( k*frequency < sampleRate/2 )
        {
          double a = -2.0 / (k*PI);
          for(int n=0; n<N; n++)
            x[n] += a * sin(k*(w*n+PI) + phase);
          k++;
        }
      }
    }
    break;
  case SQUARE:
    {
      if( antiAlias == false )
      {
        for(int n=0; n<N; n++)
          x[n] = RAPT::rsSqrWave(w*n + phase);
      }
      else
      {
        int k = 1;
        while( k*frequency < sampleRate/2 )
        {
          double a = -4.0 / (k*PI);
          for(int n=0; n<N; n++)
            x[n] += a * sin(k*(w*n+PI) + phase);
          k+=2;
        }
      }
    }
    break;
  case TRIANGLE:
    {
      if( antiAlias == false )
      {
        for(int n=0; n<N; n++)
          x[n] = RAPT::rsTriWave(w*n + phase);
      }
      else
      {
        int    k = 1;
        double s = 1.0; // sign
        while( k*frequency < sampleRate/2 )
        {
          double a = 8.0 / (k*k*PI*PI);
          for(int n=0; n<N; n++)
            x[n] += s * a * sin(k*w*n + phase);
          k +=  2;
          s *= -1.0;
        }
      }
    }
    break;
  }
}

void rosic::synthesizePulseWave(double *x, int N, double frequency, double dutyCycle,
                                double sampleRate, double phase, bool antiAlias)
{
  double w = 2*PI*frequency/sampleRate;
  RAPT::rsArrayTools::fillWithZeros(x, N);
  if( antiAlias == false )
  {
    for(int n=0; n<N; n++)
      x[n] = RAPT::rsPulseWave(w*n + phase, dutyCycle);
  }
  else
  {
    int k = 1;
    while( k*frequency < sampleRate/2 )
    {
      double a = 4.0 * sin(k*PI*dutyCycle) / (k*PI);
      for(int n=0; n<N; n++)
        x[n] += a * cos(k*w*n + phase);
      k++;
    }
  }
}

void rosic::synthesizeDecayingSine(double *x, int N, double frequency, double amplitude,
                                   double decayTime, double startPhase, double sampleRate)
{
  //ModalFilter filter;
  RAPT::rsModalFilter<double, double> filter;
  filter.setModalParameters(frequency, amplitude, decayTime, startPhase, sampleRate);
  filter.reset();

  x[0] = filter.getSample(1.0);
  for(int n=1; n<N; n++)
    x[n] = filter.getSample(0.0);
}

void rosic::synthesizeModal(double *x, int N, Vector frequencies, Vector amplitudes,
                            Vector decayTimes, Vector startPhases, double sampleRate)
{
  //ModalSynthesizer synth;
  //synth.setModalParameters(frequencies, amplitudes, decayTimes, startPhases);

  RAPT::rsAssertFalse("Needs to be updated to work with new RAPT version of the modal synth");

  /*
  Vector attackTimes(frequencies.dim);
  attackTimes.initWithZeros();
  RAPT::rsModalFilterBank<double, double> synth;
  synth.setModalParameters(frequencies, amplitudes, attackTimes, decayTimes, startPhases);

  synth.setSampleRate(sampleRate);
  x[0] = synth.getSample(1.0);
  for(int n=1; n<N; n++)
    x[n] = synth.getSample(0.0);
  RAPT::rsArrayTools::normalize(x, N, 1.0);
  */
}

void rosic::synthesizeModalPluckedString(double *x, int N, double frequency, double sampleRate,
    double decayTime, double decayExponent, double amplitudeExponent, double inharmonicity,
    double phase, double evenAmplitudeScaler)
{
  double f0 = frequency;     // fundamental frequency
  double d0 = decayTime;     // fundamental decay-time
  double h;                  // partial number
  int numModes = (int) floor(sampleRate/(2*frequency));
  Vector f(numModes);
  Vector a(numModes);
  Vector d(numModes);
  Vector p(numModes);
  for(int m = 0; m < numModes; m++)
  {
    h      = m + 1;
    f.v[m] = f0 * h * sqrt(1+inharmonicity*h*h);
    d.v[m] = d0 / pow(h, decayExponent);
    p.v[m] = phase;
    a.v[m] = 1  / pow(h, amplitudeExponent);
    if( RAPT::rsIsEven((int)h) )
      a.v[m] *= evenAmplitudeScaler;
    if( f.v[m] > 0.5*sampleRate ) // could happen due to inharmonicity
      a.v[m] = 0.0;
  }
  synthesizeModal(x, N, f, a, d, p, sampleRate);
}

void rosic::synthesizeModalRodFreeFree(double *x, int N, double frequency, double sampleRate,
                                       double decayTime, double decayExponent,
                                       double amplitudeExponent, double phase)
{
  // calculate frequency ratios (following Dave Benson's Math and Music):
  std::vector<double> lambdas;
  std::vector<double> frequencies;
  double c, r, lambda;
  double fTmp = frequency;
  //double fScale;

  lambdas.push_back(4.7300407448627040260240481);
  lambdas.push_back(7.8532046240958375564770667);
  lambdas.push_back(10.9956078380016709066690325);
  lambdas.push_back(14.1371654912574641771059179);
  lambdas.push_back(17.2787596573994814380910740);
  lambdas.push_back(20.4203522456260610909364112);
  c       = (1+0.5)*PI;
  lambda  = c - pow(-1.0,(double)1)*2*exp(-c) - 4*exp(-2*c);
  c       = pow(lambdas[0], 2);
  r       = (lambda*lambda)/c;
  fTmp    = r * frequency;
  int m = 7;
  while( fTmp <= 0.5*sampleRate )
  {
    c       = (m+0.5)*PI;
    lambda  = c - pow(-1.0,(double)m)*2*exp(-c) - 4*exp(-2*c);
    lambdas.push_back(lambda);
    c       = pow(lambdas[0], 2);
    r       = (lambda*lambda)/c;
    fTmp    = r * frequency;
    m++;
  }
  lambdas.pop_back(); // the last one will give a frequency above fs/2 - dont want

  for(m = 0; m < (int) lambdas.size(); m++)
  {
    lambda  = lambdas[m];
    c       = pow(lambdas[0], 2);
    r       = (lambda*lambda)/c;
    fTmp    = r * frequency;
    frequencies.push_back(fTmp);
  }

  int numModes = (int) frequencies.size();
  //double f0    = frequency;     // fundamental frequency
  double d0    = decayTime;     // fundamental decay-time
  double h;
  Vector f(numModes);
  Vector a(numModes);
  Vector d(numModes);
  Vector p(numModes);
  for(int m=0; m<numModes; m++)
  {
    f.v[m] = frequencies[m];
    h      = f.v[m] / f.v[0];
    d.v[m] = d0 / pow(h, decayExponent);
    p.v[m] = phase;
    a.v[m] = 1  / pow(h, amplitudeExponent);
  }
  synthesizeModal(x, N, f, a, d, p, sampleRate);
}

// Interpolation:

void rosic::upsampleLinear(double *in, int inLength, double *out, int upsamplingFactor)
{
  int offset = 0;
  for(int n = 1; n < inLength; n++)  // should start at n = 0 -> for the case n == 0, we need special treatment
  {
    for(int m = 0; m < upsamplingFactor; m++)
    {
      double d      = 1.0 - (double) m / (double) upsamplingFactor;
      out[offset+m] = RAPT::getDelayedSampleLinear(d, &in[n]);
    }
    offset += upsamplingFactor;
  }

  // tail:
  double tmpBuffer[2];
  tmpBuffer[0] = in[inLength-1];
  tmpBuffer[1] = 0.0;
  for(int m = 0; m < upsamplingFactor; m++)
  {
    double d      = 1.0 - (double) m / (double) upsamplingFactor;
    out[offset+m] = RAPT::getDelayedSampleLinear(d, &tmpBuffer[1]);
  }
}

void rosic::upsampleHermiteAsymmetric1(double *in, int inLength, double *out, int upsamplingFactor,
                                       double shape)
{
  int offset = 0;
  for(int n = 1; n < inLength; n++)
  {
    // special handling for the input sample at index 1 - it has only one predecessor,
    // but the cubic interpolator needs two:
    double *tmpPointer;
    double tmpBuffer[3];
    if( n == 1 )
    {
      tmpBuffer[2] = in[n];
      tmpBuffer[1] = in[n-1];
      tmpBuffer[0] = 0.0;
      tmpPointer = &tmpBuffer[2];
    }
    else
      tmpPointer = &in[n];

    for(int m = 0; m < upsamplingFactor; m++)
    {
      double d      = 1.0 - (double) m / (double) upsamplingFactor;
      out[offset+m] = RAPT::getDelayedSampleAsymmetricHermite1(d, tmpPointer, shape);
    }
    offset += upsamplingFactor;
  }
}

void rosic::upsampleHermiteAsymmetricM(double *in, int inLength, double *out,
                                       int upsamplingFactor, int M, double shape)
{
  int n, m, i;
  int offset        = 0;
  int bufferSize    = M+2;
  double *tmpBuffer = new double[bufferSize];
  double *tmpPointer;

  for(n = 1; n < inLength; n++)
  {
    if( n < M+2 )
    {
      RAPT::rsArrayTools::fillWithZeros(tmpBuffer, M+2);
      for(i = 0; i <= n; i++)
        tmpBuffer[M+1-i] = in[n-i];
      tmpPointer = &tmpBuffer[M+1];
    }
    else
      tmpPointer = &in[n];
    for(m = 0; m < upsamplingFactor; m++)
    {
      double d      = 1.0 - (double) m / (double) upsamplingFactor;
      out[offset+m] = RAPT::getDelayedSampleAsymmetricHermiteM(d, tmpPointer, M, shape);
    }
    offset += upsamplingFactor;
  }

  // tail:
  RAPT::rsArrayTools::fillWithZeros(tmpBuffer, M+2);
  for(i = 1; i <= RAPT::rsMin(M+1, inLength); i++)
    tmpBuffer[M+1-i] = in[inLength-i];
  tmpBuffer[M+1] = 0.0;
  for(m = 0; m < upsamplingFactor; m++)
  {
    double d      = 1.0 - (double) m / (double) upsamplingFactor;
    out[offset+m] = RAPT::getDelayedSampleAsymmetricHermiteM(d, &tmpBuffer[M+1], M, shape);
  }

  delete[] tmpBuffer;
}

// Filtering:

void rosic::filterButterworth(double *x, double *y, int N, double frequency, double sampleRate,
                              int mode, int prototypeOrder, double gain, bool forwardBackward)
{
  rsEngineersFilterMono filter;
  filter.setFrequency(frequency);
  filter.setSampleRate(sampleRate);
  filter.setMode(mode);
  filter.setPrototypeOrder(prototypeOrder);
  filter.setGain(gain);
  for(int n=0; n<N; n++)
    y[n] = filter.getSampleDirect1(x[n]);
  if( forwardBackward == true )
  {
    for(int n=N-1; n>=0; n--)
      y[n] = filter.getSampleDirect1(y[n]);
  }
}


// Others:

/*
void rosic::estimateEnvelope(double *x, double *y, int N, double sampleRate, double attackTime,
                             double releaseTime, int mode, bool forwardBackward)
{

}
*/

void rosic::fft(double *signalBlock, int blockSize, Complex *spectrum, int fftSize)
{
  rassert(fftSize >= blockSize);
  //static FourierTransformerBluestein transformer;  // produces memleak in test
  FourierTransformerBluestein transformer;
  transformer.setBlockSize(fftSize);
  transformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);
  for(int n=0; n<blockSize; n++)
    spectrum[n] = Complex(signalBlock[n]);
  for(int n=blockSize; n<fftSize; n++)
    spectrum[n] = Complex(0.0);
  transformer.transformComplexBufferInPlace(spectrum);
}

void rosic::ifft(Complex *spectrum, int fftSize, Complex *signalBlock)
{
  //static FourierTransformerBluestein transformer;
  FourierTransformerBluestein transformer;
  transformer.setBlockSize(fftSize);
  transformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);
  transformer.setDirection(FourierTransformerRadix2::INVERSE);
  transformer.transformComplexBuffer(spectrum, signalBlock);
}

void rosic::ifftReal(Complex *spectrum, int fftSize, double *signalBlock)
{
  Complex *tmp = new Complex[fftSize];
  ifft(spectrum, fftSize, tmp);
  for(int n=0; n<fftSize; n++)
    signalBlock[n] = tmp[n].re;
  delete[] tmp;
}

void rosic::fftMagnitudesAndPhases(double *signalBlock, int blockSize, double *magnitudes,
                                   double *phases, int fftSize, bool scale)
{
  double s = 1.0;
  if(scale)
    s = fftSize; // or should it be blockSize?

  Complex *spectrum = new Complex[fftSize];
  fft(signalBlock, blockSize, spectrum, fftSize);
  int kMax;

  if( RAPT::rsIsEven(fftSize) ) // superfluous - integer division will take care of this anyway
    kMax = fftSize/2;
  else
    kMax = (fftSize-1)/2;

  for(int k=0; k<=kMax; k++)
    magnitudes[k] = s * spectrum[k].getRadius();
  if( phases != NULL )
  {
    for(int k=0; k<=kMax; k++)
      phases[k] = spectrum[k].getAngle();
  }
  delete[] spectrum;
}


void rosic::signalToRealCepstrum(double *signal, int numSamples, double *cepstrum)
{
  Complex *spectrum = new Complex[numSamples];
  fft(signal, numSamples, spectrum, numSamples);
  for(int n = 0; n < numSamples; n++)
    spectrum[n] = log(spectrum[n].getRadius());
  ifftReal(spectrum, numSamples, cepstrum);
  delete[] spectrum;
}

void rosic::realCepstrumToSignal(double *cepstrum, int numSamples, double *signal)
{
  Complex *spectrum = new Complex[numSamples];
  fft(cepstrum, numSamples, spectrum, numSamples);
  for(int n = 0; n < numSamples; n++)
    spectrum[n] = expC(spectrum[n]);
  ifftReal(spectrum, numSamples, signal);
  delete[] spectrum;
}

void rosic::minimumPhaseReconstruction(double *input, int numSamples, double *output)
{
  double *c = output; // use output buffer for intermediate results in cepstarl domain
  int N = numSamples;

  signalToRealCepstrum(input, numSamples, c);

  // apply cepstral window to make the cepstrum causal (left-to-zero indices are reflected into the
  // second half of the buffer):
  int n;
  if( RAPT::rsIsEven(N) )
  {
    for(n = 1; n < N/2; n++)
      c[n] *= 2;
    for(n = N/2+1; n < N; n++)
      c[n] = 0;
  }
  else
  {
    for(n = 1; n < (N+1)/2; n++)
      c[n] *= 2;
    for(n = (N+1)/2; n < N; n++)
      c[n] = 0;
  }

  realCepstrumToSignal(c, numSamples, output);
}

void rosic::crossCorrelation(double *x, int xLength, double *y, int yLength, double *result)
{
  int length  = xLength + yLength - 1;
  int fftSize = RAPT::rsNextPowerOfTwo(length);
  Complex *X  = new Complex[fftSize];
  Complex *Y  = new Complex[fftSize];

  fft(x, xLength, X, fftSize);
  fft(y, yLength, Y, fftSize);
  for(int k=0; k<fftSize; k++)
    X[k] *= Y[k].getConjugate();  // X is now the cross-spectrum

  ifftReal(X, fftSize, result);

  delete X;  // triggers debug break
  delete Y;

  //delete[] X;  // hangs
  //delete[] Y;
}

double rosic::estimateFundamental(double* /*x*/, int /*N*/, double /*sampleRate*/,
                                  double /*minExpected*/, double /*maxExpected*/)
{
//  int minLag = (int) floor(sampleRate / maxExpected);
//  int maxLag = (int) ceil(sampleRate  / minExpected);
//
//  // dummy instructions to avoid compiler warnings - delete, when function is properly implemented:
//  minLag = maxLag;
//  N = 1;
//  x = 0;

  return 0.0;  // preliminary
}

/*
void rosic::estimateModalParameters(double *x, int N, Vector *frequencies, Vector *amplitudes,
                                    Vector *decayTimes, double sampleRate)
{

}
*/
