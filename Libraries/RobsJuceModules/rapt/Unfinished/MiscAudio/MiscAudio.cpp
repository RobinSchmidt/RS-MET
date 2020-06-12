template<class T>
void synthesizeWaveform(T *x, int N, int shape, T frequency, T sampleRate, T phase, bool antiAlias)
{
  T w = 2*PI*frequency/sampleRate;
  rsFillWithZeros(x, N);
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
          x[n] = rsSawWave(w*n + phase);
      }
      else
      {
        int k = 1;
        while( k*frequency < sampleRate/2 )
        {
          T a = -2.0 / (k*PI);
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
          x[n] = rsSqrWave(w*n + phase);
      }
      else
      {
        int k = 1;
        while( k*frequency < sampleRate/2 )
        {
          T a = -4.0 / (k*PI);
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
          x[n] = rsTriWave(w*n + phase);
      }
      else
      {
        int k = 1;
        T s = 1.0; // sign
        while( k*frequency < sampleRate/2 )
        {
          T a = 8.0 / (k*k*PI*PI);
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

template<class T>
void synthesizePulseWave(T *x, int N, T frequency, T dutyCycle, T sampleRate, T phase,
  bool antiAlias)
{
  T w = 2*PI*frequency/sampleRate;
  rsFillWithZeros(x, N);
  if( antiAlias == false )
  {
    for(int n=0; n<N; n++)
      x[n] = rsPulseWave(w*n + phase, dutyCycle);
  }
  else
  {
    int k = 1;
    while( k*frequency < sampleRate/2 )
    {
      T a = 4.0 * sin(k*PI*dutyCycle) / (k*PI);
      for(int n=0; n<N; n++)
        x[n] += a * cos(k*w*n + phase);
      k++;
    }
  }
}

template<class T>
void synthesizeDecayingSine(T *x, int N, T frequency, T amplitude, T decayTime, T startPhase,
  T sampleRate)
{
  rsModalFilter<T,T> filter;
  filter.setModalParameters(frequency, amplitude, decayTime, startPhase, sampleRate);
  filter.reset();
  x[0] = filter.getSample(1.0);
  for(int n=1; n<N; n++)
    x[n] = filter.getSample(0.0);
}

template<class T>
void synthesizeModal(T *x, int N, std::vector<T> frequencies, std::vector<T> amplitudes,
  std::vector<T> decayTimes, std::vector<T> startPhases, T sampleRate)
{
  //rsModalSynthesizer synth;
  rsModalFilterBank<T,T> synth;

  // preliminary: fixed (to zero) attack times - make this a user parameter later:
  std::vector<T> attackTimes(frequencies.dim);

  synth.setModalParameters(frequencies, amplitudes, attackTimes, decayTimes, startPhases);
  synth.setSampleRate(sampleRate);
  x[0] = synth.getSample(1.0);
  for(int n=1; n<N; n++)
    x[n] = synth.getSample(0.0);
  rsNormalize(x, N, 1.0);
}

template<class T>
void synthesizeModalPluckedString(T *x, int N, T frequency, T sampleRate,
  T decayTime, T decayExponent, T amplitudeExponent, T inharmonicity,
  T phase, T evenAmplitudeScaler)
{
  T f0 = frequency;     // fundamental frequency
  T d0 = decayTime;     // fundamental decay-time
  T h;                  // partial number
  int numModes = (int) floor(sampleRate/(2*frequency));
  std::vector<T> f(numModes);
  std::vector<T> a(numModes);
  std::vector<T> d(numModes);
  std::vector<T> p(numModes);
  for(int m = 0; m < numModes; m++)
  {
    h      = m + 1;
    f.v[m] = f0 * h * sqrt(1+inharmonicity*h*h);
    d.v[m] = d0 / pow(h, decayExponent);
    p.v[m] = phase;
    a.v[m] = 1  / pow(h, amplitudeExponent);
    if( rsIsEven((int)h) )
      a.v[m] *= evenAmplitudeScaler;
    if( f.v[m] > 0.5*sampleRate ) // could happen due to inharmonicity
      a.v[m] = 0.0;
  }
  synthesizeModal(x, N, f, a, d, p, sampleRate);
}

template<class T>
void synthesizeModalRodFreeFree(T *x, int N, T frequency, T sampleRate,
  T decayTime, T decayExponent, T amplitudeExponent, T phase)
{
  // calculate frequency ratios (following Dave Benson's Math and Music):
  std::vector<T> lambdas;
  std::vector<T> frequencies;
  T c, r, lambda;
  T fTmp = frequency;
  //T fScale;

  // numbers taken from???
  lambdas.push_back(4.7300407448627040260240481);
  lambdas.push_back(7.8532046240958375564770667);
  lambdas.push_back(10.9956078380016709066690325);
  lambdas.push_back(14.1371654912574641771059179);
  lambdas.push_back(17.2787596573994814380910740);
  lambdas.push_back(20.4203522456260610909364112);
  c       = (1+0.5)*PI;
  lambda  = c - pow(-1.0,(T)1)*2*exp(-c) - 4*exp(-2*c);
  c       = pow(lambdas[0], 2);
  r       = (lambda*lambda)/c;
  fTmp    = r * frequency;
  int m = 7;
  while( fTmp <= 0.5*sampleRate )
  {
    c       = (m+0.5)*PI;
    lambda  = c - pow(-1.0,(T)m)*2*exp(-c) - 4*exp(-2*c);
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

  int numModes = (int) frequencies.size();  // maybe use rsUint32
  //T f0    = frequency;     // fundamental frequency
  T d0    = decayTime;     // fundamental decay-time
  T h;
  std::vector<T> f(numModes);
  std::vector<T> a(numModes);
  std::vector<T> d(numModes);
  std::vector<T> p(numModes);
  for(int m = 0; m < numModes; m++)
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

template<class T>
void upsampleLinear(T *in, int inLength, T *out, int upsamplingFactor)
{
  int offset = 0;
  for(int n = 1; n < inLength; n++)  // should start at n = 0 -> for the case n == 0, we need special treatment
  {
    for(int m = 0; m < upsamplingFactor; m++)
    {
      T d      = 1.0 - (T) m / (T) upsamplingFactor;
      out[offset+m] = getDelayedSampleLinear(d, &in[n]);
    }
    offset += upsamplingFactor;
  }

  // tail:
  T tmpBuffer[2];
  tmpBuffer[0] = in[inLength-1];
  tmpBuffer[1] = 0.0;
  for(int m = 0; m < upsamplingFactor; m++)
  {
    T d      = 1.0 - (T) m / (T) upsamplingFactor;
    out[offset+m] = getDelayedSampleLinear(d, &tmpBuffer[1]);
  }
}

template<class T>
void upsampleHermiteAsymmetric1(T *in, int inLength, T *out, int upsamplingFactor,
  T shape)
{
  int offset = 0;
  for(int n = 1; n < inLength; n++)
  {
    // special handling for the input sample at index 1 - it has only one predecessor,
    // but the cubic interpolator needs two:
    T *tmpPointer;
    T tmpBuffer[3];
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
      T d      = 1.0 - (T) m / (T) upsamplingFactor;
      out[offset+m] = getDelayedSampleAsymmetricHermite1(d, tmpPointer, shape);
    }
    offset += upsamplingFactor;
  }
}

template<class T>
void upsampleHermiteAsymmetricM(T *in, int inLength, T *out,
                                       int upsamplingFactor, int M, T shape)
{
  int n, m, i;
  int offset        = 0;
  int bufferSize    = M+2;
  T *tmpBuffer = new T[bufferSize];
  T *tmpPointer;

  for(n = 1; n < inLength; n++)
  {
    if( n < M+2 )
    {
      rsFillWithZeros(tmpBuffer, M+2);
      for(i = 0; i <= n; i++)
        tmpBuffer[M+1-i] = in[n-i];
      tmpPointer = &tmpBuffer[M+1];
    }
    else
      tmpPointer = &in[n];
    for(m = 0; m < upsamplingFactor; m++)
    {
      T d      = 1.0 - (T) m / (T) upsamplingFactor;
      out[offset+m] = getDelayedSampleAsymmetricHermiteM(d, tmpPointer, M, shape);
    }
    offset += upsamplingFactor;
  }

  // tail:
  rsFillWithZeros(tmpBuffer, M+2);
  for(i = 1; i <= rsMin(M+1, inLength); i++)
    tmpBuffer[M+1-i] = in[inLength-i];
  tmpBuffer[M+1] = 0.0;
  for(m = 0; m < upsamplingFactor; m++)
  {
    T d      = 1.0 - (T) m / (T) upsamplingFactor;
    out[offset+m] = getDelayedSampleAsymmetricHermiteM(d, &tmpBuffer[M+1], M, shape);
  }

  delete[] tmpBuffer;
}

// Filtering:

template<class T>
void filterButterworth(T *x, T *y, int N, T frequency, T sampleRate,
                              int mode, int prototypeOrder, T gain, bool forwardBackward)
{
  rsError("We need to get the high-order IIR filter code into RSLib to make this work again");
  /*
  EngineersFilter filter;
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
  */
}


// Others:

/*
void rosic::estimateEnvelope(T *x, T *y, int N, T sampleRate, T attackTime,
                             T releaseTime, int mode, bool forwardBackward)
{

}
*/

template<class T>
void fft(T *signalBlock, int blockSize, std::complex<T> *spectrum, int fftSize)
{
  rsAssert(fftSize >= blockSize);
  static rsFourierTransformerBluestein<T> transformer;
  transformer.setBlockSize(fftSize);
  transformer.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_FORWARD_TRAFO);
  for(int n=0; n<blockSize; n++)
    spectrum[n] = std::complex<T>(signalBlock[n]);
  for(int n=blockSize; n<fftSize; n++)
    spectrum[n] = std::complex<T>(0.0);
  transformer.transformComplexBufferInPlace(spectrum);
}

template<class T>
void ifft(std::complex<T> *spectrum, int fftSize, std::complex<T> *signalBlock)
{
  static rsFourierTransformerBluestein<T> transformer;
  transformer.setBlockSize(fftSize);
  transformer.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_FORWARD_TRAFO);
  transformer.setDirection(rsFourierTransformerRadix2<T>::INVERSE);
  transformer.transformComplexBuffer(spectrum, signalBlock);
}

template<class T>
void ifftReal(std::complex<T> *spectrum, int fftSize, T *signalBlock)
{
  std::complex<T> *tmp = new std::complex<T>[fftSize];
  ifft(spectrum, fftSize, tmp);
  for(int n=0; n<fftSize; n++)
    signalBlock[n] = tmp[n].re;
  delete[] tmp;
}

template<class T>
void fftMagnitudesAndPhases(T *signalBlock, int blockSize, T *magnitudes,
                                   T *phases, int fftSize)
{
  std::complex<T> *spectrum = new std::complex<T>[fftSize];
  fft(signalBlock, blockSize, spectrum, fftSize);
  int kMax;
  if( rsIsEven(fftSize) )
    kMax = fftSize/2;
  else
    kMax = (fftSize-1)/2;
  for(int k=0; k<=kMax; k++)
    magnitudes[k] = spectrum[k].getRadius();
  if( phases != NULL )
  {
    for(int k=0; k<=kMax; k++)
      phases[k]     = spectrum[k].getAngle();
  }
  delete[] spectrum;
}

template<class T>
void signalToRealCepstrum(T *signal, int numSamples, T *cepstrum)
{
  std::complex<T> *spectrum = new std::complex<T>[numSamples];
  fft(signal, numSamples, spectrum, numSamples);
  for(int n = 0; n < numSamples; n++)
    spectrum[n] = log(spectrum[n].getRadius());
  ifftReal(spectrum, numSamples, cepstrum);
  delete[] spectrum;
}

template<class T>
void realCepstrumToSignal(T *cepstrum, int numSamples, T *signal)
{
  std::complex<T> *spectrum = new std::complex<T>[numSamples];
  fft(cepstrum, numSamples, spectrum, numSamples);
  for(int n = 0; n < numSamples; n++)
    spectrum[n] = rsExpC(spectrum[n]);
  ifftReal(spectrum, numSamples, signal);
  delete[] spectrum;
}

template<class T>
void minimumPhaseReconstruction(T *input, int numSamples, T *output)
{
  T *c = output; // use output buffer for intermediate results in cepstral domain
  int N = numSamples;

  signalToRealCepstrum(input, numSamples, c);

  // apply cepstral window to make the cepstrum causal (left-to-zero indices are reflected into the
  // second half of the buffer):
  int n;
  if( rsIsEven(N) )
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

template<class T>
void crossCorrelation(T *x, int xLength, T *y, int yLength, T *result)
{
  int length  = xLength + yLength - 1;
  int fftSize = rsNextPowerOfTwo(length);
  std::complex<T> *X  = new std::complex<T>[fftSize];
  std::complex<T> *Y  = new std::complex<T>[fftSize];

  fft(x, xLength, X, fftSize);
  fft(y, yLength, Y, fftSize);
  for(int k=0; k<fftSize; k++)
    X[k] *= Y[k].getConjugate();  // X is now the cross-spectrum

  ifftReal(X, fftSize, result);

  delete[] X;
  delete[] Y;
}

// maybe move this into rsArrayTools - it should be useful for finding the maximum of *any* array 
// with subsample precision - there's nothing particular about correlation here - but it makes use
// of rsPolynomial which is in the "Math" folder which is at a higher level of the library's 
// dependency structure than the "Data" folder (to which rsArrayTools belongs)...hmmm
template<class T>
T rsMaxCorrelationLag(T *r, int N)
{
  int nMax = rsArrayTools::maxIndex(r, N);
  if( nMax == 0 || nMax == N-1 )
    return nMax; // no subsample precision possible at ends of the array

  // find exact location of maximum by fitting a parabola through 3 successive correlation values
  // and using the maximum of the parabola:
  T a[3];
  rsPolynomial<T>::fitQuadratic_0_1_2(a, &r[nMax-1]);
  T offset = 0.0;
  if( a[2] != 0.0 )
    offset = -0.5*a[1]/a[2];
  return nMax - 1 + offset; // -1, because fitQuadratic_0_1_2 assumes x-coordinates 0,1,2
}

template<class T>
T rsMaxCorrelationLag(T *x1, T *x2, int N, bool deBias)
{
  // obtain windowed signals:
  T *y1 = new T[N];  // x1, windowed
  T *y2 = new T[N];  // x2, windowed
  rsWindowFunction::hamming(y1, N);  
    // now, y1 contains the window. Hamming seems to be better than Hann and Blackman...hmm...but 
    // not for the cycle marks - maybe the choice of the window should be up to client code - for 
    // whole samples, we may want to use no window at all (i.e. rectangular) and if only initial 
    // sections are used, a one-sided window that tapers off only to the right may be most suitable

  for(int n = 0; n < N; n++)
  {
    y2[n] = y1[n] * x2[n];
    y1[n] = y1[n] * x1[n];
  }

  // obtain cross-correlation sequence of windowed signals:
  T *r = new T[N];
  rsCrossCorrelation(y1, y2, N, r);
  if( deBias == true )
    rsRemoveCorrelationBias(r, N, r);

  // find the maximum of the cross-correlation sequence with subsample precision:
  T lag = rsMaxCorrelationLag(r, N);

  // clean up and return lag:
  delete[] y1;
  delete[] y2;
  delete[] r;
  return lag;
}

template<class T>
T rsGetShiftForBestMatch(T *x1, T *x2, int N, bool deBias)
{
  T lag1 = rsMaxCorrelationLag(x2, x1, N, deBias); // how much x1 lags behind x2
  T lag2 = rsMaxCorrelationLag(x1, x2, N, deBias); // how much x2 lags behind x1

  // If only one of the lags for best match is > 0.0, return it (with the correct sign), if both
  // are > 0.0, return the one with smaller absolute value (also with the correct sign):
  if(lag1 > 0.0 && lag2 > 0.0)
  {
    if(lag2 < lag1)
      return lag2;
    else
      return -lag1;
  }
  if(lag1 > 0.0)
    return -lag1;
  if(lag2 > 0.0)
    return lag2;
  return 0.0;

  // the old logic was buggy - may be deleted soon:
  //if(lag2 < lag1 || lag1 == 0.0) // sometimes, lag1 == 0.0, even when we need to shift (?)
  //  return lag2;
  //else
  //  return -lag1;

  // todo: Maybe the criterion which of the lags should be returned should not be based only on
  // the amount of shift that would be necessary for a best match with shifting in either
  // direction, but on the actual "goodness" of both matches. Instead of choosing the lag with
  // smaller absolute value for the time-shift variable, we could choose the lag with larger
  // actual value of the cross-correlation function evaluated at the respective lag. In
  // rsMaxCorrelationLag, we would not only determine the subsample-precision lag value, but also
  // the value of the cross-correlation function *at* this lag and use this y-value in the decision
  // logic. Maybe rsMaxCorrelationLag should return an rsPoint2D object with the lag value as
  // x-coordinate and the cross-correlation function value as y coordinate at this x. In the
  // (unlikely) event, that these y-values are indeed exactly equal (or maybe within a threshold),
  // the logic above could be used. Or maybe, we could establish a "soft" logic, using the smaller
  // lag whenever some measure of the difference in the lags is greater than some measure in the
  // difference in the function values - it's a mess....

  // todo: for finding the best match for full signals, we may want to do it without windowing or
  // perhaps with one-sided windowing
}

/*
T RSLib::estimateFundamental(T *x, int N, T sampleRate,
                                  T minExpected, T maxExpected)
{
  int minLag = (int) floor(sampleRate / maxExpected);
  int maxLag = (int) ceil(sampleRate  / minExpected);


  // dummy instructions to avoid compiler warnings - delete, when function is properly implemented:
  minLag = maxLag;
  N = 1;
  x = 0;


  return 0.0;  // preliminary
}
*/

/*
void rosic::estimateModalParameters(T *x, int N, Vector *frequencies, Vector *amplitudes,
                                    Vector *decayTimes, T sampleRate)
{

}
*/

