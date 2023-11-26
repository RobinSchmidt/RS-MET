

SpectralShifter::SpectralShifter(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)
  : SpectralProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  int N = getMaxSpectrumSize();
  tmpSpectrum = new Complex[N];
  mag.resize(N);
  phs.resize(N);
  phsOld.resize(N);
}

SpectralShifter::~SpectralShifter()
{
  delete[] tmpSpectrum;
}


void SpectralShifter::processSpectrum(Complex* spec, int size)
{
  using Algo = Algorithm;
  switch(algo)
  {
  case     Algo::RobSchm1: shiftViaRS1(spec, size); break;
  case     Algo::RobSchm2: shiftViaRS2(spec, size); break;
  case     Algo::LaroDols: shiftViaLD( spec, size); break;
  case     Algo::JuilHirs: shiftViaJH( spec, size); break;
  default: RAPT::rsError("Unknown algorithm in SpectralShifter::processSpectrum");
  }
  frameIndex++;
}

void SpectralShifter::shiftViaLD(Complex* spectrum, int spectrumSize)
{
  RAPT::rsError("Not yet implemented");

  // Implements the algorithm explained here:
  //  https://www.ee.columbia.edu/~dpwe/papers/LaroD99-pvoc.pdf
  //
  // Notation used in the paper and in this implementation:
  //   N    : FFT size
  //   R    : hop size
  //   K    : overlap factor (= N/R)
  //   j    : imaginary unit
  //   w    : "omega": normalized radian frequency of input sine (= 2*pi*frequency/sampleRate)
  //   beta : desired frequency scale factor such that         wNew = w * beta
  //   dw   : "Delta omega": desired frequency shift such that wNew = w + dw
  //   Z_u  : twiddle factor Z_u = e^(j*dw*R). Should accumulate: Z_{u+1} = Z_u * dw_{u+1} * R
  //          ...that's a bit unclear. These two formulas seem to conflict.
  //
  // Notation used in the paper that we don't need here:
  //   n           : sample index
  //   u           : frame index (I think)
  //   a           : Not explained but isn't used in the formulas anyway
  //   A           : input sinusoid's ampltiude
  //   phi         : input sinusoid's start phase at n = 0
  //   t_ua        : physical time corresponding to frame index u (I think)
  //   h(n)        : analysis window (typically Hann)
  //   x(n)        : input sine x(n) = A * e^(j(w*n + phi))
  //   H(Om)       : Fourier trafo of h(n)
  //   X(Om, t_ua) : STFT of x(n) at time t_ua




  // Notes:
  // -The two formulas to compute the phase twiddle factors:
  //    Z_u = e^(j*dw*R)  and  Z_{u+1} = Z_u * dw_{u+1} * R
  //  seem to be conflicting and the second doesn't seem to make a whole lot of sense anyway. The
  //  dw_{u+1} * R factor isn't even of unit magnitude. Perhaps they mean something like:
  //    Z_0 = 1, Z_{u+1} = Z_u * e^(j*dw*R)
  //  That looks kinda more plausible for an accumulating phase twiddler. I think, it is this phase
  //  twiddle factor that can be (soft) resetted to 1 on transients.
}

void SpectralShifter::shiftViaJH(Complex* spectrum, int spectrumSize)
{
  // Under Construction - Does not yet work

  // Implements the algorithm explained here:
  // https://www.researchgate.net/publication/261078164_Low_latency_audio_pitch_shifting_in_the_frequency_domain
  //
  // Notation used in the paper and in this implementation:
  //   a    : Source bin index
  //   b    : Destination bin index
  //   Om_x : Complex STFT value "Omega_x" at bin with index x, x is placeholder for a or b
  //   O    : Overlap factor (typically 2,4,8)
  //   N    : FFT size (typically 512..2048)
  //   k    : Frequency scaling factor (typically 0.5..2.0)
  //   p    : STFT frame index


  using AT = RAPT::rsArrayTools;

  // Define algorithm variables using the notation from the paper:
  Complex* Om    = spectrum;       // Omega array, the array of the "phasors"
  Complex* OmTmp = tmpSpectrum;    // Temporary storage for the Omega array
  int      N     = spectrumSize;   // Maybe it should be multiplied by zeroPaddingFactor?
  int      O     = overlapFactor;  
  int      p     = frameIndex;
  int      m     = paddingFactor;
  Complex  i     = Complex(0, 1);  // Imaginary unit

  // Prepare input and output buffers:
  AT::copy(Om, OmTmp, N);          // Prepare temp buffer to read from
  AT::fillWithZeros(Om, N);        // Clear buffer to write into.


  // Main loop over the bins of the current spectrum:
  double k = freqScale;
  for(int a = 1; a < N; a++)       // we start at 1 because we leave DC as is
  {
    int b = (int) (k*a + 0.5);     // Eq. 1
    if(b >= N) break;              // Avoid reading beyond the end
    Om[b] = OmTmp[a];              // Copy value as explained in section 3.2

    // Optionally apply phase correction according to modified Eq. 2 on top-right on page 3:
    if(phaseFormula == PhaseFormula::useMultiplier)
    {
      Om[b] *= expC(-i * ((double(b-a)*p)/(O)) * (2*PI/N));  // for test
      //Om[b] *= expC(-i * ((double(b-a)*p)/(m*O)) * (2*PI/N)); ..i think, it's wrong
      // In the paper in section 3.4, the factor m is introduced because the synthesis FFT size is 
      // chosen to be m times the analysis FFT size. But this factor also affects b like so:
      // b = (int) (m*k*a + 0.5). Yes - in our implementation here, the "m" is already absorbed 
      // into the FFT-size N, so I think, we should not use m in the denominator
    }
  }

  // For debugging - plot (partial) spectrum of 8th STFT frame:
  //if(frameIndex == 8) 
  //  rsPlotComplexArrays(spectrumSize/2, (double*)tmpSpectrum, (double*)spectrum);

  if(p == 8) rsPlotComplexArrays(N/2, (double*)Om);
  //int dummy = 0;


  // ToDo:
  // -Gibe the variable b-a a name: dk (for "delta k"). It's the difference between the new bin b 
  //  and the old bin a, i.e. a frequency difference measured in "number of bins". That's the 
  //  value that LaroDols calls Delta-omega.
  // -The LD paper says in section 3.5: "the phase rotations hosuld be cumulated from one frame to
  //  the next". Maybe such a consideration applies to this algo also? If so, then the phase 
  //  formula needs to be modifeid accordingly. Maybe introduce another field in the enum that 
  //  switches to such an accumulating formula and try it.
  // -Check if skipping the DC bin is the right thing to do in general. Maybe zero out the 
  //  Nyquist freq in case of upshifting (imag part of bin zero).
  // -Check, if zeroing Om in the preparation step is really the right thing to do. What happens
  //  if we just don't do it?
  // -Check what happens, if we do  Om[a] = OmTmp[b];  instead of  Om[b] = OmTmp[a];  I'm not so
  //  sure, which way around they mean it. Hmm - when doing it that way, it seems to apply the 
  //  reciprocal freq scaling, which makes sense. But interestingly, when k < 1 then 
  //  Om[a] = OmTmp[b] does the upward shift by 1/k without amp-modulation artifacts while
  //  Om[b] = OmTmp[a] does the downward 
  //  shift without amp-modulation. Maybe it has to do with overwriting?
  // -Verify the formula for the phase-twiddle factor. Try to rederive it to learn where it's 
  //  coming from. I think, we assume a zero phase at sample zero or at the center of the frame 
  //  with index zero. This is supposed to be a phase difference that results after p frames when 
  //  the two frequencies are W_a and W_b respectively, I think.
  // -Plot the spectra for each frame for inspection
  // -Maybe use double for a and b to avoid the conversions inside the loop
  // -Use interpolation instead of rounding
  // -Try using Om[b] += OmTmp[b] * w. Rationale: if we write into the same bin multiple times, the
  //  values accumulate rather than letting the last value overwrite whatever was there before.
  // -Precompute the phase twiddle factors to save the call to expC
  // -Maybe get rid of the conditional for the phaseFormula. Maybe always use the twiddle. But 
  //  maybe if we keep the conditional, drag it out of the loop. Will that make a difference
}

void SpectralShifter::shiftViaRS1(Complex* spectrum, int spectrumSize)
{
  using AT = RAPT::rsArrayTools;

  // Prepare:
  AT::copy(spectrum, tmpSpectrum, spectrumSize);

  // Experimental - for energy normalization:
  double inEnergy = AT::sumOfSquares((double*) spectrum, 2*spectrumSize);

  double inSumAbs = AT::sumOfAbsoluteValues((double*) spectrum, 2*spectrumSize);
  // But wait - this is |re| + |im|. Maybe we should sum the magnitudes of the complex values
  // instead!


  int w;    // write index    ( maybe use double to avoid type conversion in loop)
  for(w = 1; w < spectrumSize; w++)    // we start at 1 because we leave DC as is
  {
    double r = w / freqScale;              // read position - todo: precompute k = 1/shift

    // Linear interpolation:
    double rFloor = floor(r);
    double rFrac  = r - rFloor;
    int    rInt   = (int) rFloor;
    if(rInt >= spectrumSize - 1)
      break;                           // We are done - quit the loop
    spectrum[w] = (1-rFrac) * tmpSpectrum[rInt] + rFrac * tmpSpectrum[rInt+1];
    // ToDo: factor out, try cubic and maybe quintic interpolation
    // Maybe AT::valueAt?


    // Test: try to use the phase formula from the JH algo:
    if(phaseFormula == PhaseFormula::useMultiplier)
    {
      double a = r;
      int    b = w;
      double p = frameIndex;
      double O = overlapFactor;
      int    m = paddingFactor;
      int    N = spectrumSize;
      Complex i(0,1);

      spectrum[b] *= expC(-i * ((double(b-a)*p)/(O)) * (2*PI/N));  // for test

      //spectrum[b] *= expC(-i * ((double(b-a)*p)/(m*O)) * (2*PI/N));  // normal


      //spectrum[b] *= expC(-i * PI/4);  // test
      //spectrum[b] *= expC(-i * 2.0 * p);  // test

      // I think, we should not include the m in the denominator.
      // But shouldn't the phase-shifts accumulate?
     
    }
  }


  // Experimental:
  // A trick to keep the overall spectral energy the same:
  double gain1 = 1.0;
  double gain2 = 1.0;
  double eps   = std::numeric_limits<double>::epsilon();

  double outSumAbs = AT::sumOfAbsoluteValues((double*) spectrum, 2*spectrumSize);
  if(inSumAbs > eps && outSumAbs > eps)
    gain1 = inSumAbs / outSumAbs;

  double outEnergy = AT::sumOfSquares((double*) spectrum, 2*spectrumSize);
  if(inEnergy > eps && outEnergy > eps)
    gain2 = sqrt(inEnergy / outEnergy);


  AT::scale((double*) spectrum, 2*spectrumSize, 0.47*gain1 + 0.53*gain2);

  // Normalizing the spectrum either with respect to RMS or with respect to sum-of-abs values
  // does not seem to be the solution to the problem that somtimes the amplitude is too low.
  // When shifting an octave down, the signal is actually perfect without this normalization and
  // in this case, the normalization makes the signal too loud. It was just an ad hoc idea anyway.
  // But maybe a similar apporach in the time domain could be more successful

  // The sum-of-abs normalization makes it a bit too loud, the energy-based is too quite. Maybe
  // something in between can be used? OK  - done - I'm using a weighted average that works well 
  // for the sine with period 128. We need to figure out, if this generalizes or is just working
  // well for thsi particular signal
  // See also comments under double inSumAbs = ...

  //// For debug:
  //if(frameIndex == 8) 
  //  rsPlotComplexArrays(spectrumSize/2, (double*)tmpSpectrum, (double*)spectrum);

  int dummy = 0;


  // ToDo:
  // -Multiply spectrum[w] by an appropriate root-of-unity factor as explained in (JH)
  // -When shifting upward, zero out Nyquist freq in tmpSpectrum (imag part of bin 0) before 
  //  reading from it. Maybe we need to zero it out in spectrum even before copying it over into
  //  tmpSpectrum
  // -When shifting downward, maybe the DC bin should be zeroed ou in tmpSpectrum before reading 
  //  from it.
}


void SpectralShifter::shiftViaRS2(Complex* spectrum, int spectrumSize)
{
  // UNDER CONSTRUCTION - is still wrong

  using AT = RAPT::rsArrayTools;

  // Try a magnitude interpolation with a free-running phase that (soft) resets on transients.
  // Maybe use magnitude squared or maybe even dB values (log-magnitude). Let's see what works 
  // best. ....

  int N = spectrumSize;
  int H = getHopSize();
  int B = getBlockSize();
  int P = getZeroPaddingFactor();
  Complex i(0, 1); 

  // Compute magnitudes and phases of current input spectrum:
  for(int k = 1; k < N; k++)
  {
    mag[k] = spectrum[k].getRadius();
    phs[k] = spectrum[k].getAngle();
  }


  // Do a linear interpolation of the magnitudes and use a free-running phase:
  for(int kw = 1; kw < N; kw++)
  {
    double kr = kw / freqScale;              // read position - todo: precompute 1/freqScale

    // Linear interpolation:
    double krFloor = floor(kr);
    double krFrac  = kr - krFloor;
    int    krInt   = (int) krFloor;
    if(krInt >= N - 1)
      break;                           // We are done - quit the loop
    double kMag = (1-krFrac) * mag[krInt] + krFrac * mag[krInt+1];

    // Compute free-running phase:
    double kPhs = phsOld[kw] + (2*PI*kw*H) / N;
    //double kPhs = phsOld[kw] + (2*PI*kw*H) / (P*N);
    phs[kw] = kPhs;
    // VERIFY the formula! I'm not sure about it.
    // Hmm - if we do it like this, we actually do not need the phsOld buffer. The phs buffer would
    // be enough - we could update the value directly there like:
    // phs += (2*PI*H) / (N);


    // Compute an additional desired phase-shift that has the effect of circularly shifting the
    // output buffer in such a way that the energy blip gets concentrated in the inital segement of
    // the zero-padded buffer:

    //int    sampleShift = (frameIndex * blockSize + blockSize/2) % (paddingFactor * blockSize);
    // Works only for overlap = 2

    //int    sampleShift = (frameIndex * blockSize + hopSize) % (paddingFactor * blockSize);

    int    sampleShift = (2*frameIndex * hopSize + hopSize) % (paddingFactor * blockSize);


    double phaseShift  = (-PI * kw * sampleShift) / N; 
    // ToDo: explain this better - it's a bit messy. It has to do with the phase-reference point of 
    // the (padded) buffer. See:
    // https://dsp.stackexchange.com/questions/70909/is-there-a-fft-algorithm-with-the-circular-buffering
    // https://de.mathworks.com/help/matlab/ref/fftshift.html
    // https://numpy.org/doc/stable/reference/generated/numpy.fft.fftshift.html


    // ToDo: include some sort of reset strategy here based on (per bin) transients 

    // Write the new complex value into the complex output:
    spectrum[kw] = kMag * expC(-i * (kPhs + phaseShift));  // Verify the minus!
    //spectrum[kw] = -kMag * expC(-i * (kPhs + phaseShift));
    //spectrum[kw] = kMag * expC(i * (kPhs + phaseShift));
    //spectrum[kw] = kMag * expC(i * (kPhs - phaseShift));




    int dummy = 0;
  }

  // Update phase buffer:
  AT::copy(&phs[0], &phsOld[0], spectrumSize);


  int dummy = 0;
}


/*

Resources:



*/
