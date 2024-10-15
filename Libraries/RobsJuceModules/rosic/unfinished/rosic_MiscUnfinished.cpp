

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

  int B = getBlockSize();
  int P = getZeroPaddingFactor();
  int O = getOverlapFactor();
  int N = spectrumSize;              // == P*B/2
  int H = getHopSize();              // == B/O
  Complex i(0, 1); 

  // Compute magnitudes and phases of current input spectrum:
  for(int k = 1; k < N; k++)
  {
    mag[k] = spectrum[k].getRadius();
    phs[k] = spectrum[k].getAngle();
  }
  // What about bin index 0? It has this special meaning: the real part encodes DC, the imag part 
  // encodes the Nyquist freq. How should we deal with that? Maybe we should just leave it 
  // untouched? Maybe set the Nyquist bin to zero in upshifting?


  // Test - when uncommented, it should give us identity resynthesis:
  //for(int kw = 0; kw < N; kw++)
  //  spectrum[kw] = mag[kw] * expC(i * phs[kw]);
  //return;
  // -The identity resynthesis still works (at least visually from looking at plots) when kw starts
  //  at 1. It seems to still work even for starting at e.g. kw = 10 - as long as all freq 
  //  components are above that index, I guess. But this happens only exactly when an integer 
  //  number of cycles fits into the block. If this isn't the case, we may have leakage into the
  //  DC bin and may need to include it for identity resynthesis.
  // -Replacing the formula by spectrum[kw] = mag[kw] * expC(-i * phs[kw]); has the effect to 
  //  shifting the padded block for resynthesis to the end of the padded buffer. That means, the 
  //  minus would be wrong. The exponent should *not* use the minus.

  //rsPlotVectors(phs, phsOld);

  // OK - now let's try to use only the magnitude from the input spectrum and generate a phase from
  // the previous phase and a phase increment per hop and per bin:
  int sampleShift = rsMod(-frameIndex * H - B/2, P*B);
  // Try to use only positive values: 
  //   sampleShift = -rsMod(frameIndex * H + B/2, P*B);
  //   sampleShift = -(frameIndex*H + B/2) % (P*B);
  //   sampleShift = -((frameIndex % O) * H + B/2) % (P*B);  // avoid overflow


  for(int k = 0; k < N; k++)
  {
    // Grab phase and magnitude of current bin k:
    double pk = phsOld[k];
    double mk = mag[k];      
    // Maybe we should either use magOld[k] or grab the phase after the phase update. Otherwise it 
    // would seem that we use magnitudes of frame M together with phases of frame M-1, right? I'm 
    // not totally sure. Maybe if we do a conditional phase-reset here like:
    //
    //  if( isTransient(k) )
    //    phsOld[k] = phs[k];
    //
    // then it would be appropriate to do it like that?
    //
    // ToDo:
    // -replace the "grab magnitude" code by implementing linear interpolation of mag array


    // Update the free running phases:
    double wk = (PI * k) / N;           // Normalized radian frequency of bin k, N = fftSize/2
    double phsDelta = H * wk;           // Phase change per hop at bin k
    phsOld[k] += phsDelta;              // Update the phase
    while(phsOld[k] >= PI)              // Keep the phase in -pi...+pi
      phsOld[k] -= 2*PI; 
    // ToDo: 
    // -Maybe use fmod for the phase wrapping. Maybe set up a benchmark and then meausre which one
    //  is better.

    // Compute and apply a phase-twiddle factor that shifts the center of energy of the padded 
    // block into its first section:
    double phaseShift = (PI * k * sampleShift) / N;
    spectrum[k] = mk * expC(i * (pk + phaseShift));
    // I think, this may sometimes lead to phase-cancellations between blocks. We adjust the 
    // envelope nicely - but what about the phase of the sine itself?
    //
    // ToDo:
    // -Explain the phaseShift. That was tricky to figure out! It shifts the enveloped wave-packet 
    //  into the right place (namely to the front) of the padded output buffer. Without it, it will
    //  be centered in different places in the output buffer (and therefore be ignored) at 
    //  different frames/blocks
  }
  return;




  // Do a linear interpolation of the magnitudes and use a free-running phase:
  double readScale = 1.0 / freqScale;        // Scaler for read-position wrt write-position 


  for(int kw = 1; kw < N; kw++)
  {
    double kr = kw * readScale;              // Read position

    // Linear interpolation:
    double krFloor = floor(kr);
    double krFrac  = kr - krFloor;
    int    krInt   = (int) krFloor;
    if(krInt >= N - 1)
      break;                           // We are done - quit the loop
    double kMag = (1-krFrac) * mag[krInt] + krFrac * mag[krInt+1];

    // Compute free-running phase:
    //double kPhs = phsOld[kw] + (2*PI*kw*H) / N;
    double kPhs = phsOld[kw] + (0.5*PI*kw*H) / N;      // test
    phs[kw] = kPhs;
    // VERIFY the formula! I'm not sure about it.
    // Hmm - if we do it like this, we actually do not need the phsOld buffer. The phs buffer would
    // be enough - we could update the value directly there like:
    // phs += (2*PI*H) / (N);
    // -I think, the factor 2 might be 2 much because our n here already is blocSize/2
    // -Looking at the output spectral plot, the phase looks like it's moving too fast between the
    //  bins


    // Compute an additional desired phase-shift that has the effect of circularly shifting the
    // output buffer in such a way that the energy blip gets concentrated in the inital segement of
    // the zero-padded buffer:


    int sampleShift = 0;

    //int    sampleShift = (2*frameIndex * hopSize + hopSize) % (paddingFactor * blockSize);

    //int    sampleShift = (2*frameIndex * H + H) % (P*B); // Worked with B=1024, H=512, P=4
    //int    sampleShift = (2*frameIndex * H    ) % (P*B); // Worked with B=1024, H=256, P=4

    //int    sampleShift = (2*frameIndex * H - H) % (P*B); 

    //double phaseShift  = (-PI * kw * sampleShift) / N;
    double phaseShift  = (PI * kw * sampleShift) / N;



    // ToDo: explain this better - it's a bit messy. It has to do with the phase-reference point of 
    // the (padded) buffer. See:
    // https://dsp.stackexchange.com/questions/70909/is-there-a-fft-algorithm-with-the-circular-buffering
    // https://de.mathworks.com/help/matlab/ref/fftshift.html
    // https://numpy.org/doc/stable/reference/generated/numpy.fft.fftshift.html


    // ToDo: include some sort of reset strategy here based on (per bin) transients 

    // Write the new complex value into the complex output:
    //spectrum[kw] = kMag * expC(-i * (kPhs + phaseShift));  // Verify the minus!
    //spectrum[kw] = -kMag * expC(-i * (kPhs + phaseShift));
    //spectrum[kw] = kMag * expC(i * (kPhs + phaseShift));
    //spectrum[kw] = kMag * expC(i * (kPhs - phaseShift));

    spectrum[kw] = kMag * expC(i * (kPhs + phaseShift));



    // The minus is wrong - see the commented test with the identity resynthesis


    // Compensate for bandwidth compression/expansion:
    spectrum[kw] *= readScale;
    // The rationale is that spectral peaks have some bandwidth and if we scale the peaks in 
    // frequency, their bandwidths (and therfore their energy content) will get scaled as well.
    // If a spectral peak it only half as wide after scaling, we should make it twice as tall.
    // Maybe we should raise readScale to some power, like 1/2, to make it an energy-preserving 
    // factor? Maybe experiment with that - maybe give the user a parameter like 
    // setBandwidthCompensation or something like that.
    // Abother idea would be to compute the total spectral input and output energies and compute
    // an appropriate scale factor from this ratio. That strategy would be more adaptive to the 
    // signal. Maybe try both or let the user decide.

    int dummy = 0;
  }



  // Update phase buffer:
  AT::copy(&phs[0], &phsOld[0], spectrumSize);


  int dummy = 0;
}


/*

Resources:
https://github.com/kupix/bungee  Phase vocoder based pitch-shifter/time-stretcher



*/

//=================================================================================================

rsFlatZapper::rsFlatZapper() : allpassChain(maxNumStages)
{
  initSettings(true);
}

void rsFlatZapper::initSettings(bool initAlsoSampleRate)
{
  allpassChain.setNumStages(50);

  freqLo    = 20.0;
  freqHi    = 20000.0;
  freqShape = 0.0;
  qLo       = 1.0;
  qHi       = 1.0;
  qShape    = 0.0;
  mode      = Mode::biquad;

  if(initAlsoSampleRate)
    sampleRate = 44100.0;

  setDirty();
}

void rsFlatZapper::updateCoeffs()
{
  double fsR = 1.0 / sampleRate;
  double *b0 = allpassChain.getAddressB0();
  double *b1 = allpassChain.getAddressB1();
  double *b2 = allpassChain.getAddressB2();
  double *a1 = allpassChain.getAddressA1();
  double *a2 = allpassChain.getAddressA2();

  // Helper function to map the unit interval 0..1 to istelf via a curve determined by our shape
  // parameter. This is used in the computation of the stage-index dependent tunign frequency and
  // Q for the allpass stage at the given index:
  auto shape = [](double x, double shapeParam) 
  { 
    double s = RAPT::rsPow(2.0, shapeParam);  // Slope at x = 0
    double a = (s-1)/(s+1);                   // Function parameter for rational map in -1..+1
    return RAPT::rsRationalMap_01(x, a);
  }; // For convenience
  // ToDo: 
  // -Avoid conversion from s to a. Using s directly leads ot a simpler formula.
  // -Use the more flexible 3-parametric shape from rsLinearFractionalInterpolator. Have 3 
  //  parameters 
  //  -slope: -inf...+inf - this is the current one
  //  -sigmoidVsSpikey: -inf...+inf - controls the center portion of the shape
  //  -asymmetry: compute s2 = s1^(asymmetry-1) - rationale: when asymmetry == 0, the top-right 
  //   slope should be the reciprocal of the left slope. s1 is the "s" in the function above, i.e. 
  //   the slope of the curve at bottom-left.

  // More helper functions:
  using  BQD = BiquadDesigner;
  auto setupBiquadAllpassStage = [&](int i, double f, double q)
  {
    BQD::calculateCookbookAllpassCoeffs(b0[i], b1[i], b2[i], a1[i], a2[i], fsR, f, q);
    a1[i] = -a1[i];  // The design routine uses a different convention for the sign of the
    a2[i] = -a2[i];  // a-coeffs than the class rsBiquadCascade
  };
  auto setupOnePoleAllpassStage = [&](int i, double f)
  {
    BQD::calculateFirstOrderAllpassCoeffs(b0[i], b1[i], b2[i], a1[i], a2[i], fsR, f);
    a1[i] = -a1[i];
  };

  // Error handler:
  auto handleUnknownMode = [&]()
  {
    RAPT::rsError("Unknown mode in rsWhiteZapper::updateCoeffs");
    allpassChain.resetAllCoeffs();
  };


  // The case for one allpass stage must be handled separately to avoid a division by zero. In the
  // case of one stage, we just use the settings for the lowest stage:
  int numStages = allpassChain.getNumStages();
  if(numStages == 1)
  {
    switch(mode)
    {
    case Mode::onePole: setupOnePoleAllpassStage(0, freqLo);       break;
    case Mode::biquad:  setupBiquadAllpassStage( 0, freqLo, qLo);  break;
    default:            handleUnknownMode();
    }
    dirty = false;
    return;
  }

  // All other cases (i.e. numStages != 1), can be handled by the code below. This includes the
  // numStages == 0 case in which case the loops are just not entered at all:
  double scaler = 1.0 / double(numStages-1); 
  switch(mode)
  {

  case Mode::onePole:
  {
    for(int i = 0; i < numStages; i++)
    {
      double p = scaler * i;   // Goes from 0 to 1
      double f = RAPT::rsLinToExp(shape(p, freqShape), 0.0, 1.0, freqLo, freqHi);
      setupOnePoleAllpassStage(i, f);
    }
  } break;

  case Mode::biquad:
  {
    for(int i = 0; i < numStages; i++)
    {
      double p = scaler * i;   // Goes from 0 to 1
      double f = RAPT::rsLinToExp(shape(p, freqShape), 0.0, 1.0, freqLo, freqHi);
      double q = RAPT::rsLinToExp(shape(p, qShape),    0.0, 1.0, qLo,    qHi);
      setupBiquadAllpassStage(i, f, q);
    }
  } break;

  default:
  {
    handleUnknownMode();
  }

  }

  dirty = false;

  // ToDo:
  // -Implement first order mode
  // -Optimize the calls to rsLinToExp. There are certain things can be precomputed which are
  //  recomputed there multiple times, I think. Maybe create a class rsLinToExpMapper that has
  //  a setup(T inMin, T inMax, T outMin, T outMax) function and a map(T x) function. See 
  //  rsMapper in RAPT (Math/Functions/Mappers.h/cpp)
  // -precompute 1 / (numStages-1)
  // -avoid the typecast to double per iteration - maybe use a double variable as second loop 
  //  counter that just counts up a double in parallel. We still need the int i, though because of
  //  the call to setup...AllpassStage which takes an int as first parameter.
  // -Use simd
}

// See also:
// https://kilohearts.com/products/disperser
// https://github.com/robbert-vdh/diopser
// They seem to be similar in what they do. Diopser also seems to have a mor where the allpass 
// tunings are distributed linearly. Try that, too. I think, we just need to replace the calls to
// rsLinToExp with calls to rsLinToLin. But even better would be a continuous parameter to morph 
// between liner and exponential and maybe go beyond. Check the functions used the 
// sineSweepBassdrum() experiment. There's this adjustable exp mapping that is also used for the 
// shape of the envelope segments in the sampler.

//=================================================================================================

const double rsFreqSweeper::refFreq = 50.0;

rsFreqSweeper::rsFreqSweeper()
{
  initSettings(true);
  reset();
}

void rsFreqSweeper::initSettings(bool initAlsoSampleRate)
{
  freqHi      = 10000.0;
  freqLo      =     0.0;
  chirpAmount =     0.0;
  chirpShape  =     0.0;
  sweepTime   =     0.2; 
  phase       =     0.0;
  phaseStereo =     0.0;
  //waveShape   =     0.0;

  if(initAlsoSampleRate)
    sampleRate = 44100.0;

  setDirty();
}

void rsFreqSweeper::updateCoeffs()
{
  double ep = 0.5*(chirpAmount + chirpShape);
  double eq = 0.5*(chirpAmount - chirpShape);
  a = freqHi;
  p = pow(2, ep);
  q = pow(2, eq);
  c = (pow(a/refFreq, 1/q) - 1) / pow(sweepTime, p); // refFreq = a / (1 + c * sweepTime^p)^q
  b = freqLo * pow(c, q);
}

//=================================================================================================

double rsPhaseShaper::powerLaw(double p, double a) 
{
  p = RAPT::rsLinToLin(p, 0.0, 1.0, -1.0, +1.0);   // 0..1 -> -1..+1
  if(p >= 0)
    p =  pow( p, a);
  else
    p = -pow(-p, a);
  p = RAPT::rsLinToLin(p, -1.0, +1.0, 0.0, 1.0);   // -1..+1 -> 0..1
  return p;

  // ToDo:
  // -Optimize: Replace the (moderately expensive) rsLinToLin calls by simpler formulas.
}

//=================================================================================================

rsSweepKicker::rsSweepKicker()
{
  //freqSweeper.setWaveForm(RAPT::rsSin<double>);


  //auto waveFunc = [](double p) { return rsSin<double> };


  initSettings(true);

  //waveParam = -0.5; // test

  // This is a phase-shaped sine wave:
  auto waveFunc = [this](double p) 
  { 
    p = fmod(p, 1.0);
    p = rsPhaseShaper::powerLaw(p, pow(2.0, 2.0 * waveParam));
    // The scaler 2.0 is rather ad hoc. The goal is that the user gets a parameter in -1..+1 where
    // the ends correspond to bright waves

    return sin(2*PI*p);
    // ToDo:
    // -Try to optimize. I think, we need the fmod because due to the phase-offset parameters. 
    //  But maybe a simple if statement is good enough? But maybe not when we use phase-modulation.
    //  In this case, all bets are off. But maybe have an optimized path when waveParam == 0. In 
    //  this case, it's just a sine wave.
  };
  // Maybe try using feedback-FM to turn the wwaveshape from sin to saw. It doesn't need to be ZDF
  // feedback. UDF is good enough

  freqSweeper.setWaveForm(waveFunc);
  reset();
}

void rsSweepKicker::initSettings(bool initAlsoSampleRate)
{
  frqLo       =     0;
  frqLoByKey  =   100;
  frqLoByVel  =     0;
  frqHi       = 10000;
  frqHiByKey  =   100;
  frqHiByVel  =     0;
  swpTm       =     0.2;
  swpTmByKey  =     0;
  swpTmByVel  =     0;
  fadeOutTime =     0;
  waveParam   =     0;
  freqSweeper.initSettings(initAlsoSampleRate);
  fadeOutEnv.setNumFadeSamples(RAPT::rsRoundToInt(fadeOutTime * getSampleRate()));
}

void rsSweepKicker::noteOn(int key, int vel)
{
  double factor;

  // Set up frequencies:
  double maxFreq = 0.5 * getSampleRate(); // We clip the freq at the Nyquist limit
  double freq;
  factor = RAPT::rsMidiKeyAndVelToFreqFactor(key, vel, frqHiByKey, frqHiByVel);
  freq   = RAPT::rsMin(factor * frqHi, maxFreq);
  freqSweeper.setHighFreq(freq);
  factor = RAPT::rsMidiKeyAndVelToFreqFactor(key, vel, frqLoByKey, frqLoByVel);
  freq   = RAPT::rsMin(factor * frqLo, maxFreq);
  freqSweeper.setLowFreq(freq);

  // Set up speed, shape, etc:
  double time = swpTm;
  freqSweeper.setSweepTime(time);


  freqSweeper.reset();
  // I'm not yet sure if always doing a hard reset is the right thing here. Maybe we should reset 
  // only under certain conditions. Maybe only if the fadeOutEnv has reached its end?


  fadeOutEnv.noteOn();
  currentNote = key;
}

void rsSweepKicker::noteOff(int key)
{
  if(key == currentNote)
    fadeOutEnv.noteOff();
  // Without the conditional, it would behave strangely when the player plays a note, then a second
  // note, then releases the first. It would go into realease even though the most recent note is 
  // still held.
}
