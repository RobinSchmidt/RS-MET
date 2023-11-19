

SpectralShifter::SpectralShifter(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)
  : SpectralProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  tmpSpectrum = new Complex[getMaxSpectrumSize()];
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
  case     Algo::RobSchmt: shiftViaRS(spec, size); break;
  case     Algo::LaroDols: shiftViaLD(spec, size); break;
  case     Algo::JuilHirs: shiftViaJH(spec, size); break;
  default: RAPT::rsError("Unknown algorithm in SpectralShifter::processSpectrum");
  }
  frameIndex++;
}

void SpectralShifter::shiftViaLD(Complex* spectrum, int spectrumSize)
{
  RAPT::rsError("Not yet implemented");
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
  double k = shift;
  for(int a = 1; a < N; a++)       // we start at 1 because we leave DC as is
  {
    int b = (int) (k*a + 0.5);     // Eq. 1
    if(b >= N) break;              // Avoid reading beyond the end
    Om[b] = OmTmp[a];              // Copy value as explained in section 3.2

    // Optionally apply phase correction according to modified Eq. 2 on top-right on page 3:
    if(phaseFormula == PhaseFormula::useMultiplier)
      Om[b] *= expC(-i * ((double(b-a)*p)/(m*O)) * (2*PI/N));
  }

  // For debugging - plot (partial) spectrum of 8th STFT frame:
  //if(frameIndex == 8) 
  //  rsPlotComplexArrays(spectrumSize/2, (double*)tmpSpectrum, (double*)spectrum);

  //if(p == 8) rsPlotComplexArrays(N/2, (double*)Om);
  //int dummy = 0;


  // ToDo:
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


void SpectralShifter::shiftViaRS(Complex* spectrum, int spectrumSize)
{
  using AT = RAPT::rsArrayTools;

  // Prepare:
  AT::copy(spectrum, tmpSpectrum, spectrumSize);

  // Experimental - for energy normalization:
  double inEnergy = AT::sumOfSquares((double*) spectrum, 2*spectrumSize);
  double inSumAbs = AT::sumOfAbsoluteValues((double*) spectrum, 2*spectrumSize);


  int w;    // write index    ( maybe use double to avoid type conversion in loop)
  for(w = 1; w < spectrumSize; w++)    // we start at 1 because we leave DC as is
  {
    double r = w / shift;              // read position - todo: precompute k = 1/shift

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
      spectrum[b] *= expC(-i * ((double(b-a)*p)/(m*O)) * (2*PI/N));
    }
  }


  // Experimental:
  // A trick to keep the overall spectral energy the same:
  double gain = 1.0;
  double eps  = std::numeric_limits<double>::epsilon();

  //double outSumAbs = AT::sumOfAbsoluteValues((double*) spectrum, 2*spectrumSize);
  //if(inSumAbs > eps && outSumAbs > eps)
  //  gain = inSumAbs / outSumAbs;

  //double outEnergy = AT::sumOfSquares((double*) spectrum, 2*spectrumSize);
  //if(inEnergy > eps && outEnergy > eps)
  //  gain = sqrt(inEnergy / outEnergy);

  //AT::scale((double*) spectrum, 2*spectrumSize, gain);

  // Normalizing the spectrum either with respect to RMS or with respect to sum-of-abs values
  // does not seem to be the solution to the problem that somtimes the amplitude is too low.
  // When shifting an octave down, the signal is actually perfect without this normalization and
  // in this case, the normalization makes the signal too loud. It was just an ad hoc idea anyway.
  // But maybe a similar apporach in the time domain could be more successful

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



/*

Resources:



*/
