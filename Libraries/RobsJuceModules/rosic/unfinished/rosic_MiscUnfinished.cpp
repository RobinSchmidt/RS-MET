

SpectralShifter::SpectralShifter(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)
  : SpectralProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  tmpSpectrum = new Complex[getMaxSpectrumSize()];
}

SpectralShifter::~SpectralShifter()
{
  delete[] tmpSpectrum;
}


void SpectralShifter::processSpectrum(Complex* spectrum, int spectrumSize)
{
  shiftViaRS(spectrum, spectrumSize);

  // ToDo: use switch(algorithm) to switch between the different algorithms
}




void SpectralShifter::shiftViaLD(Complex* spectrum, int spectrumSize)
{

}

void SpectralShifter::shiftViaJH(Complex* spectrum, int spectrumSize)
{
  // Implements the algorithm explained here:
  // https://www.researchgate.net/publication/261078164_Low_latency_audio_pitch_shifting_in_the_frequency_domain
  //
  // Notation used in the paper and in this implementation:
  //   a        : Source bin index
  //   b        : Destination bin index
  //   Om_x     : Complex STFT value "Omega_x" at bin with index x, x is placeholder for a or b
  //   O        : Overlap factor (typically 2,4,8)
  //   N        : FFT size (typically 512..2048)
  //   k        : Frequency scaling factor (typically 0.5..2.0)
  //
  // Other notation from paper not used or needed here:
  //   f_s      : Sample rate (typically 44100)
  //   f_a      : Frequency of input sine
  //   B        : Bandwidth of an FFT bin in Hz (== f_s / N)
  //   phi      : Phase of input sine in first STFT frame
  //   m        : Multiplier for synthesis FFT size (typically 2 or 4)
  //   p        : STFT frame index
  //   E:       : Error ratio between actually synthesized freq and desired freq of output sine
  //   s1,s2,s3 : 3 input sinusoids
  //   eps      : freq difference between the 3 test input sines



}


void SpectralShifter::shiftViaRS(Complex* spectrum, int spectrumSize)
{
  using AT = RAPT::rsArrayTools;

  AT::copy(spectrum, tmpSpectrum, spectrumSize);



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
  }

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
