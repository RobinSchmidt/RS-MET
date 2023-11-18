

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
  using AT = RAPT::rsArrayTools;

  AT::copy(spectrum, tmpSpectrum, spectrumSize);



  int w;    // write index    ( maybe use double to avoid type conversion in loop)
  for(w = 1; w < spectrumSize; w++)    // we start at 1 because we leave DC as is
  {
    double r = w * shift;              // read position

    // Linear interpolation:
    double rFloor = floor(r);
    double rFrac  = r - rFloor;
    int    rInt   = (int) rFloor;
    if(rInt >= spectrumSize - 1)
      break;                           // We are done - quite the loop
    spectrum[w] = (1-rFrac) * tmpSpectrum[rInt] + rFrac * tmpSpectrum[rInt+1];
    // ToDo: factor out, try cubic and maybe quintic interpolation
    // Maybe AT::valueAt?
  }





  int dummy = 0;

  // ToDo: 
  // -When shifting upward, zero out Nyquist freq in tmpSpectrum (imag part of bin 0) before 
  //  reading from it. Maybe we need to zero it out in spectrum even before copying it over into
  //  tmpSpectrum
  // -When shifting downward, maybe the DC bin should be zeroed ou in tmpSpectrum before reading 
  //  from it.
}

