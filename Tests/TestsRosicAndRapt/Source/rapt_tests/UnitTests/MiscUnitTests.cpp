bool testSpectrogramResynthesis(int blockSize, int hopSize, int signalLength, int fftSize, 
  int windowType = RAPT::rsWindowFunction::HANNING_WINDOW_ZN)
{
  bool r = true;  


  return r;
}


bool spectrogramUnitTest()
{
  bool r = true;      // test result

  r &= testSpectrogramResynthesis(128, 64, 500, 128);

  return r;
}