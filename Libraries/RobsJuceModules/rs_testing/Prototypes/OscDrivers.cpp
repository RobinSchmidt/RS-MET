template<class T, class TOsc, class TBlep>
rsSuperBlepOsc<T, TOsc, TBlep>::rsSuperBlepOsc()
{
  setMaxNumOscillators(8);
}

template<class T, class TOsc, class TBlep>
void rsSuperBlepOsc<T, TOsc, TBlep>::updateIncrements()
{
  // this is the place, where we should experiment with different strategies to compute the 
  // ratios ...maybe have functions like getMetallicRatios, getIntegerSquareRoots, etc., call them
  // here and after that is done, 

  // compute ratios of the increments:
  //rsArray::fillWithValue(&incs[0], numOscs, T(1)); // preliminary

  if(numOscs == 1)
    incs[0] = 1;
  else
  {
    // factor out into fillIncrementRatios:

    // linear progression of increments:
    T minRatio = T(1) - T(0.5) * detune;
    T maxRatio = T(1) + T(0.5) * detune;
    T ratioInc = (maxRatio-minRatio) / T(numOscs-1);
    for(int i = 0; i < numOscs; i++)
      incs[i] = minRatio + i * ratioInc;
  }


  // todo: other options: linear progression of frequencies (increment reciprocals), geometric
  // spacing...hmm - can we somehow "invert" the generalized mean formula to obtain a generalized
  // spreding function - we can do it for arithmetic mean, harmonic mean and geometric mean - but
  // what should be in between?
  // ...



  // see SuperOscillator::setFreq for other spacing strtegies





  

  // todo: maybe compute the mean ratio and use it to center the mean frequency/increment - but 
  // what mean should we use? probably the arithmetic, but the geometric or harmonic seems 
  // plausible as well - maybe make it a user option (use a generalized mean and let the user set 
  // the exponent). actually i think, the harmonic mean of the increments (corresponding to the
  // arithmetic mean of the frequencies) makes more sense - experimentation needed....

  // multiply in the reference increment:
  for(int i = 0; i < numOscs; i++)
    incs[i] *= refInc;
}




template<class T, class TBlep>
void rsSyncOsc<T, TBlep>::reset()
{
  master.reset();
  slave.reset();
  blep.reset();
}



template<class T, class TBlep>
void rsDualBlepOsc<T, TBlep>::reset()
{
  osc1.reset();
  osc2.reset();
  blep1.reset();
  blep2.reset();
}