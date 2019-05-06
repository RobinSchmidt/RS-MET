




template<class T, class TOsc, class TBlep>
void rsSuperBlepOsc<T, TOsc, TBlep>::updateIncrements()
{
  // this is the place, where we should experiment with different strategies to compute the 
  // ratios ...maybe have functions like getMetallicRatios, getIntegerSquareRoots, etc., call them
  // here and after that is done, multiply in the reference increment

  // compute ratios of the increments:
  rsArray::fillWithValue(&incs[0], numOscs, T(1)); // preliminary

  // todo: maybe compute the mean ratio and use it to center the mean frequency/increment - but 
  // what mean should we use? probably the arithmetic, but the geometric or harmonic seems 
  // plausible as well - maybe make it a user option (use a generalized mean and let the user set 
  // the exponent). actually i think, the harmonic mean of the increments (corresponding to the
  // arithmetic mean of the frequencies) makes more sense - experimentation needed....


  for(int i = 0; i < numOscs; i++)
    incs[i] *= inc;
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