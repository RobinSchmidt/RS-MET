




template<class T, class TOsc, class TBlep>
void rsSuperBlepOsc<T, TOsc, TBlep>::updateIncrements()
{
  // this is the place, where we should experiment with different strategies to compute the 
  // ratios ...maybe have functions like getMetallicRatios, getIntegerSquareRoots, etc., call them
  // here and after that is done, multiply in the reference increment
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