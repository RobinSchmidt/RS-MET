
template<class T>
std::vector<T> rsMultiplicativeSynth<T>::renderOutput(int numSamples)
{
  std::vector<T> y(numSamples);

  // ToDo:
  // -Create an array of sine-oscillators and set up the frequencies
  // -Synthesize sound according to selected algorithm
  // -post-process (1st order lowpass filter for the 1/f spectrum, 1st order highpass to block DC,
  //  maybe normalize)

  int numPartials = getNumPartials();

  using SineIt = RAPT::rsSineIterator<T>;;
  std::vector<SineIt> sines;


  return y;
}

template<class T>
int rsMultiplicativeSynth<T>::getNumPartials()
{
  size_t n = opFreqFactors.size();
  n = std::min(n, cmWeightsA.size());
  n = std::min(n, cmWeightsB.size());
  n = std::min(n, cmWeightsP.size());
  return (int) n;

  // ToDo: Write a (variadic template) function that takes an arbitrary number of stdd:vector and
  // returns the minimum of all of the lengths. maybe rsMinSize(...vectors...). That could be often 
  // convenient - for example, here. Sdd it to RAPT StandardContainerTools.h
}

/*

ToDo:

Ideas:
-Let operators be a bit more complex: let them have also a decay-time, maybe an attack, too. Maybe
 use a modal filter with white noise as input. Maybe the frqeuency could be time-varying, too. 
 Maybe use an LFO and/or pitch envelope.

*/