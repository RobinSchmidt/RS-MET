
template<class T>
std::vector<T> rsMultiplicativeSynth<T>::renderOutput(int numSamples)
{
  // Create an array of sine-oscillators and set up the frequencies:
  int numStages = getNumPartials();
  using SineIt = RAPT::rsSineIterator<T>; // maybe keep that as member
  std::vector<SineIt> sines(numStages);
  for(int i = 0; i < numStages; i++)
  {
    T wi = 2*PI*baseFreq*opFreqFactors[i] / sampleRate;
    sines[i].setup(wi, 0.0, 1.0);
  }

  // -Synthesize sound according to selected algorithm
  std::vector<T> y(numSamples);

  // Algo 1: each stage i takes as input A the output of the previous stage and as input B the 
  // output of operator i 
  T A, B, P;
  for(int n = 0; n < numSamples; n++)
  {
    P = 0;
    for(int i = 0; i < numStages; i++)
    {
      A = P;
      B = sines[i].getValue();
      P = cmWeightsA[i]*A + cmWeightsB[i]*B + cmWeightsP[i]*A*B;
    }
    y[n] = P;
  }
  // wrap inner loop into function: getSampleCascadedLoToHi


  // ToDo:

  // -post-process (1st order lowpass filter for the 1/f spectrum, 1st order highpass to block DC,
  //  maybe normalize)


  return y;
}


// rename to getNumStages - a stage is a pair of an operator and a combinator
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