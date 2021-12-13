
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

  // Apply lowpass to obtain 1/f spectrum:
  RAPT::rsOnePoleFilter<T, T> flt;
  flt.setMode(flt.LOWPASS_IIT);
  flt.setCutoff(baseFreq);
  flt.reset();
  for(int n = 0; n < numSamples; n++)
    y[n] = flt.getSample(y[n]);

  // Apply highpass to block DC (but what about subharmonics? maybe the highpass freq should be 
  // adjustable relative to baseFreq):
  flt.setMode(flt.HIGHPASS_BLT);
  flt.setCutoff(baseFreq);
  flt.reset();
  for(int n = 0; n < numSamples; n++)
    y[n] = flt.getSample(y[n]);
  // maybe the highpass should be optional - instead we could just subtract DC in the normalization
  // step

  // Apply gain (this formula is still a bit ad-hoc - todo: figure out, if it's good):
  T p = getSumOfSquaresOfWeights();   // power...i think
  T g = 1 / sqrt(p);
  RAPT::rsArrayTools::scale(&y[0], numSamples, g*gain);

  // ToDo: normalize or apply a gain factor derived from the weight-arrays
  //RAPT::rsArrayTools::normalize(&y[0], numSamples);


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

  // and/or: write a function, that checks, if all vectors have equal length
}

template<class T>
T rsMultiplicativeSynth<T>::getSumOfSquaresOfWeights()
{
  using AT = RAPT::rsArrayTools;
  int N = getNumPartials();
  T sum(0);
  sum += AT::sumOfSquares(&cmWeightsA[0], N);
  sum += AT::sumOfSquares(&cmWeightsB[0], N);
  sum += AT::sumOfSquares(&cmWeightsP[0], N);
  return sum;
}

/*

Algorithms:

LinearForward: Input A is output of previous stage and B is always the next oscillator:
P = 0; A = P;
B = osc1; P = cmb1(A,B);
B = osc2; P = cmb2(A,B);
B = osc3; P = cmb3(A,B);
B = osc4; P = cmb4(A,B);
...
out = P

BinaryTree (shown for 8 oscs):
P_1_1 = cmb1(osc1,osc2)          1st level
P_1_2 = cmb2(osc3,osc4)
P_1_3 = cmb3(osc5,osc6)
P_1_4 = cmb4(osc7,osc8)
P_2_1 = cmb5(P_1_1, P_1_2)       2nd level
P_2_2 = cmb6(P_1_3, P_1_4)
P_3_1 = cmb7(P_2_1, P_2_2)       3rd level
out   = P_3_1

Algos to do:
-LinearBackward: like LinearForward but traversing array in reverse order
-maybe have two lanes of LinearForward and combine them
-what about feedback configurations? can we do zdf? maybe because the equations are actually all
 linear






ToDo:

Ideas:
-Let operators be a bit more complex: let them have also a decay-time, maybe an attack, too. Maybe
 use a modal filter with white noise as input. The oscillator behavior can be recovered by setting
 the decay-time to infinity (or better: decay-rate to 0) and feeding a unit impulse. Decay rates 
 should be specified in dB/sec.
-The weights need modulation, especially envelopes. Maybe use the filter-based AD-envelopes 
 (important, if we want to make the whole thing filter-based instead of oscillator based)
-Maybe the frequencies could be time-varying, too.  Maybe use an LFO and/or pitch envelope.





*/