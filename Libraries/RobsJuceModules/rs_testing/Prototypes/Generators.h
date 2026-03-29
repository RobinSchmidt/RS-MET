#pragma once

//=================================================================================================

/** This class contains some prototypical implementations of the production of pitch-dithered 
waveforms. By "pitch dithering" I mean a technique where the user wants to produce a periodic 
waveform at some given non-integer period length P but we restrict ourselves to producing cycles of
integer lengths.The dithering comes into play when we produce different integer lengths that 
straddle the actually desired length. For example, if the desired period is P = i.f where i is the
integer part of P and f is its fractional part, we could produce cycles of length i with 
probability (1-f) and cycles of length i+1 with probability f. This would be the simplest case of a
probabilistic pithc dithering algorithm. In this class we implement different algorithms to produce
the integer lengths of the cycles to be produced and their probabilities of production. We also 
implement various functions that can be used actually produce various waveforms such as saw, sine, 
etc. using various probabilistic and deterministic pitch dithering algorithms. Additonally, we 
implement some supporting functionality that an be used to statistically measure the quality of the
algorithms.

...TBC... Class is still under construction.  */

template<class T> 
class rsPitchDitherProto
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Cycle length and probability computation

  /** Structure to represent a probability distribution of cycle lengths. Cycles of lengths 
  L1,L2,L3 will be produced with probabilities p1,p2,p3 respectively. */
  struct CycleDistribution
  {
    T   p1, p2, p3;   // p3 = 1 - (p1 + p2)
    int L1, L2, L3;   // L2 = L1 + 1, L3 = L2 + 1 = L1 + 2
    // We accept the redundancies of p3, L2, L3 because this is just a proof of concept prototype 
    // and they may add some clarity.
  };

  /** Produces a cycle distribution that minimizes the variance of the distribution of cycles for 
  the given period. If the period is an exact integer, the distribution that minimizes the variance
  is the one that always produces the cycle with that given period, i.e. with probability 1. The 
  problem with this distribution is that waveforms with integer periods will sound very different
  from those with non-integer periods. The latter one will sound noisier. The noise variance is 
  greatest for half-integer periods. This distribution is not supposed to be useful in practice. 
  It's just nice to have for reference in experiments. In practice, we want the noisiness of the 
  waveform to be invariant with respect to the fractional part of the period. This is what the 
  distributionEqualVariance() is made for. */
  static void distributionMinVariance(T period, CycleDistribution* cd);

  /** Produces a cycle distribution based on a geometrical overlap consideration. We imagine a 
  ruler that has a segment of unit length for each integer length that can be produced and a slider
  of length 2 that we can slide along the ruler. The probabilities of the 3 lengths that straddle 
  our desired period length are proportional to the overlap of the slider with the segments (the 
  proportionality factor is 0.5). When computing the probabilities like this, the middle length L2
  will always receive a probability of 0.5. ...TBC...  */
  static void distributionViaOverlap(T period, CycleDistribution* cd);
  // Maybe rename "ViaOverlap" to something like equalCenter or equalMiddle or equalMidProb. Some 
  // name that reflects that the probability of the middle length is equalized, i.e. invariant to
  // changes in the fractional part of the period. Or maybe for shorter function names use
  // distribEqMid, ..EqDev, ..EqVar, ..MinVar,  or just
  // distEqMid, distEqDev, distEqVar, distMinVar

  /** Produces a cycle distribution that ensures that the expected absolute error (which we call 
  the deviation) of the length of the cycles is independent from the fractional part of the desired
  period. */
  static void distributionEqualDeviation(T period, CycleDistribution* cd);

  /** Produces a cycle distribution that ensures that the expected squared error (i.e. the variance 
  of the error) of the length of the cycles is independent from the fractional part of the desired 
  period. This distribution turned out experimentally to the right one if the goal is to make the
  noisiness of the waveform invariant with respect to the fractional part of the period. That is: 
  With this distributions, exact integer periods will sound the same as the worst case of 
  half-integer periods. A period P = 100.0 will sound equally noisy as one of P = 100.5. The 
  intermediate cases like 100.3 will of course also sound the same. */
  static void distributionEqualVariance(T period, CycleDistribution* cd);
  // This distribution is the one that produces the best results in the sense that the spectra are
  // most consistent as function of the frcational part of the cycle length. In a production 
  // implementation, we may only need to implement this one. Maybe write this into the 
  // documentation as well.

  // Make a function distributionMinVariance that has always 0 for p1 (or maybe 0 for p3 can also
  // occur in an edge case? But I don't think so.)

  //-----------------------------------------------------------------------------------------------
  // \name Cycle production

  /** Fills the buffer x of length N with one cycle of a sawtooth wave. */
  static void fillSawCycle(T* x, int N, T amp = T(1));

  /** Fills a section of a given "length" of the vector "x" with a saw cycle. The section starts
  at the given "start" value...TBC...  */
  static void fillSawCycle(std::vector<T>& x, int* start, int length, 
    T amp = T(1), int* counter = nullptr);
  // Maybe put the counter last and make it optional - done
  // Maybe take start by value and return the new start. 

  static void fillDitherSawMinVariance(
    std::vector<T>& x, T period, unsigned long seed = 0, T amp = T(1));


  static std::vector<T> getSawCycle(int N, T amp = T(1));


  static std::vector<T> getSawOld(int N, const CycleDistribution& cd, 
    unsigned long seed = 0, T amp = T(1));

  static std::vector<T> getSaw(int N, const CycleDistribution& cd, 
    unsigned long seed = 0, T amp = T(1));

  //static void fillDitherSaw(
  //  std::vector<T>& x, T period, unsigned long seed = 0, T amp = T(1));
  // Rename to fillDitherSawMinError or ..MinAccum which stands for "minimum accumulated error"

  //-----------------------------------------------------------------------------------------------
  // \name Algorithm assessment

  static bool isCycleDistributionValid(T period, const CycleDistribution& cd);

  /** A struct to store various error measures. */
  struct CycleErrorMeasures
  {
    T e1, e2, e3;   // Individual errors for the 3 lengths
    T mae, var;     // Mean absolute error, Variance
  };

  /** Computes the various error measures for a desired noninteger "period" length when we actually
  produce integer period lengths L1,L2,L3 with probabilities p1,p2,p3 respectively. */
  static CycleErrorMeasures getErrorMeasures(T period, int L1, T p1, int L2, T p2, int L3, T p3);
  // Instead of taking L1,p1,L2,p2,L3,p3 take a struct CycleDistribution by const reference.
  // ToDo: Try to move implementation out of the class. But I get compilation errors when trying to
  // do so. There is a commented out-of-class implementation below.

  static CycleErrorMeasures getErrorMeasures(T period, const CycleDistribution& cd)
  { return getErrorMeasures(period, cd.L1, cd.p1, cd.L2, cd.p2, cd.L3, cd.p3); }

};
