﻿#ifndef RAPT_NOISEGENERATOR_H_INCLUDED
#define RAPT_NOISEGENERATOR_H_INCLUDED

/** A simple noise generator based on the linear congruential method. It generates uniformly 
distributed random number in the range that you can set up via setRange. By default, the range is 
between -1 and +1 (not 0 and 1 because this is meant for audio). The underlying integer pseudo 
random number generator is a linear congruential with period length of 2^32. It is based on 
Numerical Recipies in C (2nd edition), page 284. 

\todo 
-make subclasses that produce random numbers with different distributions, for example by adding
 outputs of the underlying basic generator and/or using waveshaping (maybe atanh, sinh could be
 useful shaping functions - something with low slope around the origina would contract values
 near the origin - high slope far away from the origin spreads them out - or maybe the rational
 mapping could be nice - try with FuncShaper - maybe we need a histogram analyzer for that)
-make a colored noise generator by using the SlopeFilter (in rosic - needs to be dragged to rapt)

*/

template<class T>
class rsNoiseGenerator
{

public:

  //rsNoiseGenerator() = default;
	//~rsNoiseGenerator() = default;

  /** Sets the seed (initial state) of the PRNG and sets the current state to the seed value. */
  inline void setSeed(unsigned long newSeed) { state = seed = newSeed; }

  /** Sets the seed without resetting the state. */
  inline void setSeedWithoutReset(unsigned long newSeed) { seed = newSeed; }

  /** Sets the range for the numbers to be produced. */
  inline void setRange(T min, T max)
  {
    scale = T((max-min)/4294967296.0);
    shift = min;
  }

  /** Resets the internal state to the seed value. */
  inline void reset() { state = seed; }

  /** Updates the internal state of the integer PRNG */
  inline void updateState()
  {
    state = (1664525*state + 1013904223) & 4294967295; 
      // the bitmask performs the modulo operation. when unsigned long is 32 bit, it's not 
      // necesarry because then the mod occurs implicitly due to overflow, but when it's 64 bit
      // we need to do it explicitly (on mac it is required);
  }

  /** Produces one output sample at a time */
  inline T getSample()
  {
    updateState();
    return scale * state + shift;
  }

protected:

  T scale = T(2.0/4294967296.0);
  T shift = T(-1);

  unsigned long seed  = 0;
	unsigned long state = 0;
};

//=================================================================================================

/** Subclass of rsNoiseGenerator that creates the noise by adding up several noise samples in order
to approach a gaussian distribution. The order parameter determines how many noise samples are 
added - with only 1: you get the uniform distribution, 2: triangular, 3: sort of parabolic, 
4: looks already rather gaussianish

not yet tested - probably doesn't work yet
*/

template<class T>
class rsNoiseGenerator2 : public rsNoiseGenerator<T>
{

public:

  inline void setOrder(unsigned long newOrder) { order = newOrder; }


  inline T getSample()
  {
    unsigned long accu = 0;
    for(unsigned long i = 1; i <= order; i++) {
      this->updateState();
      accu += this->state; }
    return this->scale * accu + this->shift;
    // todo: scale and shift should depend on order - need to be updated in setOrder, setRange
    // hmm...maybe state needs to be an rsUint64 here to avoid overflow in the accumulation
    // ...maybe make performance tests to figure out, if that makes a difference
  }

protected:

  unsigned long order = 1;

};

// todo: allow for bi-, tri- and multimodal distributions: for example to get 3 bells at -1, 0, 1,
// first select (according to some probability), which bell is 

// getSampleTriModal
//   double selector = selectorGenerator.getSample(); // in 0..1
//   double randomVal = randomGenerator.getSample();
//   if(selector < thresh1)
//     return randomVal + offset1;
//   if(selector > thresh2)
//     return randomVal + offset2; 
//   return randomVal

// thresh1/2 would determine the weights of the 3 modes, for example with thresh1 = 0.3, 
// thresh2 = 0.7, we would have a 30% chance to get a sample of the low mode a 40% chance for the
// middle mode and again a 30% chance for the high mode - we could give the user parameters
// modeCenter, modeSpread, modeSkew



#endif