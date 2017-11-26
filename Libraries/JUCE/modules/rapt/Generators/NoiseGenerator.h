#ifndef RAPT_NOISEGENERATOR_H_INCLUDED
#define RAPT_NOISEGENERATOR_H_INCLUDED

/** A simple noise generator based on the linear congruential method. It generates uniformly 
distributed random number between in the range that you can set up via setRange, by default,
the range is between -1 and +1. The underlying integer pseudo random number generator is a linear 
congruential with period length of 2^32. It is based on Numerical Recipies in C (2nd edition), 
page 284. 

\todo 
-make subclasses that produce random numbers with different distributions, for example by adding
 outputs of the underlying basic generator
-make a colored noise generator by using the SlopeFilter (in rosic - needs to be dragged to rapt)

*/

template<class T>
class rsNoiseGenerator
{

public:

  rsNoiseGenerator() = default;
	~rsNoiseGenerator() = default;

  /** Sets the seed (initial state) of the PRNG. Note that this function doesn't reset the 
  state. */
  inline void setSeed(unsigned long newSeed) { seed = newSeed; }

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

#endif