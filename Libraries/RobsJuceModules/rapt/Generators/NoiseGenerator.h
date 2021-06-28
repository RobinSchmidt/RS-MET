#ifndef RAPT_NOISEGENERATOR_H_INCLUDED
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
-maybe make implementations that use different values for the factor and offset (but only such
 values that guarantee the maximum possible period - look up the conditions that must be sasisfied)

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
    // These numbers are taken from Numerical Recipies in C, 2nd Ed, page 284. The bitmask performs
    // the modulo operation. When unsigned long is 32 bit, it's not necesarry because then the mod
    // occurs implicitly due to overflow, but when it's 64 bit we need to do it explicitly (on mac
    // it is required). ToDo: either figure out at compile time, if it is required and use 
    // conditional compilation, or (better): make sure that it uses a 32 bit integer type (i.e. use 
    // rsUint32 instead of unsigned long for the state).
  }

  /** Produces one output sample at a time */
  inline T getSample()
  {
    updateState();
    return scale * state + shift;
  }

  inline unsigned long getSampleRaw()
  {
    updateState();
    return state;
  }

protected:

  T scale = T(2.0/4294967296.0);
  T shift = T(-1);

  unsigned long seed  = 0;
	unsigned long state = 0;
};

//=================================================================================================

/** Subclass of rsNoiseGenerator that creates the noise by adding up several noise samples in order
to approach a Gaussian distribution. The order parameter determines how many noise samples are
added - with only 1: you get the uniform distribution, 2: triangular (piecewise linear), 3: sort of
parabolic spline, 4: cubic spline - looks already rather Gaussian'ish. Note that with higher 
orders, the samples concentrate more toward the center, such that the overall variance goes down 
with the order.

In general, we get the Irwin-Hall distribution, see:
https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution
https://en.wikipedia.org/wiki/Bates_distribution
https://www.youtube.com/watch?v=-2PA7SbWoJ0&t=17m50s (german)

...i have also a sympy notebook somewhere, that computes the convolutions and gives the
distributions as piecewise polynomials

maybe rename to rsNoiseGeneratorIrwinHall  */

template<class T>
class rsNoiseGenerator2 : public rsNoiseGenerator<T>
{

public:

  inline void setOrder(unsigned long newOrder) { order = newOrder;  updateCoeffs();}

  inline void setRange(T newMin, T newMax) { min = newMin; max = newMax; updateCoeffs(); }

  inline T getSample()
  {
    unsigned long long accu = 0;
    for(unsigned long i = 1; i <= order; i++) {
      this->updateState();
      accu += this->state; }
    return this->scale * accu + this->shift;
  }

protected:

  void updateCoeffs()
  {
    this->scale = T( (max-min) / (order*4294967296.0) );
    this->shift = min;
  }

  unsigned long order = 1;
  T min = T(-1.0), max = T(+1.0);

};

// ToDo: maybe allow to create correlated noise by doing only one state-update per sample and
// doing the sum over the past N states, like:
//   accu -= state;   // subtract old state
//   updateState();   // compute new state
//   accu += state;   // add new state
//   return scale * accu + shift;
// where accu is remembered between calls to getSample. But that's actually a differencing filter, 
// i.e. a highpass, so it will change the spectrum. Maybe such things as filtering should better be
// left to client code.

// ToDo: allow for bi-, tri- and multimodal distributions: for example to get 3 bells at -1, 0, 1,
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
// modeCenter, modeSpread, modeSkew...there is some prototype code for this in the experiments
//
// Idea:
// In the modal synthesizer, we could make these chances dependent on the output signal to
// establish a nonlinear, probabilistic feedback loop interaction between exciter and resonator.
// When the output signal value is strongly negative, we should have a high chance of getting a
// positive excitation impulse value (i.e. choose the generator with positive center) and vice
// versa - this sort of simulates the probability for slip/slide events in a bowed string - if the
// string is under tension in one direction it has a higher chance to slip into the other 
// direction. Or Maybe we need a trimodal distribution with thresholds like 0.01, 0.09, and also 
// have amplitude weights for the 3 possible outputs, like 1.0, 0.01, 1.0. Then 99% of the time, we 
// would select the middle mode and only output a very quiet noise and the remaining 1 % we would 
// see positive or negative spikes.....experimentation needed



#endif
