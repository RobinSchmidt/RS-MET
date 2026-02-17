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

  //-----------------------------------------------------------------------------------------------
  // \name Setup

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

  // ToDo: setState()

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the current state of the linear congruential generator. */
  inline unsigned long getState() const { return state; }


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Produces one output sample at a time */
  inline T getSample()
  {
    updateState();
    return scale * state + shift;
  }

  /** Resets the internal state to the seed value. */
  inline void reset() { state = seed; }

  /** Updates the internal state of the integer PRNG */
  inline void updateState()
  {
    state = (1664525*state + 1013904223) & 4294967295;
    // These numbers are taken from Numerical Recipies in C, 2nd Ed, page 284. The bitmask performs
    // the modulo operation. When unsigned long is 32 bit, it's not necesarry because then the mod
    // occurs implicitly due to overflow, but when it's 64 bit we need to do it explicitly (on Mac
    // it is required). ToDo: either figure out at compile time, if it is required and use 
    // conditional compilation, or (better): make sure that it uses a 32 bit integer type (i.e. use 
    // rsUint32 or uint32_t instead of unsigned long for the state). But when doing that, maybe 
    // benchmark both versions. Maybe using 64 bit integers turns out to be faster because that's
    // the natural word length of the CPU? But maybe that could depend on how the object is aligned
    // in memory and what the type T is? Do the required tests and document the results! Create 
    // also unit tests that test sizeof(rsNoiseGenerator<T>) for T = float and T = double and check
    // there if the sizes are as expected. I think, for T = double, it should be 24 bytes (2 * 8 
    // for each of the double values and 2 * 4 for each of the uint32 values) and for T = float, it
    // should be 16 bytes (for 2 float and 2 uint32 values).

  }

  inline unsigned long getSampleRaw()
  {
    updateState();
    return state;
  }



protected:

  // By default, we produce numbers in the interval -1..+1:
  T scale = T(2.0/4294967296.0);
  T shift = T(-1);

  unsigned long seed  = 0;
	unsigned long state = 0;
  // ToDo: Use uint32_t

  // ToDo: 
  // 
  // - Maybe use:
  // 
  //     static const T modulus = T(4294967296); 
  // 
  //   and replace the occurences of the magic number by that constant.
  // 
  // - I think, with these default values for scale and shift, the interval of the random numbers 
  //   that are produced is left closed and right open, i.e. the number is in the interval [0,1). 
  //   Verify and document that! Maybe try to make it such that the default interval is closed to 
  //   both sides, i.e. [-1,+1].
  //
  // - Use uint32_t instead of unsigned long for the state and seed. Then get rid of the manual
  //   bitmasking in updateState(). Maybe leave the old code as comment for reference.
  //
  // - Factor out a class rsRandomGenerator that has only the state as member variable. To seed 
  //   it, the user can use a function like setState(). It could have a member function 
  //   getSampleInt() which returns the raw int value and getSampleFloat() which returns a numbe 
  //   in the interval [0,1) or maybe [0,1]. Or maybe have both versions. The required scale and 
  //   shift coeffs should be hardcoded such that they take up no space when creating objects.
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

  inline void setOrder(unsigned long newOrder) 
  { 
    rsAssert(newOrder > 0);
    order = newOrder;
    updateCoeffs();
  }

  inline void setRange(T newMin, T newMax) 
  { 
    min = newMin; 
    max = newMax; 
    updateCoeffs(); 
  }

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
    // The formula for  scale  is different than in the baseclass. Here, we divide by 
    // order * modulus  rather than just by the  modulus.
  }

  unsigned long order = 1;
  T min = T(-1.0), max = T(+1.0);

};

//=================================================================================================
/*

Ideas:

- In analogy to updateState(), implement a function "downdateState" of "updateStateReverse" or 
  something like that. The idea is that from the current state x, we cannot only compute the *next*
  state y by the formula y = (a*x + b) % m but also solve the formula for the *previous* state 
  like: x = ((y-b) / a) % m where the division by a is realized as multiplication by the modular 
  inverse of a for the given modulus m. I think, the condition for such a modular inverse to exist 
  is that a and m are coprime which is the case here because a is odd and m is a power of 2. I 
  think, the modular inverse can be found by the extended Euclidean algorithm. Implement this stuff 
  in the experiments and if it works, drag over the code into the library. We could use a PRNG of 
  that kind for the randomization features in the acid sequencer. There, we'd interpret the bit 
  patterns for the accents and slides as k-bit numbers where where the modulus is given by m = 2^k
  with k being the number of steps in the sequence (typically 16). We can give the user next/prev 
  buttons to skip forward/backward through random patterns. Being able to skip back to the previous
  "random" pattern would be a rather unique (and useful) feature. The prev button would be like an 
  "undo" button for a "randomize" function (realized by the "next" button) but we could do it 
  without actually implementing an expensive undo/redo infrastructure. See experiment 
  noiseReverseMode() for some first tests - it seems to work!

- Maybe allow to create correlated noise by doing only one state-update per sample and
  doing the sum over the past N states, like:
    accu -= state;   // subtract old state
    updateState();   // compute new state
    accu += state;   // add new state
    return scale * accu + shift;
  where accu is remembered between calls to getSample. But that's actually a differencing filter, 
  i.e. a highpass, so it will change the spectrum. Maybe such things as filtering should better be
  left to client code.

- Allow for bi-, tri- and multimodal distributions: for example to get 3 bells at -1, 0, 1,
  first select (according to some probability), which bell is used, then generate a sample from 
  that bell's distribution. Maybe like:
   getSampleTriModal
     double selector = selectorGenerator.getSample(); // in 0..1
     double randomVal = randomGenerator.getSample();
     if(selector < thresh1)
       return randomVal + offset1;
     if(selector > thresh2)
       return randomVal + offset2;
     return randomVal;
   thresh1,2 would determine the weights of the 3 modes, for example with thresh1 = 0.3,
   thresh2 = 0.7, we would have a 30% chance to get a sample of the low mode a 40% chance for the
   middle mode and again a 30% chance for the high mode - we could give the user parameters
   modeCenter, modeSpread, modeSkew...there is some prototype code for this in the experiments. 
   See noiseTriModal(). By the way: The term "mode" here refers to a maximum in the probability
   distribution (as in unimodal, bimodal etc.). It does not mean a mode in the sense of modal 
   synthesis. Maybe to realize different variances for the 3 modes, we could do:
    getSampleTriModal
     double selector = selectorGenerator.getSample(); // in 0..1
     double randomVal = randomGenerator.getSample();
     if(selector < thresh1)
       return varLo * randomVal + meanLo;
     if(selector > thresh2)
       return varHi * randomVal + meanHi;
     return varMid * randomVal + meanMid;


- In the modal synthesizer, we could make these chances dependent on the output signal to
  establish a nonlinear, probabilistic feedback loop interaction between exciter and resonator.
  When the output signal value is strongly negative, we should have a high chance of getting a
  positive excitation impulse value (i.e. choose the generator with positive center) and vice
  versa. I think, this kinda simulates the probability for slip/slide events in a bowed string: If
  the string is under tension in one direction it has a higher chance to slip into the other 
  direction. Or maybe we need a trimodal distribution with thresholds like 0.01, 0.09, and also 
  have amplitude weights for the 3 possible outputs, like 1.0, 0.01, 1.0. Then 99% of the time, we
  would select the middle mode and only output a very quiet noise and the remaining 1 % we would 
  see positive or negative spikes.....experimentation needed

- Try using various types of allpass filters (chirp-allpass, comb-allpass, reverb-allpass, etc.).
  These should change the noise characteristic without altering the magnitude spectrum. They will 
  probably have an impact on the distribution. Try to characterize these effects mathematically and
  perceptually.

*/



#endif
