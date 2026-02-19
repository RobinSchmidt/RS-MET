#ifndef RAPT_NOISEGENERATOR_H_INCLUDED
#define RAPT_NOISEGENERATOR_H_INCLUDED


/** Implements a simple linear congruential pseudo random number generator. The coefficients are
taken from Numerical Recipies in C, 2nd Ed, page 284.  */

class rsRandomGenerator
{

public:

  /** Sets the current state of the generator. */
  inline void setState(uint32_t newState) { state = newState; }

  /** Returns the current state of the generator. */
  inline uint32_t getState() const { return state; }

  /** Updates the internal state of the integer generator. */
  inline void updateState() { state = 1664525 * state + 1013904223; }


protected:

  static const uint64_t modulus = 4294967296ull;  // Too big for uint32_t so we need uin64_t.

	uint32_t state = 0;

};


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
class rsNoiseGenerator : public rsRandomGenerator
{

public:

  using Base = rsRandomGenerator;

  //static_assert(std::is_floating_point_v<T>, "rsNoiseGenerator requires floating-point T");
  // We could use this to trigger compile-time errors when someone tries to instantiate the class
  // with an integer type T. But I think, we currently may actually do have such an instantiation
  // somewhere (probably using it with getSampleRaw()), so we can't do that at the moment. 
  // Hmmm...but the TestsRosicAndRapt project seems to compile even when the line is not commented
  // out. However, just in case, I leave it commented out for now. Maybe it can be uncommented 
  // later. Or maybe someday we can switch to using a C++20 "float" concept for T. We'll see.

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the seed (initial state) of the PRNG and sets the current state to the seed value. */
  inline void setSeed(uint32_t newSeed) { state = seed = newSeed; }

  /** Sets the seed without resetting the state. */
  inline void setSeedWithoutReset(uint32_t newSeed) { seed = newSeed; }

  /** Sets the range for the numbers to be produced. */
  inline void setRange(T min, T max)
  {
    //scale = T((max-min)/4294967296.0);     // Old
    scale = (max-min) / T(modulus);          // New
    //scale = (max-min) * (T(1)/T(modulus)); // Maybe this could work too and be more efficient?
    shift = min;
  }
  // Maybe have an optional boolean parameter "maxIncluded" that controls if the max is included or
  // not. For example with min = 0 and max = 1, true would mean that we produce numbers in the 
  // closed interval [0,1] and false would mean that the numbers are in the half-open interval 
  // [0,1). I think, to achieve both behaviors, we should divide either by the modulus as is or by 
  // modulus-1. But: If T is single precision float, it may not work correctly due to rounding, I 
  // guess. Worse: whether or not it works may depend on the actual values of min and max. Verify 
  // that! Currently we divide by the modulus itself, so we should get the half open interval 
  // because the maximum possible value for the state is modulus-1. Document all of these behaviors
  // for T = float and T = double and write unit tests that verify these behaviors! Maybe try also
  // to replace the division by a multiplication. Maybe some DSP algorithm wants to modulate the 
  // range of a PRNG at sample rate so it may be worth to optimize this. But this change may also
  // modify the behavior with respect to rounding and hence, half-open or closed interval. By the 
  // way, the C++ standard library produces random numbers in an half-open interval. See:
  // 
  // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution.html
  // https://en.cppreference.com/w/cpp/numeric/random.html  (Not relevant here. Just for info.)
 

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  inline T getMappedState() const { return scale * state + shift; }
  // Maybe we need to use this->state or Base::state. MSVC accepts it as is but GCC may not. Not 
  // sure.

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Produces one output sample at a time. */
  inline T getSample()
  {
    Base::updateState();
    return getMappedState();
  }

  /** Resets the internal state to the seed value. */
  inline void reset() { state = seed; }

  /** Returns a raw integer random sample from the underlying integer linear congruential 
  generator. Here, "raw" means that the mapping function is not yet applied such that the range is 
  from 0 to the maximum value of uint32_t which is 2^32-1. */
  inline uint32_t getSampleRaw()
  {
    Base::updateState();
    return Base::getState();
  }
  // Maybe move this to the baseclass. Add there also function getSample01() or 
  // getSampleInUnitInterval() or getSampleUnitRange() that produces value in the range 0..1. Maybe
  // we should have two versions of this function for the closed and half-open unit interval


protected:

  // By default, we produce numbers in the interval -1..+1:
  T scale = T(2.0 / double(modulus));
  T shift = T(-1);

  uint32_t seed = 0;

  // ToDo: 
  // 
  // - I think, with these default values for scale and shift, the interval of the random numbers 
  //   that are produced is left closed and right open, i.e. the number is in the interval [-1,1). 
  //   Verify and document that! Maybe try to make it such that the default interval is closed to 
  //   both sides, i.e. [-1,+1].
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

  inline void setOrder(uint32_t newOrder) 
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
    //unsigned long long accu = 0;
    uint64_t accu = 0;                           // Accumulator needs extended range
    for(uint32_t i = 1; i <= order; i++) {
      this->updateState();
      accu += this->state; }
    return this->scale * accu + this->shift;
  }

protected:

  void updateCoeffs()
  {
    this->scale = T( (max-min) / (order * 4294967296.0) );
    this->shift = min;
    // The formula for  scale  is different than in the baseclass. Here, we divide by 
    // order * modulus  rather than just by the modulus.
  }

  uint32_t order = 1;
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
