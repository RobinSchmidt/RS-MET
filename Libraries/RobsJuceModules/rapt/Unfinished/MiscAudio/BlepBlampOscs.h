#pragma once

/** An oscillator that can generate naive sawtooth-, square waves and along with the naive waves, 
it produces the information that is required to apply anti-alising via the bandlmited step ("BLEP") 
technique. This additional information consists of the exact sub-sample position and the size of 
the step discontinuities. With this information, a driver object can correct them via bleps as a 
post-processing step. Separating the application of the blep-correction from the waveform 
generation as post-processing step allows for greater flexibility. The driver may choose to use 
table-based or polynomial bleps, linear- or minimum-phase ones and may also correct hard-sync 
between two oscillators. For example, the driver could use it in conjunction with 
rsStepBandLimiter. The oscillator can also generate naive triangle waves along with the information
required to anti-alias them via bandlimited ramps ("BLAMPS").

hmm...maybe for hard-sync, it needs more information from the osc object(s) to figure out the
step position/height?

*/

template<class T>
class rsBlepReadyOscBase
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */



  //inline void setAmplitude(T newAmplitude) { amp = newAmplitude; }

  //inline void setStartPosition(T newPosition) { start = newPosition; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline T getPhase()          const { return pos; }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /* Produces one sawtooth wave output sample at a time using a phase-increment given by inc and 
  assigns the values stepDelay, stepAmp according to any step discontinuity that may have occured 
  between the previous and the current sample. When the stepAmp is nonzero, it means that a step 
  discontinuity has occured and stepDelay will contain the corresponding fractional delay. If no 
  step has occurred, stepAmp will be zero and stepDelay will we left untouched. An outlying driver
  object should check if stepAmp is nonzero after a call to getSample and if a discontinuity did 
  occur, the driver should prepareForStep in a rsStepBandLimiter object. */

  inline T getSampleSaw(T inc, T* stepDelay, T* stepAmp)
  {
    *stepAmp = T(0);
    updatePhase(inc);
    // actually, it's not really nice to do the phase update before computing the sample - i think
    // i can do it the other way around, if we make so that the osc will set stepAmp nonzero 
    // *before* getSample produces the step - the driver would have to take this into account and
    // may change the order of calls to the osc and the blepper (call blepper's getSample first)
    // alternatively, we could set the pos to 1-inc in reset - but that's unelegant
    // or maybe we should use 0.5 as start-phase - start at the zero crossing - that may also 
    // counteract tarnsient artifacts (the very first step goes from 0 to -1, so it has size -1
    // and not -2) ..but maybe we should base the step height on the actual difference between the
    // two samples around the step and not the ideal theoretical step height? then we should 
    // remember the previous sample ...but the blepper already does this - maybe the blepper should
    // figure out the height itself by taking the difference between current and previous sample?
    // ...try both! ..maybe try a squarewave first - in this case, both sizes would be the same
    // but no - the blepper doesn't know the uncorrected previous input sample anymore - the 
    // delayline has stored a corrected sample - i think that would be wrong to use - maybe we
    // should additionally store uncorrected samples for such stepsize estimation tasks


    // Figure out, if (and when) a wrap-around occured. In such case, we produce a step 
    // discontinuity of size -2
    if(pos < inc) {                    // or should it be <= ?
      *stepAmp   = T(-2);               // downward step by -2
      *stepDelay = pos/inc;
    }
    // what if the increment is negative - then we should check if pos > (1-inc) and if so, we have
    // an upward step - maybe avoid division by keeping incInv = 1./inc as member

    return sawValue(pos);
  }

  inline T getSampleSquare(T inc, T* stepDelay, T* stepAmp)
  {
    *stepAmp = T(0);
    updatePhase(inc);

    if(pos < inc) {                    // or should it be <= ?
      *stepAmp   = T(-2);               // downward step by -2
      *stepDelay = pos/inc;
    }
    else if(pos >= T(0.5) && pos - T(0.5) < inc) 
    {
      *stepAmp   = T(2);                // upward step by 2
      *stepDelay = (T(pos)-0.5)/inc;
    }

    return squareValue(pos);
  }

  inline T getSampleTriangle(T inc, T* cornerDelay, T* cornerAmp)  
  {
    *cornerAmp = T(0);
    updatePhase(inc);

    // produce info for blamp, if corner has occurred:
    if(pos < inc) {
      *cornerAmp   = T(8)*inc;
      *cornerDelay = pos/inc;
    }
    else if(pos >= T(0.5) && pos - T(0.5) < inc) {
      *cornerAmp   = -T(8)*inc;
      *cornerDelay = (T(pos)-0.5)/inc;
    }

    return triangleValue(pos);
  }

  /** Not yet finished!
  A waveform consisting of two linear pieces:
  x1(t) = a1*t + b1   for t <  h
  x2(t) = a2*t + b2   for t >= h
  this allows for sawtooth (up and down), triangle and square waves:
  saw-up:    x1(t) = x2(t) =  2*t - 1
  saw-down:  x1(t) = x2(t) = -2*t + 1
  square:    x1(t) = -1, x2(t) = +1
  triangle:  x1(t) = 4*t - 1, x2(t) = -4*t + 3  */
  inline T getSampleTwoPiece(T inc, T* stepCornerDelay, T* stepAmp, T* cornerAmp, 
    T h, T a1, T b1, T a2, T b2)
  {
    *stepAmp   = T(0);
    *cornerAmp = T(0);
    updatePhase(inc);

    // prepare blep and blamp info:
    //...

    return twoPieceValue(pos, h, a1, b1, a2, b2);
  }

  /** Resets the phase. You should pass a start-phase and a phase-increment. The increment is 
  needed because we increment the phase before producing a sample in getSample which implies that 
  we need to reset the phase to the desired start-phase minus the increment in an reset in order to
  have the very first output sample that getSample... produces at the desired start phase. */
  inline void resetPhase(T start, T inc)
  {
    pos = start - inc; // -inc, because we increment pos before producing a sample
    wrapPhase(pos);
  }

  /** Returns the value of a sawtooth wave at given position in [0,1). */
  static inline T sawValue(T pos)
  {
    return T(2) * pos - T(1);
  }
  // rename to sawUpValue

  /** Returns the value of a square wave at given position in [0,1). */
  static inline T squareValue(T pos)
  {
    if(pos < T(0.5))  // or <=? but i think, < is fine
      return T(-1);
    else
      return T(1);
  }
  // todo: introduce duty-cycle parameter and use if pos < dutyCycle

  static inline T triangleValue(T pos)
  {
    if(pos < T(0.5)) 
      return (T(4) * pos - T(1));
    else
      return (T(1) - T(4) * (pos-T(0.5)));
  }
  // maybe invoke the twoPieceValue function from sawValue/squareValue/triangleValue
  // add a sawDown function

  static inline T twoPieceValue(T pos, T h, T a1, T b1, T a2, T b2)
  {
    if(pos < h)
      return a1 * pos + b1;
    else
      return a2 * pos + b2;
  }

  inline void updatePhase(T inc)
  {
    pos += inc;
    wrapPhase(pos);
  }

  static inline void wrapPhase(T& phase)
  {
    while( phase < T(0) )
      phase += T(1);
    while( phase >= T(1) )
      phase -= T(1);
  }


protected:

  T pos   = 0.5;   // position/phase in the range [0,1)

};


//=================================================================================================

/** Extends rsBlepReadyOscBase by having a member variable for the phase increment. */

template<class T>
class rsBlepReadyOsc : public rsBlepReadyOscBase<T>
{

public:

  /** Sets the phase-increment per sample. Should be set to frequency/sampleRate. */
  inline void setPhaseIncrement(T newIncrement)
  {
    inc = newIncrement;
  }

  inline T getPhaseIncrement() const { return inc; }



  inline T getSampleSaw(T* stepDelay, T* stepAmp) 
  { 
    return rsBlepReadyOscBase<T>::getSampleSaw(inc, stepDelay, stepAmp);
  }

  inline T getSampleSquare(T* stepDelay, T* stepAmp) 
  { 
    return rsBlepReadyOscBase<T>::getSampleSquare(inc, stepDelay, stepAmp);
  }

  inline T getSampleTriangle(T* cornerDelay, T* cornerAmp) 
  { 
    return rsBlepReadyOscBase<T>::getSampleTriangle(inc, cornerDelay, cornerAmp);
  }

  inline T getSampleTwoPiece(T* stepCornerDelay, T* stepAmp, T* cornerAmp,
    T h, T a1, T b1, T a2, T b2)
  {
    return rsBlepReadyOscBase<T>::getSampleTwoPiece(
      inc, stepCornerDelay, stepAmp, cornerAmp, h, a1, b1, a2, b2);
  }

  inline void resetPhase(T start = T(0))
  {
    rsBlepReadyOscBase<T>::resetPhase(start, inc);
  }

  inline void reset(T start = T(0))
  {
    resetPhase(start);
  }


protected:

  T inc   = 0;     // phase increment per sample

  // maybe the driver should be an oscillator pair that also allows hardsync
  // maybe the TriSawOsc would be ideally suited for the belp/blamp technique...at leatst, until
  // the shape-function is applied - but even with the shape function, we may be able to figure out
  // target values for the discontinuities of various orders (we just need to evaluate the 
  // derivatives to both sides of the discontinuity - it will include a step, derivative-change,
  // curvature-change and so on)....but maybe that's so expensive that oversampling would be
  // more efficient? we'll see...

};



//=================================================================================================

/*
todo: allow for waveforms that have discontinuities of both types (steps and corners/ramps) - maybe
a 6 segment waveform consisting of 6 linear segments that can morph between square/trapez/triangle,
triangle/saw (like in TrisawOsc) and maybe more, like this:
 
        2   
      *---*             
    1/     \3          
    *       *        *
             \4     /6
               *---*
                 5

the user parameters adjust the relative lengths of the segments in some meaningful way.
...or maybe 1,6 and 3,4 should be both be combined into 1 segment usch thatwe have 4 segments 
overall

and/or: make a PulseSaw oscillator that uses two sawtooth waves (1 up, 1 down) to produce 
pulse-waves - user may continuously adjust phase-shift (determines duty-cycle) and mix (morphs 
between saw-up and saw-down with pulse-wave in between) - maybe with an adjustable integrator 
filter to allow for triangle, too as integrated square-wave...this filter could actually perhaps
also include a highpass and/or allpass for further shaping options
maybe the two saws could be also detuned a bit for further flexibility



other idea: splice together the waveform from just two linear segments:

x1(t) = a1 * t + b1    for t in 0.0...0.5
x2(t) = a2 * t + b2    for t in 0.5...1.0

this allows for sawtooth (up and down), triangle and square waves:
saw-up:    x1(t) = x2(t) =  2*t - 1
saw-down:  x1(t) = x2(t) = -2*t + 1
square:    x1(t) = -1, x2(t) = +1
triangle:  x1(t) = 4*t - 1, x2(t) = -4*t + 3
see: https://www.desmos.com/calculator/weda63qpdr for saws and triangle

-this can be generalized to move the midpoint point from 0.5 to an arbitrary position
-the coefficents (a1,b1,a2,b2) can be interpolated to give intermediate waveforms
 -maybe they can be suitably normalized (with respect to peak value or energy)
-maybe one user parameter should be saw-vs-square (maybe sawup->square->sawdown) and another
 square-or-saw vs triangle (some sort of "softness" parameter)
-between square and triangle, we should see trapezoidal waves
-non-square pulse-waves are obtained by moving the midpoint
-can be perfectly anti-aliased via blep + blamp
-or: do it like in the trisaw (sawup->triangle->sawdown) but then use amplify+clip to obtain the 
 squareish-stuff via a 2nd parameter...but that would effectively create 2 more segments (the 2 
 flat sections) which would make anti-aliasing about two times more expensive (it doubles the 
 number required of bleps/blamps)
-but from a perceptual perspective, having saw-vs-square in one user parameter and another softness
 parameter would be more desirable - the softness is less important because similar effects can be 
 obtained by a 1st order lowpass - but with sync, it's a different story - we want to sync triangle 
 waves



*/