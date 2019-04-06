#pragma once

/** A naive implementation of an oscillator that can generate sawtooth-, square- and triangle 
waves. It's meant to be used in conjunction with rsStepBandLimiter to turn the initially naive
wavefrom generation into an anti-aliased one via the blep/blamp technique.

maybe rename to rsBlepReadyOsc

*/

template<class T>
class rsNaiveOsc
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the phase-increment per sample. Should be set to frequency/sampleRate. (times 2?) */
  inline void setPhaseIncrement(T newIncrement)
  {
    inc = newIncrement;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline T getStepAmplitude()
  {
    return stepAmp;
  }

  inline T getStepDelay()
  {
    return stepDelay;
  }



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline T getSampleSaw()
  {
    stepAmp = 0;

    //cornerAmp = 0;

    //T y = T(2) * pos - T(1);

    updatePhase();
    if(pos < inc) {      // or should it be <= ?
      stepAmp = T(-2);   // downward step by -2
      //stepDelay = ???;
    }
    // actually, it's not really nice to do the phase update before computing the sample - i think
    // i can do it the other way around, if we make so that the osc will set stepAmp nonzero 
    // *before* getSample produces the step - the driver would have to take this into account and
    // may change the order of calls to the osc and the blepper (call blepper's getSample first)
    // alternative, we could se the pos to 1-inc in reset - but that's unelegant


    // what if the increment is negative - then we should check if pos > (1-inc) and if so, we have
    // an upward step

    return T(2) * pos - T(1);
  }

  inline T getSampleSquare()
  {
    return T(0);
  }

  /*
  inline T getSampleTriangle()
  {
    return T(0);
  }
  */

  inline void updatePhase()
  {
    pos += inc;
    while( pos < T(0) )
      pos += T(1);
    while( pos >= T(1) )
      pos -= T(1);
    // maybe, here we should update the phase wrap-around flag
  }

protected:

  T pos = 0;     // position/phase in the range [0,1)
  T inc = 0;     // phase increment per sample

  // Values for discontinuities - when the amplitudes are nonzero, it means that a discontinuity of 
  // the respective kind has occured. These should be read out by an outlying driver object after a 
  // call to getSample and if a discontinuity did occur, the driver should 
  // prepareForStep/prepareForCorner in a rsStepBandLimiter object:
  T stepDelay   = 0;
  T stepAmp     = 0;
  //T cornerDelay = 0;
  //T cornerAmp   = 0;
  // maybe we don't need separate delay values for steps and corners? the two types of 
  // discontinuities always occur simultaneously

  // maybe the driver should be an oscillator pair that also allows hardsync
  // maybe the TriSawOsc would be ideally suited for the belp/blamp technique...at leatst, until
  // the shape-function is applied - but even with the shape function, we may be able to figure out
  // target values for the discontinuities of various orders (we just need to evaluate the 
  // derivatives to both sides of the discontinuity - it will include a step, derivative-change,
  // curvature-change and so on)....but maybe that's so expensive that oversampling would be
  // more efficient? we'll see...

  //bool wrapOccurred = false;




};