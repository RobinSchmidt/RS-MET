#pragma once

/** A naive implementation of an oscillator that can generate sawtooth-, square- and triangle 
waves. It's meant to be used in conjunction with rsStepBandLimiter to turn the initially naive
waveform generation into an anti-aliased one via the blep/blamp technique. The idea is that the 
oscillator produces the naive waveform and along with it the information that is required for 
constructing appropriate correction signals - this information consists of the exact sub-sample 
position and the size of the (step- or corner-) discontinuities. With this information, a driver 
object can correct them via bleps/blamps as a post-processing step. Separating the application of 
the blep-correction out as post-processing step allows for great flexibility. The driver may choose
to use table-based or polynomial bleps, linear- or minimum-phase ones and may also correct 
hard-sync between two oscillators.

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

  inline void setAmplitude(T newAmplitude)
  {
    amp = newAmplitude;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline T getStepAmplitude() { return stepAmp; }

  inline T getStepDelay()     { return stepDelay; }



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline T getSampleSaw()
  {
    stepAmp = T(0);

    //cornerAmp = 0;

    //T y = T(2) * pos - T(1);

    updatePhase();
    // actually, it's not really nice to do the phase update before computing the sample - i think
    // i can do it the other way around, if we make so that the osc will set stepAmp nonzero 
    // *before* getSample produces the step - the driver would have to take this into account and
    // may change the order of calls to the osc and the blepper (call blepper's getSample first)
    // alternative, we could se the pos to 1-inc in reset - but that's unelegant
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
    if(pos < inc) {                // or should it be <= ?
      stepAmp   = T(-2)*amp;       // downward step by -2a
      stepDelay = T(1) - pos/inc;  // is this correct?
    }

    // what if the increment is negative - then we should check if pos > (1-inc) and if so, we have
    // an upward step

    return amp * (T(2) * pos - T(1));
  }

  inline T getSampleSquare()
  {
    stepAmp = T(0);
    updatePhase();


    if(pos < inc) {                // or should it be <= ?
      stepAmp   = T(-2)*amp;       // downward step by -2a
      //stepAmp   = T(-2);       // test
      stepDelay = T(1) - pos/inc;  // is this correct?
    }
    else if(pos > T(0.5) && pos - T(0.5) < inc) 
    {
      stepAmp   = T(2)*amp;                 // upward step by -2a
      //stepAmp   = T(2); 
      stepDelay = T(1) - (T(pos)-0.5)/inc;  // is this correct?
    }


    if(pos <= T(0.5))  // or <?
      return T(-amp);
    else
      return T(amp);


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
  T amp = 1;     // amplitude

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