#pragma once

/** An oscillator that can generate naive sawtooth-, square waves and along with the naive waves, 
it produces the information that is required to apply anti-alising via the bandlmited step ("BLEP") 
technique. This additional information consists of the exact sub-sample position and the size of 
the step discontinuities. With this information, a driver object can correct them via bleps as a 
post-processing step. Separating the application of the blep-correction from the waveform 
generation as post-processing step allows for greater flexibility. The driver may choose to use 
table-based or polynomial bleps, linear- or minimum-phase ones and may also correct hard-sync 
between two oscillators. For example, the driver could use it in conjunction with 
rsStepBandLimiter. 

hmm...maybe for hard-sync, it needs more information from the osc object(s) to figure out the
step position/height?

*/

template<class T>
class rsBlepReadyOsc
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the phase-increment per sample. Should be set to frequency/sampleRate. */
  inline void setPhaseIncrement(T newIncrement)
  {
    inc = newIncrement;
  }

  //inline void setAmplitude(T newAmplitude) { amp = newAmplitude; }

  //inline void setStartPosition(T newPosition) { start = newPosition; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline T getStepAmplitude()  const { return stepAmp; }

  inline T getStepDelay()      const { return stepDelay; }

  inline T getPhaseIncrement() const { return inc; }

  inline T getPhase()          const { return pos; }





  //-----------------------------------------------------------------------------------------------
  /** \name Processing */



  inline T getSampleSaw()
  {
    stepAmp = T(0);
    updatePhase();
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
      stepAmp   = T(-2);               // downward step by -2
      stepDelay = pos/inc;
    }
    // what if the increment is negative - then we should check if pos > (1-inc) and if so, we have
    // an upward step - maybe avoid division by keeping incInv = 1./inc as member

    return sawValue(pos);
    //return (T(2) * pos - T(1));
  }

  inline T getSampleSquare()
  {
    stepAmp = T(0);
    updatePhase();

    if(pos < inc) {                    // or should it be <= ?
      stepAmp   = T(-2);               // downward step by -2
      stepDelay = pos/inc;
    }
    else if(pos >= T(0.5) && pos - T(0.5) < inc) 
    {
      stepAmp   = T(2);                // upward step by 2
      stepDelay = (T(pos)-0.5)/inc;
    }

    return squareValue(pos);

    /*
    if(pos < T(0.5))  // or <?
      return T(-1);
    else
      return T(1);
    */
  }

  inline void reset(T start = T(0))
  {
    pos = start - inc; // -inc, because we increment pos before producing a sample
    wrapPhase();
  }

  /** Returns the value of a sawtooth wave at given position in [0,1). */
  static inline T sawValue(T pos)
  {
    return T(2) * pos - T(1);
  }

  /** Returns the value of a square wave at given position in [0,1). */
  static inline T squareValue(T pos)
  {
    if(pos < T(0.5))  // or <=? but i think, < is fine
      return T(-1);
    else
      return T(1);
  }
  // todo: introduce duty-cycle parameter and use if pos < dutyCycle


  inline void updatePhase()
  {
    pos += inc;
    wrapPhase();
  }

  inline void wrapPhase()
  {
    while( pos < T(0) )
      pos += T(1);
    while( pos >= T(1) )
      pos -= T(1);
  }



protected:

  T pos   = 0.5;   // position/phase in the range [0,1)
  T inc   = 0;     // phase increment per sample



  //T start = 0.5;   // start positon
  //T amp   = 1;     // amplitude
  // maybe get rid of these two variables, handle start by passing it as parameter to reset - makes
  // the classe's memory footprint smaller - relevant when we use arrays of oscs later




  // Values for discontinuities - when the stepAmp is nonzero, it means that a step discontinuity 
  // has occured and stepDelay will contain the corresponding fractional delay. An outlying driver 
  // object should check if getStepAmplitude is nonzeor after a call to getSample and if a 
  // discontinuity did occur, the driver should prepareForStep in a rsStepBandLimiter object. The 
  // required fractionla delay can be inquired via getStepDelay
  T stepDelay = 0;
  T stepAmp   = 0;


  // maybe the driver should be an oscillator pair that also allows hardsync
  // maybe the TriSawOsc would be ideally suited for the belp/blamp technique...at leatst, until
  // the shape-function is applied - but even with the shape function, we may be able to figure out
  // target values for the discontinuities of various orders (we just need to evaluate the 
  // derivatives to both sides of the discontinuity - it will include a step, derivative-change,
  // curvature-change and so on)....but maybe that's so expensive that oversampling would be
  // more efficient? we'll see...
};

//=================================================================================================

/** Oscillator that can generate naive triangle waves along with the information required to 
anti-alias them via bandlimited ramps ("BLAMPS"). Works analogously to the rsBlepReadyOsc. */

template<class T>
class rsBlampReadyOsc : public rsBlepReadyOsc<T>
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline T getCornerAmplitude() { return cornerAmp; }


  inline T getCornerDelay() { return stepDelay; }
  // we don't need a new variable for the corner delay - we can re-use the inherited stepDelay 
  // member for this new purpose - but for client code consistency, we define a separate function 
  // to inquire the corner delay

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline T getSampleTriangle()
  {
    stepAmp = cornerAmp = T(0);  // reset info for blamp (also for blep, but that's only relevant,
    updatePhase();               // if the driver looks for the blep info, which it may not for triangle)

    // produce info for blamp, if corner has occurred:
    if(pos < inc) {
      cornerAmp = T(8)*inc;
      stepDelay = pos/inc;
    }
    else if(pos >= T(0.5) && pos - T(0.5) < inc) {
      cornerAmp = -T(8)*inc;
      stepDelay = (T(pos)-0.5)/inc;
    }

    // produce output signal:
    if(pos < T(0.5)) 
      return (T(4) * pos - T(1));
    else
      return (T(1) - T(4) * (pos-T(0.5)));
  }

protected:


  T cornerAmp = 0;

};


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

*/