#pragma once


/** Implements a pair of two sawtooth/pulse oscillators whose waveforms get anti-aliased via a
blep object. It also implements hardsync via bleps... */

template<class T, class TBlep> // T: type for signal and parameter, TBlep: class for BLEP object
class rsDualBlepOsc
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */


  void setPhaseIncrement1(T newIncrement) { osc1.setPhaseIncrement(newIncrement); }


  void setPhaseIncrement2(T newIncrement) { osc2.setPhaseIncrement(newIncrement); }

  /** Lets osc1 act as sync master for osc2. Whenever osc1 has completed a full cycle, osc2 will 
  reset its phase (in addition to its own phase resets). */
  void setSync12(bool shouldSync) { sync12 = shouldSync; }

  //void setSync21(bool shouldSync) { sync21 = shouldSync; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline void getSamplePair(T* x1, T* x2)
  {
    // preliminary (later switch waveforms)
    *x1 = osc1.getSampleSaw();

    *x2 = osc2.getSampleSaw(); // wait - we should call this *after* a potential reset

    // apply blep corrections to waveform discontinuities:
    if(osc1.getStepAmplitude() != 0.0)
    {
      blep1.prepareForStep(osc1.getStepDelay(), osc1.getStepAmplitude());
      if(sync12)
      {
        //osc2.reset(); 
        // preliminary - it should actually reset to its start-phase plus something where that 
        // something is determined by the current phasor of osc1 (the reset has occured some time 
        // before the "now" time instant)

        T oldPhase = osc2.getPhase();

        // i think, the correct phase for osc2 is osc1.pos * osc2.inc / osc1.inc
        T newPhase = osc1.getPhase() * osc2.getPhaseIncrement() / osc1.getPhaseIncrement();
        osc2.reset(newPhase);
        // maybe, we should have a special function osc2.syncReset(phs)

        T stepAmp = osc2.sawValue(oldPhase) - osc2.sawValue(newPhase); // or the other way around?

        blep2.prepareForStep(newPhase, stepAmp);


        // todo: determine the precise time instant of the reset and the associated amplitude jump
        // and prepare blep2 accordingly ...or maybe the call to reset() should set up the step
        // variables in osc2? ...maybe that's better because otherwise we may sometimes prepare the
        // blep2 twice?

      }
    }
    if(osc2.getStepAmplitude() != 0.0)
      blep2.prepareForStep(osc2.getStepDelay(), osc2.getStepAmplitude());

    // commented out for test:
    //*x1 = blep1.getSample(*x1);
    //*x2 = blep2.getSample(*x2);

    // todo: apply sync:
    // ...
  }


  inline T getSample()
  {
    T x1, x2;
    getSamplePair(&x1, &x2);
    return x1 + x2;  // maybe have mixing coefficients here? ...or maybe custom mixing be handled
                     // by outlying object?i
  } 



  void reset();


protected:

  rsBlepReadyOsc<T> osc1, osc2;
  TBlep blep1, blep2;  
  // maybe make embedded objects public

  bool sync12 = false;
  bool sync21 = false;

  //rsPolyBlep2<T, T> blep1, blep2;
  // todo: maybe make the class of the blep a template parameter such that we can select the
  // type of blep at instantiation time - done

};