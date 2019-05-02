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

  inline void getSamplePairNaive(T* x1, T* x2)
  {
    // preliminary (later switch waveforms)
    *x1 = osc1.getSampleSaw();
    *x2 = osc2.getSampleSaw(); // wait - we should call this *after* a potential reset - or not?

    // apply blep corrections to waveform discontinuities:
    if(osc1.getStepAmplitude() != 0.0) {
      //blep1.prepareForStep(osc1.getStepDelay(), osc1.getStepAmplitude());//commented for debug
      if(sync12) {
        T oldPhase = osc2.getPhase();

        // test:
        oldPhase += osc2.getPhaseIncrement();
        osc2.wrapPhase(oldPhase);

        T newPhase = osc1.getPhase() * osc2.getPhaseIncrement() / osc1.getPhaseIncrement();

        // todo: maybe have a continuous sync12 parameter between 0 and 1 and compute 
        // newPhase = (1-sync12)*oldPhase + sync12*newPhase

        //osc2.reset(newPhase);   // maybe we should add the increment?
        //osc2.reset(newPhase + osc2.getPhaseIncrement()); // like this?
        osc2.resetPhase(newPhase + osc2.getPhaseIncrement()); // like this?
        //osc2.reset(newPhase - osc2.getPhaseIncrement()); // like this?
        // because in reset, the increment is subtracted

          

        //T stepAmp = osc2.sawValue(oldPhase) - osc2.sawValue(newPhase);
        T stepAmp = osc2.sawValue(newPhase) - osc2.sawValue(oldPhase);
        //T stepAmp = osc2.sawValue(newPhase + osc2.getPhaseIncrement()) - osc2.sawValue(oldPhase);
        // or the other way around? also - later we need to switch between sawValue/squareValue
        // i think, new-old: when new = -1, old = +1, stepAmp = -2 - a downward step by 2

        T stepDly = newPhase / osc2.getPhaseIncrement(); 
        // is this correct? shouldn't this be osc1.phase / osc1.increment ...but that's probably
        // the same value
        //T stepDly2 = osc1.getPhase() / osc1.getPhaseIncrement(); // yep - same value

        blep2.prepareForStep(stepDly, stepAmp);
      }
    }




    if(osc2.getStepAmplitude() != 0.0)
    {
      blep2.prepareForStep(osc2.getStepDelay(), osc2.getStepAmplitude());
      if(sync21)
      {
        // ...
      }
    }
  }


  inline void applyBleps(T* x1, T* x2, T* y1, T* y2)
  {
    *y1 = blep1.getSample(*x1);
    *y2 = blep2.getSample(*x2);
  }

  inline void getSamplePair(T* x1, T* x2)
  {
    getSamplePairNaive(x1, x2);
    applyBleps(x1, x2, x1, x2);
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
  // maybe make embedded objects public - if we would not need to be able produce both osc-signals
  // separately, one single blep object would actually be sufficient - maybe make a variant of the
  // osc that does it like this - but then, getSamplePair would not be possible anymore...we'll see

  bool sync12 = false;
  bool sync21 = false;

  //rsPolyBlep2<T, T> blep1, blep2;
  // todo: maybe make the class of the blep a template parameter such that we can select the
  // type of blep at instantiation time - done

};