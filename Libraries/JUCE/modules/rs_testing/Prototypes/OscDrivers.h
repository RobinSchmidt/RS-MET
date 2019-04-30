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



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline void getSamplePair(T* x1, T* x2)
  {
    // prelimiary (later switch waveforms)
    *x1 = osc1.getSampleSaw();
    *x2 = osc2.getSampleSaw();

    // todo: apply blep corrections...

  }


  inline T getSample()
  {
    T x1, x2;
    getSamplePair(&x1, &x2);
    return x1 + x2;  // maybe have mixing coefficients here? ...or maybe custom mixing be handled
                     // by outlying object?
  } 



  void reset();


protected:

  rsBlepReadyOsc<T> osc1, osc2;
  TBlep blep1, blep2;

  //rsPolyBlep2<T, T> blep1, blep2;
  // todo: maybe make the class of the blep a template parameter such that we can select the
  // type of blep at instantiation time - done

};