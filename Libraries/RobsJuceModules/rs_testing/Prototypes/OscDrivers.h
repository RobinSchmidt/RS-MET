#pragma once


// template parameters: 
// T: type for signals and parameters, TOsc: class for the (blep-ready) osc, TBlep: class for blep
template<class T, class TOsc, class TBlep> 

/** A class for producing waveforms like supersaw, supersquare, etc. using a blep-ready oscillator 
class as basis and applies a blep for anti-aliasing. */ 

class rsSuperBlepOsc
{

public:


  //rsSuperBlepOsc();

  /** Sets the reference phase increment which determines the center frequency of the osc stack 
  around which all other frequencies are arranged. */
  inline void setReferenceIncrement(T newIncrement) 
  { 
    inc = newIncrement; 
    updateIncrements();
  }

  /** Sets the number of oscillators to use. */
  inline void setNumOscillators(int newNumber)
  {
    numOscs = rsMin(newNumber, getMaxNumOscillators());
  }

  /** Sets the maximum number of oscillators that can be used. May cause memory (re)allocation and
  should probably be called once at startup time. */
  inline void setMaxNumOscillators(int newMaximum)
  {
    incs.resize(newMaximum);
    oscs.resize(newMaximum);
    numOscs = rsMin(numOscs, newMaximum);
    updateIncrements();
  }

  /** Returns the maximum number of oscillators that can be used. */
  inline int getMaxNumOscillators() const { return (int) oscs.size(); } 


protected:

  /** Updates our array of phase increments according to the desired reference increment and 
  settings of detune, etc.. */
  void updateIncrements();

  int numOscs = 1;        // current number of oscillators
  T inc = 0;              // reference increment
  std::vector<TOsc> oscs; // oscillator array
  std::vector<T> incs;    // array of phase increments
  TBlep blep; // a single blep is shared among all oscs

  // todo: have an array of pan positions - if we do stereo later, we will need two blep objects 
  // - one for each channel
};





//=================================================================================================

/** A class that produces a phasor (i.e. a sawtooth wave from 0 to 1) running at some "slave" 
frequency but also syncing to some (typically lower) "master" frequency. For preliminary 
investigations for oscillator sync (it's simpler to consider the phasor's first). */


template<class T, class TBlep> // T: type for signal and parameter, TBlep: class for BLEP object
class rsSyncPhasor
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  inline void setMasterIncrement(T newIncrement) { masterInc = newIncrement; }

  inline void setSlaveIncrement(T newIncrement) { slaveInc  = newIncrement; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline T getMasterIncrement() const { return masterInc; }

  inline T getSlaveIncrement()  const { return slaveInc;  }

  inline T getMasterPosition()  const { return masterPos; }

  inline T getSlavePosition()   const { return slavePos;  }


  //-----------------------------------------------------------------------------------------------
  /** \name Embedded object access */

  /** Returns a reference to the embedded blep object */
  inline const TBlep& getBlep() const { return blep; };


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline T getSample()
  {
    return applyBlep(getSampleNaive());
  }

  inline T getSampleNaive()
  {
    // increment phase variables (we use increment-before-output):
    masterPos += masterInc;
    slavePos  += slaveInc;
    handleWrapArounds();
    return slavePos;
  }

  inline T applyBlep(T x)
  {
    return blep.getSample(x);
  }



  inline void handleMasterWrapAround(T newMasterPos)
  {
    T masterStepDelay = newMasterPos / masterInc;
    T oldPos  = slavePos;
    T newPos  = masterStepDelay * slaveInc;
    slavePos  = newPos;
    //T stepAmp = slavePos - oldPos;   // or 0 - oldPhase ?
    T stepAmp = T(0) - oldPos;   // ..like this?
    T stepDly = newPos / slaveInc;
    blep.prepareForStep(stepDly, stepAmp);
    masterPos = newMasterPos;
  }
  // rename to handleMasterReset

  inline void handleSlaveWrapAround(T newSlavePos)
  {
    T slaveStepDelay = newSlavePos  / slaveInc;

    T slaveStepAmp   = T(-1);
    // i think, this is wrong in general - if a master wraparound occurs before the slave 
    // wraparound, the amplitude is much less - i think, the case can be identified by
    // if(slavePos < 1) or if(slavePos <= 1)

    if(slavePos < T(1))
    {
      //slaveStepAmp = T(0);  // preliminary
      //slaveStepAmp = T(0) - slavePos;
      slaveStepAmp = newSlavePos - slavePos;
    }

    blep.prepareForStep(slaveStepDelay, slaveStepAmp);
    slavePos = newSlavePos;
  }

  inline bool isMasterResetFirst(T wrappedMasterPos, T wrappedSlavePos)
  {
    return wrappedMasterPos * slaveInc > wrappedSlavePos * masterInc;
    // Conceptually, we do this:
    // T masterStepDelay = wrappedMasterPos / masterInc;
    // T slaveStepDelay  = wrappedSlavePos  / slaveInc;
    // return masterStepDelay > slaveStepDelay;
    // ...but without divisions
  }

  /** Figures out if the master and/or slave phasor needs a reset/wraparound and if so, handles it
  and prepares the blep object accordingly. There are 4 cases: master-wrap-only, slave-wrap-only, 
  master-then-slave-wrap, slave-then-master-wrap. */
  inline void handleWrapArounds()  // rename to handleResets
  {
    T wrappedMasterPos = T(-1);
    T wrappedSlavePos  = T(-1);
    if(masterPos >= T(1)) wrappedMasterPos = masterPos - T(1);
    if(slavePos  >= T(1)) wrappedSlavePos  = slavePos  - T(1);
    if(wrappedMasterPos >= T(0) || wrappedSlavePos >= T(0)) {
      if(wrappedSlavePos == T(-1))                         // master reset only
        handleMasterWrapAround(wrappedMasterPos);          
      else if(wrappedMasterPos == T(-1))                   // slave reset only
        handleSlaveWrapAround(wrappedSlavePos);            
      else {                                               // master and slave reset
        if( isMasterResetFirst(wrappedMasterPos, wrappedSlavePos) ) {
          handleMasterWrapAround(wrappedMasterPos);        // master reset first...
          handleSlaveWrapAround(wrappedSlavePos);          // ...slave reset second
        } else {
          handleSlaveWrapAround(wrappedSlavePos);          // slave reset first...
          handleMasterWrapAround(wrappedMasterPos);        // ...master reset second
        }
      }
    }
  }
  // can the logic be simplified? what happens if we do the wraps in the wrong order? does that
  // actually matter? if not, the logic can indeed be simplified a lot
  
  // Here is a great discussion of how to process the reset logic:
  // https://www.kvraudio.com/forum/viewtopic.php?p=5944009
  // mystran suggests partial sample steps - which seems to a good alternative to the logic above


  inline void resetPhase()
  {
    masterPos = T(0) - masterInc; wrapPhase(masterPos);
    slavePos  = T(0) - slaveInc;  wrapPhase(slavePos);
    // subtract increments because in getSample, we increment before producing output
  }

  inline void reset()
  {
    resetPhase();
    blep.reset();
  }

  static inline void wrapPhase(T& phase)
  {
    while( phase < T(0) )
      phase += T(1);
    while( phase >= T(1) )
      phase -= T(1);
  }



protected:

  TBlep blep;
  // maybe a pointer can be used - allows to share the blep object among several oscs - they can 
  // then all create naive signals and after each naive signal generation the blep is prepared
  // but not yet applied - it is only applied to the mix of all oscs - optimizes memory and cpu use
  // maybe make protected and use an access function that returns a reference ...like 
  // getBlepReference/getBlepRef or just getBlep - maybe this idiom of function returning 
  // references to embedded objects should be used consistently throughout the library


  T masterInc = T(0);
  T masterPos = T(0);

  T slaveInc = T(0);
  T slavePos = T(0);

};






//=================================================================================================

/** Implements sawtooth/pulse oscillator which can be additionally "synced" to another 
"master" frequency. The step discontinuities of the waveform itself and also those due to sync are
anti-aliased via a BLEP. */

template<class T, class TBlep> // T: type for signal and parameter, TBlep: class for BLEP object
class rsSyncOsc
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */


  void setMasterIncrement(T newIncrement) { master.setPhaseIncrement(newIncrement); }

  void setSlaveIncrement(T newIncrement) { slave.setPhaseIncrement(newIncrement); }

  // void setSyncAmount


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */




  inline T getSample()
  {
    return applyBlep(getSampleNaive());
  }


  inline T getSampleNaive()
  {
    T stepDelay, stepAmp;

    T dummy = master.getSampleSaw(&stepDelay, &stepAmp); // later use just a phasor - output not actually used


    if(stepAmp != 0.0)
    {
      T oldPhase = slave.getPhase();
      T newPhase = master.getPhase() * slave.getPhaseIncrement() / master.getPhaseIncrement();


      slave.resetPhase(newPhase + slave.getPhaseIncrement()); // verify the +inc
      //slave.resetPhase(newPhase);
      // ..i think, we may need to get rid of the +inc but also take it into when computing the
      // step amplitude ..ther we may need the +inc or something...it's tricky because of
      // pre-increment in slave.getSampleSaw

      T stepAmp = slave.sawValue(newPhase) - slave.sawValue(oldPhase);
      T stepDly = newPhase / slave.getPhaseIncrement(); 

      // i think, maybe the step-amp should alway be -1 - sawValue(oldPhase) because it always 
      // jumps down to -1 *at* the continuous step instant - we don't want to read out the osc at
      // the position where it advances to within the sample
      //T stepAmp = T(-1) - slave.sawValue(oldPhase);

      blep.prepareForStep(stepDly, stepAmp);
    }

    T out = slave.getSampleSaw(&stepDelay, &stepAmp);
    if(stepAmp != 0.0)
      blep.prepareForStep(stepDelay, stepAmp);
    return out;
  }

  inline T applyBlep(T x)
  {
    return blep.getSample(x);
  }


  void reset();



  rsBlepReadyOsc<T> master, slave; // use simple phasor for slave later
  TBlep blep;


protected:



};




//=================================================================================================


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

    T stepDelay1, stepAmp1, stepDelay2, stepAmp2;

    // preliminary (later switch waveforms)
    *x1 = osc1.getSampleSaw(&stepDelay1, &stepAmp1);
    *x2 = osc2.getSampleSaw(&stepDelay2, &stepAmp2); // wait - we should call this *after* a potential reset - or not?

    // apply blep corrections to waveform discontinuities:
    if(stepAmp1 != 0.0) {
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




    if(stepAmp2 != 0.0)
    {
      blep2.prepareForStep(stepDelay2, stepAmp2);
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