#ifndef RAPT_SCOPESCREENSCANNER_H_INCLUDED
#define RAPT_SCOPESCREENSCANNER_H_INCLUDED

/** This is a class for generating the sawtooth-shaped waveform used for scanning over the screen
horizontally in an oscilloscope in 1D mode. It provides synchronization with the incoming waveform 
by using a zero crossing detector. 

\todo:
-add setNumCyclesShown, setPhaseOffset, setFreqMultiplier (Elan requested that - what should it 
 do?)
-how to implement the PhaseOffset - maybe with a delayline?
*/

template<class T>
class rsScopeScreenScanner
{

public:

  //---------------------------------------------------------------------------------------------
  // \name Construction/Destruction

  /** Constructor. */
  rsScopeScreenScanner();

  //---------------------------------------------------------------------------------------------
  // \name Setup:

  /** Sets the sample-rate. */
  void setSampleRate(T newSampleRate);

  /** Sets the frequency that should be used when we are not in sync mode. */
  void setScanFreqNoSync(T newFrequency);

  /** Switches synchronization of the sawtooth to the input signal on/off */
  void setSync(bool shouldSync);

  /** Sets the number of cycles shown in sync mode. */
  void setNumCyclesShown(int numCycles);

  ///** Sets a zoom factor. */
  //void setZoom(double newZoom);


  //---------------------------------------------------------------------------------------------
  // \name Audio Processing

  /** Generates one sawtooth output sample at the time. You must pass the input signal value that
  is used for the pitch analysis. The output value is between 0 and 1. */
  inline T getSample(T in);

  //---------------------------------------------------------------------------------------------
  // \name Misc:

  /** Resets the internal state. */
  void reset();


protected:

  T sampleRate, scanFreq, zoom;
  T sawPhase, sawInc;  // sawtooth phase and increment
  bool sync;

  T xOld;

  int samplesSinceReset, samplesSinceLastZero; 
  int zeroCrossingCount, minZeroDistance, numZerosToReset;

  rsLadderFilter<T, T> lowpass;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class T>
inline T rsScopeScreenScanner<T>::getSample(T in)
{
  T result = sawPhase;
  sawPhase += sawInc;

  if(!sync)
  {
    if(sawPhase > 1)
      sawPhase -= 1;
    return result;
  }

  // If we are here, we are in sync mode:
  in = lowpass.getSample(in); 
  samplesSinceReset++;
  samplesSinceLastZero++;
  if(samplesSinceLastZero >= minZeroDistance)
  {
    if(xOld < 0 && in >= 0) // new zero detected
    {
      zeroCrossingCount++;
      samplesSinceLastZero = 0;
      if(zeroCrossingCount >= numZerosToReset)
        reset();
    }
  }
  xOld = in;
  return result;
}

#endif
