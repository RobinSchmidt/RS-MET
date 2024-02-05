#ifndef RS_REVERBSTUFF_H
#define RS_REVERBSTUFF_H


//=================================================================================================

/** An allpass delay that realizes the transfer function and difference equation:

          a +     z^(-M)
  H(z) = ----------------,    y[n] = a * x[n] + x[n-M] - a * y[n-M]
          1 + a * z^(-M)

so it's like a first order allpass filter with coefficient a in which the unit delay was replaced
by a delay line of length M. This is also known as a Schroeder allpass section. Such allpass delays 
can be used as building blocks for reverbs, for example. 


See:
https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html


...TBC...  */


template<class T>
class rsAllpassDelay
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsAllpassDelay() {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */


  void setMaximumDelayInSamples(int newMaxDelay);

  void setDelayInSamples(int newDelay);

  void setAllpassCoeff(T newCoeff) { allpassCoeff = newCoeff; }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline T getSample(T in);


  void reset();



protected:

  T allpassCoeff = 0.0;  // TPar

  RAPT::rsBasicDelayLine<T> inputDelayLine, outputDelayLine;  // TSig

};

// ToDo:
// -Optimize this to use only a single delayline by switching from a DF1-like implementation to a 
//  DF2-like implementation. But keep this implementation as prototype for unit tests of the 
//  optimized one.
// -use TSig and TPar


template<class T>
void rsAllpassDelay<T>::setMaximumDelayInSamples(int newMaxDelay)
{
  inputDelayLine.setMaximumDelayInSamples(newMaxDelay);
  outputDelayLine.setMaximumDelayInSamples(newMaxDelay);
}

template<class T>
void rsAllpassDelay<T>::setDelayInSamples(int newDelay)
{
  inputDelayLine.setDelayInSamples(newDelay);
  outputDelayLine.setDelayInSamples(newDelay);
}

template<class T>
T rsAllpassDelay<T>::getSample(T x)
{
  T xM = inputDelayLine.getSample(x);                            // x[n-M]
  T yM = outputDelayLine.getSampleSuppressTapIncrements(T(0));   // y[n-M]
  T y  = allpassCoeff * x + xM - allpassCoeff * yM;              // y[n], our current output
  outputDelayLine.addToInput(y);
  outputDelayLine.incrementTapPointers();
  return y;
  // ToDo: verify that this does the right thing with respect to the order of reading, writing and
  // incrementing the taps of the outputDelayLine. Maybe write a unit test that uses a delay of 
  // M = 1 and compare output to a regular first order allpass filter.
  //
  // We want to realize:
  //
  //          a +     z^(-M)
  //  H(z) = ----------------,    y[n] = a * x[n] + x[n-M] - a * y[n-M]
  //          1 + a * z^(-M)
}

template<class T>
void rsAllpassDelay<T>::reset()
{
  inputDelayLine.reset();
  outputDelayLine.reset();
}

// ToDo:
// -Build a nested allpass in which the z^(-M) term has been replaced by another allpass filter.
// -When building a series connection of those, as is done for a diffuser stage for a reverb, the 
//  output delay line of one stage becomes the input delayline for the next stage, so we can share
//  some delaylines compared toa naive implementation.



//=================================================================================================

/** Under construction....TBC...

This is supposed to be a naive prototype implementation. We just use a std::vector of 
rsAllpassDelay. This is suboptimal because in such a chain, the output delayline of stage i can at 
the same time serve as input delayline for stage i+1 such that we can get rid of almost half of the
delaylines in an optimized implementation. But we don't do that here. Well - when the 
rsAllpassDelay is implemented so as to use only one delayline (e.g. switch from DF1 to DF2 or TDF2, 
see https://www.dsprelated.com/freebooks/filters/Four_Direct_Forms.html), then this optimization 
here won't be needed anymore.


See:
https://ccrma.stanford.edu/~jos/pasp/Schroeder_Allpass_Sections.html
https://www.dsprelated.com/freebooks/pasp/Schroeder_Allpass_Sections.html

*/


template<class T>
class rsAllpassDelayChain  // maybe rename to rsAllpassDelayChainNaive
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxNumStages(int newMaxNumStages);

  void setNumStages(int newNumStages);

  void setMaxDelayInSamples(int stageIndex, int newMaxDelay);

  void setDelayInSamples(int stageIndex, int newDelay);

  void setAllpassCoeff(int stageIndex, T newCoeff);


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  int getMaxNumStages() const { return (int) allpassDelays.size(); }



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline T getSample(T in);

  void reset();



protected:

  std::vector<rsAllpassDelay<T>> allpassDelays;
  int numStages = 0;

};


template<class T>
void rsAllpassDelayChain<T>::setMaxNumStages(int newMaxNumStages)
{
  allpassDelays.resize(newMaxNumStages);
}

template<class T>
void rsAllpassDelayChain<T>::setNumStages(int newNumStages)
{
  RAPT::rsAssert(newNumStages <= getMaxNumStages());
  numStages = newNumStages;
}

template<class T>
void rsAllpassDelayChain<T>::setMaxDelayInSamples(int stageIndex, int newMaxDelay)
{
  RAPT::rsAssert(stageIndex < getMaxNumStages());
  allpassDelays[stageIndex].setMaximumDelayInSamples(newMaxDelay);
}

template<class T>
void rsAllpassDelayChain<T>::setDelayInSamples(int stageIndex, int newDelay)
{
  RAPT::rsAssert(stageIndex < getMaxNumStages());
  allpassDelays[stageIndex].setDelayInSamples(newDelay);
}

template<class T>
void rsAllpassDelayChain<T>::setAllpassCoeff(int stageIndex, T newCoeff)
{
  RAPT::rsAssert(stageIndex < getMaxNumStages());
  allpassDelays[stageIndex].setAllpassCoeff(newCoeff);
}


template<class T>
T rsAllpassDelayChain<T>::getSample(T in)
{
  T tmp = in;
  for(int i = 0; i < numStages; i++)
    tmp = allpassDelays[i].getSample(tmp);
  return tmp;
}

template<class T>
void rsAllpassDelayChain<T>::reset()
{
  for(int i = 0; i < getMaxNumStages(); i++)
    allpassDelays[i].reset();
}



#endif