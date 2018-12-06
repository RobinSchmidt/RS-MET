#ifndef RAPT_SINUSOIDALMODEL_H
#define RAPT_SINUSOIDALMODEL_H

//=================================================================================================

template<class T>
class rsInstantaneousSineParams
{

public:


  /** Constructor. */
  //rsInstantaneousSineParams();


  /** \name Inquiry */

  inline T getTime() const { return time; }

  inline T getFrequency() const { return freq; }

  inline T getAmplitude() const { return gain; }

  inline T getWrappedPhase() const { return phase; }

  inline T getUnwrappedPhase() const { return phase + 2*PI*cycles; }


protected:

  T time = 0;      // in seconds
  T freq = 0;      // in Hz
  T gain = 0;      // as raw factor
  T phase = 0;     // radians in [0,2*pi) ...really - or is it [-pi, pi)? ...check atan2
  int cycles = 0;  // number of cycles passed: unwrapped phase = phase + 2*pi*cycles

};


//=================================================================================================

template<class T>
class rsSinusoidalPartial
{

public:





  /** \name Setup */

  inline void appendInstantaneousParameters(const rsInstantaneousSineParams<T>& params)
  {
    rsAssert(params.time() >= getEndTime()); // can only append to the end for insertion at random
                                             // position use insert...
    rsAppend(instParams, params);
  }

  inline void prependInstantaneousParameters(const rsInstantaneousSineParams<T>& params)
  {
    rsAssert(params.time() <= getStartTime()); // can only prepend at the start
    rsPrepend(instParams, params);
  }





  /** \name Inquiry */

  inline T getStartTime() const { return instParams[0].getTime(); }

  inline T getEndTime() const { return instParams[instParams.size()-1].getTime(); }


  //rsInstantaneousSineParams<T> getInstantaneousParameters() const;

  //std::vector<T> renderPartial(); // maybe it should take a (sub-sample) time-offset?


protected:

  std::vector<rsInstantaneousSineParams<T>> instParams;

};


//=================================================================================================

template<class T>
class rsSinusoidalModel
{

public:




  /** \name Setup */

  //inline void setBlockSize(int newBlockSize) const;

  //inline void setHopSize(int newHopSize) const;

  //setRootKey/setFundamentalFrequency, 

  /** Adds a new partial to the model. */
  void addPartial(const rsSinusoidalPartial<T>& partialToAdd) { partials.push_back(partialToAdd); }

  /** Removes the partial with given index from the model. */
  void removePartial(size_t index) { rsRemove(partials, index); }



  /** \name Processing */

  void analyze(T* sampleData, int numSamples);

  std::vector<T> synthesize() const;

protected:

  std::vector<rsSinusoidalPartial<T>> partials;


  rsPhaseVocoder<T> phsVoc; // or maybe that should be left to a rsSinusoidalModeler class that
  // spits out an rsSinusoidalModel? would have the advantage that we could make different 
  // implementations of the modeler that all spit out the same model data-structure

};



#endif