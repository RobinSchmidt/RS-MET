#ifndef RAPT_SINUSOIDALMODEL_H
#define RAPT_SINUSOIDALMODEL_H

//=================================================================================================

template<class T>
class rsInstantaneousSineParams
{

public:


  /** Constructor. */
  rsInstantaneousSineParams(T time, T frequency, T amplitude, T phase/*, int numCycles*/)
  {
    this->time   = time;
    this->freq   = frequency;
    this->gain   = amplitude;
    this->phase  = phase;
    //this->cycles = numCycles;
  }


  /** \name Inquiry */

  inline T getTime() const { return time; }

  inline T getFrequency() const { return freq; }

  inline T getAmplitude() const { return gain; }

  inline T getWrappedPhase() const { return phase; }

  //inline T getUnwrappedPhase() const { return phase + 2*PI*cycles; }


protected:

  T time = 0;      // in seconds
  T freq = 0;      // in Hz
  T gain = 0;      // as raw factor
  T phase = 0;     // radians in [0,2*pi) ...really - or is it [-pi, pi)? ...check atan2
  //int cycles = 0;  // number of cycles passed: unwrapped phase = phase + 2*pi*cycles

};

//=================================================================================================

template<class T>
class rsSinusoidalPartial
{

public:


  /** \name Setup */

  inline void appendDataPoint(const rsInstantaneousSineParams<T>& params)
  {
    rsAssert(params.getTime() >= getEndTime()); // can only append to the end for insertion at random
                                             // position use insert...
    rsAppend(instParams, params);
  }

  inline void prependDataPoint(const rsInstantaneousSineParams<T>& params)
  {
    rsAssert(params.getTime() <= getStartTime()); // can only prepend at the start
    rsPrepend(instParams, params);
  }

  // maybe have insertDataPoint function


  /** \name Inquiry */

  inline T getStartTime() const 
  { 
    if(instParams.size() == 0) 
      return T(0);
    return instParams[0].getTime(); 
  }

  inline T getEndTime() const 
  { 
    if(instParams.size() == 0) 
      return T(0);
    return instParams[instParams.size()-1].getTime(); 
  }


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

  /** Adds a new partial to the model. */
  void addPartial(const rsSinusoidalPartial<T>& partialToAdd) { partials.push_back(partialToAdd); }

  /** Removes the partial with given index from the model. */
  void removePartial(size_t index) { rsRemove(partials, index); }


  /** \name Inquiry */

  /** Returns the end time of the whole sound, i.e. the end time of the partial that ends last. */
  T getEndTime() const;

  /** Returns the number of partials in this model. */
  size_t getNumPartials() const { return partials.size(); }

  /** Returns a reference to the partial with given index. */
  const rsSinusoidalPartial<T>& getPartial(size_t index) const { return partials[index]; }

protected:

  std::vector<rsSinusoidalPartial<T>> partials;

};

//=================================================================================================

template<class T>
class rsSinusoidalSynthesizer
{

public:

  /** Given a sinusoidal model data-structure, this function synthesizes the audio signal from that 
  model. */
  std::vector<T> synthesize(const rsSinusoidalModel<T>& model, T sampleRate) const;

};

//=================================================================================================

template<class T>
class rsSinusoidalAnalyzer
{

public:


  /** \name Setup */

  //inline void setBlockSize(int newBlockSize) const;

  //inline void setHopSize(int newHopSize) const;

  //setRootKey/setFundamentalFrequency, 

  /** \name Processing */

  /** Takes in an array of audio samples and returns the sinusoidal model that approximates the
  sample data. */
  rsSinusoidalModel<T> analyze(T* sampleData, int numSamples, T sampleRate);

protected:

  rsPhaseVocoder<T> phsVoc; 

};


#endif