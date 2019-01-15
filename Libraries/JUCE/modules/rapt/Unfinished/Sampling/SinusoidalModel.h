#ifndef RAPT_SINUSOIDALMODEL_H
#define RAPT_SINUSOIDALMODEL_H

//=================================================================================================

template<class T>
class rsInstantaneousSineParams
{

public:


  /** Constructor. */
  rsInstantaneousSineParams(T time = 0, T frequency = 0, T amplitude = 0, T phase = 0)
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


//protected:

  T time = 0;      // in seconds
  T freq = 0;      // in Hz
  T gain = 0;      // as raw factor
  T phase = 0;     // radians in [-pi, pi]
  //int cycles = 0;  // number of cycles passed: unwrapped phase = phase + 2*pi*cycles

};

//=================================================================================================

template<class T>
class rsSinusoidalPartial  // maybe rename to rsSinusoidalTrack
{

public:


  /** \name Setup */

  /** Appends the given datapoint at the end. */
  inline void appendDataPoint(const rsInstantaneousSineParams<T>& params)
  {
    rsAssert(params.getTime() >= getEndTime()); // can only append to the end for insertion at random
                                             // position use insert...
    rsAppend(instParams, params);
  }

  /** Prepends given datapoint at the front. */
  inline void prependDataPoint(const rsInstantaneousSineParams<T>& params)
  {
    rsAssert(params.getTime() <= getStartTime()); // can only prepend at the start
    rsPrepend(instParams, params);
  }

  /** Prepends a zero amplitude datapoint at the front with frequency equal to the start 
  frequency, a time-stamp given by the start time minus the fadeTime and appropriately chosen 
  phase. */
  void applyFadeIn(T fadeTime);


  /** Appends a zero amplitude datapoint at the front with frequency equal to the end
  frequency, a time-stamp given by the end time plus the fadeTime and appropriately chosen 
  phase. */
  void applyFadeOut(T fadeTime);

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

  /** Returns the total length in seconds. */
  inline T getLength() const
  {
    return getEndTime() - getStartTime();
  }



  inline T getStartFreq() const 
  { 
    if(instParams.size() == 0) 
      return T(0);
    return instParams[0].getFrequency(); 
  }

  inline T getEndFreq() const 
  { 
    if(instParams.size() == 0) 
      return T(0);
    return instParams[instParams.size()-1].getFrequency(); 
  }


  inline size_t getNumDataPoints() const { return instParams.size(); }


  /** Returns the data point at given index. */
  inline rsInstantaneousSineParams<T> getDataPoint(size_t index) const 
  { 
    return instParams[index]; 
  }

  /** Returns the first data point. */
  inline rsInstantaneousSineParams<T> getFirstDataPoint() const 
  { 
    return instParams[0]; 
  }

  /** Returns the last data point. */
  inline rsInstantaneousSineParams<T> getLastDataPoint() const 
  { 
    return instParams[instParams.size()-1]; 
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

  /** Appends multiple partials to the model at once. */
  void addPartials(const std::vector<rsSinusoidalPartial<T>>& partialsToAdd)
  {
    for(size_t i = 0; i < partialsToAdd.size(); i++)
      addPartial(partialsToAdd[i]);
    // we need to make a deep copy
  }


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
  rsSinusoidalModel<T> analyze(T* sampleData, int numSamples, T sampleRate) const;

protected:

  rsSpectrogram<T> phsVoc; 

};


#endif