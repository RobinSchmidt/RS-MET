#ifndef RAPT_SINUSOIDALMODEL_H
#define RAPT_SINUSOIDALMODEL_H

//=================================================================================================

/** Data structure to hold and edit the instantaneous parameters of a sinusoidal partial at one 
instant of time. */

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
    //this->phase  = phase;
    this->phase  = rsWrapToInterval(phase, -PI, PI);
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
  T gain = 0;      // as raw factor - rename to amp
  T phase = 0;     // radians in [-pi, pi]
  //int cycles = 0;  // number of cycles passed: unwrapped phase = phase + 2*pi*cycles

};

//=================================================================================================

/** Data structure to hold and edit the information about one single sinusoidal partial. This is 
effectively an array of rsInstantaneousSineParams objects ordered by time-stamp.  */

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

  /** Sets the frequency for the datapoint at index i. */
  void setFrequency(int i, T newFreq) { instParams[i].freq = newFreq; }

  /** Sets the data for the instantaneous parameters for the datapoint with given index. */
  void setData(int i, T time, T freq, T gain, T phase)
  {
    rsAssert(i >= 0 && i < (int)instParams.size(), "Invalid index");
    instParams[i].time  = time;
    instParams[i].freq  = freq;
    instParams[i].gain  = gain;
    instParams[i].phase = phase;
  }

  /** Initializes our array of data-points. If you pass a nonzero number, then memory will be 
  allocated for that number of datapoints. You can then fill in the actual data via setData. */
  void init(int numDataPoints = 0)
  {
    instParams.clear();
    if(numDataPoints > 0)
      instParams.resize(numDataPoints);
  }

  /** Modifies the frequency values in this partial, such that when they are numerically 
  integrated (via a trapezoidal rule), the resulting phase values end up at values that are 
  consistent with the stored phase values, i.e. differ from the stored values only by a multiple of
  2*pi. This is helpful to remove a bias in the estimated frequency values that may have occured 
  during analysis, so it is a recommended post-processing step to refine the frequency estimates. 
  Such a bias can result in temporary phase desynchronization issues of a resynthesized signal with
  respect to the original signal, leading to sinusoidal bursts in the resiudal - which are clearly 
  undesirable. Such desync bursts would occur whenever the accumulated frequency bias crosses a 
  multiple of pi (i think - verify). */
  //void makeFreqsConsistentWithPhases();
  // or maybe call it deBiasFreqEstimates ...or maybe the name should somehow reflect, that the 
  // frequencies are modified, because we could also make it consistent by adjusting the phases - 
  // which is no good idea (it is, in fact, actually exactly what the synthesizer does in case of
  // inconsistency), because phase errors don't accumulate, so we can assume that they are more or 
  // less correct at each datapoint (which is not true for accumulated freq, if there's a bias in
  // the freq estimate) ..makeFreqsConsistentWithPhase
  // maybe return the maximum phase difference that occured between the integrated freq and adjusted 
  // (by k*2*pi) stored phase from one point to the next - this value should be much less than pi. 
  // if it gets close to pi, it may mean that the hop-size was too small, so we may use this as 
  // feedback for the user, if the hopsize parameter was good enough, so they may decide to analyze 
  // again with smaller hopSize
  // maybe this function should be part of the rsSinusoidalPartial class - it may be useful for
  // transformation/effect algorithms that mangle the datapoints, too
  // after a couple of tests, it actually seeems to disimprove the freq estimates - sometimes they 
  // tends to alternate between two wrong values below and above the correct value...maybe the 
  // end-condition / additional equation is a bad choice? or maybe the whole thing is a bad idea
  // anyway? -> more experiments needed


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

  /** Computes and returns the mean frequency of this partial. */
  T getMeanFreq() const;

  /** Returns the minimum frequency of this partial. */
  T getMinFreq() const;

  /** Returns the minimum frequency of this partial. */
  T getMaxFreq() const;

  /** Returns the maximum absolute difference between phases that would be computed by numerically
  integrating frequencies and actually stored phase values (always taking appropriate wrapping into
  account). */
  T getMaxFreqPhaseInconsistency() const;

  /** Returns true, when this partial will alias at a given sample-rate. If "allTheTime" is true, 
  it will return true only if it will alias all the time, i.e. if its minimum frequency is higher
  than sampleRate/2. If "allTheTime" is false (the default), it will return true, even if the 
  partial will alias only temporarily, i.e. if its maximum frequency is higher than 
  sampleRate/2. */
  bool willAlias(T sampleRate, bool allTheTime = false) const;

  /** Returns true, if the sample-instants (i.e. time-stamps of datapoints) of this partial are all
  the same as the sample-instants of the given other partial. Currently, this checks for exact 
  equality of all time-stamps - maybe allow for a tolerance later (optional, zero by default).  */
  bool isSampledInSyncWith(const rsSinusoidalPartial<T>& otherPartial) const;

  /** Returns the number of data points in this partial */
  inline size_t getNumDataPoints() const { return instParams.size(); }


  /** Returns a copy of the data point at given index. */
  inline rsInstantaneousSineParams<T> getDataPoint(size_t index) const 
  { 
    return instParams[index]; 
  }

  /** Returns a reference to the data point at given index. */
  inline rsInstantaneousSineParams<T>& getDataPointRef(size_t index)
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

  /** Fills the 4 given vectors with the data from our data points - useful for plotting and 
  creating interpolated data arrays. */
  void getDataArrays(std::vector<T>& time, std::vector<T>& freq, std::vector<T>& gain, 
    std::vector<T>& phase) const;


  std::vector<T> getTimeArray() const;

  std::vector<T> getFrequencyArray() const;

  std::vector<T> getAmplitudeArray() const;

  std::vector<T> getPhaseArray() const;





  /** Returns reference to given datapoint for manipulation (use with care) */
  rsInstantaneousSineParams<T>& getDataRef(int dataIndex)
  { 
    return instParams[dataIndex];
    // todo: range-check index value
  }

  /** Returns a constant reference to given datapoint for inquiry. */
  const rsInstantaneousSineParams<T>& getConstDataRef(int dataIndex) const
  { 
    return instParams[dataIndex];
  }



  //rsInstantaneousSineParams<T> getInstantaneousParameters() const;

  //std::vector<T> renderPartial(); // maybe it should take a (sub-sample) time-offset?


protected:

  std::vector<rsInstantaneousSineParams<T>> instParams;

  //friend class rsSinusoidalModel<T>; // doesn't compile - why?

};

//=================================================================================================

/** Data structure to hold the information of a sinusoidal model for a complete sound. This is 
effectively an array of rsSinusoidalPartial objects. 

todo: maybe allow the time to be represented in samples or seconds - samples may have numerical
advantages in computations, when the time-stamps are at integers or half-integers - maybe have a
sampleRate member - if it is zero (or nan or whatever special value), it means that time-stamps 
are in seconds, otherwise in samples at the given sample rate

todo: currently, the partials are in no particular order - maybe order them somehow (perhaps by 
start-time and then by frequency - or maybe by amplitude or even a psychoacoustic "importance"
measure?) . yes - it would be interesting to have a getPartialImportance() function that analyzes
according to amplitude, masking, etc. - can be used to "simplify" models - but maybe that should be
done in a dedicated class. we could also "sparsify" models by downsampling the partials datapoint 
density - we could remove datapoints whose information is well reconstructed by the interpolation 
process - we could measure the importance of a datapoint by reconstructing with and without it and
measure the amplitude and phase difference at the instant of the datapoint in the interpolated
trajectories - but that will depend on the selected interpolation - maybe a model should also 
contain data about its prefered interpolation mode


*/

template<class T>
class rsSinusoidalModel
{

public:

  /** \name Setup */

  /** Adds a new partial to the model. */
  void addPartial(const rsSinusoidalPartial<T>& partialToAdd) { partials.push_back(partialToAdd); }

  /** Appends multiple partials to the model at once. */
  void addPartials(const std::vector<rsSinusoidalPartial<T>>& partialsToAdd)
  {
    for(size_t i = 0; i < partialsToAdd.size(); i++)
      addPartial(partialsToAdd[i]);
    // we need to make a deep copy
  }

  /** Removes the partial with given index from the model. */
  void removePartial(size_t index) { rsRemove(partials, index); }

  /** Removes all partials except those whose indices are given in the passed vector. */ 
  void keepOnly(std::vector<size_t> partialsToKeep) 
  { 
    partials = rsSelect(partials, partialsToKeep); 
  }




  /** Removes the partials above given index from the model. */
  void removePartialsAbove(size_t maxIndexToRetain) 
  { 
    rsRemoveRange(partials, maxIndexToRetain+1, partials.size()-1); 
  }


  /** Sets the data for the given partial- and data-point index. */
  void setData(int partialIndex, int dataIndex, T time, T freq, T gain, T phase)
  {
    rsAssert(partialIndex >= 0 && partialIndex < (int)partials.size(), "Invalid index");
    partials[partialIndex].setData(dataIndex, time, freq, gain, phase);
  }

  /** Initializes the model. Creates the given number of partials where each partial has the given 
  number of (initially empty, i.e. zero-valued) datapoints. This is meant to be used to 
  pre-allocate memory for the model data to be subsequently filled in by calls to setData. */
  void init(int numPartials = 0, int numFrames = 0);

  // maybe have a replacePartial function

  /** Calls the function of the same name for each partial. See comment in rsSinusoidalPartial for
  details what it does. */
  void makeFreqsConsistentWithPhases();
  // move to rsSinusoidalProcessor


  /** \name Inquiry */

  /** Returns the start time of the whole sound, i.e. the start time of the partial that starts 
  first. */
  T getStartTime() const;

  /** Returns the end time of the whole sound, i.e. the end time of the partial that ends last. */
  T getEndTime() const;

  /** Returns the start-sample index at a given sample-rate as determined by the start time. */
  int getStartSampleIndex(T sampleRate) const { return (int)floor(getStartTime() * sampleRate); }

  /** Returns the end-sample index at a given sample-rate as determined by the end time. */
  int getEndSampleIndex(T sampleRate) const { return (int)ceil(getEndTime() * sampleRate); }

  /** Returns the total length in samples for the sound at a given sample rate. */
  int getLengthInSamples(T sampleRate) const 
  { 
    return  getEndSampleIndex(sampleRate) - getStartSampleIndex(sampleRate) + 1;
  }

  /** Returns the number of partials in this model. */
  size_t getNumPartials() const { return partials.size(); }

  /** Returns a constant reference to the partial with given index. Use this when you need inquire
  information about this partial and its datapoints but don't intend to modify the data .*/
  const rsSinusoidalPartial<T>& getPartial(size_t index) const { return partials[index]; }
  // maybe rename to getConstPartialRef

  /** Returns a non-constant, i.e. modifiable, reference to the partial with given index. Use this, 
  when you intend to manipulate the data of the partial. */
  rsSinusoidalPartial<T>& getModifiablePartialRef(size_t index) { return partials[index]; }
  // maybe rename to getPartialRef

  /** Returns a non-constant, i.e. modifiable, reference to a particular datapoint inside a 
  particular the partial with given index. Use this, when you intend to manipulate the datapoint of
  the partial. */
  rsInstantaneousSineParams<T>& getDataRef(int partialIndex, int dataIndex)
  { 
    return partials[partialIndex].getDataRef(dataIndex);
    // todo: range-check index values
  }

  /** Returns true, if all partials are sampled at the same time-instants. */
  bool isSampledSynchronously();


protected:

  std::vector<rsSinusoidalPartial<T>> partials;


  //friend class ::SinusoidalAnalyzer; // needs dirct access to the partials array
  // todo: remove the namespace-level up ::, hen the analyzer is moved to rapt

};

//=================================================================================================

// Analyzer:





//=================================================================================================

// Synthesizer:





//=================================================================================================

// Transformations:

// ideas: removeNoisyPartials, denoisePartials, harmonifyPartials, phaseLockPartials, ...
// -maybe it would make sense to consider the de-trended, unwrapped phase signal and apply signal
//  processing algorithsm to that
// -in general, this seems to be a good oppottunity to apply non-uniform filters (although the 
//  analyzer will produce uniformly sampled data, a general model may be non-uniformly sampled)
//  -maybe to test such algorithms, take the analysis data and "poke holes" in it (i.e. delete some 
//   datapoints in the middle)
//  -maybe a sinusoid with vibrato and tremolo is a good prototypical test case, maybe the waveform
//   of the vibrato could be initially squarish and we apply a lowpass to get a more sinusoidal 
//   vibrato (or maybe use TriSaw vibrato)






#endif