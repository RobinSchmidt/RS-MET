#pragma once


/** Baseclass for processing data of a sinusoidal model. Subclasses of that can be injected into
(re)synthesis pipeline to modify the data in order to achieve various effects. Two sorts of
processes can be injected at two points in the synthesis pipeline: first, you can modify the
model data itself and second, you can modify the (interpolated) trajectories for instantaneous
amplitude and phase. The first process operates at whatever (in general non-uniform) rate the
data of the model happens to be defined and the second process operates at the target sample-rate.

todo: maybe it's better to make two different baseclasses for these two types of processes - maybe
we should have a sort of facade class that encapsulates the analyzer(s), the synthesizer(s), a
chain of data processors and a chain of trajectory processors

The baseclass also defines a couple of static functions that fullfill basic processing tasks that
are required already in the analysis stage but don't really fit there very well because they may
be useful in other contexts as well. */

template<class T>
class rsSinusoidalProcessor
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Static utility functions

  /** Given an array of time-stamps and corresponding frequency and wrapped phase values, this
  function computes the corresponding array of unwrapped phase values by numerically integrating
  the frequency array and then re-adjusting the resulting (unwrapped) phases to have a value that
  is consistent with the wrappedPhase values (i.e. a suitable multiple of 2*pi shifted from the
  stored value). To integrate the frequency data, we also need the time axis because the data may
  be nonuniformly sampled. */
  static std::vector<T> unwrapPhase(const std::vector<T>& time,
    const std::vector<T>& freq, const std::vector<T>& wrappedPhase);



  /** Modifies the frequency values in the given partial, such that when they are numerically
  integrated (via a trapezoidal rule), the resulting phase values end up at values that are
  consistent with the stored phase values, i.e. differ from the stored values only by a multiple of
  2*pi. This is helpful to remove a bias in the estimated frequency values that may have occured
  during analysis, so it is a recommended post-processing step to refine the frequency estimates.
  Such a bias can result in temporary phase desynchronization issues of a resynthesized signal with
  respect to the original signal, leading to sinusoidal bursts in the resiudal - which are clearly
  undesirable. Such desync bursts would occur whenever the accumulated frequency bias crosses a
  multiple of pi (i think - verify). */
  static void makeFreqsConsistentWithPhases(rsSinusoidalPartial<T>& partial);
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

  /** Modifies the instantaneous phase values such that they are all equal to what the numrical
  integration (via trapezoidal rule) of the instantaneous frequency values would give. You must
  also pass an index at which the phase value should remain unchanged - all other phases will be
  adjusted with respect to the anchor point. It's typically a good idea to choose and index from
  the middle of the sample, where there's no transient stuff going on or the sample has faded out
  already. */
  static void makePhasesConsistentWithFreqs(rsSinusoidalPartial<T>& partial, int referenceIndex);



  static void refineFreqsViaPhaseDerivative(rsSinusoidalPartial<T>& partial);

  static void refineFreqsViaPhaseDerivative(rsSinusoidalModel<T>& model);





  /** Sets all instantaneous frequency values of the given partial to a fixed value and also
  updates the instantaneous phase values to be consistent with the new frequencies.
  todo: maybe make the phase re-adjustment optional...or maybe not - maybe it makes no sense to
  leave the phases as is anyway... */
  static void fixPartialFrequency(rsSinusoidalPartial<T>& partial, T freq, T pickPhaseAt = T(0.5));

  /** Sets the frequencies of all partials in the given model to exact integer multiples of the
  given fundamental frequency (at all time-instants) and adjusts the phase values to be consistent
  with the integrated frequency values. */
  static void makeStrictlyHarmonic(rsSinusoidalModel<T>& model, T fundamentalFreq,
    T inharmonicity, T pickPhaseAt = T(0.5));

  /** Extracts the partials whose mean frequency is less than or equal to the given splitFreq. */
  static rsSinusoidalModel<T> extractLowpassPart(rsSinusoidalModel<T>& model, T splitFreq);

  /** Extracts the partials whose mean frequency is greater than the given splitFreq. */
  static rsSinusoidalModel<T> extractHighpassPart(rsSinusoidalModel<T>& model, T splitFreq);




  //-----------------------------------------------------------------------------------------------
  // \name Processing functions to override


  /** This is the function that subclasses must override to implement their actual processing
  algorithm. The function gets a reference to the sinusoidal model object and is supposed to
  manipulate its data directly. */
  virtual void processModel(rsSinusoidalModel<T>& model) = 0;


  // virtual void processPartialEnvelopes
  // ..hmm...this should really be in another baseclass - maybe rsInstantaneousEnvelopeProcessor

};

//=================================================================================================

template<class T>
class rsPartialBeatingRemover : public rsSinusoidalProcessor<T>
{

public:

  /** Constuctor. Sets up the embedded envExtractor object to its default settings. */
  rsPartialBeatingRemover();


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up the parameters for the phase-smoothing filter. This is the filter that is used to
  smooth away the discontinuities in the (de-trended) phase data. These phase discontinuities are a
  feature of the beating - for an explanation, why these occur and why we need to get rid of them,
  see here: https://github.com/RobinSchmidt/RS-MET/issues/280 */
  virtual void setPhaseSmoothingParameters(T cutoff, int filterOrder, int numFilterPasses)
  {
    phaseSmootherCutoff    = cutoff;
    phaseSmootherOrder     = filterOrder;
    phaseSmootherNumPasses = numFilterPasses;
  }

  
  // obsolete
  //virtual void setMaxEnvelopeSampleSpacing(T newSpacing) 
  //{ 
  //  envExtractor.setMaxSampleSpacing(newSpacing);
  //}
  // this is a bad name - choose a better ...use minPrecision or just Precision...or Resolution
  // get rid of this parameter - automatically choose an appropriate value - given by the maximum
  // distance between actual peaks in the meta envelope

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Applies beating-removal to all partials. */
  virtual void processModel(rsSinusoidalModel<T>& model) override;

  /** Applies the beating removal to a selection of partials */
  virtual void processModel(rsSinusoidalModel<T>& model, const std::vector<int>& partialIndices);

  /** Applies beating-removal to the given partial. */
  virtual void removeBeating(rsSinusoidalPartial<T>& partial);

  /** De-beats the amplitude envelope of the given partial. */
  virtual void removeAmplitudeBeating(rsSinusoidalPartial<T>& partial);

  /** De-beats the phase trajectory of the givne partial. */
  virtual void removePhaseBeating(rsSinusoidalPartial<T>& partial);

  /** Smoothes the given phase-array... */
  static std::vector<T> smoothPhases(
    const std::vector<T>& time, const std::vector<T>& freqs, const std::vector<T>& phases,
    T cutoff, int order, int numPasses);
  // maybe factor this out - it could be useful for other processors, too ...maybe have a class
  // rsPhaseSmoother that lets the user select cutoff-freq, filter-type, order, numPasses, etc.
  // todo: make input vectors const (-> make GNUPlotter const-correct) - but this gives weird 
  // compiler errors
  // ...hmm - or maybe the filter-type should be fixed - to Bessel or Gaussian?



  /** Embedded envelope extrator object - made public, so client code can set up its options
  directly by accessing it via dot-syntax. */
  rsEnvelopeExtractor<T> envExtractor;

protected:

  // settings for the phase smoothing filter:
  T phaseSmootherCutoff = T(5);
  int phaseSmootherOrder = 1;
  int phaseSmootherNumPasses = 4;

};
