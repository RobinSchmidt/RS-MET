#pragma once


/** Baseclass for processing data of a sinusoidal model. Subclasses of that can be injected into
(re)synthesis pipeline to modify the data in order to achieve various effects. Two sorts of 
processes can be injected at two points in the synthesis pipeline: first, you can modify the 
model data itself and second, you can modify the (interpolated) trajectories for instantaneous
amplitude and phase. The first process operates at whatever (in general non-uniform) rate the
data of the model happens to be defined and the second process operates at the target sample-rate. 

The baseclass also defines a couple of static functions that fullfill basic processing tasks that 
are required already in the analysis stage but don't really fit there very well because they may
be useful in other contexts as well. */

template<class T>
class rsSinusoidalProcessor
{

public:


  /** Given an array of time-stamps and corresponding frequency and wrapped phase values, this 
  function computes the corresponding array of unwrapped phase values by numerically integrating
  the frequency array and then re-adjusting the resulting (unwrapped) phases to have a value that 
  is consistent with the wrappedPhase values (i.e. a suitable multiple of 2*pi shifted from the 
  stored value). To integrate the frequecy data, we also need the time axis because the data may
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


  static void refineFreqsViaPhaseDerivative(rsSinusoidalPartial<T>& partial);

  static void refineFreqsViaPhaseDerivative(rsSinusoidalModel<T>& model);



  // virtual void processModelData
  // virtual void processPartialEnvelopes



};