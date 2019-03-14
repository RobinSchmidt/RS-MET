template<class T>
void rsSinusoidalPartial<T>::applyFadeIn(T fadeTime)
{
  rsInstantaneousSineParams<T> params = getFirstDataPoint();
  params.time  -= fadeTime;
  params.gain   = 0.0;
  params.phase -= 2*PI*fadeTime*params.freq;               // is this correct? 
  params.phase  = rsWrapToInterval(params.phase, -PI, PI); // or should it be 0..2*PI? check synthesis
  prependDataPoint(params);
}

template<class T>
void rsSinusoidalPartial<T>::applyFadeOut(T fadeTime)
{
  rsInstantaneousSineParams<T> params = getLastDataPoint();
  params.time  += fadeTime;
  params.gain   = 0.0;
  params.phase += 2*PI*fadeTime*params.freq;
  params.phase  = rsWrapToInterval(params.phase, -PI, PI);
  appendDataPoint(params);
}
// idea: maybe the fade-in/out should not be a fixed time but somehow extend over a time between 
// one half cycle and one full cycle - or some number of cycles - so the sinusoid starts at a time
// instant that coincides with a zero crossing. that alone would avoid discontinuties, even if 
// there's no amplitude fade. (i think) with an additional linear fade-in, it would avoid 
// discontinuities in the derivative, too
// the fade-in/out, if too short, can also lead to overshooting artifacts in the amplitude 
// envelope, when cubic interpolation is used
// maybe, instead of doing fade-in/out by appending additional data points, we should just 
// extrapolate the envelope beyond the end points and apply an *additional* envelope
// ->experimentation is required




template<class T>
void rsMinSqrDiffWithGivnSum(T* x, T* s, int N)
{
  // x: values to be computed, s: desired sum-values (constraints)


  // todo: this is very preliminary, just for proof of concept and uses a standard Gaussian 
  // elimination algorithm to solve the linear system, which has a complexity of O(N^3). But the 
  // system is pentadigonal, so an O(N) algorithm is possible...so we need an algorithm for 
  // pentadiagonal systems in the library - or even better: an algorithm for general band-diagonal
  // matrices

  // this thesis: https://web.stanford.edu/group/SOL/dissertations/bradley-thesis.pdf
  // says that for symmetric, positive definite matrices, scaling to unit diagonal is effective for
  // making the problem better conditioned


  // maybe this is suitable:
  // https://www.boost.org/doc/libs/1_69_0/libs/numeric/ublas/doc/index.html

  // https://www.boost.org/doc/libs/1_69_0/libs/numeric/ublas/doc/banded.html
}

template<class T>
void rsSinusoidalPartial<T>::makeFreqsConsistentWithPhases()
{
  std::vector<T> t = getTimeArray();
  std::vector<T> f = getFrequencyArray();
  std::vector<T> p = getPhaseArray();
  int M = (int) t.size(), m;
  RAPT::rsAssert(M >= 2);  // valid partials should have at least two datapoints
  // for optimization, we could do away with obtaining these arrays, working on them and then 
  // writing the frequency back. Instead, we could operate directly on the datapoints- 
  // but the algorithm is clearer that way - maybe optimize later

  // wrap the phases from -pi..+pi to 0..2pi - it's more convenient to deal with this interval and
  // we don't change the original data anyway, we just use it here internally for our computations:
  for(m = 0; m < M; m++)
    p[m] = RAPT::rsWrapToInterval(p[m], 0, 2*PI);
  //plotVector(p); 

  std::vector<T> a(M-1);   // average frequencies (optimized code could avoid this array, too)
  for(m = 0; m < M-1; m++) {
    T dt = t[m+1] - t[m];                   // length of time interval t[m]...t[m+1] "delta-t"
    a[m] = T(0.5) * (f[m] + f[m+1]);        // "old" average freq in interval t[m]...t[m+1]
    T q  = p[m] + a[m] * dt * 2*PI;         // computed phase at end of interval
    T ps = p[m+1];                          // stored phase at end of current interval


    T qp = rsFindCosistentUnwrappedValue(q, p[m+1], 0.0, 2*PI);  // new
    //T qp = rsFindCosistentPhase(p[m+1], q); // q' - adjusted phase ...old



    a[m] = (qp-p[m])/(dt*2*PI);             // "new" average freq, consistent with p[m] and p[m+1]

    T dq = q - qp;  // |dq| should be (much) less than pi - otherwise the hopSize is too small for
    // correctly estimating frequencies from phase-differences - maybe return the maximum dp as
    // feedback and/or maybe have a function getMaximumPhaseDeviation (where the deviation is 
    // measured with respect to integrated frequency)

    // at m=31 and m = 32, we get a large dq (around 1.7) - this leads to wildly alternating new
    // frequencies - something is wrong - and decreasing the hop-size tends to make the problem
    // even worse - check findConsistentPhase - could id have to with range 0..2pi vs. -pi..+pi?

    // check, if new a[m] is indeed consistent (for debug):
    q = p[m] + a[m] * dt * 2*PI;                   // same computation as above - should now give
    RAPT::rsAssert(rsArePhasesConsistent(q, p[m+1])); // a phase q that is consistent with p[m+1]
  }

  // OK - we have our new desired average frequencies for the segments - from these, we now compute
  // the new frequencies at the datapoints (we are actually one equation short of determining all 
  // frequencies, so we make the choice that the last two datapoint should have the same frequency
  // as our additional equation/condition):
  f[M-1] = f[M-2] = a[M-2];
  for(m = M-3; m >= 0; m--)
    f[m] = 2*a[m] - f[m+1];
  // maybe we should have a switch, if we run the loop over the datapoints forward or backward - in
  // the backward case, we set the two last frequencies equal and in the forward case the two first
  // frequencies...can a more symmetric way be found and one that doens't enforce two equal 
  // frequencies - think about, how we could provide the "missing equation" in other ways - maybe 
  // an equation that involves all datapoints on the same footing? maybe a condition that minimizes
  // the tendency to produce alternating corrections in successive datapoints? ...maybe keep the
  // old freqs data available, compute the new freqs, obtain their difference, 
  // lowpass this difference, apply the lowpassed corrections to the old data
  // or: use the f array as computed above as preliminary and then apply an "equalize" function
  // that looks at a pair of neighbours at a time and adjusts their frequencies in a way to make 
  // them as close to equal as possible while maintaining the phase constraint, choose first 
  // pairs (0,1),(2,3),(4,5) etc. and then do it again with (1,2),(3,4),(5,6) to couple 1 to 2, too
  // etc maybe iterate until it converges to something
  // see the textfile MinDiffGivenSum.txt - the additional equation should be to minimize the 
  // sum-of-squared-differences - this leads to a constrained optimization problem soluble via 
  // Lagrange mulitpliers
  
  // preliminary:
  std::vector<T> sum = 2.0 * a;
  rsMinSqrDiffWithGivnSum(&f[0], &sum[0], M);

  // finally, write the new frequencies into the datapoints of the partial:
  for(m = 0; m < M; m++)
    setFrequency(m, f[m]);

  // this actually increases the freq-estimate bias in sinusoidalAnalysis2 - do i have a 
  // theory bug? the implementation seems good..the assert doesn't trigger - or maybe the hopsize
  // is indeed too small and we get an adjustment by more than pi?
}




//-------------------------------------------------------------------------------------------------
// inquiry

template<class T>
T rsSinusoidalPartial<T>::getMeanFreq() const
{
  size_t M = getNumDataPoints();
  if(M == 0) return T(0);
  if(M == 1) return getStartFreq();
  T freqSum   = T(0);
  T weightSum = T(0);
  for(size_t m = 0; m < M-1; m++)
  {
    T fm = 0.5 * (instParams[m].freq + instParams[m+1].freq); // mean freq between point m and m+1
    T dt = instParams[m+1].time - instParams[m].time;         // duration of segment
    freqSum   += dt*fm;
    weightSum += dt;
  }
  return freqSum / weightSum;
}

template<class T>
T rsSinusoidalPartial<T>::getMinFreq() const
{
  T minFreq = RS_INF(T);
  for(size_t i = 0; i < instParams.size(); i++)
    if(instParams[i].freq < minFreq)
      minFreq = instParams[i].freq;
  return minFreq;
}

template<class T>
T rsSinusoidalPartial<T>::getMaxFreq() const
{
  T maxFreq = -RS_INF(T);
  for(size_t i = 0; i < instParams.size(); i++)
    if(instParams[i].freq > maxFreq)
      maxFreq = instParams[i].freq;
  return maxFreq;
}

template<class T>
bool rsSinusoidalPartial<T>::willAlias(T sampleRate, bool allTheTime) const
{
  if(allTheTime == false)
    return getMaxFreq() > T(0.5) * sampleRate;
  else
    return getMinFreq() > T(0.5) * sampleRate;
}

template<class T>
void rsSinusoidalPartial<T>::getDataArrays(
  std::vector<T>& t, std::vector<T>& f, std::vector<T>& a, std::vector<T>& p) const
{
  size_t M = getNumDataPoints();
  t.resize(M);
  f.resize(M);
  a.resize(M);
  p.resize(M);
  for(size_t m = 0; m < M; m++) {
    rsInstantaneousSineParams<T> dp = getDataPoint(m);
    t[m] = dp.getTime();         // time data
    f[m] = dp.getFrequency();    // frequency data
    a[m] = dp.getAmplitude();    // amplitude data
    p[m] = dp.getWrappedPhase(); // (wrapped) phase data
  }
}

template<class T>
std::vector<T> rsSinusoidalPartial<T>::getTimeArray() const
{
  size_t M = getNumDataPoints();
  std::vector<T> t(M);
  for(size_t m = 0; m < M; m++)
    t[m] = getDataPoint(m).getTime();
  return t;
}

template<class T>
std::vector<T> rsSinusoidalPartial<T>::getFrequencyArray() const
{
  size_t M = getNumDataPoints();
  std::vector<T> f(M);
  for(size_t m = 0; m < M; m++)
    f[m] = getDataPoint(m).getFrequency();
  return f;
}

template<class T>
std::vector<T> rsSinusoidalPartial<T>::getAmplitudeArray() const
{
  size_t M = getNumDataPoints();
  std::vector<T> a(M);
  for(size_t m = 0; m < M; m++)
    a[m] = getDataPoint(m).getAmplitude();
  return a;
}

template<class T>
std::vector<T> rsSinusoidalPartial<T>::getPhaseArray() const
{
  size_t M = getNumDataPoints();
  std::vector<T> p(M);
  for(size_t m = 0; m < M; m++)
    p[m] = getDataPoint(m).getWrappedPhase();
  return p;
}



//=================================================================================================

template<class T>
void rsSinusoidalModel<T>::init(int numPartials, int numFrames)
{
  partials.clear();
  if(numPartials > 0) {
    partials.resize(numPartials);
    for(int i = 0; i < numPartials; i++)
      partials[i].init(numFrames);
  }
}

template<class T>
void rsSinusoidalModel<T>::makeFreqsConsistentWithPhases()
{
  for(auto p : partials)
    p.makeFreqsConsistentWithPhases();
}

template<class T>
T rsSinusoidalModel<T>::getStartTime() const
{
  if(partials.size() == 0)
    return T(0);
  T start = RS_INF(T);
  for(size_t i = 0; i < partials.size(); i++)
    start = rsMin(start, partials[i].getStartTime());
  return start;
}

template<class T>
T rsSinusoidalModel<T>::getEndTime() const
{
  if(partials.size() == 0)
    return T(0);
  T end = -RS_INF(T);
  for(size_t i = 0; i < partials.size(); i++)
    end = rsMax(end, partials[i].getEndTime());
  return end;
}


//=================================================================================================
// Analyzer:



//=================================================================================================
// Synthesizer:




//=================================================================================================
// Transformations:




/*

Ideas:

-frequency estimation:
-to find the actual frequency of a partial in an FFT frame, often parabolic interpolation of the
 log magnitude is used (the frequency is to be taken the maximum of the parabola that goes through
 3 adjacent bins)
-idea for refining this estimate: the transform a single sinusoid is given by an appropriately 
 shifted (and scaled) transform of the window function - try to find the optimum cross-correlation 
 lag, such that maximizes the cross-correlation between the frame's spectrum and the 
 window-transform...restricted to a certain correlation length (maybe 3...20)...maybe it could
 be useful to do it in the complex domain? maybe that would automatically account for phase
 information?
 -or: use a gaussian window - in this case the parabolic interpolation is actually exact:
  https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
  https://ccrma.stanford.edu/~jos/sasp/Gaussian_Window_Transform_I.html#fig:gaussianWindow

-peak picking/continuation:
-to match a peak in the k-th frame with exisitng active sinusoids, find the one that is closest in 
 frequency to a linear continuation of (k-1)-th and (k-2)-th frame 
-maybe special treatment has to be given to sine tracks that currently have only 1 frame - or maybe 
 rule them out by requiring that every stable partial must be at least 2 frames long to be 
 considered a partial at all...but that would require the peak picking algorithm not at one frame
 at a time but (at least) at 2 - but that may be a good idea anyway
-test the algorithm by using two sines that cross each other in the spectrogram (one ascending, one
 descending)

-create different variations of analysis/synthesis algorithms

*/