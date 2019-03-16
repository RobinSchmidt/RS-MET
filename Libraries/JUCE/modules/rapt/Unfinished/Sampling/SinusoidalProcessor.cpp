
template<class T>
std::vector<T> rsSinusoidalProcessor<T>::unwrapPhase(const std::vector<T>& t,
  const std::vector<T>& f, const std::vector<T>& wp)
{
  size_t M = t.size();
  RAPT::rsAssert(f.size()  == M);
  RAPT::rsAssert(wp.size() == M);
  std::vector<T> up(M);  // unwrapped phase

  // obtain preliminary uwrapped phase data points by numerically integrating the frequency:
  rsNumericIntegral(&t[0], &f[0], &up[0], (int)M, wp[0]);
  up = 2*PI*up; // convert from "number of cycles passed" to radians

  // incorporate the target phase values into the unwrapped phase:
  bool accumulatePhaseDeltas = true; // not accumulating makes no sense, i think
  for(size_t m = 0; m < M; m++) {
    T wp1 = RAPT::rsWrapToInterval(up[m], 0, 2*PI);  // 0..2*pi
    T wp2 = RAPT::rsWrapToInterval(wp[m], 0, 2*PI);  // 0..2*pi
    T d   = wp2-wp1;            // -2*pi..2*pi, delta between target phase and integrated frequency

    // these adjustments here are related to the phase-desync bursts in resynthesized signals (the
    // resynthesized phases get temporarily desynchronized from the original phases, leading to
    // short sinusodial bursts in the residual) - but commenting them out leads to even more bursts
    if(d < 0) d += 2*PI;        // 0..2*PI
    if(d > PI)                  // choose adjustment direction of smaller phase difference
      d -= 2*PI;                // -pi..pi
    up[m] += d;                 // re-adjust final unwrapped phase
    if(accumulatePhaseDeltas)
      for(size_t k = m+1; k < M; k++) // re-adjustment at m should also affect m+1, m+2, ...
        up[k] += d;
  }
  // maybe provide an alternative implementation that uses the measured (unwrapped) phases y and 
  // the measured instantaneous frequencies as y'' in a Hermite interpolation scheme (this is how
  // it's described in the literature). to unwrap the phase of datapoint m, take the phase of m-1
  // and add the integral over t[m-1]..t[m] of the average frequency (f[m-1]+f[m])/2 and choose
  // p[m] + 2*pi*k as unwrapped where k is chosen such that p[m] is closest to value obtained from
  // integration - use this function here as dispatcher between the different algorithms

  //rsPlotVector(up);

  return up;
}




// x: values to be computed, s: desired sum-values (constraints)
template<class T>
void rsMinSqrDiffWithGivnSum(T* x, T* s, int N)
{
  T *w = nullptr;
  rsMinSqrDifFixSum(x, N, s, w);
  // todo: when N is even (or odd? look at experiment) create a weight-vector with values 
  // of 1/2 for the first and last value (and 1 everywhere else);


  // this thesis: https://web.stanford.edu/group/SOL/dissertations/bradley-thesis.pdf
  // says that for symmetric, positive definite matrices, scaling to unit diagonal is effective for
  // making the problem better conditioned
}
template<class T>
void rsSinusoidalProcessor<T>::makeFreqsConsistentWithPhases(rsSinusoidalPartial<T>& partial)
{
  std::vector<T> t = partial.getTimeArray();
  std::vector<T> f = partial.getFrequencyArray();
  std::vector<T> p = partial.getPhaseArray();
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
    //T qp = rsConsistentUnwrappedValue(q, p[m+1], 0.0, 2*PI);  // q' - adjusted phase new - delete
    T qp = rsConsistentUnwrappedValue0(q, p[m+1], 2*PI);        // q' - adjusted phase
    a[m] = (qp-p[m])/(dt*2*PI);             // "new" average freq, consistent with p[m] and p[m+1]
    // ...wait - what if this becomes negative?
    rsAssert(a[m] > T(0));

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

  rsPlotVector(a); 

  // OK - we have our new desired average frequencies for the segments - from these, we now compute
  // the new frequencies at the datapoints (we are actually one equation short of determining all 
  // frequencies, so we make the choice that the last two datapoint should have the same frequency
  // as our additional equation/condition):
  f[M-1] = f[M-2] = a[M-2];
  //f[M-1] = f[M-2] = T(0.5) * (a[M-2] + a[M-3]);  // test
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
    partial.setFrequency(m, f[m]);

  // this actually increases the freq-estimate bias in sinusoidalAnalysis2 - do i have a 
  // theory bug? the implementation seems good..the assert doesn't trigger - or maybe the hopsize
  // is indeed too small and we get an adjustment by more than pi?
}
