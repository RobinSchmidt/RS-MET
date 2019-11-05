
// helper functions (maybe  move into class rsSinusoidalPartial )
/*
// useful for freq de-biasing and for phase interpolation in synthesis
template<class T>
T rsFindCosistentPhase(T storedPhase, T computedPhase) 
{
  T p = storedPhase;     // 0..2pi
  T q = computedPhase;   // 0..inf
  T tmp = p;
  T k   = 0;
  T d;
  while(true) {
    tmp = p + 2*k*PI;
    d   = q - tmp;
    if(rsAbs(d) <= PI + 1.e-13)
      break;
    k  += T(1);
  }
  return tmp;
  // this may lead to an endless loop when 
  // storedPhase = 6.2804880025844536, computedPhase = 0.0044408067227543532
  // -> make unit test - i think, the condition rsAbs(d) < PI is not enough ..maybe
  // it needs a tolerance? nope: using rsAbs(d) <= PI + 1.e-13 also doesn't work
}

// whoa - this is very tricky - isn't there a simpler way for this?
template<class T>
T rsPhaseError(T p1, T p2)
{
  p1  = RAPT::rsWrapToInterval(p1, 0, 2*PI);  // 0..2pi  // rename to rsWrapToRange...or just rsWrap
  p2  = RAPT::rsWrapToInterval(p2, 0, 2*PI);  // 0..2pi
  T d = p2-p1;                                // -2pi..2pi
  d   = RAPT::rsWrapToInterval(d,  0, 2*PI);  // 0..2pi
  return d;
}

template<class T>
bool rsArePhasesConsistent(T p1, T p2, T tol = 1.e-13)
{
  T d = rsPhaseError(p1, p2);
  if(d < tol)
    return true;
  if(d-2*PI < tol)
    return true;
  return false;
}
*/


//=================================================================================================

template<class T>
std::vector<T> rsSinusoidalSynthesizer<T>::synthesize(const rsSinusoidalModel<T>& model) const
{
  int N = model.getLengthInSamples(sampleRate);
  std::vector<T> x(N);
  RAPT::rsArray::fillWithZeros(&x[0], N);
  for(size_t i = 0; i < model.getNumPartials(); i++)
    synthesizePartial(model.getPartial(i), &x[0], N, -model.getStartTime()); 
  return x;
}

template<class T>
void rsSinusoidalSynthesizer<T>::synthesizePartial(
  const rsSinusoidalPartial<T>& partial, T* x, int xLength, T timeShift) const 
{
  // figure out number of samples to produce:
  int nStart = (int) floor(sampleRate * (partial.getStartTime() + timeShift));
  int nEnd   = (int) ceil( sampleRate * (partial.getEndTime()   + timeShift)) + 1;
  nStart = rsClip(nStart, 0, xLength);
  nEnd   = rsClip(nEnd,   0, xLength);
  int N = nEnd - nStart;  // number of samples to generate

  // create time axis and interpolate instantaneous amplitude and phase up to sample rate:
  std::vector<T> td = partial.getTimeArray();
  td = td + timeShift; // todo: implement vector += scalar operator and use: td += timeShift;
  T Ts = T(1) / sampleRate;  // sampling interval
  std::vector<T> t(N);
  for(size_t n = 0; n < N; n++)          // fill time-array
    t[n] = (nStart + n) * Ts;

  std::vector<T> a = getInterpolatedAmplitudes(partial, td, t);
  std::vector<T> p = getInterpolatedPhases(    partial, td, t);

  //rsPlotVector(a);
  //rsPlotVector(p);
  //rsPlotVectors(a, p);
  //rsPlotVector(rsDifference(p));

  // synthesize the sinusoid and add it to what's already there:
  std::vector<T> s(N); // needed here only for plotting, remove for production code
  for(size_t n = 0; n < N; n++)
    s[n] = x[nStart+n] += a[n] * cos(p[n]);
  // we use the cosine (not the sine) because that's what's used in the literature - probably 
  // because it's consistent with representing real sinusoids as the real part of complex
  // sinusoids (using the sine, we would have to take the imaginary part instead)


  // it seems cubic interpolation for the phase and linear for the amplitude is most suitable,
  // although, for the amplitude, we may also use cubic - but linear for the phase leads to audible
  // artifacts (sort of clicks) at the segement junctions
  // -maybe linear interpolation of frequency with subsequent integration would work (due to the
  //  smoothing effect of integration) - but that would make it difficult to incorporate the 
  //  target phases (i think, we would have to produce an interpolated phase-delta array and add 
  //  that to an interpolated preliminary unwrapped-phase array)
  // -maybe using cubic interpolation for frequency, than integrating and then adding an 
  //  interpolated phase-delta array could give and even smoother freq-trajectory? maybe try it
  // -or maybe use higher order numeric integration on the non-interpolated freq-data?
  // -but when smoothing comes into play later, linear interpolation of the phase might be not so
  //  bad as it is without smoothing
  // -however - when the data comes from an analysis, it will be typically much denser in time than
  //  in this test (like one datapoint per cycle), so the difference between the interpolation
  //  methods may be not that important anymore - but for sparse data, it is crucial
  // -maybe the synthesizer should have "presets" for most useful combinations of synthesis 
  //  parameters

  // maybe let the select to apply a bidirectional smoothing filter to the amplitude and/or phase
  // ...and/or maybe to the frequency data before integration (but then it would have to opreate on
  // non-uniformly sampled data). when one is to be applied to the phase, we should first subtract
  // out the linear trend, then apply it and add the trend back

  // maybe provide hook-functions that can manipulate the amplitude/phase data after interpolation
  // maybe have subclasses that also work with real and imaginare parts instead of magnitude/phase
  // (and apply filters to those) ...but that is already something for the transformations


  //GNUPlotter plt;
  ////plt.addDataArrays(M, &td[0], &fd[0]);
  ////plt.addDataArrays((int)M, &td[0], &upd[0]);
  //plt.addDataArrays((int)N, &t[0],  &p[0]); // interpolated phase
  //plt.addDataArrays((int)N, &t[0],  &a[0]);   // interpolated amplitude
  ////plt.addDataArrays((int)N, &t[0],  &s[0]);   // produced sinusoid
  //plt.plot();
}

template<class T>
std::vector<T> rsSinusoidalSynthesizer<T>::getInterpolatedAmplitudes(
  const RAPT::rsSinusoidalPartial<T>& partial, 
  const std::vector<T>& td, 
  const std::vector<T>& t) const
{
  int M = (int) td.size(); // td: time axis data
  int N = (int) t.size();  // t: interpolated time axis (to sample-rate)
  std::vector<T> a(N);     // a: interpolated amplitude values
  std::vector<T> ad  = partial.getAmplitudeArray();
  RAPT::rsAssert(ad.size() == M);
  if(cubicAmplitude) 
    rsNaturalCubicSpline(&td[0], &ad[0], (int)M, &t[0], &a[0], (int)N);
  else               
    rsInterpolateLinear( &td[0], &ad[0], (int)M, &t[0], &a[0], (int)N);
  return a;
}

template<class T>
std::vector<T> rsSinusoidalSynthesizer<T>::getInterpolatedPhases(
  const RAPT::rsSinusoidalPartial<T>& partial, 
  const std::vector<T>& td, 
  const std::vector<T>& t) const
{
  typedef PhaseInterpolationMethod PIM;
  switch(phaseInterpolation)
  {
  case PIM:: tweakedFreqIntegral: return phasesViaTweakedIntegral(partial, td, t);
    //case PIM:: cubicHermite:        return phasesCubicHermite(      partial, td, t);
    //case PIM:: quinticHermite:      return phasesQuinticHermite( partial, td, t);
  default: return phasesHermite(partial, td, t, 
    phaseInterpolation == PhaseInterpolationMethod::quinticHermite);
  }
}

template<class T>
std::vector<T> rsSinusoidalSynthesizer<T>::phasesViaTweakedIntegral(
  const RAPT::rsSinusoidalPartial<T>& partial,
  const std::vector<T>& td,
  const std::vector<T>& t) const
{
  int M = (int) td.size(); // td: time axis data
  int N = (int) t.size();  // t: interpolated time axis (to sample-rate)
  std::vector<T> p(N);     // p: interpolated phase values
  std::vector<T> fd  = partial.getFrequencyArray();
  std::vector<T> wpd = partial.getPhaseArray();       // rename to pd
  std::vector<T> upd = rsSinusoidalProcessor<T>::unwrapPhase(td, fd, wpd);

  bool cubicPhase = true; // was user parameter - but linear phase interpolation makes no sense
  //cubicPhase = false;     // test
  if(cubicPhase) 
    rsNaturalCubicSpline(&td[0], &upd[0], (int)M, &t[0], &p[0], (int)N);
  else           
    rsInterpolateLinear( &td[0], &upd[0], (int)M, &t[0], &p[0], (int)N);

  //rsPlotVector(rsDifference(p));

  return p;
}

template<class T>
std::vector<T> rsSinusoidalSynthesizer<T>::phasesHermite(
  const RAPT::rsSinusoidalPartial<T>& partial,
  const std::vector<T>& td,
  const std::vector<T>& t, bool quintic) const
{
  // the quintic hermite would have to duplicate most of the code her, maybe instead let a boolean
  // "qunitic" parameter be passed here and do if(quitic) here at appropriate places

  //bool quintic = false;
  //bool quintic = true; 
  // does not yet work - maybe the computed k works only with cubic  interpolation and with
  // quitic phase inetrpolation we should use also another formula for k? ...but why should
  // that be the case...but it seems indeed that the computed k messes up the quintic
  // interpolation method - so maybe we must be consistent - for quintic phase use cubic 
  // frequency estimation of fa

  // init data arrays:
  std::vector<T> pd = partial.getPhaseArray();     // pd: phase data
  std::vector<T> fd = partial.getFrequencyArray(); // fd: freq data
  int M = (int) td.size();      // td:  time axis data
  int N = (int) t.size();       // t:   interpolated time axis (to sample-rate)
  std::vector<T> p(N);          // p:   interpolated phase values
  std::vector<T> fdd;           // fdd: frequency derivative data
  if(quintic) {
    fdd.resize(M);
    rsNumericDerivative(&td[0], &fd[0], &fdd[0], M, false);
    //plotVector(fdd);
  }

  // declare some variables and loop over datapoints:
  T f2w = T(2*PI);         // factor to convert from frequency f to radian frequency omega w
  T Ts  = T(1)/sampleRate; // sampling interval
  T y0[3];                 // p0,w0,w0' phase, omega, omega' at start of current segment
  T y1[3];                 // p1,w1,w1' phase, omega, omega' at end of current segment
  T a[6];                  // polynomial coefficients (3rd or 5th order)
  int n = 0;               // sample index
  for(int m = 0; m < M-1; m++) {  // loop over the datapoints
    T t0 = td[m],  t1 = td[m+1];  // start and end time of segment
    T p0 = pd[m],  p1 = pd[m+1];  // start and preliminary end phase of segment
    T f0 = fd[m],  f1 = fd[m+1];  // start and end frequency of segment
    T w0 = f2w*f0, w1 = f2w*f1;   // start and end omega of segment
    T dt = t1 - t0;               // length of segment (in seconds)



    T s = 1.0;  // make parameter, 1: regular quintic, 0: quitic with 0 2nd derivative
    T f0d = 0, f1d = 0; // frequency derivative at start and end
    T w0d = 0, w1d = 0; // for the quintic interpolation
    if(quintic)  {
      f0d = s*fdd[m]; f1d = s*fdd[m+1]; // setting them to zero gives an interpolant similar to the
      w0d = f2w*f0d;  w1d = f2w*f1d;    // cubic, but not exactly the same
    }

    T fa  = 0.5 * (f0 + f1);      // average frequency of segment (in Hz) - todo: maybe try geometric mean, or generalized mean
    //T fa  = sqrt(f0 * f1);
    //T k   = round(dt * fa);       // number of cycles passed between t0 and t1
    //T k   = floor(dt * fa); 
    //k -= 1;
    //p1   += 2*k*PI;               // adjust end-phase of segment
    // ...phase adjustment is still questionable - is the formula for k really the most reasonable?
    // and/or maybe computation of fa should also take p0 and p1 into account? for quintic 
    // interpolation use the integral over a cubic hermite frequency interpolant
    // make a function, meanFreq(f0, f1, dt, f0p, f1p)

    //p1 = rsFindCosistentPhase(p1, p0 + 2*PI*fa*dt); // old
    p1 = rsConsistentUnwrappedValue(p0 + 2*PI*fa*dt, p1, 0.0, 2*PI); // new - needs test
    // todo: use rsConsistentUnwrappedValue0



    // compute cubic or quintic Hermite coefficients for the current segment:
    T scl = dt;
    y0[0] = p0;            // initial phase
    y1[0] = p1;            // final phase
    y0[1] = w0*scl;        // initial phase derivative
    y1[1] = w1*scl;        // final phase derivative
    if(quintic) {
      y0[2] = w0d*scl;     // initial omega derivative
      y1[2] = w1d*scl;     // final omega derivative
      getHermiteCoeffs2(y0, y1, a);
    }
    else
      getHermiteCoeffs1(y0, y1, a);

    // evaluate polynomial at the sample-points:
    scl = T(1) / scl;
    int order = 3 + 2*quintic;     // 3 or 5
    while(n*Ts < t1) {
      p[n] = rsPolynomial<T>::evaluate(scl*(t[n]-t0), a, order);
      n++;
      if(n == N)
        break;    // we are done, interpolated phases are ready
    }
  }

  return p;
}
// idea: maybe the frequency interpolation should be done in the log-domain, i.e. interpolate
// pitches instead of frequencies


/*
template<class T>
std::vector<T> rsSinusoidalSynthesizer<T>::unwrapPhase(const std::vector<T>& t,
  const std::vector<T>& f, const std::vector<T>& wp) const
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
*/


// template instantiation:
//template class SinusoidalSynthesizer<double>;

/*
Ideas:
-let the user select between different phase-interpolation methods:
-integrate-freq and re-adjust (unwrapped) phase (with or without accumulation - or maybe only 
with), this can be done with linear interpolation or natural cubic splines
-hermite interpolation where the target-derivative values are obtained from the frequency data, 
like (x0,y0,y0') = (t0,p0,w0) and (x1,y1,y1') = (t1,p1,w1) where w0 = 2*pi*f0/fs, etc.
...but instead of p1 use p1 + k*2*pi for a suitably chosen k, maybe k = round(dt*fa) where
dt = t1-t0, fa = (f0+f1)/2 (time-delta and average frequency over the interval)
-the formula for fa is simple but maybe a better estimate for the average frequency can be 
obtained by taking into account the frequency derivative (which can be computed numerically) and
using the integral of hermite polynomial: fa = 1/(t1-t0) * integral_t0^t1 f(t) dt where f(t)
is the cubic hermite polynomial
-if we have a numeric frequency derivative available anyway, we could also use quintic hermite 
interpolation for the phase
-for amplitude interpolation, we could also let the user choose hermite interpolation - maybe 
that's less susceptible to overshooting artifacts when we add small fade-in/out sections... but
maybe it's similar to natural cubic splines -> try it, be sure to try numeric differentiation with
and without extrapolation
*/
