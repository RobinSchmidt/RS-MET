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