template<class T>
void rsSinusoidalPartial<T>::applyFadeIn(T fadeTime)
{
  rsInstantaneousSineParams<T> params = getFirstDataPoint();
  params.time  -= fadeTime;
  params.gain   = 0.0;
  //params.gain   = 0.5;  // for test
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
  //params.gain   = 0.5;   // for test
  params.phase += 2*PI*fadeTime*params.freq;
  params.phase  = rsWrapToInterval(params.phase, -PI, PI);
  appendDataPoint(params);
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

template<class T>
std::vector<T> rsSinusoidalSynthesizer<T>::synthesize(
  const rsSinusoidalModel<T>& model, T sampleRate) const
{
  std::vector<T> x;


  return x;
}



template<class T>
rsSinusoidalModel<T> rsSinusoidalAnalyzer<T>::analyze(
  T* sampleData, int numSamples, T sampleRate) const
{
  rsSinusoidalModel<T> model;


  return model;
}


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