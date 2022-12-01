std::vector<double> synthesizeModal(
  const rsModalFilterParameters<double>& p, double fs, int N)
{
  std::vector<double> x(N);
  rosic::rsModalFilterWithAttackDD flt;
  flt.setModalParameters(p.freq, p.amp, p.att, p.dec, p.phase, fs);
  flt.reset();
  x[0] = flt.getSample(1);
  for(int n = 1; n < N; n++)
    x[n] = flt.getSample(0);
  return x;
}

std::vector<double> synthesizeModal(
  const std::vector<rsModalFilterParameters<double>>& p, double fs, int N)
{
  std::vector<double> x(N);
  rsArrayTools::fillWithZeros(&x[0], N);
  rosic::rsModalFilterWithAttackDD flt;
  for(size_t i = 0; i < p.size(); i++) {
    flt.setModalParameters(p[i].freq, p[i].amp, p[i].att, p[i].dec, p[i].phase, fs);
    flt.reset();
    x[0] += flt.getSample(1);
    for(int n = 1; n < N; n++)
      x[n] += flt.getSample(0); }
  return x;
}

template<class T>
std::vector<rsModalFilterParameters<T>> rsModalAnalyzer<T>::getModalModel(
  const RAPT::rsSinusoidalModel<T>& model)
{
  std::vector<rsModalFilterParameters<T>> p(model.getNumPartials());
  for(size_t i = 0; i < model.getNumPartials(); i++)
  {
    p[i] = getModalModel(model.getPartial(i));


    //// test - just for the rhodes:
    //if(i == 8) {  
    //  p[i].dec *= 30;
    //  //p[i].amp *= 10;    
    //}
  }
  return p;
}

template<class T>
rsModalFilterParameters<T> rsModalAnalyzer<T>::getModalModel(
  const RAPT::rsSinusoidalPartial<T>& inputPartial)
{
  // todo: estimate the actual parameters:
  // frequency: take average freq of partial
  // amplitude: take peak/max amplitude fo partial
  // attack:    time-instant of peak
  // decay:     compare average amplitudes of 1st and 2nd half after the peak
  // phase:     start phase in model


  rsPartialBeatingRemover<T> beatRemover;
  RAPT::rsSinusoidalPartial<T> partial = inputPartial;

  //beatRemover.removeAmplitudeBeating(partial);
  // this doesn't seem to work well with the rhodes sample



  rsModalFilterParameters<T> params;

  int peakIndex = partial.getMaxAmpIndex();
  params.amp = partial.getDataPoint(peakIndex).getAmplitude();
  params.att = partial.getDataPoint(peakIndex).getTime();

  T thresh = 1.e-05;
  if(params.amp < thresh)
  {
    // below the noise threshold, we may get confusing data...
    params.amp   = 0.0;
    params.freq  = 1000.0;
    params.dec   = 0.1;
    params.att   = 0.01;
    params.phase = 0.0;
    return params;
    //return rsModalFilterParameters<T>(); 
  }


  // maybe use parabolic interpolation for more accurate estimation


  int M = (int) partial.getNumDataPoints();
  int refIndex = peakIndex;
  //refIndex = M / 2;  // test 

  //params.freq  = partial.getMeanFreq();
  //params.freq  = partial.getMeanFreq(0, refIndex); // maybe have member function estimateFreq
  //params.freq  = partial.getMeanFreq(refIndex, M-1);

  //params.freq = estimateFrequency(partial, 0, M-1);
  params.freq = estimateFrequency(partial, 4, M-5); // 4, M-5 ad-hoc
  // maybe we should cut off the transient before taking the mean?
  // i think, the freq-estimation is still inaccurate for the faster decaying partials

  //params.phase = partial.getFirstDataPoint().getWrappedPhase();
  params.phase = estimatePhaseAt(partial, refIndex, params.freq, T(0));
  // phase seems to be still wrong - but wait - isn't the measured phase supposed to occur halfway
  // between a pair of datapoints?



  // maybe this is not the best strategy either - maybe take the phase at the amplitude peak and
  // compute the start-phase from that and the mean frequency ...it seems, it would be a good idea
  // to implement a class where the



  // estimate decay (we take the average of the log-amplitudes of the two halfs of the remaining
  // signal after the peak)
  params.dec = 0.5;  // preliminary

  /*
  int numDataPoints = (int) partial.getNumDataPoints();
  int start1 = peakIndex;
  int length = numDataPoints - peakIndex;
  int start2 = start1 + length/2;
  T mean1 = 0, mean2 = 0;
  size_t i;

  int count = 0;   // get rid
  for(i = start1; i < start2; i++) {
    mean1 += log(partial.getDataPoint(i).getAmplitude());
    count++; }
  mean1 /= count;

  count = 0;
  for(i = start2; i < numDataPoints; i++) {
    mean2 += log(partial.getDataPoint(i).getAmplitude());
    count++; }
  mean2 /= count;
  // that doesn't work because some amplitudes may be zero
  */

  /*
  int searchStart = peakIndex + (partial.getNumDataPoints()-peakIndex)/2;
  int peakIndex2 = partial.getMaxAmpIndex(searchStart);
  T t1 = partial.getDataPoint(peakIndex).getTime();
  T a1 = partial.getDataPoint(peakIndex).getAmplitude();
  T t2 = partial.getDataPoint(peakIndex2).getTime();
  T a2 = partial.getDataPoint(peakIndex2).getAmplitude();
  T dt = t2 - t1; // time difference
  T ra = a1 / a2; // amplitude ratio todo: catch a2 == 0 as special case
  // from dt and ra, we can compute the decay time tau...-> look up formula....
  */

  // ...hmm - maybe this is not so good - maybe it would be better to search through the 
  // amplitude array for the index/time, where the amplitude is peakAmp/e - but for this, we need 
  // to assume a monotonically decreasing amplitude envelope after the peak - maybe if it's not
  // obtain "meta-envelopes" repeatedly until it is monotonically decreasing
  T targetAmp = params.amp / EULER;  // peakAmp / e
  //T tau = 0;  // or should we init with inf?
  int i;
  for(i = peakIndex+1; i < M; i++)
  {
    if(partial.getDataPoint(i).gain < targetAmp)
    {
      params.dec = partial.getDataPoint(i).time - params.att;
      // this is very coarse - todo: interpolate (linearly on the dB-scale)

      break;
    }
  }
  if(i == M)
  {
    //rsError(); 
    params.dec = 1.1*params.att;
    // preliminary - todo repeat search with c = 0.5 * c, see comment below
    // ..doesn't happen with the rhodes sample
  }

  // maybe sometimes, the sample isn't long enough to contain the data, where it has decayed to
  // peak/e (because of an early fade-out or palm-muting, whatever) - then we should use 
  // c*peak/e for some c < 1 and multiply the decay-time by that same c - maybe have a loop
  // with exponentially decreasing c, i.e. c = 1,0.5,0.25,0.125,...
  // what about samples with no decay at all, i.e. sustained sounds?

  // maybe instead of averaging, take the maximum value of the section that starts halfway
  // after the peak - maybe the function getMaxAmpIndex can take a start-search index as parameter

  // this is just for development - later, we have to do something more appropriate
  if(params.att >= 0.99*params.dec)
  {
    params.att = 0.99*params.dec;
    //params.dec = params.att/0.99;
  }
  // maybe instead of shortening the attack, we should lengthen the decay? ..doesn't seem to make
  // a big difference in the rhodes sample - it seems, the underestimated amplitude is to blame 
  // there
  // oh - but lengthening artificially increases the amplitudes of high harmonics - shorten seems
  // better overall
  params.att = rsMax(params.att, 0.0001); // attack of 0 is not allowed

  // maybe use gradient descent for more accurate estimates


  //plotModeVsSineAmpEnv(params, partial);  // for development
  // ok - the rhodes sample really has a problem with the 9th partial which causes the estimated
  // decay to be way too short - i think, that's the main reason for the "wiggle" in the 
  // resynthesized signal: an (amost) missing partial (due to short decay) which is supposed to be
  // there - the amp envelopes of some of the partial just do not fit the attack/decay model
  // very well

  return params;
}

template<class T>
T rsModalAnalyzer<T>::estimatePhaseAt(
  const RAPT::rsSinusoidalPartial<T>& partial, int i, T f, T t)
{
  rsAssert(i >= 0 && i < (int) partial.getNumDataPoints()-1);
  T ti = partial.getTime(i);   // time stamp at datapoint i

  // test - phase data is estimated halfway between two datapoints (i think):
  //ti += (partial.getTime(i+1) - ti)/2;
  // that doesn't seem to help - try it with a simpler sound containing only one single mode
  // ...hmm...it seems more accurate without -> figure out why

  T pi = partial.getPhase(i);  // phase at time ti
  T dt = t - ti;               // time difference
  T pt = pi + 2*PI*f*dt;       // extrapolated phase at time t (assuming const freq in t..ti)
  pt += PI/2;                  // because analysis computes cosine phase
  pt  = rsWrapToInterval(pt, -PI, PI);
  return rsRadiantToDegree(pt);
}

template<class T>
T rsModalAnalyzer<T>::estimateFrequency(
  const RAPT::rsSinusoidalPartial<T>& partial, int start, int end)
{
  // move this code to rsSinusoidalPartial (getUnwrappedPhase or something)
  std::vector<T> t = partial.getTimeArray();
  std::vector<T> f = partial.getFrequencyArray();
  std::vector<T> p = partial.getPhaseArray();
  std::vector<T> u = rsSinusoidalProcessor<T>::unwrapPhase(t, f, p); // unwrapped phase

  T freq = (u[end]-u[start]) / (2*PI*(t[end]-t[start])); // more accurate

  // and/or have a function getMeanFreqAccurate ...or let the getMeanFreq function have a boolean
  // flag accountForPhase or something

  //freq = partial.getMeanFreq(start, end); // simple, coarse

  return freq;

  //return partial.getMeanFreq(start, end);

  // i think, for more accurate freq estimates, we need to take also the instantaneous phases into 
  // account - maybe dismiss the first and last K cycles in the frequency estimation (their data
  // may contain transient artifacts), obtain an unwrapped phase array for that middle section and
  // take the average freq as (endPhase-startPhase)/(2*PI*(endTime-startTime))
}

//=================================================================================================

template<class T>
std::vector<rsModalFilterParameters<T>> rsModalAnalyzer2<T>::analyze(T* x, int N)
{
  //...............................................................................................
  // Step 1: 
  // Figure out the mode frequencies using a big FFT on the whole signal and find the  peak freqs:
  using Vec = std::vector<T>;
  int N2 = rsNextPowerOfTwo(N);
  Vec x2(N2), mags(N2);
  rsZero(x2);                              // May not be needed
  for(int n = 0; n < N; n++)
    x2[n] = x[n];
  rsFourierTransformerRadix2<T> ft;
  ft.setDirection(ft.FORWARD);
  ft.setBlockSize(N2);
  ft.setNormalizationMode(ft.NORMALIZE_ON_FORWARD_TRAFO);
  ft.getRealSignalMagnitudes(&x2[0], &mags[0]);
  // ToDo:
  // -Maybe restrict the maximum FFT lenth: maybe use a freq-resolution parameter, by deafult 
  //  1 Hz, and take its reciprocal as max length in seconds for the FFT to avoid excessively 
  //  large FFTs. 1 second should be enough for the preliminary analysis which is only meant for 
  //  tuning the filter-bank anyway...or maybe have a max-pre-analysis length in seconds.

  //rsPlotSpectrum(mags, sampleRate, -100.0, false);  // for development
  // When peak-finding is implemented, maybe plot the specttrum with markers at the found peak
  // frequencies.

  // Apply peak-masking:
  Vec magsMasked = mags;
  rsPeakMasker<T> pm;
  T freqSeparation = 10;         // make user parameter, find better name
  //T binDecay  = freqSeparation * N2 / sampleRate;
  pm.setDecaySamples(freqSeparation * N2 / sampleRate);
  pm.applyForward( &magsMasked[0], &magsMasked[0], N2);
  pm.applyBackward(&magsMasked[0], &magsMasked[0], N2);

  // Find relevant peaks:
  T threshRatio = 0.0005;  // make use parameter (in dB)
  using PF = rsPeakFinder<T>;
  using AT = rsArrayTools;
  const int precision = 1;         // 1: use a parabolic fit
  Vec peakPositions, peakHeights;  // todo: maybe reserve some memory here, maybe maxNumModes
  T pos, height, maxHeight;
  int kMax = AT::maxIndex(&mags[0], N2);
  PF::exactPeakPositionAndHeight(&mags[0], N2, kMax, precision, &pos, &maxHeight); // global max
  for(int n = 1; n < N2-1; n++) {
    if( AT::isPeak(&magsMasked[0], n) ) {
      PF::exactPeakPositionAndHeight(&mags[0], N2, n, precision, &pos, &height);
      if(height >= threshRatio * maxHeight) {
        peakPositions.push_back(pos);
        peakHeights.push_back(height); }}}

  // ToDo:
  // -Keep only the maxNumModes modes with the highest heights. Maybe for that, we need to create
  //  a struct that contains height and position 
  // -Maybe use rsvector2D for that and define a < operator that compares based on x first then on
  //  y. Instead of having parallel arrays peakPoisitions/Heights, we'd use one array of 2D vectors
  //  storing the height in x and the position iny, such that sorting works as intended.

  // Plot results of the pre-analysis (for development):
  GNUPlotter plt;
  Vec freqs(N2);
  ft.binFrequencies(&freqs[0], N2, sampleRate);
  peakPositions = peakPositions * sampleRate / T(N2);
  plt.addDataArrays(N2, &freqs[0], &mags[0]);
  plt.addDataArrays(N2, &freqs[0], &magsMasked[0]);
  plt.addDataArrays((int)peakHeights.size(), &peakPositions[0], &peakHeights[0]);
  plt.setGraphStyles("lines", "lines", "points");
  plt.setPixelSize(1200, 400);
  plt.plot();
  // OK - this looks good. Maybe later we could use rsPeakPicker which includes masking plus some
  // more sophisticated ideas. For the time being, the masking works well

  //...............................................................................................
  // Step 2:
  // Analyze each mode one at a time by bandpassing the signal with a bandpass tuned to the 
  // respective modal frequency and then using an envelope follower on the bandpassed signal

  using ModalParams = rsModalFilterParameters<T>;
  std::vector<ModalParams> mp;

  
  // Notes:
  // -The bandwidth of the bandpass should sufficiently suppress adajacent partials (calling for 
  //   smaller bandwidth) but without introducing too much time-domain smoothing on the estimated 
  //   envelope (calling for a larger bandwidth). 
  // -Maybe we can strike an optimal compromise by (somehow) making it dependent on some 
  //  preliminary mode-bandwidth measurement from the FFT spectrum? 
  // -Maybe to improve the separation of partials in the analysis, use a bi-notch in addition to
  //  the bandpass. the notches should be placed at the two adjacent partials. The bandwidths of 
  //  the notches should be tweaked such that there's a local maximum at the mode's freq in the 
  //  bi-notch response, too. -> Some filter design work necessary: derive equation for the freq
  //  of the local maximum of the bi-notch. Maybe make a class rsBiNotchFilter or 
  //  rsDoubleNotchFilter for this purpose. For the bottommost and topmost modes use only one 
  //  notch - or maybe place the other notch at DC and fs/2. The ringing time of the notch should 
  //  be such that it doesn't increase the overall ringing time too much.
  // -Maybe it would also be good, if the bandpasses feature a flat-top, especially when the 
  //  partials have a time-varying frequency?
  // -Or maybe we should try to approximate a Gaussian filter for an optimal compromise between
  //  bandwidth and time-domain smoothing/smearing of the envelope?
  // -To the bandpassed signal, apply an envelope follower and from the envelope, find the 
  //  peak and estimate the decay by fitting an exponential. Somewhere, I already have some code
  //  to fit a sum of exponentials - maybe use that. Apply it to the section from the peak onwards.
  //  maybe we can use error-weigths of zero for the section before? It will result a decay and 
  //  amplitude for the decay. With that in hand, maybe we can compute the remaining parameters
  //  that are responsible for the attack...
  //  

  return mp;
}




// instantiations (maybe move elsewhere):
//template class rsModalAnalyzer<double>;

//template std::vector<rsModalFilterParameters<double>> 
//  getModalModel(const RAPT::rsSinusoidalModel<double>& model);

