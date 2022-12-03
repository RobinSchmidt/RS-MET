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


  // Re/allocate memory for our temoprary processing buffers
  int N2 = rsNextPowerOfTwo(N);
  buf1.resize(N2); rsZero(buf1);
  buf2.resize(N2); rsZero(buf2);
  // We play a bit of buffer ping pong with our two buffers buf1, buf2 here. St various steps of 
  // the pre-analysis algo, they contain the following content:
  //   (1)  buf1: zero-padded input signal
  //   (2)  buf2: FFT magnitudes
  //   (3)  buf1: FFT magnitudes (copied from buf2)
  //   (4)  buf2: masked FFT magnitudes
  // After the pre-analysis has finished, buf1 and buf2 can be re-used for other purposes in the 
  // main analysis step. See comments there for how they will be used there...

  // Find magnitude spectrum of zero-padded input x (result goes into buf1):
  for(int n = 0; n < N; n++)
    buf1[n] = x[n];
  ft.setDirection(ft.FORWARD);
  ft.setBlockSize(N2);
  ft.setNormalizationMode(ft.NORMALIZE_ON_FORWARD_TRAFO);
  ft.getRealSignalMagnitudes(&buf1[0], &buf2[0]);  // Why can't this be used in place?
  rsCopy(buf2, buf1);
  // ToDo:
  // -Maybe restrict the maximum FFT lenth: maybe use a freq-resolution parameter, by default 
  //  1 Hz, and take its reciprocal as max length in seconds for the FFT to avoid excessively 
  //  large FFTs. 1 second should be enough for the preliminary analysis which is only meant for 
  //  tuning the filter-bank anyway...or maybe have a max-pre-analysis length in seconds.
  // -If we do this, the logic how long the buffers need to be becomes more complex. I think
  //  we need N2 = max(N, min(rsNextPowerOfTwo(N), maxLengthFFT))  ...but this needs verifciation.

  //rsPlotSpectrum(buf2, sampleRate, -100.0, false);  // for development
  // When peak-finding is implemented, maybe plot the specttrum with markers at the found peak
  // frequencies.

  // Apply peak-masking to FFT magnitudes (result goes into buf2):
  rsPeakMasker<T> pm;
  pm.setDecaySamples(maskWidth * N2 / sampleRate);
  pm.applyForward( &buf2[0], &buf2[0], N2);
  pm.applyBackward(&buf2[0], &buf2[0], N2);

  // Find relevant peaks:
  peakPositions.clear();
  peakHeights.clear();
  T threshRatio = rsDbToAmp(threshDb);
  using PF = rsPeakFinder<T>;
  using AT = rsArrayTools;
  const int precision = 1;         // 1: use a parabolic fit
  T pos, height, maxHeight;
  int kMax = AT::maxIndex(&buf1[0], N2);
  PF::exactPeakPositionAndHeight(&buf1[0], N2, kMax, precision, &pos, &maxHeight); // global max
  for(int n = 1; n < N2-1; n++) {
    if( AT::isPeak(&buf2[0], n) ) {
      PF::exactPeakPositionAndHeight(&buf1[0], N2, n, precision, &pos, &height);
      if(height >= threshRatio * maxHeight) {
        peakPositions.push_back(pos);
        peakHeights.push_back(height); }}}
  peakPositions = peakPositions * sampleRate / T(N2);
  int numModes = (int) peakHeights.size();

  // Extract the maxNumModes loudest modes. The peakPositions/Heights arrays are parallel arrays 
  // sorted by ascending frequency. We need a sorting by descending height/amplitude, so we create
  // a helper array of 2D vectors with tze same data where the x-coordinate stores the height and
  // y-coordinate the freq (because sorting on 2D vectors using the < operator copares based on x 
  // first:
  //std::vector<rsVector2D<T>> peaks(numModes); // x: height, y: freq,  todo: use member
  peaks.resize(numModes); // x: height, y: freq
  for(int m = 0; m < numModes; m++)
    peaks[m] = rsVector2D<T>(peakHeights[m], peakPositions[m]);
  rsHeapSort(&peaks[0], numModes);
  AT::reverse(&peaks[0], numModes);
  // Maybe we later should use a 3D vector and store as 3rd element the index m because we may need
  // it later to figure out the adjacent/neighbor modes (for tuning the bi-notch) where the 
  // adjacency is again defined on the freq-axis

 
  /*
  // Plot results of the pre-analysis (for development):
  GNUPlotter plt;
  std::vector<T> freqs(N2);
  ft.binFrequencies(&freqs[0], N2, sampleRate);
  plt.addDataArrays(N2, &freqs[0], &buf1[0]);
  plt.addDataArrays(N2, &freqs[0], &buf2[0]);
  plt.addDataArrays((int)peakHeights.size(), &peakPositions[0], &peakHeights[0]);
  plt.setGraphStyles("lines", "lines", "points");
  plt.setPixelSize(1200, 400);
  plt.plot();
  */
  // OK - this looks good. Maybe later we could use rsPeakPicker which includes masking plus some
  // more sophisticated ideas. For the time being, the masking works well

  //...............................................................................................
  // Step 2:
  // Analyze each mode one at a time by bandpassing the signal with a bandpass tuned to the 
  // respective modal frequency and then using an envelope follower on the bandpassed signal. We 
  // re-use buf1 for storing the filtered signal, buf2 as ringout-buffer for the filtering process
  // and then to store the envelope...tbc...

  using ModalParams = rsModalFilterParameters<T>;
  numModes = rsMin(numModes, maxNumModes);
  std::vector<ModalParams> mp(numModes);  // maybe avoid the alloc - let the use pass a reference
  T decayTimeScale =  T(1) / (sampleRate * rsLog(decayAmp1/decayAmp2));
  for(int m = 0; m < numModes; m++)       // to such a vector and just resize it here, thenn fill it
  {
    T f = peaks[m].y;
    T a = peaks[m].x;

    mp[m].freq = f;    // preliminary
    //mp[m].amp  = a;    // preliminary

    // ToDo: 
    // -set up a bandpass tuned to the (preliminary) mode frequency f
    // -filter the signal x bidirectionally using the bandpass
    // -apply envelope-follower
    // -find out peak-amplitude and position - use these values as estimates for mode-amplitude
    //  and attack
    // -refine the frequency estimate by measuring a distance between a given number of 
    //  zero-crossings (maybe start at the peak and use a number of reliable cycles - maybe until
    //  the amplitude is down to 0.25 or soemthing - experiemntation necessary

    T bw = 2 * maskWidth;  // bandwidth for bandpass
    // preliminary - ToDo: use something based on measuring the width of the actual peak in the 
    // FFT magnitudes

    extractMode(x, &buf1[0], N, f, bw); 
    // ToDo: pass buf2 for workspace, perhaps use not the full length but rather just enough
    // to let the filter ring out properly.

    // Exctract the envelope of the mode - result goes to buf2:
    extractModeEnvelope(&buf1[0], &buf2[0], N);
    
    //rsPlotArrays(10000, &buf1[0], &buf2[0]);
    // Some of the extracted modes show an amplitude modulation that's not supposed to be there. 
    // That must be caused by some nearby mode leaking into the signal. We could try to remedy this
    // by using more aggressive filtering - ideas:  
    // -Use a smaller bandwidth and/or multiple passes
    // -Notch out neighboring modes
    // -Maybe the modulation itself can be detected, modeled and be translated to an influence
    //  of a nearby mode. Maybe we can model the beating of two modes via:
    //    x(t)  = (A + (m/A) * sin(wm*t + pm)) * sin(wc*t + pc)
    //         ?= A1 * sin(w1*t + p1) + A2 * sin(w2 + p2)
    //  where A is the overall amplitude, m is a modulation amount. Multiply that out to see what
    //  superposition of sines corresponds to that signal. That should give us the conversion 
    //  formulas between the measured A, m, wm, pm, wc, pt and the desired A1, w1, p1, A2, w2, p2.
    //  This may not be so useful here because we may see more than one mode leaking into out 
    //  extracted mode. But maybe, it we see some sinusoidal modulation of a mode, we can conclude
    //  that this should be modeled by two modes and the above idea can be used to figure out the
    //  parameters of these two modes from the "single" (modulated/beating) mode.

    // Find position and value of maximum of envelope and compute the estimates for out attack 
    // time and mode amplitude:
    int nMax = AT::maxIndex(&buf2[0], N);
    PF::exactPeakPositionAndHeight(&buf2[0], N, nMax, 1, &pos, &height);
    mp[m].att = pos / sampleRate;
    mp[m].amp = height;


    // Find time instants, where the amp-env goes through two specified reference amplitudes for 
    // the first time after the peak. From these instants compute the decay parameter:
    int n = nMax;
    while(n < N) {
      if(buf2[n] <= height * decayAmp1)
        break;
      n++; }
    T n1 = n;       // todo: be more precise using linear interpolation
    while(n < N) {
      if(buf2[n] <= height * decayAmp2)
        break;
      n++; }
    T n2 = n;       // todo: be more precise using linear interpolation
    mp[m].dec = decayTimeScale * (n2 - n1);  // verify!
    // ToDo: 
    // Compute n1, n2 more accurately, i.e. with subsample precision using linear interpolation:
    // To compute n1, Look at buf2[n-1] and buf2[n]. buf2[n-1] should be > height*decayAmp1 and
    // buf2[n] <= height*decayAmp1. Somewhere in between n-1 and n is the actual location where
    // the nev goes through height*decayAmp1 ...

    // Refine frequency estimate by counting cycles. We again start at the peak of the envelope and
    // count cycles to the right:
    using ZCF = rsZeroCrossingFinder;
    n = nMax;
    int count =  0;  // number of zero-crossings counted so far
    int nL    = -1;  // index of leftmost zero-crossing - maybe re-use n1
    int nR    = -1;  // index of rightmost zero-crossing - maybe re-use n2
    while(n < N) {
      if(ZCF::isUpwardCrossing(&buf1[0], n)) {
        count++;
        nR = n;
        if(nL == -1)
          nL = n;    }

      // We stop counting when the signal amplitude has fallen below the lower decay measurement
      // threshold - this should be a good place to stop:
      if(buf2[n] <= height * decayAmp2)
        break;
      n++;  }
    n = nMax-1;
    while(n >= 0)  {    // We also look for zeros leftward from the peak
      if(ZCF::isUpwardCrossing(&buf1[0], n)) {
        count++;
        nL = n;   }
      if(buf2[n] <= height * decayAmp2)
        break;
      n--;      }
    if(nL != -1 && nR != -1 && nR > nL)
    {
      int precision = 2;
      T tL = nL + ZCF::upwardCrossingFrac(&buf1[0], N, nL, precision);
      T tR = nR + ZCF::upwardCrossingFrac(&buf1[0], N, nR, precision);
      T p = (tR - tL) / T(count-1);  // period in samples
      p /= sampleRate;               // period in seconds
      f  = 1/p;                      // frequency in Hz
      mp[m].freq = f;   
      // OK, we get values in the right ballpark but they may not yet as precise as possible. The
      // improvement of the accuracy is a little bit disappointing. Or is that because the 
      // preliminary accuracy is already so good?

      // ToDo:
      // -Maybe try a higher precision in upwardCrossingFrac - doesn't seem to help
      // -Maybe try to use more or less cycles...more should be better, I think - but try it
      // -Try to use a more aggressive bandpass to isolate the mode
      // -Simplify the formula - use one division
    }
    else
    {
      rsWarning("Frequency estimation refinement failed");
      // This should really not happen for any halfway sane input signal
    }




    // ToDo:
    // -To estimate the decay, we can then perhaps just compare average amplitudes of some section
    //  of a length L after the peak and another section of the same length immediately right to 
    //  the first section...but yeah - for that, an actual envelope might actually be better for 
    //  computing accurate averages
    //  OR: find instants where env falls below height/2 and height/2 for the first time and use
    //  the difference as estimate for the time required to decay to 1/2. We don't start scanning
    //  from the peak because we actually want the asymptotic decay and the rounding around th peak
    //  may mess with our measurement.
    // -To ge a refined frequency estimate, find the next positive zero-crossing after the peak, 
    //  then find more zero crossings, counting them, going up to a pointwh where the envelope has
    //  decayed to 1/2 or 1/4 or something and then use the distance and count to compute a more 
    //  exact frequency estimate. Use rsZeroCrossingFinder
    // -To estimate the start phase, use the time-instant of the left zeor-crossing used above
    //  and extrapolate the phase back to time 0, using the refined frequency estimate. Maybe let 
    //  the user decide if the phase should be matched exactly at the zero-crossing near the peak
    //  or near the start. Both ways could make sense: match near the peak: better time-domain 
    //  subtraction results with original, match near start: may model the transient better.












    int dummy = 0;
  }


  
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
  // -Maybe before returning the mp array, it should be sorted by frequency...not sure about that
  //  though. Sorting by decreasing amplitude can also make sense

  // Maybe use:
  // rsExponentialEnvelopeMatcher, rsEnvelopeExtractor, rsPeakFinder::connectPeaks 
  // (MiscUnfinished.h), rsCurveFitter::fitExponentialSum (not yet suitable as is)


  return mp;
}

template<class T>
void rsModalAnalyzer2<T>::extractMode(const T* x, T* y, int N, T centerFreqHz, T bandwidthHz,
  T* wrk = nullptr, int N_wrk = 0)
{
  // This code is a bit ugly - re-work the API for biquad filters to make it more convenient:
  //using BD = rsBiquadDesigner<T>;
  //rsBiquadDF1<T> filter;

  // Set up the filter:
  using SVF = rsStateVariableFilter<T, T>;  // todo: maybe use a simpler direct-form biquad
  T bwOct = rsBandwidthConverter::absoluteBandwidthToOctaves(bandwidthHz, centerFreqHz);
  SVF filter;
  filter.setSampleRate(sampleRate);
  filter.setMode(SVF::modes::BANDPASS_PEAK);
  filter.setFrequency(centerFreqHz);
  filter.setBandwidth(bwOct);

  // Forward pass:
  for(int n = 0; n < N; n++)
    y[n] = filter.getSample(x[n]);

  // Let the filter ring out and warm up if the user has passed a buffer for that purpose:
  if(wrk != nullptr)
  {
    // This code needs to be tested...
    for(int n = 0; n < N_wrk; n++)
      wrk[n] = filter.getSample(T(0));
    for(int n = N_wrk-1; n >= 0; n--)
      wrk[n] = filter.getSample(wrk[n]);
  }
  else
    filter.reset();
    // I'm not so sure, if it's better to reset or not reset in case of no ringout/warm-up. Maybe
    // do some experiments.

  // Backward pass:
  for(int n = N-1; n >= 0; n--)
    y[n] = filter.getSample(y[n]);

  // ToDo:
  // -Maybe allow multiple passes of the filter. But then we need to scale (widen) the bandwidth
  //  for each single pass. But doing that right would also require some sort of pre-padding. But
  //  maybe that can be avoided by just using a single pass with a higher-order filter. Maybe
  //  EngineersFilter with the Gaussian response type could be the most suitable choice.
}

template<class T>
void rsModalAnalyzer2<T>::extractModeEnvelope(const T* x, T* y, int N)
{
  for(int n = 0; n < N; n++)
    y[n] = rsAbs(x[n]);
  rsPeakFinder<T>::connectPeaks(y, N, y, false); // try "true" for "useParabola"
}
