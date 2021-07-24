/** Returns the smallest multiple of N that is greater or equal to x. */
int rsCeilMultiple(int x, int N)
{
  rsAssert(x >= 0 && N > 0); 
  return ((x+N-1)/N)*N;
}
// helper function, needs tests, move to library



template<class T, int N>
void rsSineSweeperBankIterative<T, N>::setMaxNumOscillators(int newLimit)
{
  int numSimdGroups = rsCeilMultiple(newLimit, N) / N;
  // ToDo: rsCeilMultiple multiplies by N and here we divide by N again - that can be optimized.


  simdGroups.reserve(numSimdGroups);
}

template<class T, int N>
void rsSineSweeperBankIterative<T, N>::setNumActiveGroups(int newNumber)
{
  rsAssert(newNumber <= (int) simdGroups.capacity());
  simdGroups.resize(newNumber);
}

template<class T, int N>
void rsSineSweeperBankIterative<T, N>::processFrame(float* left, float* right)
{
  rsSimdVector<float, N> v(0.f);
  for(size_t i = 0; i < simdGroups.size(); i++)
    v += simdGroups[i].getSine();
  *left = *right = v.sum();
}
//
// https://rust-lang.github.io/packed_simd/perf-guide/vert-hor-ops.html

template<class T, int N>
void rsSineSweeperBankIterative<T, N>::reset()
{

}

//=================================================================================================

/*
template<int N>
void rsAdditiveSynthVoice<N>::EditablePatch::convertUserToAlgoUnits(
  float sampleRate, bool logAmplitude)
{
  float freqToOmega = 2*PI/sampleRate;
  for(size_t i = 0; i < breakpoints.size(); i++)
  {
    Breakpoint* bp = &breakpoints[i];
    bp->time = round(sampleRate * bp->time);
    for(size_t j = 0; j < bp->params.size(); j++)
    {
      SineParams* sp = &(bp->params[j]);
      sp->freq *= freqToOmega;
      sp->phase = RAPT::rsDegreeToRadiant(sp->phase);
      if(logAmplitude)
        sp->gain = log(sp->gain);
    }
  }
}
*/
// superfluous - directly use PlayablePatch::setupFrom

template<int N>
void rsAdditiveSynthVoice<N>::EditablePatch::createArtificialPhases()
{ 
  for(int i = 1; i < getNumBreakpoints(); i++)
  {
    Breakpoint* bpL = getBreakpoint(i-1);
    Breakpoint* bpR = getBreakpoint(i);
    double dt = bpR->time - bpL->time;
    for(int j = 0; j < bpL->getNumPartials(); j++)
    {
      SineParams* spL = &(bpL->params[j]);
      SineParams* spR = &(bpR->params[j]);
      double fM  = 0.5 * (spL->freq + spR->freq);        // mean frequency
      double wM  = fM * 360;                             // omega in degrees per second
      double pT  = spL->phase + wM * dt;                 // target phase
      spR->phase = rsWrapToInterval(pT, -180.0, +180.0); // wrapped into -180..+180
      int dummy = 0;
    }
  }
}
// needs tests...

// todo: implement a function that generates amplitude derivatives to make the amp-env smoother


template<int N>
bool rsAdditiveSynthVoice<N>::EditablePatch::isWellFormed() const
{
  size_t numBreakpoints = breakpoints.size();
  if(numBreakpoints == 0)        return true;  // Patch is empty. This is allowed.
  if(numBreakpoints == 1)        return false; // One breakpoint does not define a segment.
  if(breakpoints[0].time != 0.0) return false; // Time must start at zero.

  // Check, if time-stamps are strictly increasing with a time-delta of at least 1/20000 seconds.
  // Rationale: We assume that we may have one breakpoint per cycle and 20 kHz is more than high 
  // enough for a fundamental frequency. In fact, it's already unrealistically high, but
  // from a technical point of view, we must just make sure that the spacing between breakpoints is
  // at least one sample (right?) and we assume that the engine runs at least at 44.1 kHz 
  // sample-rate.
  double minTimeDelta = 1.0 / 20000.0;
  for(int i = 0; i < numBreakpoints-1; i++)
    if(breakpoints[i+1].time - breakpoints[i].time < minTimeDelta)
      return false;

  // Check, if all breakpoints have the same number of partials. This is required because...tbc...
  int numPartials0 = breakpoints[0].getNumPartials();
  for(int i = 1; i < numBreakpoints; i++)
    if(breakpoints[i].getNumPartials() != numPartials0)
      return false;

  return true;

  // ToDo: 
  // -Try to restrict the API of EditablePatch such that client code can't even produce malformed 
  //  patches. Then, the utility of this function should be merely for internal sanity checks for 
  //  debugging and should probably not be compiled into release versions.
}

template<int N>
int rsAdditiveSynthVoice<N>::EditablePatch::getNumPartials() const
{
  int num = 0;
  for(size_t i = 0; i < breakpoints.size(); i++)
    num = rsMax(num, breakpoints[i].getNumPartials());
  return num;

  // ToDo:
  // -If we can assume that all breakpoints have the same number of partials, then it's sufficient
  //  to just return breakpoints[0].getNumPartials() ...or zero if the size of breakpoints is zero.
}




template<int N>
void rsAdditiveSynthVoice<N>::PlayablePatch::setupFrom(
  const rsAdditiveSynthVoice<N>::EditablePatch& patch, double sampleRate, bool applyLog)
{
  rsAssert(patch.isWellFormed());

  //bool usePhase = true;  
  // make parameter...or maybe remove - artificial phase should be created by the EditablePatch

  // Memory (re-)allocation:
  int numPartials = patch.getNumPartials();
  numSimdGroups   = rsCeilMultiple(numPartials, N) / N;
  numBreakpoints  = patch.getNumBreakpoints();
  params.setShape(numBreakpoints, numSimdGroups);
  //params.memsetAllZero();  // superfluous
  timeStamps.resize(numBreakpoints);

  // Conversion factors:
  double freqToOmega = 2.0*PI/sampleRate;
  double degToRad    = PI/180.0;

  // Shorthands:
  using Breakpoint = rsAdditiveSynthVoice<N>::EditablePatch::Breakpoint;
  using SineParams = rsAdditiveSynthVoice<N>::EditablePatch::SineParams;

  // Temporaries:
  double tL, tR, nL, nR, pL, pR, wL, wR, aL, aR, rL, rR; // params for one partial
  rsSweepParameters<rsSimdVector<float, N>> p;           // params for one simd group
  memset(&p, 0, sizeof(p));

  // Loop over the breakpoints:
  for(int breakpointIndex = 0; breakpointIndex < numBreakpoints-1; breakpointIndex++)
  {
    int i = breakpointIndex;
    const Breakpoint* bpL = patch.getBreakpoint(i);
    const Breakpoint* bpR = patch.getBreakpoint(i+1);

    tL = sampleRate * bpL->time;
    tR = sampleRate * bpR->time;
    nL = round(tL);
    nR = round(tR);
    timeStamps[i] = (int)nL;
    double dt = tR - tL;
    double dn = nR - nL;
    // Or should it be nR-nL? That may also have to depend on how we handle phase correction etc.

    // Loop over the partials in current breakpoint:
    int partialIndex = 0; // index of current partial within the patch
    int groupIndex   = 0; // index of current simd-group (to which the current partial belongs)
    int indexInGroup = 0; // index of current partial within the current simd group
    //while(partialIndex < numPartials)
    while(partialIndex < N*numSimdGroups)
    {
      int j = partialIndex;  // shorthand


      if(j < numPartials)
      {
        const SineParams* spL = &(bpL->params[j]);
        const SineParams* spR = &(bpR->params[j]);

        // Convert parameters from user- to algo-units
        pL = degToRad * spL->phase;
        pR = degToRad * spR->phase;
        wL = freqToOmega * spL->freq;
        wR = freqToOmega * spR->freq;

        // unwrap phase (factor out):
        double wm  = 0.5 * (wL + wR);     // mean omega during segment

        double tmp = pL + (nR-nL-1)*wm;   // unwrapped end phase should be near this value
                                          // is the -1 correct?
        pR = rsConsistentUnwrappedValue(tmp, pR, 0.0, 2.0*PI);

        // nope:
        //if(usePhase)
        //  pR = rsConsistentUnwrappedValue(tmp, pR, 0.0, 2.0*PI);
        //else
        //{
        //  pR = tmp;
        //  // i think, this is still wrong - we also need to use the old pR for pL
        //}


        aL = spL->gain;  // maybe use gL for "gain"
        aR = spR->gain;
        if(applyLog) {
          aL = log(aL);
          aR = log(aR);
        }
        rL = (aR - aL) / dt; // verify this! i think, we need dn?
        //rL = (aR - aL) / dn; // test ..doesn't seem to make a difference
        rR = rL;             // we currently only use a linear (log-)amp env with constant slope
      }
      else
      {
        pL = pR = wL = wR = 0;  // phase and omega is 0
        aL = aR = rL = rR = 0;  // log-amp and rise is 0
        nL = 0;                 // left sample time is 0
        nR = 1;                 // right sample time is taken to be 1 to avoid nans
      }

      // Write the converted parameters into the appropriate positions in the simd vectors:
      j = indexInGroup;
      p.t0[j] = float(0 ); p.t1[j] = float(nR-nL);
      p.p0[j] = float(pL); p.p1[j] = float(pR);
      p.w0[j] = float(wL); p.w1[j] = float(wR);
      p.l0[j] = float(aL); p.l1[j] = float(aR);
      p.r0[j] = float(rL); p.r1[j] = float(rR);

      // Do loop variable increments and and handle assignment of simd-vector and wraparound, if 
      // we have reahced the end of a simd group:
      partialIndex++;
      indexInGroup++;
      if(indexInGroup == N) 
      {
        params(breakpointIndex, groupIndex) = p;
        groupIndex++;
        indexInGroup = 0;  
      }
    }


  }

  // ToDo:
  // -For those partials that do not have valid data, we need to use some dummy values other than 
  //  zero because zero values will give rise to nans... i think, maybe it's p.t1 variable that is to 
  //  blame - maybe just setting that to 1 could fix it?
  // -There's some redundancy in the computations: we may use the old converted pR,wR,aR,rR as
  //  pL,wL,aL,rL in the next iteration...oh - but no:...that would work only, if the inner loop
  //  would be over the breakpoints and the outer over the partials...but that may have a less 
  //  desirable access pattern and introduce other complications
  // -We may account for the rounding of the time stamps by adjusting the phase and amplitude 
  //  values a bit: take the difference d between exact and rounded value and advance phases by
  //  d*w and amplitude by d*r. ..can we also account for changes in w and r themselves? that would
  //  require an estimate of the time derivative of w and r which can be computed by (wR-wL)/dt
  //  and (rR-rL)/dt...if we do this, we should correct w,r before using them to correct p,a 
  //  (i think)...or maybe use the average of original and corrected values - that's and estimate
  //  of w,r at dt/2...i think...
}

template<int N>
void rsAdditiveSynthVoice<N>::startPlaying()
{
  if(sweeperBank == nullptr || patch == nullptr)
    return;
  sweeperBank->setNumActiveGroups(patch->getNumSimdGroups());
  playing = true;
  handleBreakpoint(0, true);
}

template<int N>
void rsAdditiveSynthVoice<N>::handleBreakpoint(int i, bool reInitAmpAndPhase)
{
  if(patch->getNumBreakpoints() > i+1)
  {
    // When there is another breakpoint after the breakpoint with index i, we set the sweepers up
    // to play the segment between breakpoints i and i+1:
    initSweepers(i, i+1, reInitAmpAndPhase);
    nextBreakpointIndex = i+1;
    samplesToNextBreakpoint = patch->getBreakpointTime(i+1) - patch->getBreakpointTime(i);
  }
  else
  {
    // When there is no more breakpoint after the breakpoint with index i, i.e. when we have now
    // reached the final breakpoint, we stop the playback (and do some cleanup):
    playing = false;
    nextBreakpointIndex = 0;     // i think, it's not really relevant, but it's just cleaner
    samplesToNextBreakpoint = 0; // ditto
  }
}

template<int N>
void rsAdditiveSynthVoice<N>::processFrame(float* left, float *right)
{
  if(samplesToNextBreakpoint == 0)
    handleBreakpoint(nextBreakpointIndex, alwaysReInit);
  if(!playing)
    return;
  sweeperBank->processFrame(left, right);
  samplesToNextBreakpoint--;
}

template<int N>
void rsAdditiveSynthVoice<N>::reset()
{
  sweeperBank->reset();
  nextBreakpointIndex = 0;
  samplesToNextBreakpoint = 0;
  playing = false;
}

template<int N>
void rsAdditiveSynthVoice<N>::initSweepers(int i0, int i1, bool reInitAmpAndPhase)
{
  rsAssert(i0 >= 0 && i1 < patch->getNumBreakpoints());
  rsAssert(i1 == i0+1); // maybe lift that restriction later to allow loops

  using SimdVector  = rsSimdVector<float, N>;
  using GroupParams = RAPT::rsSweepParameters<SimdVector>;


  if(reInitAmpAndPhase)
  {
    for(int j = 0; j < patch->getNumSimdGroups(); j++)
      sweeperBank->setup(j, patch->getParams(i0, j));
  }
  else
  {
    // We need to compute one sample more to put the bank into the desired state to retrieve the 
    // instantaneous phase and amplitude (or do we? figure out by using very dense datapoints):
    float fDummy; sweeperBank->processFrame(&fDummy, &fDummy);

    for(int j = 0; j < patch->getNumSimdGroups(); j++)
    {
      GroupParams p = patch->getParams(i0, j);

      SimdVector p0, l0, dl, dt, r0;
      p0 = sweeperBank->getPhase(j);
      l0 = rsLog(sweeperBank->getAmplitude(j));

      // for test/debug:
      SimdVector phaseErr, lgAmpErr;
      phaseErr = p0 - p.p0;
      lgAmpErr = l0 - p.l0;

      // Modify p.p0, p.l0, p.r0, (p.r1):
      p.p0 = p0;
      p.l0 = l0;
      dl   = p.l1 - p.l0;
      dt   = p.t1 - p.t0;
      r0   = dl / dt;
      p.r0 = r0;





      sweeperBank->setup(j, p);
    }
  }
}


/*

Ideas:

-The user should be able to set up the engine in terms of breakpoints where each breakpoint 
 contains a set of rsInstantaneousSineParams for each partial. Note that this enforces to have the
 breakpoints for each partial in sync. The rsSinusoidalModel is actually more flexible in this 
 regard (but the rsHarmnonicAnalyzer actually delivers models that don't make use of that 
 additional flexibility). But from a performance point of view, we really want the breakpoints of 
 all partials to be synced up, such that we can handle breakpoint arrival once and for all and not 
 for each partial seperately (which would thwart the simd parallelization scheme). But maybe not
 *all* partials should have their breakpoint arrival at the same instant. In order to spread the 
 computational load in time, it may make sense to have them in groups, e.g. partials 0..127 have
 their breakpoints at 0,1024,2048,3072, partials 128,..,255 at 0,256,1280,2304,3328, etc. At n=0
 it's unavaoidable to initialize the states for all partials, but the periodic re-inits can be 
 interleaved in groups. Such a scheme would require to pre-process the model to obtain data for
 these (artificial) breakpoint times.

-The API should allow to take an rsSinusoidalModel as input. If the model happens to have 
 asynchronous breakpoints for the partials, we may need to convert it here to obtain synchronized
 breakpoints. To read off rsInstantaneousSineParams from the model at arbitrary sample instants,
 we should interpolate the datapoints in the same way in which the rsSinusoidalSynthesizer engine
 would do it.

-We may also want to be able to take different models for different keys. In the extreme case one 
 for each key, in practice perhaps more something like one for each octave with parameter
 interpolation in between the keys. ...although, that interpolation may be problematic because
 the parameters are actually too low-level. To provide meaningful interpolation would actually need
 a more high-level description, like, instead of instantaneous amplitudes at each breakpoint, 
 something like parameters for an overall amp-envelope and maybe amp-LFO (tremolo). Maybe it makes 
 sense to define the set of high-level parameters first and then measure them from real world 
 sample data.

-Maybe split the partials conceptually into groups and let the user access certain processing 
 options per group. For example, it may make sense to have one group for even and one for odd 
 partials and let the user set their gains individually. It may also make sense to split into 
 low/mid/high frequencies or groups containing harmonic and inharmonic content.
-Maybe the sizes of the groups should be independent from the SIMD vector size, and maybe even
 variable: e.g. group 1 has 10 harmonics, group 2 has 17, group 3 has 42, etc. Maybe we can take
 inspiration from how groups are used in sfz.

-When we arrive at a new breakpoint, due to error accumulation, the instantenous amplitude and 
 phase may be slightly different from the values stored values at the breakpoint (which we were 
 aiming for when we initialized the iterators at the previous breakpoint). Using the actually 
 stored values for the re-initialization may therefore lead to (small) discontinuities. In order 
 to avoid them, we should not use the stored values (e.g. for p0,w0,l0,r0 in rsSineSweepIterator) 
 but those values, where the oscillator actually ended up with. This requires formulas to retrieve 
 instantaneous log-amp, log-amp derivative, phase and omega from the oscillator. Maybe we csn just
 read the phase and mangitude that the osc currently has (stored in y[3]), then call getValue once 
 more, look at those new/updated phases and magnitudes, too and estimate the derivatives by a 
 finite difference. Maybe we can find a more exact formula but if not, this approach shouldn't be 
 too bad. Maybe to find a formula for the derivatives that the osc ended up with, we can somehow 
 reconstruct the cubic polynomial coeffs for phase and log-amplitude from the state variables 
 y[0],..,y[3] and differentiate the polynomial?

-Numerical experiments indicate that it may be problematic to use cubic envelopes for the 
 log-amplitude, i.e. envelopes where r0 != r1 in rsSineSweepIterator. Maybe use linear ones, i.e. 
 ones where where r0 == r1 and the value is given by (l1-l0)/(N-1) where N is the length of the 
 segment, i.e. the number of samples until the next breakpoint will be reached. This will restrict
 the envelope shapes to exponential segments and the amp-enve may become less smooth, but the algo
 may becoem numerically more stabe -> experiments are needed. As said above, as l0 we don't use the
 log of the stored breakpont amplitude but the instantaneous amplitude that the osc ended with 
 (or, well, actually, that of one sample later).

-Maybe optionally apply waveshaping to each partial (maybe with Chebychev polynomials). That could
 be a cheap way add some dynamic richness to teh output. ..or maybe the waveshaping should be done
 per group. That would also create soem intermodulation components...maybe they are useful to 
 "glue" the partials together

-Maybe actually do not prescribe target values for instantaneous phase and instantaneous 
 apmlitude raise at the end. Just take whatever comes out. Maybe then we can get away with a 
 degree 2 polynomial for the iterators. we just prescribe instantaneous phase, freq, amp, raise at 
 the startpoint and at the endpoint we only prescribe freq and amp. The target value for the 
 phase at the startpoint is just the final osc-phase is and the amp-raise is determined by the 
 difference of log-amps of the target amp at the next breakpoint and the final osc-amp.

Notes:

Maybe the patch data should not contain a data field for the instantaneous phase because that makes
it difficult to transpose a patch to a different fundamental frequencies. Muliplying all 
instantaneous frequency datapoints by a factor invalidates the instantaneous phase data. We could
compute new target phases by a heuristic algorithm, but that algorithm can just as well run on the
fly.

Maybe what the suer sees as patch should here actually be a "PatchSet" because we my want to use
different patches for different keys

Currently, the largest SIMD vector size available for single precision float in standard 
hardware is 16 (for example in AVX-512). Maybe twice that value just in case that something
like AVX-1024 will become available in the future, in which case the code may directly benefit
from it without any change.
...hmm...not sure, if that's a good idea...i think, from a perceptual point of view, smaller 
groups may make more sense...but then, when we provide an API for setting up settings for
groups, their size may actually have to be different from the simd vector size anyway

// see also:
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=553578&p=7911364

*/