template<class T, int N>
void rsSineSweeperBankIterative<T, N>::setMaxNumOscillators(int newLimit)
{

}

template<class T, int N>
void rsSineSweeperBankIterative<T, N>::setNumOscillators(int newLimit)
{

}

template<class T, int N>
void rsSineSweeperBankIterative<T, N>::init(int index, float t0, float t1, float w0, 
  float w1, float a0, float a1, float p0, float p1, float r0, float r1)
{

}

template<class T, int N>
void rsSineSweeperBankIterative<T, N>::processFrame(float* left, float* right)
{

}

template<class T, int N>
void rsSineSweeperBankIterative<T, N>::reset()
{

}

//=================================================================================================

void rsAdditiveKeyPatch::convertUserToAlgoUnits(float sampleRate, bool logAmplitude)
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

//=================================================================================================

template<int N>
void rsAdditiveSynthVoice<N>::startPlaying()
{
  goToBreakpoint(0, true);
  playing = true;
}

template<int N>
void rsAdditiveSynthVoice<N>::goToBreakpoint(int index, bool reInitAmpAndPhase)
{
  if(patch->getNumBreakpoints() > index+1)
  {
    const rsAdditiveKeyPatch::Breakpoint* start = patch->getBreakpoint(index);
    const rsAdditiveKeyPatch::Breakpoint* end   = patch->getBreakpoint(index+1);
    initSweepers(start, end, reInitAmpAndPhase);
    nextBreakpointIndex = index+1;
    samplesToNextBreakpoint = (int) end->time - (int) start->time;
  }
  else
  {
    // When there is no breakpoint after bpIndex, we have reached the end and stop the playback:
    nextBreakpointIndex = 0;
    samplesToNextBreakpoint = 0;
    playing = false;
  }
}

template<int N>
void rsAdditiveSynthVoice<N>::processFrame(float* left, float *right)
{
  if(samplesToNextBreakpoint == 0)
    goToBreakpoint(nextBreakpointIndex, alwaysReInit);
  if(!playing)
    return;
  sweeperBank->processFrame(left, right);
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
void rsAdditiveSynthVoice<N>::initSweepers(const rsAdditiveKeyPatch::Breakpoint* bpStart,
  const rsAdditiveKeyPatch::Breakpoint* bpEnd, bool reInitAmpAndPhase)
{
  rsAssert(bpStart->params.size() == bpEnd->params.size());
  // number of partials must remain the same during the playback...maybe lift that restriction 
  // later

  // ToDo: We actually need the breakpoint-data in simd-goups, too...so yeah, i think, we really
  // need two vesrions of the patch class: one for editing and one for playback and a conversion 
  // routine

  /*
  //RAPT::rsSineSweepIterator::Parameters tmpParams;

  size_t numPartials = bpStart->params.size();  
  rsAdditiveKeyPatch::SineParams paramsStart, paramsEnd;
  for(size_t i = 0; i < numPartials; i++)
  {
    paramsStart = bpStart->params[i];
    paramsEnd   = bpEnd->params[i];

    //tmpParams.t0 = 0.f;
    //tmpPara


    int dummy = 0;

  }
  */


  int dummy = 0;
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

*/