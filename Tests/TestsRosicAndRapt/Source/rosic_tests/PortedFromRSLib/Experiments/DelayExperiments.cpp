#include "DelayExperiments.h"

void algoVerb()
{
  //int N = 5000;  // number of samples to generate

  int N = 88200;  // number of samples to generate

  using Real = double;
  using Vec  = std::vector<Real>;
  //using FDN  = rosic::FeedbackDelayNetwork;
  using FDN  = rosic::FeedbackDelayNetwork16;
  using DD   = FDN::delayDistributions;

  FDN* fdn = new FDN; 
  // Trying to allocate rosic::FeedbackDelayNetwork16 on the stack gives a stack overflow. 
  // Apparently, it's too data-heavy. ToDo: fix that by allocating the memory for the delaylines on
  // the heap within the class by using std::vector instead of raw arrays. Maybe rsMatrix, if 
  // needed for the 2D array. This will also allow to change the maxDelay at runtime, if necessary.
  // This may be needed when the sample rate changes at runtime.

  fdn->setAllpassMode(false);
  fdn->setFeedbackMatrix(FDN::HADAMARD);

  //fdn->assignRelativeDelayTimesAlgorithmically(DD::GEOMETRIC_MEANS, 1.0, GOLDEN_RATIO);
  //fdn->assignRelativeDelayTimesAlgorithmically(DD::GEOMETRIC_MEANS, 1.0, SQRT2);
  //fdn->assignRelativeDelayTimesAlgorithmically(DD::GEOMETRIC_MEANS, 1.0, 2.0);
  //fdn->assignRelativeDelayTimesAlgorithmically(DD::LINEAR, 1.0, 2.0); // impulse-train!
  //fdn->assignRelativeDelayTimesAlgorithmically(DD::LINEAR, 1.0, GOLDEN_RATIO);
  //fdn->assignRelativeDelayTimesAlgorithmically(DD::DISTANCE_DECAY, 1.0, 0.867231);
  //fdn->assignRelativeDelayTimesAlgorithmically(FDN::delayDistributions::PRIME_ALGO_1, 0.5, 0.5);
  // seems to make no difference - ah - the code is commented out
  // maybe try something base on the golden ratio - or maybe use rsRatioGenerator




  Vec t(N), hL(N), hR(N);
  RAPT::rsArrayTools::fillWithIndex(&t[0], N);
  hL[0] = hR[0] = 1.0;
  for(int n = 0; n < N; n++)
    fdn->processFrame(&hL[n], &hR[n]);

  //rsPlotVectorsXY(t, hL, hR);

  rosic::writeToStereoWaveFile("ImpRespFDN.wav", &hL[0], &hR[0], N, 44100);

  delete fdn;

  // Observations:
  //
  // - With the current settings, it sounds rather metallic. I guess this is due to the settings of 
  //   the delayline lengths. In the testFeedbackDelayNetwork, we get a far better result also with
  //   16 delaylines. -> Figure that out!
  //
  // - After some nice delayline setting has been figured out, implement a diffusor and put it in 
  //   front. Maybe that diffusor should be optional so save CPU.
  //
  //
  // ToDo:
  //
  // - Set up an APE project, where we can manually enter the relative delay times
  //
  // - Start with 2 delaylines and tweak the 2nd delay-time until it sounds least tonal
  //
  // - Then add in a 3rd and tweak its delaytime also until it sounds least tonal
  //   and so on
}

void basicIntegerDelayLine()
{
  static const int N = 20;
  double t[N], h[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  rsBasicDelayLineD dl;
  dl.setMaximumDelayInSamples(5);
  dl.setDelayInSamples(5);
  RAPT::getImpulseResponse(dl, h, N);
  plotData(N, t, h);

  // Observations:
  //
  // - The impulse response is an impulse shifted by 5 samples, i.e. centered at sample index 5
  //   rather than 0. That's the expected result.
}

void twoPoleAllpassDelay()
{
  // We experiment with a chain of rsTwoPoleAllpassDelay filters to see what sort of impulse 
  // responses we can achieve with it.

  using Real = double;
  using VecI = std::vector<int>;
  using VecR = std::vector<Real>;
  using APF  = rsTwoPoleAllpassDelay<Real, Real>;

  // Helper function to generate the output of a chain of rsTwoPoleAllpassDelay filters with the
  // given delays, normalized radian frequencies and Qs:
  auto create = [](const VecI& delays, const VecR& omegas, const VecR& Qs, int N)
  {
    size_t numStages = delays.size();
    rsAssert(omegas.size() == numStages);
    rsAssert(Qs.size()     == numStages);

    // Create and set up the filters:
    using APF = rsTwoPoleAllpassDelay<Real, Real>;
    std::vector<APF> filters(numStages);
    rsStateVariableFilter<Real, Real> svf;
    Real b0, b1, b2, a1, a2;
    for(size_t i = 0; i < numStages; i++)
    {
      svf.setupAllpass(omegas[i], Qs[i]);
      svf.convertToBiquad(&b0, &b1, &b2, &a1, &a2);

      filters[i].setMaxDelayInSamples(delays[i]);
      filters[i].setDelayInSamples(   delays[i]);
      filters[i].setAllpassCoeffs(a1, a2);
    }

    // Create helper function to apply the filters:
    auto applyFilters = [&](Real x)
    {
      for(size_t i = 0; i < numStages; i++)
        x = filters[i].getSample(x);
      return x;
    };

    // Generate impulse response:
    std::vector<Real> h(N);
    h[0] = applyFilters(1.0);
    for(size_t n = 1; n < N; n++)
      h[n] = applyFilters(0.0);

    return h;

    // ToDo:
    //
    // - It's a bit silly to use the state variable filter for the purpose of designing biquad 
    //   coeffs so maybe replace this code later with the direct RBJ biquad design formulas.
  };
  // Maybe this should become a library function someday. 


  // Helper function to plot impulse response of filter with given settings:
  auto plot = [&](const VecI& delays, const VecR& omegas, const VecR& Qs, int N)
  {
    std::vector<Real> h = create(delays, omegas, Qs, N);
    rsPlotVectors(h);
  };

  // Helper function to write impulse response of filter with given settings to a wave file
  auto write = [&](const VecI& delays, const VecR& omegas, const VecR& Qs, int N, 
                   const std::string& path)
  {
    std::vector<Real> h = create(delays, omegas, Qs, N);
    rosic::writeToMonoWaveFile(path, &h[0], N, 44100);
  };


  // Show some plots:
  plot({ 1   }, 
       { 0.2 }, 
       { 2.5 }, 100);

  plot({ 5   }, 
       { 0.2 }, 
       { 2.5 }, 500);

  plot({ 5,   9   }, 
       { 0.2, 0.3 }, 
       { 2.5, 2.5 }, 500);

  plot({ 5,   9,   14  }, 
       { 0.2, 0.3, 0.5 }, 
       { 2.5, 9.5, 2.5 }, 1000);


  // Render some wave files:
  write({ 7   }, 
        { 0.2 }, 
        { 25  }, 4096, "TwoPoleAllpass_7_0.2_25.wav");

  write({ 7,   5   }, 
        { 0.2, 0.3 }, 
        { 25,  25  }, 16384, "TwoPoleAllpass_7_0.2_25__5_0.3_25.wav");


  // Observations:
  //
  // - For a single stage, the output looks like an initial spike followed by an undulating spike 
  //   train. The distance between the spikes is given by the chosen delay. The undulation 
  //   frequency is probably the filter's resonance frequency - maybe divided by the delay-length 
  //   factor.
  //
  //
  // Questions:
  //
  // - I think, we do not really need the delay lengths to ba all prime numbers. It should be 
  //   sufficient if they are mutually coprime. This is a less restrictive condition but the effecr
  //   having the first coincidence of spikes at the lowes common multiple of the delay lengths 
  //   should still be satisfied.
  //
  // - Can we use the undulation frequencies to deliberately colorize the sound, i.e. give it some
  //   deliberate tonal character? That seems plausible. Maybe we could use slightly detuned
  //   omegas for left and right channel.
  //
  //
  // ToDo:
  //
  // - Try to achieve sign flipping of the output by fliiping the signs of the omegas. This doesn't
  //   fall out of the math. Trying to just use negative omegas (or negative Qs) just leads to 
  //   unstable filters. It would just be a convention that we apply ourselves manually to allow 
  //   the user to conveniently select the sign of a particular filter output. It would work by 
  //   using abs(w) for the filter design and sign(w) to scale the output.
  //
  // - Try to parametrize the filter not in terms of omega but in terms of physical frequency. We 
  //   really want to figure out the relationship between the omega parameter and the undulation 
  //   frequency. Maybe  omega*delay = 2*pi*f/fs  or  omega/delay = 2*pi*f/fs?
}

void dampedAllpassComb1()
{
  // This experiment is basically my initial step by step derivation of what later became the class
  // rsDampedAllpassCombNaive. It's based on the idea that I explain in the dcoument:
  // Notes/DSP/DampedAllpassComb.txt
  //
  // I implement an idea for starting with an arbitrary given allpass filter A(z) and arbitrary 
  // given feedback filter F(z) that sits in a feedback loop with unit delay around that allpass. I 
  // design a compensation filter that can be applied in series to this setup such that the overall
  // transfer function is allpass in nature. Without the compensation filter, this setup has the 
  // uncompensated transfer funcion:
  //
  //                   A(z)
  //  U(z) = ----------------------------
  //          1 + k * z^-1 * F(z) * A(z)  
  //
  // I have introduced also a feedback gain k for convenient parametrization although for 
  // derivations, we may want to just absorb that into F(z) to have one variable less to carry 
  // around. To completely cancel the effect of the feedback loop, we could use the following 
  // compensation filter:
  //
  //   C(z) = 1 + k * z^-1 * F(z) * A(z)
  //
  // In a first experiment, we try, if that indeed works. As A(z), we use simple delay of 
  // length M, i.e. A(z) = z^-M. 
  //
  // Later, we want to try to modify the compensation filter C(z) into one that compensates only
  // for the effect of the feedback on the magnitude response but retains features of the phase
  // response. Maybe we should try to reflect (some of) the zeros of C(z) about the unit circle. 
  // That will retain magnitude response and stability of C(z). 
  //
  // Consider our special case here with  A(z) = z^-M  and a one pole filter in the feedback loop,
  // i.e.  F(z) = (b0 + b1*z^-1) / (1 + a1*z^-1). Sustituting that into the formula for C(z) we 
  // get:
  //
  //               k * z^-1 * z^-M * (b0 + b1*z^-1)     1 + a1*d + b0*k*d^(M+1) + b1*k*d^(M+2)
  //   C(z) = 1 + ---------------------------------- = -----------------------------------------
  //                      1 + a1*z^-1                                 1 + a1*d
  //
  // With d := z^-1 for convenience. In a second step, we try to implement C(z) in this form using
  // a delayline to implement the denominator the one-pole filter 1 / (1 + a1*z^-1) for the 
  // denominator. Having the filter in this form makes it amenable for manipulating the zeros by 
  // reflecting them in the unit circle. We just need to reverse the array of numerator coeffs.
  //
  // Then, finally we actually do the reversal of the FIR part of the correction filter. Using that
  // instead of the original correction filter, we obtain an overall allpass filter but one that 
  // does not cancel out the effect of the feedback path.


  using Real    = double;
  using Vec     = std::vector<Real>;
  using Delay   = RAPT::rsBasicDelayLine<Real>;      // We use a simple delay as allpass
  using OnePole = RAPT::rsOnePoleFilter<Real, Real>; // We use a one pole as feedback filter

  int  M          =   100;     // Delay
  int  N          =  8192;     // Number of samples to generate
  Real sampleRate = 44100;     // Sample rate for writing the wavefiles
  Real dampFreq   =  1000;     // Frequency of the low shelf for feedback damping
  Real dampGain   =     0.7;   // Linear high freq damping gain
  Real k          =     0.99;  // Feedback gain factor


  // Create and set up the two given filters for A(z) and F(z):
  Delay apf;                                // Allpass filter
  apf.setMaximumDelayInSamples(M);
  apf.setDelayInSamples(M);

  OnePole fbf;                              // Feedback filter
  fbf.setMode(OnePole::modes::HIGHSHELV_BLT);
  fbf.setCutoff(dampFreq);
  fbf.setShelvingGain(dampGain);


  // Helper variable and function to implement the feedback loop filter:
  Real s = 0;                               // Output and state of uncompensated filter
  auto getSampleU = [&](Real x)
  {
    // This implements the feedback loop with unit delay without feedback filter:
    s = apf.getSample(x + k * fbf.getSample(s));
    return s;
    // The current value of s will used in the next call. This is the unit-delay feedback loop.
  };

  // Produce the uncompensated impulse response:
  Vec hu(N);
  hu[0] = getSampleU(1.0);
  for(int n = 1; n < N; n++)
    hu[n] = getSampleU(0.0);
  //rsPlotVectors(hu);                      // Decaying spike train with progressive tail smear


  // Cancel the effect of the feedback loop. For this, we re-use the existing filter objects. We 
  // can do this because they are not needed anymore for other purposes because the uncompensated 
  // output as been generated already. In a realtime implementation, we would have to use another
  // pair of filters that is identical to apf and fbf

  // Helper function that implements the compensation filter:
  rsUnitDelay<Real> ud;                     // A unit delay object for convenience.
  auto getSampleC = [&](Real x)
  {
    return x - k * ud.getSample(fbf.getSample(apf.getSample(x)));
    // C(z) = 1 + k * z^-1 * F(z) * A(z). The sign flip "1 + k * z^-1.." -> "x - k * .." is because
    // in a transfer function, the denominator signs are flipeed with respect to the difference 
    // equation. The order in which we apply k, ud, fbf, apf should not matter because everything 
    // is linear.
  };

  // This should produce a single spike at M like a delayline without feedback. That is, we compute
  // U(z) * C(z) = A(z):
  apf.reset();
  fbf.reset();
  Vec hc(N);
  for(int n = 0; n < N; n++)
    hc[n] = getSampleC(hu[n]);
  //rsPlotVectors(hc);                        // Yes - that looks good! Single spike at M.
  

  // Now we try to implement the correction filter in a different way that is more amenable to 
  // reflecting the zeros. Namely, in the form:
  //
  //           1 + a1*d + k*b0*d^(M+1) + k*b1*d^(M+2)
  //   C(z) = ----------------------------------------
  //                   1 + a1*d

  Real b0 = fbf.getB0();
  Real b1 = fbf.getB1();
  Real a1 = fbf.getA1();
  Delay   dl;                          // For the numerator
  OnePole op;                          // For the denominator
  op.setCoefficients(1.0, 0.0, a1);    // Should have same poles as fbf and 1 as denominator

  // Define the FIR coefficients to be used:
  Real c0  =  1;
  Real c1  = -a1;
  Real cM1 = -k*b0;
  Real cM2 = -k*b1;
  // I think the sign inversion for c1 has to do with rsOnePoleFilter using the other sign 
  // convention for denominator coeffs. For cM1 and cM2, I'm not sure, why it's needed. It may have
  // to do with the factors involving a k which comes from denominator but is now used in a 
  // numerator, but I'm really not sure about this. I have figured this out by trial and error.

  // Naive implementation using for each desired delay its own delayline:
  ud.reset();
  Delay dlM1;
  dlM1.setMaximumDelayInSamples(M+1);
  dlM1.setDelayInSamples(M+1);
  Delay dlM2;
  dlM2.setMaximumDelayInSamples(M+2);
  dlM2.setDelayInSamples(M+2);

  // Helper function to produce correction with alterntative correction filter:
  auto getSampleC2 = [&](Real x)
  {
    // Apply 1-pole:
    Real t = op.getSample(x);

    // Apply the FIR part. It has 4 terms:
    Real y = 0;
    y += c0  * t;
    y += c1  * ud.getSample(t);
    y += cM1 * dlM1.getSample(t);
    y += cM2 * dlM2.getSample(t);
    return y;
  };

  // Produce corrected output with alternative implementation:
  Vec hc2(N);
  for(int n = 0; n < N; n++)
    hc2[n] = getSampleC2(hu[n]);
  //rsPlotVectors(hc, hc2);


  // OK. Now we have an implementation structure of C(z) that is amenable to reflecting the zeros
  // in the unit circle. When expressing the FIR part as coefficient array, we need to just reverse
  // it. The FIR part implements:
  //
  //   y[n] = c0 * x[n] + c1 * x[n-1] + cM1 * x[n-(M+1)] + cM2 * x[n-(M+2)]
  //
  // Reversing that amounts to using:
  //
  //   r0 = cM2, r1 = cM1, rM1 = c1, rM2 = c_0

  // Define the reversed FIR coefficients:
  Real r0  = cM2;
  Real r1  = cM1;
  Real rM1 = c1;
  Real rM2 = c0;

  // Helper function to apply the reversed correction filter:
  auto getSampleR = [&](Real x)
  {
    // Apply 1-pole:
    Real t = op.getSample(x);

    // Apply the FIR part. It has 4 terms:
    Real y = 0;
    y += r0  * t;
    y += r1  * ud.getSample(t);
    y += rM1 * dlM1.getSample(t);
    y += rM2 * dlM2.getSample(t);
    return y;
  };

  // Reset objects:
  op.reset();
  ud.reset();
  dlM1.reset();
  dlM2.reset();

  // Produce output with modified correction filter:
  Vec hr(N);
  for(int n = 0; n < N; n++)
    hr[n] = getSampleR(hu[n]);
  //rsPlotVectors(hr);                        // Bipolar spike that gets progressively smeared


  // Write outputs to wavefiles for listening:
  rosic::writeToMonoWaveFile("FeedbackAllpassUncorrected.wav", &hu[0], N, sampleRate);
  rosic::writeToMonoWaveFile("FeedbackAllpassCorrected.wav",   &hr[0], N, sampleRate);

  // Plot uncorrected and reverse corrected output. These are the most interesting things:
  rsPlotVectors(hu, hr);


  // Observations:
  //
  // - Without feedback damping (i.e. dampGain == 1), we see a spike train as impulse response.
  //   the first spike occurs at 100 (== M), the second at 201 (== M + (M+1)), the third at 302
  //   (== M + (M+1) + (M+1)), etc. If the feedback gain k is less than 1, the spikes decay. The 
  //   fact that the spikes after the first ares spaced out by M+1 rtaher than M is due the fact
  //   that there is this additional unit delay in the feedback loop.
  //  
  // - With dampGain < 1 (like 0.7), the spikes will progressively turn more and more into a 
  //   spike-with-tail sort of shape. With dampGain = 0.7, we can even use a feedback gain k of 1
  //   and the overall trend is still decaying. It might be unstable at DC, though. Maybe test with
  //   longer N.
  //
  // - The compensation filter as implemented does indeed undo the effect of the feedback loop. The
  //   compensated output is a single unit spike at M as we expect from a length M delayline.
  //
  // - The peak envelope of the impulse response shows a characteristic knee like in a two-stage
  //   exponential decay. I think, this is because the peaks do not only decay but also smear out 
  //   such that they are less peaky later. Try dampFreq = 3000, dampGain = 0.7, k = 1.0. At the
  //   end, we see an almost sinusoidal wave with period = M+1 (I think). With k = 1.0, I guess,
  //   we will finally settle at some nonzero DC. Only with k < 1, this DC will also decay away to
  //   zero eventually (that's my prediction - not yet tried).
  //
  // - When listening to the wavefiles, the uncorrected one has a clear tonal character. It's a bit
  //   similar to a Karplus strong algo initialized with a unit impulse. The spectrum of it has a
  //   clearly harmonic structure. In the corrected one, the tonal character is much less present 
  //   and the spectrum is indeed flat! It wöööörks!!!! YAY!!!! This filter could be very useful!
  //
  //
  //  ToDo:
  //
  // - Try to find the zeros of the filter C(z) and reflect them in the unit circle. If we have a 
  //   direct form expression of C(z), this can be achieved by reversing the array of b-coeffs, I 
  //   think. Do this under the assumption A(z) = z^-M first before trying to attempt the 
  //   general case. The idea is that we don't really gain anything by just undoing the effect of
  //   the feedback. What we actually want is to undo only the effect of the feedback on the 
  //   magnitude response. But we want to retain effects on the phase response. Let's also assume
  //   a 1-pole/1-zero feedback filter F(z) = (b0 + b1*d) / (1 + a1*d) with d = z^-1.
  //
  // - In a production implemenation, we may use a single delayline to realize d^(M+2), d^(M+1), d.
  //   but maybe d should not be realized by the delyline. Maybe using a unit delay is more 
  //   efficient for this.
  //
  // - Wrap the filter into a class for convecient use. Make more unit tests.
  //
  // - To use it for decorrelation purposes, use a series connection of a bunch of them for left 
  //   and right channel with different delays for left and right but in such a way that the sum
  //   of the delays is equal. The parametrization of each such delay unit should be in terms of 
  //   the delay and the decay time. The decay times for all the filters in the chain should be
  //   equal, I think. Maybe wite a class rsAllpassDecorrelator. Maybe it should support more than
  //   two channels. Maybe we should also support fractional delays...but maybe not. That sounds
  //   complicated. I don't think, we can jsut shove in a fractional delayline instead of z^-M 
  //   because the formulas that I derived relied on these z^-M terms. If they get replaced by
  //   something more complex, the whole derivation has to be redone for the more complex case.
  //
  // - Maybe write about the filter on GitHub or maybe make a pdf.
  //
  // - Generalize to higher order feedback filters and maybe more complicated allpass filters. But
  //   I guess, when we do this, the direct form representations won't be sparse anymore. 
  //
  // - To optimze the delaylines, use one for both correction delays and call getSample once and 
  //   use the implict tap-move within getSample to call a second readout after getSample to get
  //   the other value.
}

void dampedAllpassComb2()
{
  // Now, we use the class rsDampedAllpassComb which encapsulates the algorithm derived above into a 
  // class. We use that class here to test some different options for the sign of the feedback gain
  // and the swapping of damper and delay. This leads to 4 different modes that all have different
  // properties and sound different.

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedAllpassComb<Real, Real>;

  // User parameters:
  int  delay      =   100;     // Main delay roundtrip length in samples. Is M-1 in the algo
  int  numSamples =  8192;     // Number of samples to generate
  Real sampleRate = 44100;     // Sample rate for writing the wavefiles
  Real dampFreq   =  1000;     // Frequency (in Hz) of the shelf filter for feedback damping
  Real dampGain   =     0.7;   // Linear high freq damping gain
  Real feedback   =     0.99;  // Feedback gain factor

  // Helper function to set up the given flt object with the given settings:
  auto setupComb = [&](Allpass& flt, 
    int delay, Real feedback, Real omega, Real hiGain, bool preDelay)
  {
    flt.setMaxDelayInSamples(delay);
    //flt.setupHighDamp(delay, feedback, omega, hiGain, preDelay);
    rsSetupHighDamp(flt, delay, feedback, omega, hiGain, preDelay);
  };

  // Helper function to produce the impulse response for given settings:
  auto impResp = [&](int delay, Real feedback, Real omega, Real hiGain, bool preDelay, int mode)
  {
    Allpass comb;
    setupComb(comb, delay, feedback, omega, hiGain, preDelay); 
    Vec h = impulseResponse(comb, numSamples, 1.0);
    return h;
  };

  // Set up the 4 damped allpass comb filters:
  int  N  = numSamples;
  int  d  = delay;
  Real w  = 2*PI*dampFreq/sampleRate;
  Real g  = dampGain;
  Real k  = feedback;

  // We have 4 different modes: unipolar/bipolar (selected by sign of k) and predelay or not 
  // (selected by bool parameter):                                             
                                            // polar  predelay 
  Vec h1 = impResp(d, +k, w, g, false, 0);  //  bi      no
  Vec h2 = impResp(d, -k, w, g, false, 0);  //  uni     no
  Vec h3 = impResp(d, +k, w, g, true,  0);  //  bi      yes
  Vec h4 = impResp(d, -k, w, g, true,  0);  //  uni     yes

  // Plot them all together and then one at a time:
  rsPlotVectors(h1, h2, h3, h4);
  //rsPlotVectors(h1);
  //rsPlotVectors(h2);
  //rsPlotVectors(h3);
  //rsPlotVectors(h4);

  // Let's compare the output with predelay to the one without - like h1 and h3 - but shift h3
  // to the left for a match. We want to see if they exactly match when shifted properly
  rsPlotArrays(N-delay+1, &h1[0], &h3[delay-1]);
  rsPlotArrays(N-delay+1, &h2[0], &h4[delay-1]);


  // Observations:
  //
  // - The sign of the feedback gain switches between unipolar and bipolar spike trains. A positive
  //   k gives bipolar, a negative unipolar ones. This might seem counterintuitive and has to do 
  //   with the math of transfer functions and the sign conventions used which imply that feedback 
  //   coeffs get a negative sign in the implementation. For a user parameter, we may want to 
  //   switch that behavior but I have decided against it in the DSP code such that it can stay 
  //   closer to the math equations.
  //
  // - The output with and without predelay are indeed exactly equal, if one shifts them 
  //   accordingly for alignment. The desired shift is exactly delay-1. That is, the predelay
  //   is exactly delay-1 samples. This is also kinda counterintuitive. A user might expect it to
  //   be "delay" samples. Maybe when turning this into an end-user facing unit, we should just
  //   artificially introduce an additional sample of delay in "pre-delay" mode.
}

void dampedAllpassComb3()
{
  // Now, we create a chain of 4 such allpass comb filters and produce its impulse response.

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedAllpassComb<Real, Real>;

  // User parameters:
  int  delay1     =    53;     // 1st main delay roundtrip length in samples. Is M-1 in the algo
  int  delay2     =    67;     // Same for 2nd comb
  int  delay3     =    83;     // ..3rd
  int  delay4     =   101;     // ..and finally 4th
  int  numSamples =  8192;     // Number of samples to generate
  Real sampleRate = 44100;     // Sample rate for writing the wavefiles
  Real dampFreq   =   500;     // Frequency (in Hz) of the low shelf for feedback damping
  Real dampGain   =     0.9;   // Linear high freq damping gain
  Real feedback   =     0.99;  // Feedback gain factor

  // Helper function to set up the given flt object with the given settings:
  auto setupComb = [&](Allpass& flt, int delay, Real feedback, Real omega, Real hiGain, 
    bool preDelay)
  {
    flt.setMaxDelayInSamples(delay); 
    //flt.setupHighDamp(delay, feedback, omega, hiGain, preDelay);
    rsSetupHighDamp(flt, delay, feedback, omega, hiGain, preDelay);
  };

  // Set up the 4 damped allpass comb filters:
  Real w  = 2*PI*dampFreq/sampleRate;
  Real g  = dampGain;
  Real k  = feedback;
  int  d1 = delay1, d2 = delay2, d3 = delay3, d4 = delay4;
  int  delaySum = d1 + d2 + d3 + d4;
  Allpass apf1; setupComb(apf1, d1, k, w, g, false);
  Allpass apf2; setupComb(apf2, d2, k, w, g, false);
  Allpass apf3; setupComb(apf3, d3, k, w, g, false);
  Allpass apf4; setupComb(apf4, d4, k, w, g, false);

  // Generate impulse response of the allpass chain:
  int N  = numSamples;
  Vec h1 = impulseResponse(apf1, N, 1.0);
  Vec h2 = filterResponse( apf2, N, h1);
  Vec h3 = filterResponse( apf3, N, h2);
  Vec h4 = filterResponse( apf4, N, h3);

  // Write outputs to wavefiles and plot:
  rosic::writeToMonoWaveFile("DampedAllpasComb1.wav", &h1[0], N, sampleRate);
  rosic::writeToMonoWaveFile("DampedAllpasComb2.wav", &h2[0], N, sampleRate);
  rosic::writeToMonoWaveFile("DampedAllpasComb3.wav", &h3[0], N, sampleRate);
  rosic::writeToMonoWaveFile("DampedAllpasComb4.wav", &h4[0], N, sampleRate);
  rsPlotVectors(h1, h2, h3, h4);


  // Observations:
  //
  // - The first spike in the output seems to occur at the sum off all (delay-1) values, i.e. at 
  //   delaySum - numDelays. So, it introduces a significant predelay/latency which is not so 
  //   desirable.
  //
  // - The outputs become progressively more complex.
  //
  // - When one reduces the decay time by using feedback gains closer to one and damping gains
  //   also closer to 1, the sound becomes more metallic. The outputs of the later stages also tend
  //   to sound more metallic than thos of the previous stage. So, adding more stages does not seem 
  //   reduce metallicness - to the contrary. However, it also gets longer, so wehn one would 
  //   compensate for that by using less feedback, the effect might be canceled. Maybe it can also 
  //   be mitigated by better tuning of the delays.
  //
  // - When using feedback = 1, dampGain = 1, the result is strange. There's only one actual spike.
  //   Figure out what is going on!
  //
  //
  // ToDo:
  //
  // - Maybe try to apply it to other signals to see what it does to them. Maybe noise, sawtooth,
  //   guitar, etc.
  //
  // - It looks like the spike spacing is not equal to delay but one sample less - check that!
  //
  // - Compute the required feedback factor and dampGain from a desired decay time that the user
  //   specifies. Look up the formulas in the FDN implementation or literature.
  //
  // - Figure out what is different between ordering the combs from short to long anf long to 
  //   short. Is ther any difference? No - this can't be the case because they are all LTI!
  //
  // - Make a helper function setupFilter()
  //
  // - Some prime numbers to experiment with:
  //   2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 67, 71, 83, 89, 101, 107, 109, 
  //   113, 127, 131, 137, 139, 149, 157, 167, 179, 181, 191, 197, 199, 211, 227, 233, 239, 251, 
  //   257, 263, 269, 281, 293, 307, 311, 317, 337, 347, 353, 359, 379, 389, 401, 409
  //
  // - The setupComb function should not directly take a feedback and highFreqGain parameter but
  //   instead a desired decay time and maybe a decaytime scaler for the high frequencies. These
  //   should then be used to compute the feedback gain and the shelver gain. This is a 
  //   parametrization closer to what we might present to a user. The formula is in
  //   FeedbackDelayNetwork16::updateDampingAndCorrectionFilters(). I think the formula for the 
  //   feedback gain is:
  //
  //     gain = pow(10.0, -3.0*delay / (time * sampleRate));
  //
  //   we need this formula as is for computing the feedback gain k and to compute the damping 
  //   gain, we use the same formula but with the desired decay time for high freqs. But maybe it 
  //   needs to be scaled by the k factor because the high freqs, actually both gain factors are 
  //   applied, i.e. the product of both. We should have asetup function like
  //   setupComb = [&](Allpass& flt, int delay, Real decayTime, Real omega, Real hiDecayTimeScale)
  //
  // - Wait! I think, the sound has changed after switching from:
  //
  //     out = mainDelay.getSample(in + k * damper.getSample(out));  to
  //     out = mainDelay.getSample(in - k * damper.getSample(out));
  //
  //   We really want to offer both modes! I think, the one with the - sounds more metallic. When 
  //   using the minus, the spikes are altenrting up and down. When using the plus, they are all
  //   upward. Well - from a user perspective, it makes more sense if passing a negative k leads to
  //   alternating mode. So, maybe revert to using + k * ... This then requires to re-introduce the
  //   minusses in the r-coeffs as well.
  //   
  // - I think, we can provide 4 different modes and they all sound different: w/o predelay
  //   and negative and positive k.
}

void dampedAllpassComb4()
{
  // Here, we try to use higher order feedback damping filters. Specifically, we use a biquad that
  // realized a wide-band dip (i.e. downward bell) frequency response.

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Complex = rsComplex<Real>;
  using Allpass = rsDampedAllpassComb<Real, Real>;

  // User parameters:
  int  N        = 8192;
  int  delay    = 100;
  Real feedback = 0.9;

  // Design a wideband dip filter to be used as damping filter:
  Real b[3];
  Real a[3];
  a[0] = 1;
  rsStateVariableFilter<Real, Real> svf;  // Only used for designing the feedback filter
  svf.setupBell(0.2, 0.3, 0.5);           // A wideband dip filter
  svf.convertToBiquad(&b[0], &b[1], &b[2], &a[1], &a[2]);

  Allpass ap;
  ap.setMaxDelayInSamples(delay);
  ap.setup(delay, feedback, 2, b, a, false);
  Vec h = impulseResponse(ap, N, 1.0);
  bool ok = isAllpass(h, 1.e-3);
  rsPlotVectors(h);

  // ToDo:
  //
  // - Add more dispersion to the feedback path by using a maximum phase rather than minimum phase
  //   biquad. This can be achieved by reversing the b-array, i.e. swapping b[0] and b[2]. But then
  //   we can't use the mode without predelay. See documentation of rsDampedAllpassComb. But maybe
  //   it can actually be made to work if we don't use the inverse damping filter as is but reflect
  //   its poles.
}

void dampedAllpassComb5()
{
  // Under construction.
  //
  // We want to set up an rsDampedAllpassComb that achieves a desired overall decay time (in the 
  // RT60 sense) and also allows that decay time to be scaled at low and high frequencies via 
  // shelving filters. ...TBC...

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedAllpassComb<Real, Real>;

  // User parameters:
  int  sampleRate = 48000;     // Sampling rate.
  int  numSamples =  8192;     // Number of samples to render.
  int  delay      =    50;     // Delay in samples.
  Real decayTime  =     1.0;   // Decay time for mid frequencies in seconds.
  Real lowFreq    =   250.0;   // Crossover freq between low and mid frequencies in Hz.
  Real lowScale   =     1.5;   // Decay time scaler for low frequencies.
  Real highFreq   =  4000.0;   // Crossover freq between mid and high frequencies in Hz.
  Real highScale  =     0.5;   // Decay time scaler for high frequencies.

  // Compute intermediate values:
  Real decaySamples = decayTime     * sampleRate;
  Real lowOmega     = 2*PI*lowFreq  / sampleRate;
  Real highOmega    = 2*PI*highFreq / sampleRate;

  // Create and set up the allpass filter:
  Allpass ap;
  ap.setMaxDelayInSamples(delay);
  rsSetupDecayTimes(ap, delay, decaySamples, lowOmega, lowScale, highOmega, highScale, false);

  // Get impulse response:
  Vec h = impulseResponse(ap, numSamples, 1.0);
  rsPlotVectors(h);

  int dummy = 0;

  // Observations:
  //
  // - With longer decay times, the initial spike grows larger with respect to the decaying tail.
  //
  //
  // ToDo:
  //
  // - Plot an energy decay relief and check if it looks like expected.
  //
  // - Maybe give the user an option to switch between positive and negative feedback.
  //
  // - Make it work also when lowScale and/or highScale is greater than 1. ..I think it should 
  //   already work -> verify that!
  //
  // - Write result to wavefile and listen to it.
  //
  // - Try creating a series connection of a bunch of such allpasses with different values for the
  //   delay. Maybe try giving them different decay times. Maybe the average of the decay times 
  //   should stay the same but the decay time of the individual allpasses should scale with the 
  //   delay somehow. Figure out, how chaining such filters affects the overall decay time. Will 
  //   the decay times add up?
}

void dampedAllpassCombComplex()
{
  // We instantiate rsDampedAllpassComb with a complex datatype for the signals. The feedback gain
  // is also complex. Using a complex feedback gain further increases the space of the things that 
  // we can do with this filter. 

  // Define types to be used:
  using Real    = double;
  //using Complex = rsComplex<Real>;
  using Complex = std::complex<Real>;
  using VecR    = std::vector<Real>;
  using VecC    = std::vector<Complex>;
  using Allpass = rsDampedAllpassComb<Complex, Real>;

  // User parameters:
  int  delay      =   100;
  int  numSamples = 16384;
  Real sampleRate = 44100;
  Real dampFreq   =   500;
  Real dampGain   =     0.8;
  Real fbGain     =     0.99;
  Real fbPhase    =    PI/4;   // Try: pi/4, 0.7, 1.7

  // Create and set up the complex allpass and plot real and imaginary part of output:
  Allpass ap;
  Complex j(0,1);
  Complex fb = fbGain * rsExp(j*fbPhase);
  Real    w  = 2*PI*dampFreq/sampleRate;
  int     N  = numSamples;
  ap.setMaxDelayInSamples(delay);
  //ap.setupHighDamp(delay, fb, w, dampGain, false);
  rsSetupHighDamp(ap, delay, fb, w, dampGain, false);
  VecC h = impulseResponse(ap, N, Complex(1));
  plotComplexVectorReIm(h);                      // Doesn't accept rsComplex - fix that!

  // Extract real and imaginary parts:
  VecR hr(N), hi(N);
  for(int n = 0; n < N; n++)
  {
    hr[n] = real(h[n]);
    hi[n] = imag(h[n]);
  }

  // Plot spectra:
  using SpecPlot = SpectrumPlotter<Real>;
  SpecPlot plt;
  plt.setFftSize(N);
  plt.plotSpectra(N, &hr[0], &hi[0]);


  // Observations:
  //
  // - Using complex numbers imprints and undulation onto the impulse reponses.
  //
  // - The spectra of real and imaginary part by themselves are not allpass like. They are more
  //   of a strange comb like structure. One has zeros between the peaks, the other doesn't. 
  //
  // - Feedback phases of 0 and pi give purely real outputs. With 0, it's unipolar, with pi its
  //   bipolar. That's how it has to be - we expect to get back to behavior of positive and
  //   negative signs. pi/2 gives a bipolar spike train both real and imaginary part. Same 
  //   for 3*pi/2.
  //
  //
  // ToDo:
  //
  // - Make sure that everything works with rsComplex and std::complex for Complex. We need to
  //   implement real/imag for rsComplex or rsReal/rsImag for both. I think, the latter way is the
  //   right way. It's consistent with how the rest of the library works.
  //
  // - I think, the feedback phase should scale with the delay to achieve a uniform undulation
  //   frequency. Figure this out!
}

void dampedAllpassCombNonLin()
{
  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedAllpassCombNonLin<Real, Real>;

  // User parameters:
  int  delay      =   100;     // Main delay roundtrip length in samples. Is M-1 in the algo
  int  numSamples =  2000;     // Number of samples to generate
  Real sampleRate = 44100;     // Sample rate for writing the wavefiles
  Real dampFreq   =  1000;     // Frequency (in Hz) of the shelf filter for feedback damping
  Real dampGain   =     0.7;   // Linear high freq damping gain
  Real feedback   =     1.5;   // Feedback gain factor

  Allpass ap;
  Real dampOmega = 2*PI*dampFreq/sampleRate;
  ap.setMaxDelayInSamples(delay);
  //ap.setupHighDamp(delay, feedback, dampOmega, dampGain, false);
  rsSetupHighDamp(ap, delay, feedback, dampOmega, dampGain, false);
  int N = numSamples;
 
  Vec h1  = impulseResponse(ap, N,  1.0);
  Vec h2  = impulseResponse(ap, N,  2.0);
  Vec h4  = impulseResponse(ap, N,  4.0);
  Vec h8  = impulseResponse(ap, N,  8.0);
  Vec h16 = impulseResponse(ap, N, 16.0);

  rsPlotVectors(h1, h2, h4, h8, h16);


  // - With feedback = 1.5, it becomes unstable despite the saturation the  feedbakc path. But:
  //   the saturation applies only to the comb. I guess, the correction filter may be the cause 
  //   for the overall instability. Maybe we need to apply some nonlinearity there, too.
}

void dampedAllpassDelayContent()
{
  // An experiment that lets us look at the content of the delaylines of rsDampedAllpassComb. To 
  // make this experiment work, one needs to temporarily move the mainDelay, corrDelay members
  // into the public section. Then one can uncomment the line:
  //
  //   rsPlotDelayLineContent(ap.mainDelay, ap.corrDelay);
  //
  // in the for-loop. What it them will do is to show a plot of the contents of the delaylines at 
  // each sample.


  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedAllpassComb<Real, Real>;

  // User parameters:
  int  delay      =    10;
  int  numSamples =   100;
  Real sampleRate = 44100;
  Real dampFreq   =   100;
  Real dampGain   =     0.5;
  Real feedback   =     0.9;
  bool predelay   = false;

  Allpass ap;
  Real dampOmega = 2*PI*dampFreq/sampleRate;
  ap.setMaxDelayInSamples(delay);
  //ap.setupHighDamp(delay, feedback, dampOmega, dampGain, predelay);
  rsSetupHighDamp(ap, delay, feedback, dampOmega, dampGain, predelay);
  int N = numSamples;

  Vec h(N);
  for(int n = 0; n < N; n++)
  {
    h[n] = ap.getSample(cos(0.8*n));
    //rsPlotDelayLineContent(ap.mainDelay, ap.corrDelay); 
    // To call this, we temporarily need to move the delaylines into the public section
  }

  // Observations:
  //
  // - The delaylines do not seem to have shared content. So we do indeed seem to need both of 
  //   them.
}

void dampedSchroederAllpass()
{
  // We test the class rsDampedSchroederAllpass here. We compare its output against two naive
  // prototype implementations and one naive implementation that we implement directly here.

  // Define types to be used:
  using Real     = double;
  using Vec      = std::vector<Real>;
  using Allpass0 = rsDampedSchroederAllpass<Real, Real>;
  using Allpass1 = rsDampedSchroederAllpassNaive<Real, Real>;
  using Allpass2 = rsDampedSchroederAllpassNaive2<Real, Real>;

  // User parameters:
  int  delay      =    20;     // Delay roundtrip length in samples. 
  int  numSamples =  1024;     // Number of samples to generate
  Real sampleRate = 44100;     // Sample rate for writing the wavefiles
  Real feedback   =    +0.9;   // Feedback gain factor

  // Abbreviations for convenience:
  int  N = numSamples;
  int  M = delay;
  Real k = feedback;

  // We use a 2-point FIR filter in the feedback path:
  Real b[2] = { 0.75, 0.25 };  // Examples: 0.75,0.25; 0.6,0.4; 0.9,0.1; 0.2,0.8

  // Create impulse- and magnitude response of proper implementation:
  Allpass0 ap0;
  ap0.setMaxDelayInSamples(M);
  ap0.setup(M, k, 1, b);
  Vec h0 = impulseResponse(ap0, N, 1.0);
  Vec m0 = rsSpectralMagnitudes(h0);

  // Create impulse- and magnitude response of first naive implementation:
  Allpass1 ap1;
  ap1.setMaxDelayInSamples(M);
  ap1.setup(M, k, 1, b);
  Vec h1 = impulseResponse(ap1, N, 1.0);
  Vec m1 = rsSpectralMagnitudes(h1);

  // Create impulse- and magnitude response of second naive implementation:
  Allpass2 ap2;
  ap2.setMaxDelayInSamples(M);
  ap2.setup(M, k, 1, b);
  Vec h2 = impulseResponse(ap2, N, 1.0);
  Vec m2 = rsSpectralMagnitudes(h2);

  // Implement the difference equation of the filter directly. With a 1st order FIR filter in
  // the feedback path,  we have the following transfer function:
  //
  //            k*b0 + k*b1*d + d^M         M=5     k*b0 + k*b1*d + d^5
  //  H(z) = -----------------------------   =   --------------------------
  //          1 + k*b1*d^(M-1) + k*b0*d^M         1 + k*b1*d^4 + k*b0*d^5
  //
  Vec d(N), h3(N);  // d[n] is the unit impulse, i.e. the Diract delta function
  d[0] = 1;
  rsInfiniteDataStream<Real> x(&d[0], N), y(&h3[0], N);
  y.setZero();
  for(int n = 0; n < N; n++)
  {
    Real tmp = 0;
    tmp += k*b[0]*x[n] + k*b[1]*x[n-1] + x[n-M];
    tmp -= k*b[1]*y[n-(M-1)] + k*b[0]*y[n-M];
    y[n] = tmp;
  }
  Vec m3 = rsSpectralMagnitudes(h3);

  // Plot all 4 impulse and magnitude responses. They should all match and be allpass:
  rsPlotVectors(h0, h1, h2, h3);
  rsPlotVectors(m0, m1, m2, m3);

  // Post-process impulse response and write it to a wave file:
  Real period = delay / sampleRate;
  Real freq   = 1 / period;
  Real cutoff = 0.5 * freq;    // Tune lowpass cutoff an octave below the signal frequency
  RAPT::rsOnePoleFilter<Real, Real> lpf;
  lpf.setMode(RAPT::rsOnePoleFilter<Real, Real>::modes::LOWPASS_IIT);
  lpf.setCutoff(cutoff);
  lpf.setSampleRate(sampleRate);
  Vec z = filterResponse(lpf, N, h0);
  rsArrayTools::normalize(&z[0], N);
  rsPlotVectors(h0, z);
  rosic::writeToMonoWaveFile("DampedSchroederAllpass.wav", &z[0], N, (int) sampleRate);


  // Observations:
  //
  // - With high positive feedback and a lowpassish feedback filter, the impulse response tends
  //   to a periodic waveform with period 2*M. With negative feedback, the period is M.
  //
  // 
  // ToDo:
  //
  // - Implement a variant that uses a sparse FIR in the feedback path. But then - how would we 
  //   ensure stability? Maybe experimentally. If we want to use such filters in the context of a
  //   reverb, we'll tune these during development anyway. Such sparse FIR feedback filter settings
  //   are not something, I'll want to expose as end-user parameters anyway.
  //
  // - Before writing the results to wavefiles, post-process them by a 1st order lowpass and 
  //   normalization. That makes the non-allpass "percussions" more musically useful.
  //
  // - Maybe tune the lowpass to a cutoff related to the length of the delayline which determines 
  //   the period of the late signal. Maybe create a sample rendering function that creates a 
  //   couple of samples that can be used as transients. Such allpass transients an perhaps be used
  //   together with modal synthesis. The modal synthesis creates the body of the sound and the 
  //   allpasses create the transients.
}

void dampedAllpassComb()
{
  dampedAllpassComb5();
  //dampedAllpassComb4();
  //dampedSchroederAllpass();

  dampedAllpassComb1();
  dampedAllpassComb2();
  dampedAllpassComb3();
  dampedAllpassComb4();
  dampedAllpassComb5();
  dampedAllpassCombComplex();
  dampedAllpassCombNonLin();
  dampedAllpassDelayContent();
  dampedSchroederAllpass();

  // ToDo: 
  //
  // - Plot the phase-delay and group-delay as function of frequency for the damped allpass comb 
  //   and damped Schroeder allpass. We can do this by evaluating the complex frequency response 
  //   using the getTransferFunctionAt() functions.
}