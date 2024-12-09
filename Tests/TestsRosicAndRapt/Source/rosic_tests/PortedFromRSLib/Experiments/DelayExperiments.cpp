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

void feedbackFilterAllpass()
{
  // I try to implement an idea for starting with an arbitrary given allpass filter A(z) and 
  // arbitrary given feedback filter F(z) that sits in a feedback loop with unit delay around that 
  // allpass. I try to design a compensation filter that can be applied in series to this setup
  // such that the overall transfer function is allpass in nature. Without the compensation 
  // filter, this setup has the uncompensated transfer funcion:
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
  
  // Consider our special case here with  A(z) = z^-M  and  F(z) = (b0 + b1*z^-1) / (1 + a1*z^-1). 
  // So we have:
  //
  //               k * z^-1 * z^-M * (b0 + b1*z^-1)     1 + a1*d + b0*k*d^(M+1) + b1*k*d^(M+2)
  //   C(z) = 1 + ---------------------------------- = -----------------------------------------
  //                      1 + a1*z^-1                                 1 + a1*d
  //
  // With d := z^-1 for convenience. In a second step, we try to implement C(z) in this form using
  // a delaline to implement the denominator the one-pole filter 1 / (1 + a1*z^-1) for the 
  // denominator. Having the filter in this form makes it amenable for manipulating the zeros by 
  // reflecting them in the unit circle. We just need to reverse the array of numerator coeffs, I 
  // think.
  //
  // In the private repo, there's a filter AllpassStuff.txt where it's explained a bit more. Maybe
  // drag that into the main repo. But it's currently too messy for a public repo.

  using Real    = double;
  using Vec     = std::vector<Real>;
  using Delay   = RAPT::rsBasicDelayLine<Real>;      // We use a simple delay as allpass
  using OnePole = RAPT::rsOnePoleFilter<Real, Real>; // We use a one pole as feedback filter

  int    M          =   100;     // Delay
  int    N          = 10000;     // Number of samples to generate
  double sampleRate = 44100;
  double dampFreq   =  3000;     // Frequency of the low shelf for feedback damping
  double dampGain   =     0.7;   // Linear high freq damping gain
  double k          =     1.0;   // Feedback gain factor


  // Create and set up the two given filters for A(z) and F(z):
  Delay apf;                                // Allpass filter
  apf.setMaximumDelayInSamples(M);
  apf.setDelayInSamples(M);

  OnePole fbf;                              // Feedback filter
  fbf.setMode(OnePole::modes::HIGHSHELV_BLT);
  fbf.setCutoff(dampFreq);
  fbf.setShelvingGain(dampGain);


  // Helper variable and function to implement the feedback loop filter:
  double s = 0;                             // Output and state of uncompensated filter
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
  rsPlotVectors(hu);                        // Decaying spike train with progressive tail smear


  // Cancel the effect of the feedback loop. For this, we re-use the existing filter objects. We 
  // can do this because they are not needed anymore for other purposes because the uncompensated 
  // output as been generated already. in a realtime implementation, we would have to use another
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
  rsPlotVectors(hc);                        // Yes - that looks good! Single spike at M.
  

  // Now we try to implement the correction filter in a different way that is more amenable to 
  // reflecting the zeros. Namely, in the form:
  //
  //           1 + a1*d + k*b0*d^(M+1) + k*b1*d^(M+2)
  //   C(z) = ----------------------------------------
  //                   1 + a1*d

  Real b0 = fbf.getB0();
  Real b1 = fbf.getB1();
  Real a1 = fbf.getA1();
  Delay   dl;                               // For the numerator
  OnePole op;                               // For the denominator
  op.setCoefficients(1.0, 0.0, a1);         // Should have same poles as fbf and 1 as denominator

  // rsOnePoleFilter uses the other sign convention for the feedback coeffs, so if we 
  // want to use the a1 coeff directly, we have to negate it:
  a1 = -a1;  

  // For some reason that I don't quite understand, we need to negate the b-coeffs, too:
  b0 = -b0;
  b1 = -b1;
  // Maybe what really happens that we need to negate k when it moves into a numerator? I'm not
  // quite sure - but it works!

  // Naive implementation using for each desired delay its own delayline:
  ud.reset();
  Delay dlM1;
  dlM1.setMaximumDelayInSamples(M+1);
  dlM1.setDelayInSamples(M+1);
  Delay dlM2;
  dlM2.setMaximumDelayInSamples(M+2);
  dlM2.setDelayInSamples(M+2);

  // Helper function:
  auto getSampleC2 = [&](Real x)
  {
    // Apply 1-pole:
    Real t = op.getSample(x);

    // Apply the FIR part:
    Real y = t;                             // 1
    y +=   a1*ud.getSample(t);              // a1*d
    y += k*b0*dlM1.getSample(t);            // b0*k*d^(M+1)
    y += k*b1*dlM2.getSample(t);            // b1*k*d^(M+2)
    return y;
    // This works! But why?! The b-coeffs are not supposed to flip signs between transfer function
    // and difference equation! Ah! I think, it's because the involve the factor k, which does 
    // indeed come from a feedback path...but that seems weird anyway.
  };

  // Produce corrected output with alternative implementation:
  Vec hc2(N);
  for(int n = 0; n < N; n++)
    hc2[n] = getSampleC2(hu[n]);

  rsPlotVectors(hc, hc2); 


  // OK. Now we have an implementation structure of C(z) that is amenable to reflecting the zeros
  // in the unit circle. When expressing the FIR part as coefficient array, we need to just reverse
  // it. The FIR part implements
  //
  //   y[n] = c_0 * x[n] + c_1 * x[n-1] + c_{M+1} * x[n-(M+1)] + c_{M+2} * x[n-(M+2)]
  //
  // where
  //
  //   c_0 = 1, c_1 = a1, c_{M+1} = k*b0, c_{M+2} = k*b1
  //
  // Reversing that amounts to using:
  //
  //   r0 = c_{M+2}, r1 = c_{M+1}, r_{M+1} = c_1, r_{M+2} = c_0
  //
  //
  // ...TBC...



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
  //   efficient for this
}