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
  // Stub. I try to implement an idea for starting with and arbitrary give allpass filter and 
  // arbitrary given feedback filter that sits in a feedback loop with unit delay around that 
  // allpass. I try to design a compensation filter that can be applied in series to this setup
  // such that the overall transfer function is allpass in nature...TBC...

  using Real = double;
  using Vec  = std::vector<Real>;
  using APF  = RAPT::rsBasicDelayLine<Real>;      // We use a simple delay as allpass
  using FBF  = RAPT::rsOnePoleFilter<Real, Real>; // We use a one pole as feedback filter


  int    M          =   100;     // Delay
  int    N          = 10000;     // Number of samples to generate
  double sampleRate = 44100;
  double dampFreq   =  1000;     // Frequency of the low shelf for damping.
  double dampGain   =     0.7;   // Linear high freq damping gain
  double k          =     0.95;  // Feedback gain factor

  APF apf;
  apf.setMaximumDelayInSamples(M);
  apf.setDelayInSamples(M);

  FBF fbf;
  fbf.setMode(FBF::modes::HIGHSHELV_BLT);
  fbf.setCutoff(dampFreq);
  fbf.setShelvingGain(dampGain);

  // Helper variables and functions to implement the feedback loop filter:
  double u;            // Uncompensated filter output. State variable of the filter.
  auto reset = [&]()
  { 
    apf.reset();
    fbf.reset();
    u = 0; 
  };
  auto getSample = [&](Real in)
  {
    // This implements the feedback loop with unit delay without feedback filter:
    u = apf.getSample(in + k * fbf.getSample(u));
    return u;
    // The current value of u will used in the next call. This is the unit-delay feedback loop.
  };

  // Produce the uncompensated impulse response:
  Vec hu(N);
  reset();
  hu[0] = getSample(1.0);
  for(int n = 1; n < N; n++)
    hu[n] = getSample(0.0);
  rsPlotVectors(hu);


  // ToDo: Design the compensation filter and apply it to hu...TBC...



  int dummy = 0;

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
  //   longer N
}