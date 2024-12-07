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
  // Under construction

  using Real = double;
  using VecI = std::vector<int>;
  using VecR = std::vector<Real>;
  //using APF  = rsTwoPoleAllpassDelay<Real, Real>;
  using APF  = rsTwoPoleAllpassDelayNaive<Real, Real>;




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
      filters[i].setAllpassCoeffs(a1, a2); // Or should we include a minus? But I don't think so.

      int dummy = 0;
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

    
    rsPlotVectors(h);

  
    int dummy = 0;

    // ToDo:
    //
    // - It's a bit silly to use the state variable filter for the purpose of designing biquad 
    //   coeffs so maybe replace this code later with the direct RBJ biquad design formulas.
  };


  create({ 5 }, { 0.2 }, { 2.5 }, 500);

  int dummy = 0;
}