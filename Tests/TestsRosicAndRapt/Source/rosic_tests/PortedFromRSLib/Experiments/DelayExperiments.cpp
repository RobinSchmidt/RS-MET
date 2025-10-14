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
  //
  // - When introducing modulation of the generalized Hadamard matrix parameters a,b,c,d, maybe use
  //   different values of a,b,c,d for the different levels of the trafo and modulate them by 
  //   different signals. The rationale is that we want to avoid situations where the matrix 
  //   becomes temporarily very sparse. This may happen if the modulation of the angle phi causes
  //   the matrix to become diagonal, for example. But when each stage/level of the trafo uses 
  //   different coeffs, this is unlikely to happen. Maybe we can even ensure that it doesn't 
  //   happen by relating the angle phi of various levels by fixed phase shifts. See
  //   rsStagedKroneckerTrafo2x2()
}

void delayLineBasic()
{
  static const int N = 20;
  using Real = double;
  Real t[N], h[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  rsDelay<Real> dl; 
  dl.setMaxDelayInSamples(5);
  dl.setDelayInSamples(5);
  RAPT::getImpulseResponse(dl, h, N);
  plotData(N, t, h);

  // Observations:
  //
  // - The impulse response is an impulse shifted by 5 samples, i.e. centered at sample index 5
  //   rather than 0. That's the expected result.
}

void delayLineLinear()
{
  // We plot the impulse responses of a linearly interpolating delayline with delays of 5.0, 5.1,
  // 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0 one after another.

  static const int N = 20;
  using Real = double;
  Real t[N], h[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  rsDelayLinear<Real, Real> dl; 
  dl.setMaxDelayInSamples(6);
  Real d0 = 5;                                // Reference delay
  for(int i = 0; i <= 10; i++)
  {
    Real f =  Real(i) / Real(10);
    Real d = d0 + f;
    dl.setDelayInSamples(d);
    RAPT::getImpulseResponse(dl, h, N);
    plotData(N, t, h);
  }


  // Observations:
  //
  // - The 1st impulse response is an impulse shifted by 5 samples, like in the case of the integer
  //   delayline. The 2nd has a spike oh hight 0.9 at 5 and a tail of 0.1 at 6. The 3rd a spike of
  //   height 0.8 at 5 with a tail of 0.2 at 6 and so on. Until it reaches a spike of height 1 at 6
  //   in the 10th response. This is all as it should be.
}

void delayLineAllpass()
{
  // This basically replicates the code of delayLineLinear but with an allpass interpolated 
  // delayline. Maybe templatize this function on the delayline type such that we can use the same
  // function for both types of delayline

  static const int N = 30;
  using Real = double;
  Real t[N], h[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  rsDelayAllpass<Real, Real> dl; 
  dl.setMaxDelayInSamples(6);

  dl.setDelayInSamples(5.01);
  RAPT::getImpulseResponse(dl, h, N);
  plotData(N, t, h);
  dl.setDelayInSamples(5.99);
  RAPT::getImpulseResponse(dl, h, N);
  plotData(N, t, h);

  Real d0 = 5;                                // Reference delay
  for(int i = 0; i <= 10; i++)
  {
    Real f = Real(i) / Real(10);
    Real d = d0 + f;
    dl.setDelayInSamples(d);
    RAPT::getImpulseResponse(dl, h, N);
    plotData(N, t, h);
  }

  // Observations:
  //
  // - The 1st impulse response is the same as for integer and linear delays. Then it gets weird.
  //   ...TBC...
  //
  // - When the delay is a little bit above an integer (like 5.01), we see an oscillation at the 
  //   Nyquist freq. It seems, its amplitude goes down, the closer we get to the integer (good) but
  //   its decay time goes up (bad).
  //
  //
  // ToDo:
  //
  // - Plot frequency responses, too.
}


bool universalCombVsOnePole()
{
  // Compares outputs of a one pole filter and universal comb to show the correspondence between 
  // the two types filters. The universal comb can be seen as a generalization to a general 1st 
  // order filter where the unit delay has been replaced by an arbitrary delay of M samples. This 
  // test demonstrates this by setting up a 1-pole filter and a uninversal comb filter with M=1 and
  // the same coefficients (up to different sign conventions for a1 unfortunately - an 
  // inconsistency in the APIs of the respective classes). In this case, both filters are supposed 
  // to produce the exact same impulse response which is what we verify here.


  bool ok = true;

  // Types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Comb    = rsUniversalCombFilter<Real, Real>;
  using OnePole = rsOnePoleFilter<Real, Real>;

  // User parameters:
  int  N  =  10;        // Number of samples to produce
  Real b0 =  0.75;      // Coeff for direct input
  Real b1 =  0.25;      // Coeff for delayed input
  Real a1 = -0.5;       // Coeff for delayed output (Negated or not? -> Document sign convention)

  // Set up one pole and obtain its impulse response:
  OnePole onePole;
  onePole.setCoefficients(b0, b1, -a1);
  Vec h1 = impulseResponse(onePole, N, 1.0);

  // Set up comb and obtain its impulse response:
  Comb comb;
  comb.setMaxDelayInSamples(10);
  comb.setDelayInSamples(1);
  comb.setCoeffs(b0, b1, a1);  // Has apparently different sign convenction than the 1-pole
  Vec h2 = impulseResponse(comb, N, 1.0);

  // Compute error and check that it's within numerical tolerance:
  Vec err = h2 - h1;
  ok &= rsIsAllZeros(err, 0.0);
  rsPlotVectors(h1, h2, err);

  return ok;


  // ToDo: 
  //
  // - Maybe move into unit tests. 
  //
  // - Explain what the replacement of the unit delay by an M sample delay does to the pole-zero 
  //   pattern and frequency response. I think, it just gets repeated M times. The effect on the 
  //   impulse response is that we intersperse M-1 samples of value zero between each pair of 
  //   samples in the 1-pole response.
  //
  // - There is actually no difference in numerical roundoff error so we use use 0.0 as our error 
  //   tolerance here. Maybe try less nice coefficients and see if we still can use zero tolerance.
  //
  // - Fix the API inconsistency - but with great(!) care because this is a (silent) breaking 
  //   change for client code. Adopt the convention that the transfer function is 
  //   H(z) = (b0 + b1*z^-1) / (a0 + a1*z^-1)  such that  y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1]
  //   consistently throughout the whole library for all direct form filters. Most of them already
  //   follow this convention but not all. The 1-pole is one example and I think some of the 
  //   RBJ-based biquads also currently use the other convention. Maybe there are more.
}

void combVsAllpassPhase()
{
  // Under construction...

  // We plot the phase response, group delay and ring response for a comb and Schroeder allpass
  // with the same feedback coefficient
}

void universalCombResponses()
{
  using Real = double;
  using Vec  = std::vector<Real>;
  using Comb = rsUniversalCombFilter<Real, Real>;

  // Helper functions (ToDo: change order of params to bl, ff, fb and update all call sites):

  auto setupComb = [](Comb& comb, int delay, Real bl, Real ff, Real fb)
  {
    comb.setMaxDelayInSamples(delay);
    comb.setDelayInSamples(delay);
    comb.setCoeffs(bl, ff, fb);
  };

  auto plotImpResp = [&](int delay, Real bl, Real ff, Real fb, int numSamples)
  {
    Comb comb;
    setupComb(comb, delay, bl, ff, fb);
    Vec h = impulseResponse(comb, numSamples, 1.0);
    rsPlotVector(h);
  };

  auto plotFreqResp = [&](int delay, Real bl, Real ff, Real fb, int numBins)
  {
    Comb comb;
    setupComb(comb, delay, bl, ff, fb);
    plotFrequencyResponse(comb, numBins, 0.0, 0.5, 1.0, false);
  };

  auto plotImpAndFreqResp = [&](int delay, Real bl, Real ff, Real fb, int numSamples, int numBins)
  {
    plotImpResp( delay, bl, ff, fb, numSamples);
    plotFreqResp(delay, bl, ff, fb, numBins);
  };


  // Plot impulse and frequency responses of the universal comb filter for various settings:

  int M = 10;   // Delay to be used in all plots below.


  // Feedback combs:
  //                     BL    FF    FB 
  plotImpAndFreqResp(M, +1.0, +0.0, -0.9, 202, 2001);
  // Feedback comb with negative feedback coeff. Because the coeff is negated again in the 
  // implementation (for consistency with the sign convention used in most of the DSP literature),
  // this produces a full series of harmonics.
  // Impulse:   Decaying upward spikes starting with unit amplitude at n = 0.
  // Magnitude: Peaks at all multiples of 1/M with height +20 dB. Rounded valleys at -5.57 dB.
  // Phase:     Rounded upward saw. Start phase is middle of the (downward) edge.

  plotImpAndFreqResp(M,  0.0, +1.0, -0.9, 202, 2001);   
  // Gives same result as above but with an added predelay of M samples. This superimposes a 
  // linear phase over the phase response of the filter above.
  // Impulse:   Decaying upward spikes starting with unit amplitude at n = M.
  // Magnitude: Peaks at all multiples of 1/M with height is +20 dB.
  // Phase:     Slanted staircase with rounded steps. I think, it's the phase response from above
  //            with a linear downward trend added.
 

  // Feedforward combs:
  //                     BL    FF    FB  
  plotImpAndFreqResp(M, +1.0, +1.0, +0.0, 30, 2001);
  // Feedforward comb with unit weights for input and for delayed signal.
  // Impulse:   Two spikes at n = 0 and n = M both with height 1.0.
  // Magnitude: Notches at odd multiples 0.5/M with infinite depth.
  // Phase:     Downward saw from +90° to -90°. Start phase is middle of the ramp, i.e. 0°.

  plotImpAndFreqResp(M, +1.0, +0.9, +0.0, 30, 2001);
  // Feedforward comb with unit weight for input and weight 0.9 for delayed signal.
  // Impulse:   Two spikes at n = 0 and n = M with height 1.0 and 0.9.
  // Magnitude: Troughs at odd multiples 0.5/M with depth -20 dB.
  // Phase:     Rounded downward saw. Start phase is middle of the ramp.

  plotImpAndFreqResp(M, +0.9, +1.0, +0.0, 30, 2001);
  // Feedforward comb with weight 0.9 for input and and unit weight for delayed signal. This swaps
  // the weights of the filter above. The result is a time reversed impulse response.
  // Impulse:   Two spikes at n = 0 and n = M with height 0.9 and 1.0.
  // Magnitude: Troughs at odd multiples 0.5/M. with depth -20 dB.
  // Phase:     Slanted staircase with rounded steps. It's a linear downward trend superimposed 
  //            with a rounded upward saw, I think.

  plotImpAndFreqResp(M, +1.0, -0.9, +0.0, 30, 2001);
  // Feedforward comb with unit weight for input and weight -0.9 for delayed signal.
  // Impulse:   Two spikes at n = 0 and n = M with height 1.0 and -0.9.
  // Magnitude: Notches at all multiples 1/M with depth -20 dB. Rounded maxima at height +5.57 dB.
  // Phase:     Rounded downward saw. Start phase is middle of the edge of the ramp, i.e. +0° but
  //            quickly rising up to +90°.

  // A couple of more feedforward cases without comments:
  //plotImpAndFreqResp(M, -0.9, +1.0, +0.0, 30, 2001);
  //plotImpAndFreqResp(M, +1.0, +1.0, +0.0, 30, 2001);
  //plotImpAndFreqResp(M, -1.0, -1.0, +0.0, 30, 2001);
  //plotImpAndFreqResp(M, +1.0, -1.0, +0.0, 30, 2001);
  //plotImpAndFreqResp(M, -1.0, +1.0, +0.0, 30, 2001);
 



  // Schroeder allpasses:
  //                     BL    FF    FB  
  plotImpAndFreqResp(M, +0.9, +1.0, +0.9, 202, 2001);
  // Allpass with positive (negated) feedback.
  // Impulse:   Initial upward spike at n = 0 with height 0.9 followed by alternating spikes. The 
  //            2nd spike goes upward, too. The 3rd goes downward, etc.
  // Magnitude: Flat at 0 dB, i.e. allpass - as it should be.
  // Phase:     Stairsteps at odd multiples 0.5/M. The plateaus between the steps are at integer
  //            multiples of 360°.

  plotImpAndFreqResp(M, -0.9, +1.0, -0.9, 202, 2001);
  // Allpass with negative (negated) feedback.



  // Blesser notchpasses:
  plotImpAndFreqResp(M, +0.9, +1.0, -0.8, 202, 2001);
  plotImpAndFreqResp(M, -0.9, +1.0, +0.8, 202, 2001);
  // Looks qualitatively right but I'm not sure, if the depths of the magnitude notches has the 
  // correct relationship with the group delay peaks (i.e. phase response cliffs). Also: the 
  // magnitude goes above unity between the peaks. When such a notchpass is part of a feedback 
  // loop, we may want to scale the input or output such that the magnitude always maxes out at 
  // unity, i.e. at 0 dB. Figure out if input or output scaling behaves better under modulation.
  // Maybe try also to implement the filter in different forms (like DF1/2, TDF1/2) and check their
  // behavior under modulation.


  // Perfect cancellation between feedforward and feedback path with FF = FB:
  plotImpAndFreqResp(M, +1.0, +0.9, +0.9, 30, 2001); // is unit impulse
  plotImpAndFreqResp(M, +1.0, -0.9, -0.9, 30, 2001); // dito
  // Apparently, when FF = FB, we get a perfect cancellation between feedback and feedforward path?


  plotImpAndFreqResp(M, +1.0, +0.9, -0.9, 202, 2001);
  plotImpAndFreqResp(M, +1.0, -0.9, +0.9, 202, 2001);
  // Phase responses are "square waves", mag responses are combish



  plotImpAndFreqResp(M, +0.8, +1.0, -0.9, 202, 2001);
  plotImpAndFreqResp(M, -0.8, +1.0, +0.9, 202, 2001);


  // Compare feedback comb with Schroeder allpass with same feedback:
  plotImpAndFreqResp(M, +1.0, +0.0, +0.9, 202, 2001);  // Comb
  plotImpAndFreqResp(M, +0.9, +1.0, +0.9, 202, 2001);  // Allpass
  plotImpAndFreqResp(M, +1.0, +0.0, -0.9, 202, 2001);  // Comb - maybe try using -1 for blend
  plotImpAndFreqResp(M, -0.9, +1.0, -0.9, 202, 2001);  // Allpass
  // The combs have a sawtooth like phase response, the allpasses a staircase like one.
  // ToDo: plot the group delay and ringing of both.

  // ...

  // Combs with freely assigned coeffs:
  plotImpAndFreqResp(10, +0.5, +1.0, +0.8, 202, 2001);
  plotImpAndFreqResp(10, +0.5, +0.9, +0.8, 202, 2001);
  plotImpAndFreqResp(10, +0.5, -0.9, +0.8, 202, 2001);


  int dummy = 0;

  // ToDo:
  //
  // - Make a ToolChain module UniComb that has these 4 parameters and plots the magnitude and 
  //   phase response such that we can conveniently experiment with the effect of different 
  //   parameters
  //
  // - Maybe make a pole/zero plot. This requires to produce the transfer function as pair of 
  //   rsPolynomial and to find the roots of these polynomials. Maybe we can use class 
  //   FilterPlotter for this (after producing the filter specification in terms of B,A)
  //
  // - Plot ringing responses.
  //
  // - Make IIR combs with different gain normalizations (see DAFX pg. 70). We can apply an overall
  //   gain factor by just scaling the ff and bl coeff by that factor
  //
  // - Compare phase and group delay responses of feedback comb with allpass that has the same 
  //   feedback coeff.
  //
  // - Try to set it up in terms of the first 3 samples of the impulse response (when delay = 1)
  //   From y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1]  with  x = 1,0,0,0,0,... we should get: 
  //   h[0] = b0, 
  //   h[1] = b0*0 + b1*1 - a1*b0 = b1 - a1*b0
  //   h[2] = b0*0 + b1*0 - a1*(b1 - a1*b0) = -a1*b1 + a1*a1*b0 = a1*(a1*b0 - b1)
  //   so b0 = h[0] immediately and b1,a1 can be computed from a quadratic system of equations, I 
  //   think. Maybe that means, we have a choice? But maybe not when we impose a stability 
  //   constraint? ...work it out! Maybe being able to ste the first 3 samples of the impulse 
  //   response directly can help to fine tune the early reflections. ...but we may have to undo
  //   the frequency response of the filter in the feedback path 
}

/** Converts a reverberation time to the corresponding time constant tau of an exponential decay
function f(t) = exp(-t/tau) that decays to a certain level (levelToReach) after reverbTime seconds
or milliseconds or samples or whatever. The physical unit doesn't matter. The output will be in the
same unit as the input. The desired target level to reach is given as a raw amplitude and it 
defaults to 1/1000 which corresponds to -60 dB which is a standard value used in audio engineering 
and acoustics to describe a reverb time. This time to decay down to -60 dB is also known as RT60 
(RT for reverberation time). The time constant tau of an exponential decay, on the other hand, is 
defined to be the time to decay down to 1/e. The formula to compute tau from the revrb time is 
given by:  tau = -reverbTime / ln(levelToReach) which is what this function encapsulates. */
template<class T>
T rsReverbTimeToTau(T reverbTime, T levelToReach = T(0.001))
{
  return -reverbTime / rsLog(levelToReach);
}
// Needs tests. If it works, move it into the library into the file AudioFunctions.h near the
// function rsDecayTimeToFeedbackGain(). OK - it seems to work well in combVsModalBank().


template<class TSig, class TPar>
void rsSetModalBankToComb(rsModalFilterBank<TSig, TPar>* mfb, int M, 
  bool odd = false, int mMin = 0, int mMax = -1)
{
  // Default values must be known at compile time so we can't use something like mMax = M in the
  // function signature. Instead, we use mMax = -1 as code to indicate that mMax should be M.
  if(mMax == -1)
    mMax = M;

  // Prepare the vectors of the modal parameters:
  std::vector<TPar> frq(M+1), amp(M+1), att(M+1), dec(M+1), phs(M+1);
  for(int m = 0; m <= M; m++)
  {
    frq[m] = m;    // Frequency relative to the fundamental
    amp[m] = 1.0;  // Linear amplitude
    att[m] = 0.0;  // Relative attack time (i.e. location of the peak)
    dec[m] = 1.0;  // Relative decay time (i.e. time to decay dwon to 1/e)
    phs[m] = 90;   // Start phase in degrees. 90 seems correct when odd == false

    // Apply the brickwall filtering to the modes:
    if(m < mMin || m > mMax) 
      amp[m] = 0.0;

    // Adjust the amplitude of the DC and Nyquist modes (see comments below why):
    if(m == 0 || m == M)
      amp[m] *= 0.5;
  }
  
  mfb->setModalParameters(frq, amp, att, dec, phs);

  // Notes:
  // 
  // - The reason for multiplying the amplitudes for m = 0 (DC) and m = M (Nyquist) by one half is
  //   the observation that otherwise, the impulse response of the modal bank will feature an 
  //   oscillation at the Nyquist frequency which is not present in an actual feedback comb filter
  //   implemented the standard way using a delay line. I think it may be because in an actual 
  //   feedback comb filter, the DC component is represented to one half by the actual DC frequency
  //   and the to the other half by the Nyquist frequency and if, in the modal filter, we give both
  //   of these modes a gain of 1, we end up boosting DC and fs/2 by a factor of two compared to 
  //   what it should be. That's my hypothesis at least. ToDo: Verify it theoretically! By the way:
  //   Just setting the gain at m = M and/or m = 0 to zero doesn't fix it. We really need both 
  //   components with gain 0.5 to make it right.
}

void combVsModalBank()
{
  // We create comparative plots of the impulse- and frequency responses of a feedback comb filter
  // and a bank of modal filters. We want to adjust the modes of the modal bank in such a way as 
  // to match the modes of the comb in terms of frequency, amplitude, decay time and phase. One 
  // goal is to figure out how to correctly set up the start phases of the modal  filters. 
  // ...TBC...

  // Types:
  using Real = double;
  using Vec  = std::vector<Real>;
  using UCF  = rsUniversalCombFilter<Real, Real>;
  using MFB  = rsModalFilterBank<Real, Real>;

  // Setup:
  int  M          = 100;      // Number of modes. Determines fundamental and delay length
  Real RT60       = 4000;     // Number of samples to decay to -60 dB
  bool odd        = false;    // If true, we produce only odd harmonics
  int  mMin       = 0;        // Lowest mode to produce. 0 is DC, 1 the fundamental
  int  mMax       = M/1;      // Highest mode to produce in the modal bank
  int  numBins    = 2001;     // Number of bins for frequency response plots
  int  numSamples = 1000;     // Number of samples for impulse response plots

  // Create and set up the feedback comb filter:
  int delay = 2*M;
  Real fb = rsDecayTimeToFeedbackGain(RT60, Real(delay), 0.001);  // 0.001 is -60 dB
  if(!odd)
    fb = -fb;
  UCF comb;
  comb.setMaxDelayInSamples(delay);
  comb.setDelayInSamples(delay);
  comb.setToFeedbackComb(fb);

  /*
  // Prepare the vectors of the modal parameters:
  Vec frq(M+1), amp(M+1), att(M+1), dec(M+1), phs(M+1);
  for(int m = 0; m <= M; m++)
  {
    frq[m] = m;    // Frequency relative to the fundamental
    amp[m] = 1.0;  // Linear amplitude
    att[m] = 0.0;  // Relative attack time (i.e. location of the peak)
    dec[m] = 1.0;  // Relative decay time (i.e. time to decay dwon to 1/e)
    phs[m] = 90;   // Start phase in degrees. 90 seems correct when odd == false

    // Apply the brickwall filtering to the modes:
    if(m < mMin || m > mMax) 
      amp[m] = 0.0;

    // Adjust the amplitude of the DC and Nyquist modes (see comments below why):
    if(m == 0 || m == M)
      amp[m] *= 0.5;
  }
  */

  // Create and set up the bank of modal filters:
  Real f0  = 0.5/M;                    // Fundamental frequency
  Real tau = rsReverbTimeToTau(RT60);  // Decay time constant 
  MFB mfb;
  mfb.setSampleRate(1.0);
  mfb.setReferenceFrequency(f0);
  mfb.setReferenceDecay(tau);
  rsSetModalBankToComb(&mfb, M, odd, mMin, mMax);
  //mfb.setModalParameters(frq, amp, att, dec, phs);  // Old

  // Produce the impulse responses:
  Vec hc = impulseResponse(comb, numSamples, 1.0);
  Vec hm = impulseResponse(mfb,  numSamples, 1.0);

  // Plot impulse- and frequency responses:
  Real scl = 1.0 / (mMax-mMin+1);      // Scale factor to obtain unit amplitude
  rsPlotVectors(hc, scl*hm); 
  plotFrequencyResponse(comb, numBins, 0.0, 0.5, 1.0, false);
  plotFrequencyResponse(mfb,  numBins, 0.0, 0.5, 1.0, false);


  // Observations:
  // 
  // - With M = 50, mMin = 0, mMax = M, odd = false, the modal response looks very similar to the 
  //   comb response but it shows tiny ripples at the Nyquist frequency. It looks like they ripple 
  //   roughly between 0 and 0.02 (which is 1/M). When we set mMin = 1, i.e. start at the 
  //   fundamental rather that at DC, they ripple between 0 and -0.02. Maybe this is an artifact
  //   relating to a phase error in the topmost harmonic? The Nyquist frequency is generally 
  //   problematic with regard to phase because a sinusoid at exactly that frequency may either be
  //   captured faithfully (when the sample instant align with the minima and maxima) or it may be 
  //   not captured at all (when the sample instants align with the zero crossings) or anything in 
  //   between (in which case it may have the wrong amplitude).
  // 
  // - When setting mMax to something like M/5 or M/10, the modal output looks like a lowpassed 
  //   version of the impulse train that the comb produces. This is exactly what we expect. When
  //   using scl = 1.0 / (mMax-mMin+1), the heights of the peaks will exactly match those of the
  //   comb. Note however, that adjusting for the height of the peaks like that is actually "wrong"
  //   in the sense that we really modify the amplitude of the lower modes in order to compensate
  //   for cutting out the higher modes. It's good for the impulse response plot, though. It may be
  //   bad for the frequency response plot, though. Verify that! Soo - maybe it would actually 
  //   better to not bake the scaling factor already into the signal but instead apply it only when
  //   plotting the impulse response. ...done
  //
  // - It doesn't make a difference if we use mMax = M or mMax = M-1. The mode with m = M seems to
  //   be an all zeros signal. That can be verified by chossing mMin = mMax = M such that only that
  //   mode is produced. The result is indeed all zeros. This was with M = 5. Seems to be the same 
  //   with other values of M.
  //
  // - With M = 10, mMin = 0, mMax = M, the frequency responses look very similar but the modal 
  //   bank shows some differences near DC and sampleRate/2. ToDo: Try to multiply the amplitudes
  //   of those modes (0 and M) by 1/2. I think, these modes are just too loud. Maybe it is because
  //   in the comb, the DC component is represented to one half by the actual DC frequency and the 
  //   to the other half by the Nyquist frequency and if, in the modal filter, we give both of 
  //   these modes a gain of 1, we end up boosting DC and fs/2 by a factor of two compared to what 
  //   it should be?
  //   OK - done! And it does indeed work! We have now an exact match. The only problem is that now
  //   when we choose mMax < M, the computation of the scl factor is now not exact anymore. ToDo:
  //   fix this! I think, we should replace scl = 1.0 / (mMax-mMin+1); by 
  //   scl = 1.0 / (mMax-mMin+0.5); or something like that. ...but maybe only when mMax < M and 
  //   otherwise keep the formula as is. Figure this out by trial and error! It's not really 
  //   important (because it only affects our plots here and nothing of it relates to any DSP that
  //   we would do in a production context) but it would be nice to have anyway.
  // 
  //
  // ToDo:
  // 
  // - Plot the frequency responses of the comb and modal bank together in a single plot to spot 
  //   the differences more easily. Maybe for this, we need to write a new plotting function.
  //
  // - Try it with odd = true. Maybe we have to adjust the number of modes and maybe the phases, 
  //   too. Oh - and the computation of the modal frequencies, too. We will produce only odd 
  //   harmonics! I think, that also means that there will be no DC component.
  // 
  // - Try it with different comb/allpass settings in the universal comb. Figure out how that 
  //   affects the required modal parameters. I think, frequencies, amplitudes and decay times 
  //   must remain the same so the only thing to adjust is the start phase.
  //
  // - Document in class rsModalFilter that "DC modes" are also allowed.
  // 
  // - Maybe instead of just setting the modes below mMin and above mMax to zero amplitude, we 
  //   should actually shorten the frq, amp, etc. vectors accordingly, i.e. to lengths mMax-mMin+1.
  //
  // - Maybe instead of baking the decay times and amplitudes into the individual mode parameters,
  //   assign them to relative values and use tau in mfb.setReferenceDecay() and scl in 
  //   mfb.setReferenceAmplitude(). The latter doesn't exist yet, I think -> add it.
  //
  // - Maybe write a function (or class) that encapsulates setting up a modal bank from (universal)
  //   comb parameters to simulate the comb exactly by the modal bank. A class could be named
  //   rsModalCombSimulator and a function could be named setupModalBankForCombSimulation. Maybe
  //   it could take a const reference to the actual universal comb. Or maybe it could just work
  //   with the comb parameters.
  //
  // - But maybe it's not possible to model the zeros of the universal comb in e.g. allpass mode 
  //   or feedforward mode. This is because the comb's zeros are complex and a regular modal filter 
  //   is usually realized as two-pole/one-zero filter and the zero is therefore restricted to be 
  //   real valued. In order to cover the full range of the universal comb's possibilities, we may
  //   have to use a full biquad filter for each mode. 
}

void delayLines()
{
  // Delaylines with different interpolation methods:
  //delayLineBasic();                // No interpolation
  //delayLineLinear();               // Linear interpolation
  //delayLineAllpass();              // Allpass interpolation

  // Other delayline based stuff:
  //universalCombVsOnePole();
  //combVsAllpassPhase();              // Under construction
  //universalCombResponses();          // Imp- and freq-responses of uniCombs with various settings
  combVsModalBank();                   // Compare comb with matched modal bank
}

//=================================================================================================

template<class T>
void rsSetupTwoPoleAllpassDelay(rsTwoPoleAllpassDelay<T, T>* filter, int delay, T omega, T Q)
{
  T b0, b1, b2, a1, a2;

  rsStateVariableFilter<T, T> svf;
  svf.setupAllpass(omega, Q);
  svf.convertToBiquad(&b0, &b1, &b2, &a1, &a2);

  filter->setMaxDelayInSamples(delay);
  filter->setDelayInSamples(   delay);
  filter->setAllpassCoeffs(a1, a2);

  // ToDo:
  //
  // - It's a bit silly to use the state variable filter for the purpose of designing biquad 
  //   coeffs so maybe replace this code later with the direct RBJ biquad design formulas.
}

void twoPoleAllpassDelayResponses()
{
  // Under construction.

  // We plot various kinds of frequency responses for filters of the rsTwoPoleAllpassDelay type.

  using Real = double;
  using APF  = rsTwoPoleAllpassDelay<Real, Real>;

  // User parameters:
  int  numBins = 5000;
  int  delay   = 20;
  Real omega   = 0.5;
  //Real omega   = 1.57;
  Real Q       = 2.5; 

  // Create the filter:
  APF flt;
  rsSetupTwoPoleAllpassDelay(&flt, delay, omega, Q);

  // Plot its frequency response:
  plotFrequencyResponse( flt, numBins, 0.0, 24000.0, 48000.0, false       );
  plotMagAndRingResponse(flt, numBins, 0.0, 24000.0, 48000.0, false, false);



  // Observations:
  //
  // - There are stairsteps with a narrower mini-step in the middle of the verticalish part.
  //
  // - With lower Q, the steps get washed out more and more.
  //
  // - Omega seems to control the width of the mini-step. Actually, it was "mini" only because I 
  //   tried first with a rather low omega of 0.3. With 1.3, the picture is quite different. Looks
  //   like with omgega = pi/2, the width of the "mini-steps" is excatly the same as the width of 
  //   the primary steps.
  //
  // - With numBins = 2000, delay = 20, omega = 0.3, Q = 9.5, we see some artifacts in the plot.
  //   Maybe these are errors in the phase-unwrapping algorithm? Figure this out!
  //
  // - The ringing response shows pairwise spikes. The distance between the two spikes is 
  //   controlled by omega. If omega = pi/2, the distance between the two spikes of a pair is such
  //   that all spikes become equidistant, i.e. there's no apparent pairing visible anymore. But 
  //   what about the effect of Q? It seems to affect the spike frequencies!
}

void twoPoleAllpassDelayChains()
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
    for(size_t i = 0; i < numStages; i++)
      rsSetupTwoPoleAllpassDelay(&filters[i], delays[i], omegas[i], Qs[i]);

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


  // Show plots for some chains of 2-pole allpass delays. The 1st list gives the delays, the 2nd
  // gives the omegas and the 3rd gives the Qs:

  // With delay = 1, the rsTwoPoleAllpassDelay reduces to a normal biquad allpass:
  plot({ 1   }, 
       { 0.2 }, 
       { 2.5 }, 100);

  // Now we replace the unit delay by a delay of M = 5 samples in the biquad allpass:
  plot({ 5   }, 
       { 0.2 }, 
       { 2.5 }, 500);


  // Let's now try some series connections of two such filters:

  // This gives a nice diffusor of length around 100-150:
  plot({ 5,   9   }, 
       { 0.5, 0.7 }, 
       { 1.5, 1.5 }, 500);

  plot({ 5,   9   }, 
       { 0.2, 0.3 }, 
       { 2.5, 2.5 }, 500);
   
  plot({ 5,   9   }, 
       { 0.2, 0.3 }, 
       { 1.0, 1.5 }, 500);


  // Let's now try some series connections of three such filters:

  plot({ 5,   9,   14  }, 
       { 0.2, 0.3, 0.5 }, 
       { 2.5, 9.5, 2.5 }, 1000);

  plot({ 5,   9,   14  }, 
       { 0.2, 0.3, 0.5 }, 
       { 1.0, 1.0, 1.0 }, 1000);



  // Render some wave files for audition:

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
  //   factor. -> Figure this out!
  //
  // - Higher Q values seem to emphasize the initial spike. The decaying part gets longer (i.e. 
  //   decays more slowly) but its overall level relative to the initial spike gets smaller. For 
  //   diffusion, the sweet spot seem to be somewhere in the range 0.5...1.5. Maybe longer delays
  //   should use more or less Q? Make more tests with 2 filters in series to figure this out.
  //   I think, that - counterintuitively - the filter with longer delays should use higher Q. 
  //   Maybe test filters with different delays in isolation and look at the ratio of the initial
  //   spike to the second. I think, we wnat to keep that ratio constant as function of the delay,
  //   i.e. keep it independent from the delay.
  //
  // - If the frequency is zero, the impulse response is just a unit impulse.
  //
  //
  // Questions:
  //
  // - I think, we do not really need the delay lengths to be all prime numbers. It should be 
  //   sufficient if they are mutually coprime. This is a less restrictive condition but the effect
  //   having the first coincidence of spikes at the lowest common multiple of the delay lengths 
  //   should still be satisfied.
  //
  // - Can we use the undulation frequencies to deliberately colorize the sound, i.e. give it some
  //   deliberate tonal character? That seems plausible. Maybe we could use slightly detuned
  //   omegas for left and right channel.
  //
  // - Maybe the undulation frequencies can be set up in such a way that the sine goes through a 
  //   zero whenever the spikes meet? For example, int the chain with 2 filters with 5 and 9 
  //   samples delay, the spikes would meet at n = 45 for the first time. But there is no spike in 
  //   the impulse response. I gues it's because either the 1st or the 2nd filter goes through an 
  //   undulation zero at that point? Figure that out!
  //
  //
  // ToDo:
  //
  // - Try to achieve sign flipping of the output by flipping the signs of the omegas. This doesn't
  //   fall out of the math. Trying to just use negative omegas (or negative Qs) just leads to 
  //   unstable filters. It would just be a convention that we apply ourselves manually to allow 
  //   the user to conveniently select the sign of a particular filter output. It would work by 
  //   using abs(w) for the filter design and sign(w) to scale the output. But maybe it's clearer
  //   to just introduce another (boolean) parameter for this. Not sure...
  //
  // - Try to parametrize the filter not in terms of omega but in terms of physical frequency. We 
  //   really want to figure out the relationship between the omega parameter and the undulation 
  //   frequency. Maybe  omega*delay = 2*pi*f/fs  or  omega/delay = 2*pi*f/fs?
}

void twoPoleAllpassDelays()
{
  twoPoleAllpassDelayResponses();  // Under construction
  twoPoleAllpassDelayChains();
}


void dampedCombFilter()
{
  // Under construction. rsDampedCombFilter should factor out the comb filter from the class 
  // rsDampedCombAllpass

  using Real = double;
  using Comb = rsDampedCombFilter<Real, Real, Real, Real>;



  // Create and set up the comb filter object:
  Comb comb;


  int dummy = 0;
}

void dampedCombAllpassIdea()
{
  // This experiment is basically my initial step by step derivation of what later became the class
  // rsDampedCombAllpassNaive. It's based on the idea that I explain in the document:
  // Notes/DSP/DampedCombAllpass.txt
  //
  // I implement an idea for starting with an arbitrary given allpass filter A(z) and arbitrary 
  // given feedback filter F(z) that sits in a feedback loop with unit delay around that allpass. I 
  // design a compensation filter that can be applied in series to this setup such that the overall
  // transfer function is allpass in nature. Without the compensation filter, this setup has the 
  // uncompensated transfer function:
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
  using Delay   = RAPT::rsDelay<Real>;      // We use a simple delay as allpass
  using OnePole = RAPT::rsOnePoleFilter<Real, Real>; // We use a one pole as feedback filter

  int  M          =   100;     // Delay
  int  N          =  8192;     // Number of samples to generate
  Real sampleRate = 44100;     // Sample rate for writing the wavefiles
  Real dampFreq   =  1000;     // Frequency of the low shelf for feedback damping
  Real dampGain   =     0.7;   // Linear high freq damping gain
  Real k          =     0.99;  // Feedback gain factor


  // Create and set up the two given filters for A(z) and F(z):
  Delay apf;                                // Allpass filter
  apf.setMaxDelayInSamples(M);
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
  dlM1.setMaxDelayInSamples(M+1);
  dlM1.setDelayInSamples(M+1);
  Delay dlM2;
  dlM2.setMaxDelayInSamples(M+2);
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
  rosic::writeToMonoWaveFile("FeedbackAllpassUncorrected.wav", &hu[0], N, (int)sampleRate);
  rosic::writeToMonoWaveFile("FeedbackAllpassCorrected.wav",   &hr[0], N, (int)sampleRate);

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

void dampedCombAllpassClass()
{
  // Now, we use the class rsDampedCombAllpass which encapsulates the algorithm derived above into a 
  // class. We use that class here to test some different options for the sign of the feedback gain
  // and the swapping of damper and delay. This leads to 4 different modes that all have different
  // properties and sound different.

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedCombAllpass<Real, Real, Real, Real>;

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

/** A structure to hold settings of a frequency response plot. */
struct rsFrequencyResponsePlotSettings
{
  // This shall eventually be moved somewhere else. Maybe into Plotting.h or DSPPlotters.h.

  /** Enumeration of the possible types of plots. */
  enum class Type
  {
    magnitude,          // Linear magnitude
    decibels,           // Magnitude in decibels
    phase,              // Phase in degrees
    decibelsAndPhase,   // Decibels together with phase
    realAndImag,        // Real and imaginary part
    phaseDelay,         // Phase delay
    groupDelay,         // Group delay
    decibelsAndRing     // Decibels together with ringing response.

    // ...TBC...maybe. Maybe phase together with group delay, group delay with ringing, ...
    // Whatever combinations seems to make sense.
  };

  // The settings:
  double minFreq =    20.0;         // Minimum frequency
  double maxFreq = 20000.0;         // Maximum frequency
  double sampRat = 44100.0;         // Sample rate (if applicable)
  int    numBins =  1000;           // Number of bins (i.e. samples along freq axis)
  Type   type    = Type::decibels;  // Type of response plot
  bool   logFreq = true;            // Use logarithmic frequency axis (or not)

};

// ToDo:
//
// - Maybe make similar structs for time responses and s- or z-domain plots. Types of time 
//   responses could be: impulse response, step response, responses to various waveforms like 
//   sawtooth, square, etc. Maybe step-response minus the step iteself (not sure, if that's 
//   useful, though). Types of z-plane or s-plane plots: poles and zeros, magnitude as 3D plot or 
//   heat map (should be clipped to some adjustable maximum at the poles, I guess), color coded 
//   phase plot, Polya potential, etc.


template<class T, class TFunc>
void rsPlotThreeFilterResponses(
  const TFunc& tfComb, const TFunc& tfCorr, const TFunc& tfFull, 
  const rsFrequencyResponsePlotSettings& plotSetup, T dummy)
{
  // The dummy is need for the compiler to infer the type T. Maybe later we can replace it by some
  // actually useful parameter of type T.

  // Retrieve the relevant plot settings:
  T    fs      = plotSetup.sampRat;
  int  N       = plotSetup.numBins;
  bool logFreq = plotSetup.logFreq; 
  T    fMin    = plotSetup.minFreq;
  T    fMax    = plotSetup.maxFreq;

  // Create the desired type of plot for the comb, corrector and full filter:
  using PT = rsFrequencyResponsePlotSettings::Type;  // Plot type
  switch(plotSetup.type)
  {
  case PT::decibelsAndPhase: {
    plotFreqRespFromTransFunc(tfComb, N, fMin, fMax, fs, logFreq);
    plotFreqRespFromTransFunc(tfCorr, N, fMin, fMax, fs, logFreq);
    plotFreqRespFromTransFunc(tfFull, N, fMin, fMax, fs, logFreq); }  break;

  case PT::decibelsAndRing: {
    plotMagAndRingRespFromTransFunc(tfComb, N, fMin, fMax, fs, logFreq, false);
    plotMagAndRingRespFromTransFunc(tfCorr, N, fMin, fMax, fs, logFreq, false);
    plotMagAndRingRespFromTransFunc(tfFull, N, fMin, fMax, fs, logFreq, false); }  break;

    // ...TBC...
  }


  // ToDo:
  //
  // - Plot also group delay and phase delay response. Maybe let the caller switch between the kind
  //   of response plot. Maybe make an enum rsFilterPlotType and let this function take a parameter
  //   of this enum to switch between different types of plots.
  //
  // - Maybe let the caller also set up the logFreq setting, min/max freqs, etc. Maybe it's worth
  //   to define a struct rsFilterPlotSettings that contains all these such that the caller doesn't
  //   have to write them all out in every call.
  //
  // - Look into how plotFrequencyResponse() determines the spacing of the ticks for the phase 
  //   response this is not yet how it should be.
  //
  // - Maybe factor out a function rsPlotResponse(TFunc transFunc, 
  //   const rsFrequencyResponsePlotSettings& plotSetup, T dummy). Then we can use it also fo 
  //   ad-hoc structures made from allpasses, provided that we define the appropriate transfer 
  //   function computation function.
}

/** Plots various responses of the given comb-allpass filter. */
template<class TFlt, class T>
void rsPlotDampedCombAllpassResponses(
  TFlt& filter, const rsFrequencyResponsePlotSettings& plotSetup, T dummy)
{
  // Create function objects for the various transfer functions:
  using Complex = std::complex<T>;
  using Func    = std::function<Complex(Complex)>;
  Func tfComb = [&](Complex z) { return filter.getCombTransferFunctionAt(z);      };
  Func tfCorr = [&](Complex z) { return filter.getCorrectorTransferFunctionAt(z); };
  Func tfFull = [&](Complex z) { return filter.getTransferFunctionAt(z);          };

  // Plot the 3 transfer functions, one after another:
  rsPlotThreeFilterResponses(tfComb, tfCorr, tfFull, plotSetup, dummy);
}

// Convenience function that doesn't need a plotsetup parameter. 
template<class TFlt, class T>
void rsPlotDampedCombAllpassResponses(TFlt& filter, rsFrequencyResponsePlotSettings::Type type,
  T sampleRate)
{
  // Create and set up the plot settings object:
  rsFrequencyResponsePlotSettings plotSetup;
  plotSetup.logFreq = false;
  plotSetup.minFreq = 0.0;
  plotSetup.maxFreq = 0.5*sampleRate;
  plotSetup.sampRat = sampleRate;
  plotSetup.numBins = (int) (sampleRate/10);  // Not sure about that
  plotSetup.type    = type;

  // Create the plot:
  rsPlotDampedCombAllpassResponses(filter, plotSetup, T(0));
}

void dampedCombAllpassResponses()
{
  // We plot the responses magnitude, phase, phase delay, group delay, ringing, etc. for a damped
  // comb allpass filter and its underlying comb and corrector.

  // Define types to be used:
  using Real     = double;
  using Vec      = std::vector<Real>;
  using Allpass  = rsDampedCombAllpass<Real, Real, Real, Real>;
  using PlotType = rsFrequencyResponsePlotSettings::Type;

  // User parameters:
  int  delay      =    20;     // Main delay roundtrip length in samples. Is M-1 in the algo
  int  numSamples =  4410;     // Number of samples to generate  ...rename to numBins
  Real sampleRate = 44100;     // Sample rate for writing the wavefiles
  Real dampFreq   = 10000;     // Frequency (in Hz) of the shelf filter for feedback damping
  Real dampGain   =     0.5;   // Linear high freq damping gain
  Real feedback   =    -0.90;  // Feedback gain factor

  // Create and set up the damped comb allpass filter:
  Allpass flt;
  flt.setMaxDelayInSamples(delay);
  Real dampOmega = 2*PI*dampFreq/sampleRate;
  rsSetupHighDamp(flt, delay, feedback, dampOmega, dampGain, false);

  // Helper function to create a plot:
  auto plotResponse = [&](PlotType type)
  { 
    rsPlotDampedCombAllpassResponses(flt, type, sampleRate); 
  };

  // Create the different types of plots:
  plotResponse(PlotType::decibelsAndPhase);
  plotResponse(PlotType::decibelsAndRing);
  // ...TBC... ToDo: plot phase delay and group delay. For that, plotResponse needs to be extended.


  // Observations:
  //
  // - The phase response of the full allpass filter features steep slopes at a harmonic series 
  //   whose fundamental is determined by the delay. The phase response looks like a downward 
  //   staircase. The general linear trend seems to come mostly from the corrector. The comb 
  //   contributes a sawtooth shaped phase response. Added together, the result in the staircase.
  //
  // - The stairsteps get steeper when the (absolute value of the) feedback is higher. When the 
  //   dampGain is less than 1, the steps become smoother to the right (i.e. at higher 
  //   frequencies).
  //
  // - The ringing response shows clear spikes at the ringing frequencies, i.e. at the stair 
  //   steps of the phase response. When dampGain = 1, the spikes are all of the same width and 
  //   height. With dampGain < 1, they become wider and less tall to the right.
  //
  // - For delay = 20, feedback = -0.95, dampGain = 1, we see a ringing resonance at the Nyquist 
  //   freq and also at DC. For feedback = +0.95, there is no such thing.
  //
  //
  // ToDo:
  //
  // - Try a series of an allpass with positive and negative feedback but otherwise the same 
  //   settings. The ringing frequencies of both filters should interleave. But maybe that's just 
  //   equivalent ot a single comb-allpass with lower (fundamental) ringing frequency? Because if 
  //   we superimpose ringing at 100,200,300,... with 50,150,250,... we get a harmonic series with
  //   a fundamental of 50.
  //
  // - Try to put an allpass in series that attempts to counteract the ringing somehow. It should 
  //   smooth out the stairsteps. But maybe if such a compensation allpass is perfect (i.e. 
  //   smoothes out the straisteps perfectly), we just end up with an allpass with linear phase 
  //   response, i.e. a pure delay? That would be boring!
  //
  // - In a musical context, tune the ringing frequencies to the key of the song. Maybe use a bunch
  //   in series tuned to chords. Or maybe do that with the multicomb allpass. Or maybe tune them
  //   to all 12 notes in the 12-tone system.
  //
  // - Plot the phase-delay and group-delay as function of frequency for the damped allpass comb 
  //   and damped Schroeder allpass. We can do this by evaluating the complex frequency response 
  //   using the getTransferFunctionAt() functions.
  //
  // - Maybe also plot the group-delay response of the underlying comb filter. According to Barry
  //   Blesser's patent for the notchpass filter, it should feature peaks at frequencies related to
  //   the delayline length. That's the reason these frequencies recirculate longer. The notchpass
  //   design counteracts this by applying notches (or dips) as these frequencies in the feedback 
  //   path such the lower gain at these frequencies compensates for the longer delay thereby 
  //   equalizing the decay time for all frequencies.
  //
  // - Try to plot the decay time as function of frequency. One way to approach this is to use the
  //   information about the group delay of the comb together with the (frequency dependent) 
  //   feedback gain. But for this task, we may also try the "ringing response" functions. We could
  //   also look at the poles. Their distances from the unit circle determines the ringing time. 
  //   But it may be not so easy to compute all the poles because the filter are of very high order 
  //   such that normal polynomial root finding algorithms may not work (not sure, though - try 
  //   it!). But for that, we first need to implement root finding algorithms for 
  //   rsSparsePolynomial.
}

void dampedCombAllpassChainOf4()
{
  // Now, we create a chain of 4 such allpass comb filters and produce its impulse response.

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedCombAllpass<Real, Real, Real, Real>;

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
  rosic::writeToMonoWaveFile("DampedAllpasComb1.wav", &h1[0], N, (int)sampleRate);
  rosic::writeToMonoWaveFile("DampedAllpasComb2.wav", &h2[0], N, (int)sampleRate);
  rosic::writeToMonoWaveFile("DampedAllpasComb3.wav", &h3[0], N, (int)sampleRate);
  rosic::writeToMonoWaveFile("DampedAllpasComb4.wav", &h4[0], N, (int)sampleRate);
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

void dampedCombAllpassFeedbackBiquad()
{
  // Here, we try to use higher order feedback damping filters. Specifically, we use a biquad that
  // realized a wide-band dip (i.e. downward bell) frequency response.

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Complex = rsComplex<Real>;
  using Allpass = rsDampedCombAllpass<Real, Real, Real, Real>;

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
  //   we can't use the mode without predelay. See documentation of rsDampedCombAllpass. But maybe
  //   it can actually be made to work if we don't use the inverse damping filter as is but reflect
  //   its poles.
}

void dampedCombAllpassFreqDependentRT60()
{
  // We want to set up an rsDampedCombAllpass that achieves a desired overall decay time (in the 
  // RT60 sense) and also allows that decay time to be scaled at low and high frequencies via 
  // shelving filters. To test it, we split the impulse response into 3 bands (low, mid, high) and
  // plot the amplitude decay of the 3 bands with a logarithmic amplitude axis
  
  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedCombAllpass<Real, Real, Real, Real>;

  // User parameters:
  Real sampleRate = 48000;     // Sampling rate.
  int  numSamples = 24000;     // Number of samples to render.
  Real delay      =    50;     // Delay in samples.
  Real decayTime  =     1.0;   // Decay time for mid frequencies in seconds.
  Real lowFreq    =   250.0;   // Crossover freq between low and mid frequencies in Hz.
  Real lowScale   =     1.5;   // Decay time scaler for low frequencies.
  Real highFreq   =  4000.0;   // Crossover freq between mid and high frequencies in Hz.
  Real highScale  =     0.2;   // Decay time scaler for high frequencies.

  // Test:
  //lowScale  = 1.0; highScale = 1.0;

  // Compute intermediate values:
  Real decaySamples = decayTime     * sampleRate;
  Real lowOmega     = 2*PI*lowFreq  / sampleRate;
  Real highOmega    = 2*PI*highFreq / sampleRate;
  Real spikeFreq    = Real(sampleRate) / Real(delay);  // Frequency of the spikes
  int  N            = numSamples;


  // Create and set up the allpass filter:
  Allpass ap;
  ap.setMaxDelayInSamples(delay);
  rsSetupDecayTimes_LinViaFb(
    ap, delay, decaySamples, lowOmega, lowScale, highOmega, highScale, false);

  // Get comb impulse response:
  Vec hc(N);
  hc[0] = ap.getSample(1.0);
  for(int n = 1; n < N; n++)
    hc[n] = ap.getSample(0.0);
  //rsPlotVectors(hc);

  // Get impulse response:
  Vec h = impulseResponse(ap, N, 1.0);

  // Split the impulse response into low, mid and high parts:
  Vec hL(N), hM(N), hH(N);
  rsEngineersFilterMono engFlt;
  engFlt.setSampleRate(sampleRate);
  engFlt.setApproximationMethod(rsPrototypeDesigner<Real>::BUTTERWORTH);
  engFlt.setPrototypeOrder(4);

  // Extract lowpass part:
  engFlt.setMode(rsInfiniteImpulseResponseDesigner<Real>::LOWPASS);
  engFlt.setFrequency(lowFreq);
  engFlt.reset();
  for(int n = 0; n < N; n++)
    hL[n] = engFlt.getSample(h[n]);

  // Extract bandpass part:
  Real midFreq   = RAPT::rsBandwidthConverter::bandedgesToCenterFrequency(   lowFreq, highFreq);
  Real bandWidth = RAPT::rsBandwidthConverter::bandedgesToBandwidthInOctaves(lowFreq, highFreq);
  engFlt.setMode(rsInfiniteImpulseResponseDesigner<Real>::BANDPASS);
  engFlt.setBandwidth(bandWidth);
  engFlt.setFrequency(highFreq);
  engFlt.reset();
  for(int n = 0; n < N; n++)
    hM[n] = engFlt.getSample(h[n]);

  // Extract highpass part:
  engFlt.setMode(rsInfiniteImpulseResponseDesigner<Real>::HIGHPASS);
  engFlt.setFrequency(highFreq);
  engFlt.reset();
  for(int n = 0; n < N; n++)
    hH[n] = engFlt.getSample(h[n]);

  // Compute the dB-values of the absolute values of the signals:
  Vec db(N), dbL(N), dbM(N), dbH(N);
  Real ampFloor = rsDbToAmp(-140.0);
  for(int n = 0; n < N; n++)
  {
    db[n]  = rsAmpToDb(rsMax(rsAbs(h[n] ), ampFloor));
    dbL[n] = rsAmpToDb(rsMax(rsAbs(hL[n]), ampFloor));
    dbM[n] = rsAmpToDb(rsMax(rsAbs(hM[n]), ampFloor));
    dbH[n] = rsAmpToDb(rsMax(rsAbs(hH[n]), ampFloor));
  }

  // Extract envelope from the dB-decay signals:
  RAPT::rsSlewRateLimiter<Real, Real> envFlw;
  envFlw.setSampleRate(sampleRate);
  envFlw.setAttackTime(0.0);
  envFlw.setReleaseTime(10 * 1000 * Real(delay) / Real(sampleRate));
  Vec env  = filterResponse(envFlw, N, db );
  Vec envL = filterResponse(envFlw, N, dbL);
  Vec envM = filterResponse(envFlw, N, dbM);
  Vec envH = filterResponse(envFlw, N, dbH);

  // Apply a moving average filter of length equal to the delay to smooth out the artifacts of the
  // envelope follower:
  RAPT::rsMovingAverage<Real, Real> ma;
  ma.setMaxLengthInSamples(rsCeilInt(delay));
  ma.setLengthInSamples(rsCeilInt(delay));
  env  = filterResponse(ma, N, env );
  envL = filterResponse(ma, N, envL);
  envM = filterResponse(ma, N, envM);
  envH = filterResponse(ma, N, envH);
  // Could we implement such a filter with fractional delay?

  // ToDo:
  // 
  // - Try a processing chain of:  abs -> ma -> amp2db  instead of the current 
  //   abs -> amp2db -> envFollow -> ma  setup. Maybe that's better.





  //RAPT::rsEnvelopeExtractor<Real> envExt;

  // Plot the decaying lowpass, bandpass and highpass parts. The decay should be linear on a dB
  // scale:
  //rsPlotVectors(db, env);
  //rsPlotVectors(envM, envH);
  //rsPlotVectors(h);
  //rsPlotVectors(hL);
  //rsPlotVectors(hM);
  //rsPlotVectors(hH);
  //rsPlotVectors(h, hL, hM, hH);
  //rsPlotVectors(env, 20.0*h);
  rsPlotVectors(env, envL, envM, envH);
  //rsPlotVectors(db, dbL, dbM, dbH);
  //rsPlotVectors(db, dbL, dbM, dbH);
  //rsPlotVectors(dbL + 40.0, dbM, dbH); // dbL + 40 to lift it up because otherwise it's covered
  // ToDo: apply an envelope follower to the dB-decay profiles. That's much better for plotting. 
  // Then we may also get rid of the +40


  // Plot a spectrogram:
  //plotSpectrogram(&h[0], N, 64, 256, 256, sampleRate);
  plotSpectrogram(&h[0], N, 32, 128, 128, sampleRate);
  //plotSpectrogram(&h[0], N, 16,  64,  64, sampleRate);


  // Create an enveloped noise as input:
  Vec noise = createNoise(N, -1.0, +1.0, 0);
  rsAttackDecayFilter<Real> adEnv;
  adEnv.setAttackSamples(100.0);
  adEnv.setDecaySamples( 300.0);
  Vec nsEnv(N);
  nsEnv[0] = adEnv.getSample(1.0);
  for(int n = 0; n < N; n++)
    nsEnv[n] = adEnv.getSample(0.0);
  Vec envdNoise = nsEnv * noise;

  // Pass the enveloped noise through the allpass:
  Vec y = filterResponse(ap, N, envdNoise);
  rsPlotVectors(h, y);

  // Plot spectrogram of the allpass output::
  plotSpectrogram(&y[0], N, 256, 2048, 2048, sampleRate);
  plotSpectrogram(&y[0], N,  32,  128,  128, sampleRate);
  // numFrames does not seem to depend on the hopSize






  int dummy = 0;

  // Observations:
  //
  //
  // - With highScale == lowScale == 1, I would expect it to behave exactly like without any 
  //   feedback filter. But this is not the case! Apparently, the shelvers become nontrivial
  //   allpass filters in this case! This is not good - I think. ...or maybe it can be good from
  //   the perspective of diffusion - which might be a goal in some contexts. But we really should
  //   reconsider the shelver design and maybe try others! 
  //
  // - With longer decay times, the initial spike grows larger with respect to the decaying tail.
  //
  // - The bandpassed decay profiles look kind of as expected but they seem to feature a sort of
  //   two-stage decay. Maybe that's because of leakage from the adjacent frequency band. Maybe we
  //   can fix this by adjusting the cutoff frequencies of the filters and/or their steepness.
  //   Maybe test, if the two-stage decay becomes more pronounced when we use lower order 
  //   bandsplitting filters. ...hmm...hard to say. But  increaing the cutoff for the highpass
  //   to something like 1.5*highFreq does indeed seem to fix it. So, yeah - I think, the 
  //   explanation for the 2-stage decay is the leakage of the mid band into the high band. So, 
  //   the actually relevant slope for the high band is the inital (steeper) slope. The more 
  //   shallow slope towards the end comes from the more slowly decaying mid band.
  //
  // - But: the 2-stage decay is also present in the unfiltered signal. Maybe there, it has to do
  //   with how the spikes are smoothed out over time. At first, we see the decay of the spikes and
  //   and soon as they get buried in the smooth sine-like part, we see the decay of the smooth 
  //   part? Yes! That is indeed the case. The crossover moment is at around n = 6700. From this 
  //   point on onwards, the smooth part of the waveform peaks higher than the spikes.
  //
  // - The low frequency output looks very different from the others. It has this weird artifact at
  //   the beginning. Maybe it has to do with the fact that the low band contains the DC component
  //   whichbehaves differently? Maybe try it with a sign-flipped feedback. Will the low band still 
  //   look so different from mid and high bands?
  //
  // - Actually, it's a 3-stage decay. The initial decay goes up to sample 399, then the medium 
  //   decay goes to around sample 6700 and then the late decay begins.
  //
  // - The spectrogram looks weird. Not really what I expected. Figure out, why it looks so 
  //   strange!
  //
  //
  // ToDo:
  //
  // - Apply envelope followers to the (bandpassed) dB signals and fit a line to measure the 
  //   slopes. Compare the measured slope to what we would predict based on the settings. I don't
  //   know the formulas for what we should predict, though. We need to convert between a desired
  //   T60 and some decay slope - perhaps given in dB/sec.
  //
  // - Try using the more complicated calculations from 
  //   FeedbackDelayNetwork16::updateDampingAndCorrectionFilters() for designing the feedback 
  //   filter.
  //
  // - Plot an energy decay relief and check if it looks like expected. Maybe plot a spectrogram.
  //   ...hmm - the spectrogram of the impulse-response looks weird. Maybe it's because the hop
  //   size is synced with the delay. Maybe try it with a noise input.
  //
  // - Split the impulse response into low, mid and high parts and plot the log of the abs of these
  //   filtered signals. Maybe we can get a better temporal resolution by this.
  //
  // - We really need to check the spectrogram computation and plotting functions. I think, they 
  //   may not yet have been ready for general use. Make an experiment to test them. Maybe with a
  //   modal synthesis guitar tone.
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
  //
  // - Plot the decay time as function of frequency. Do this by taking the magnitude response of 
  //   the damping filter and using the formula that converts from decay-time to amplitude in 
  //   reverse. I expect that this may look less than desirable. If so, try tweaking the design of
  //   feedback bi-shelver. There is some example code for hwo this can be done in the 16-delayline
  //   FDN.
}


void dampedCombAllpassFractional1()
{
  // We set up an rsDampedCombAllpass with a fractional amount of delay. This is implemented by 
  // using a (linear) interpolation filter in the feedback path. In this case here, we don't really
  // have any damping filter in the feedback path (other than the interpolation filter which does
  // a high frequency damping dependning on the amount of fractional delay - it'S most severe when
  // the fraxctional delay is 0.5 in which case there's a zero at z = -1).

  // Define types to be used:
  using Real      = double;
  using Vec       = std::vector<Real>;
  using Allpass   = rsDampedCombAllpass<Real, Real, Real, Real>;

  int  numSamples = 5000;     // Number of samples to render.
  Real delay      =  100.3;   // Delay in samples - not necessarily integer, though 
  Real feedback   =    0.9;

  int  N        = numSamples;
  int  maxDelay = (int) rsCeil(delay);   // Verify!

  // Create and set up the allpass filter:
  Allpass ap;
  ap.setMaxDelayInSamples(maxDelay);
  rsSetupFractional_LinViaFb(ap, delay, feedback, false);

  // Get impulse response of comb:
  Vec hc(N);
  hc[0] = ap.getSampleComb(1.0);
  for(int n = 1; n < N; n++)
    hc[n] = ap.getSampleComb(0.0);
  rsPlotVectors(hc);

  // Get impulse response of full allpass:
  Vec h = impulseResponse(ap, N, 1.0);
  rsPlotVectors(h);


  // Observations:
  //
  // - With delay = 100.25, we get a sort of "mexican hat" shape in the second spike of the 
  //   allpass response
  //
  //
  // ToDo:
  //
  // - Try implementing allpass interpolation. I did already try but with no success, so far. The 
  //   combs became unstable. But I'm not totally convinced that it is impossible. Maybe I just did
  //   something wrong.
}

void dampedCombAllpassFractional2()
{
  // Now we wanto to create damped allpass with fractional delay and some actual damping filters.
  // They way we do this is to absorb the interpolation filter into the damping filter.

  // Define types to be used:
  using Real         = double;
  using Vec          = std::vector<Real>;
  using Allpass      = rsDampedCombAllpass<Real, Real, Real, Real>;
  using CombSettings = rsDampedCombSettings<Real, Real, Real>;
  using TransFunc    = rsSparseTransferFunction<Real, Real>;


  // User parameters:
  Real sampleRate = 48000;     // Sampling rate.
  int  numSamples = 12000;     // Number of samples to render.
  Real delay      =   100.3;   // Delay in samples - not necessarily integer, though 
  Real decayTime  =     0.5;   // Decay time for mid frequencies in seconds.
  Real lowFreq    =   250.0;   // Crossover freq between low and mid frequencies in Hz.
  Real lowScale   =     1.5;   // Decay time scaler for low frequencies.
  Real highFreq   =  4000.0;   // Crossover freq between mid and high frequencies in Hz.
  Real highScale  =     0.2;   // Decay time scaler for high frequencies.

  // Compute intermediate values and define abbreviations:
  Real decaySamples = decayTime     * sampleRate;
  Real lowOmega     = 2*PI*lowFreq  / sampleRate;
  Real highOmega    = 2*PI*highFreq / sampleRate;
  int  N            = numSamples;
  int  maxDelay     = (int) rsCeil(delay);   // Verify!


  // Compute settings for a linearly interpolating comb once by baking the interpolator filter into
  // the feedback dampling filter and once by baking it into the delayline and compare the 
  // resulting transfer functions:
  CombSettings s1, s2;
  rsSetupDecayTimes_LinViaFb(
    s1, delay, decaySamples, lowOmega, lowScale, highOmega, highScale, false);
  rsSetupDecayTimes_LinViaDly(
    s2, delay, decaySamples, lowOmega, lowScale, highOmega, highScale, false);
  TransFunc tf1, tf2;
  s1.getCombTransferFunction(&tf1);
  s2.getCombTransferFunction(&tf2);
  bool ok = tf2.isCloseTo(tf1, 1.e-14);
  rsAssert(ok);
  // Yep - they are indeed equal. That means, at least for the linear interpolator, we can bake it
  // either into the delay or into the feedback filter. It doesn't matter. I think, this may be 
  // generally true for FIR interpolators but not for IIR interpolators - at least, that's my
  // practical experience so far. ToDo: Figure this out theoretically. Derive the full transfer 
  // functions for a comb using a 1st order allpass interpolator baked into the delay and damping
  // filter. It would actually be nice, if we always could bake the interpolator into the feedback
  // filter because that would simplify the implementation of the whole comb - we wouldn't need to 
  // use a delayline that has interpolation already built in. Although, it's no big deal to use one
  // either.


  // Create and set up the allpass filter:
  Allpass ap;
  ap.setMaxDelayInSamples(maxDelay);
  rsSetupDecayTimes_LinViaFb(
    ap, delay, decaySamples, lowOmega, lowScale, highOmega, highScale, false);

  // Get output of comb:
  Vec hc(N);
  hc[0] = ap.getSampleComb(1.0);
  for(int n = 1; n < N; n++)
    hc[n] = ap.getSampleComb(0.0);
  rsPlotVectors(hc);

  Vec h = impulseResponse(ap, N, 1.0);
  rsPlotVectors(h);


  // ToDo:
  //
  // - Maybe create 3 filters: one with floor(delay), one with delay and one with ceil(delay) and 
  //   compare their outputs
}


void dampedMultiCombAllpassIdea()
{
  // We extract the transfer functions from 3 rsDampedCombAllpass objects with different delays as
  // rsSparseRationalFunction and then use a weighted sum of them as our "multi comb" transfer 
  // function. Then we apply the usual inversion/reflection mumbo jumbo to tunn the multicomb into
  // an allpass.
  

  // Define types to be used:
  using Real      = double;
  using Complex   = rsComplex<Real>;
  using Vec       = std::vector<Real>;
  using Allpass   = rsDampedCombAllpass<Real, Real, Real, Real>;
  using TransFunc = rsSparseTransferFunction<Real, Real>;
  using SparseFlt = rsSparseFilter<Real, Real, Real>;

  // User parameters:
  Real sampleRate = 48000;     // Sampling rate.
  int  numSamples = 24000;     // Number of samples to render.

  Real delay1     =    23;     // Delay of 1st delayline
  Real delay2     =    31;     //          2nd
  Real delay3     =    41;     //          3rd
  Real delayScale =     5;     // Scale factor for all delays

  Real gain1      =     1.0;   // Gain of 1st delayline
  Real gain2      =     1.0;   //         2nd
  Real gain3      =     1.0;   //         3rd

  Real decayTime  =     1.0;   // Decay time for mid frequencies in seconds.
  Real lowFreq    =   250.0;   // Crossover freq between low and mid frequencies in Hz.
  Real lowScale   =     1.5;   // Decay time scaler for low frequencies.
  Real highFreq   =  4000.0;   // Crossover freq between mid and high frequencies in Hz.
  Real highScale  =     0.2;   // Decay time scaler for high frequencies.

  decayTime = 0.2;               // Test
  //lowScale = highScale = 1.0;  // Test

  // Compute intermediate values and define abbreviations:
  Real decaySamples = decayTime     * sampleRate;
  Real lowOmega     = 2*PI*lowFreq  / sampleRate;
  Real highOmega    = 2*PI*highFreq / sampleRate;
  Real spikeFreq1   = Real(sampleRate) / Real(delay1);  // Frequency of the spikes
  int  N            = numSamples;

  delay1 *= delayScale;
  delay2 *= delayScale;
  delay3 *= delayScale;

  // Create and set up the prototype allpass filters. We are interested mostyl in the comb transfer
  // functions here:
  Allpass ap1, ap2, ap3;

  ap1.setMaxDelayInSamples(delay1);
  rsSetupDecayTimes_LinViaFb(
    ap1, delay1, decaySamples, lowOmega, lowScale, highOmega, highScale, false);

  ap2.setMaxDelayInSamples(delay2);
  rsSetupDecayTimes_LinViaFb(
    ap2, delay2, decaySamples, lowOmega, lowScale, highOmega, highScale, false);

  ap3.setMaxDelayInSamples(delay3);
  rsSetupDecayTimes_LinViaFb(
    ap3, delay3, decaySamples, lowOmega, lowScale, highOmega, highScale, false);

  // Retrieve the comb transfer functions:
  Real tol = 1.e-13;
  TransFunc U1 = ap1.getCombTransferFunction(tol);
  TransFunc U2 = ap2.getCombTransferFunction(tol);
  TransFunc U3 = ap3.getCombTransferFunction(tol);

  // Combine the comb transfer functions into one, set up a sparse filter that realizes that sum
  // of combs and retrieve its impulse response:
  TransFunc U = gain1 * U1  +  gain2 * U2  +  gain3 * U3;
  SparseFlt comb;
  comb.setMaxDelayInSamples(U.getFilterOrder());
  comb.setup(U);
  Vec hu = impulseResponse(comb, N, 1.0);
  rsPlotVector(hu);

  // Test the non-allocating transfer fucntion computations:
  TransFunc U1n;
  ap1.getCombTransferFunction(&U1n);



  // Compute different notions of density for inspection:
  double densityNum = U.getNumeratorDensity();
  double densityDen = U.getDenominatorDensity();
  double densitySep = U.getSeparatedDensity();
  double densityCmb = U.getCombinedDensity();


  // Create the corrector filter and apply it to the comb impulse response:
  TransFunc C = U;
  C.invert();
  C.reflectZeros();
  SparseFlt corr;
  corr.setMaxDelayInSamples(C.getFilterOrder());
  corr.setup(C);
  Vec h = filterResponse(corr, N, hu);
  rsPlotVector(h);

  // Transform the comb to max-phase and obtain the impulse response:
  comb.reflectZeros();
  Vec hup = impulseResponse(comb, N, 1.0);
  rsPlotVector(hup);

  // Apply the corrector to the phased comb output:
  Vec hp = filterResponse(corr, N, hup);
  rsPlotVector(hp);

  rsPlotVectors(h, hp);


  rosic::writeToMonoWaveFile("TriCombAllpass_CombSum.wav",          &hu[0],  N, (int)sampleRate);
  rosic::writeToMonoWaveFile("TriCombAllpass_Corrected.wav",        &h[0],   N, (int)sampleRate);
  rosic::writeToMonoWaveFile("TriCombAllpass_CombSum_Phased.wav",   &hup[0], N, (int)sampleRate); 
  rosic::writeToMonoWaveFile("TriCombAllpass_Corrected_Phased.wav", &hp[0],  N, (int)sampleRate); 


  // Observations:
  //
  // - With decayTime = 1.0, the impulse response of the comb-sum shows spikes at the products of 
  //   the delays, i.e. at 713 = 23*31, 943 = 23*41, 1271 = 31*41 and their multiples, i.e. 
  //   1426 = 2 * 713, etc.
  //
  // - With delayScale = 5 and decayTime = 0.2, there's actually a build-up phase in the
  //   "corrected-phased" version. It sounds like breath-noise of some flute-ish instrument 
  //   combined with a little "thump" or "pop" transient. Th just corrected one without the 
  //   "phased" sound like a pure breath-noise without the pop.
  //
  // - The impulse response of the corrected comb sum has a strong initial bipolar spike and is 
  //   rather quiet after that. But with shorter decay times, the effect becomes smaller. Also,
  //   the spikes at 713, .. become less pronounced.
  //
  // - Try it with a full series of harmonics and with odd harmonic only (this is the current case)
  //
  // - U has 49/56 terms and therefore C has 56/49 terms. The maximum delay (i.e. total filter 
  //   order) is 101 and determined by the 56th term (index 55) of the denominator of U. That's a 
  //   pretty high density actually. Maybe define the density of a sparse filter as the number of
  //   coeffs it actually has divided by the maximum number that it could have for the given order.
  //   The maximum number it could have would be 2*max(deg(num)+1, deg(den)), I think.
  //   Looks like the desnity is 0.5833333. ...but verify the definition of the formula - if it
  //   makes sense to define density that way.
  //
  //
  // Conclusions:
  //
  // - The filter should perhaps be best used with shorter decay times and with a phased comb sum.
  //
  //
  // ToDo:
  //
  //
  // - Maybe implement a musically tunable reverb algorithm based on that idea. Maybe call it 
  //   TuniVerb. It should give the user the possibility to set up a parallel connection of (up to)
  //   some number N of combs which are then turned into an allpass via our transfer function 
  //   inversion-and-reversal magic. We may want tune the combs to musical notes.
  //
  // - Maybe the combs should all have the same decay time. Or maybe the higher combs should have 
  //   a shorter decay time. Or maybe make the scaling of the decay-time with comb-freq a user 
  //   param.
  //
  // - The user should be able to switch between all and odd harmonics - maybe per comb.
  //
  // - Try to implement a fractional delay by using allpass interpolation of the delayline. This
  //   will probably also solve the spikes problem when the delays are chosen in suitable 
  //   irrational ratios. I tried to bake the allpass interpolator into the feedback filter. But 
  //   with allpass interpolation, it didn't work - the combs were unstable. It did seem to work 
  //   with linear interpolation, though. I'm not sure anymore if the idea of baking the 
  //   interpolator into the feedback filter is right to begin with. I think, at least, when the 
  //   interpolation (allpass) filter is combined with another filter (such as the damping filter)
  //   it may not work because the combination with other filter changes the meaning of y[n-1].
  //   ...I'm not really sure about all of this. More research is needed. Ahhh! I think, baking the 
  //   inteprolator into the feedback filter F(z) is wrong indeed! It needs to be baked into the 
  //   delay filter A(z) instead! But maybe with pure FIR interpolators it really doesn't matter? 
  //   It seems to work with the linear interpolator baked into the feedback path.
  //
  // - Try longer delay lengths. They are chosen pretty short in the experiment (even the longets 
  //   is less than a millisecond). Maybe try something on the order of 10-20 milliseconds. That 
  //   should be suitable for reverb.
  //
  // - Try phasing (i.e. zero-reflecting) the damping filters (per comb).
}


void dampedMultiCombAllpassClass()
{
  using Real     = double;
  using Vec      = std::vector<Real>;
  using Filter   = rsDampedMultiCombAllpass<Real, Real, Real>;
  using PlotType = rsFrequencyResponsePlotSettings::Type;


  int  numSamples =  8192;
  Real sampleRate = 48000;
  Real frequency  =   480.0;   // Reference frequency
  Real decayTime  =     0.2;   // Decay time for mid frequencies in seconds.
  Real lowFreq    =   250.0;   // Crossover freq between low and mid frequencies in Hz.
  Real lowScale   =     1.5;   // Decay time scaler for low frequencies.
  Real highFreq   =  4000.0;   // Crossover freq between mid and high frequencies in Hz.
  Real highScale  =     0.2;   // Decay time scaler for high frequencies.
  bool serial     = false;


  // Create and set up the filter:
  Filter flt;
  flt.setFilterOrderLimits(8191, 4);
  flt.setSampleRate(sampleRate);
  flt.setFrequency(frequency);
  flt.setDecayTimeInSeconds(decayTime);
  flt.setLowCrossoverFreq(lowFreq);
  flt.setLowDecayScale(lowScale);
  flt.setHighCrossoverFreq(highFreq);
  flt.setHighDecayScale(highScale);
  flt.setNumCombs(3);
  flt.setSerialCombsMode(serial);
  flt.setCombFreqScale(1, 1.3352657623556);
  flt.setCombFreqScale(2, 1.7238234123465);
  // The numbers were just randomly created via the hsdjfgasdfhgskajdfhk method


  //// Test parallel vs serial combs:
  //flt.setSerialCombsMode(false);   Vec hp = impulseResponse(flt, 500, 1.0);
  //flt.setSerialCombsMode(true);    Vec hs = impulseResponse(flt, 500, 1.0);
  //Vec diff = hp - hs;
  //flt.setSerialCombsMode(serial);  // Set it back to the user's setting
  //rsPlotVectors(hp, hs, diff);
  //// They are indeed the same! 


  bool ok = true;

  // Get an plot the impulse response:
  int N = numSamples;
  Vec h = impulseResponse(flt, N, 1.0);
  ok &= isAllpass(h, 0.02);
  rsPlotVector(h);

  // Helper function to create a frequency response plot:
  auto plotResponse = [&](PlotType type)
  { 
    rsPlotDampedCombAllpassResponses(flt, type, sampleRate); 
  };

  // Plot various frequency responses:
  plotResponse(PlotType::decibelsAndPhase);
  plotResponse(PlotType::decibelsAndRing);


  // Observations:
  //
  // - The impulse responses of parallel and serial mode are the same! There doesn't even seem to
  //   a difference in roundoff error. 
  //
  //
  // ToDo: 
  //
  // - Figure out why serial and parallel connection gives the same result. Is this a bug or is 
  //   this normal? Check the math! Well, that the denominators (i.e. poles) are the same is
  //   expected. But what about the numerators (i.e. zeros)? Let's try it with two combs. The 
  //   parallel connection has the transfer function: 
  //   Up(z) = U1(z) + U2()z   and  Us(z) = U1(z) * U2(z)
  //   where  U1(z) = z^-M1 / (1 + k * z^-1 * z^-M1), U2(z) = z^-M2 / (1 + k * z^-1 * z^-M2)
  //   Try to bring Up(z) into same the same form as Us(z)...
  //
  //
}


// This is needed for the case when we want ot have a complex TPar. It's a bit dirty to use 
// the explicit type double. Maybe it should be a third template parameter TReal or something.
template<class TSig, class TPar, class TDly, class TTol>
void rsSetupHighDampComplex(rsDampedCombAllpass<TSig, TPar, TDly, TTol>& flt,
  int delay, TPar feedback, double dampOmega, double dampGain, bool predelay)
{
  double a[2], b[2]; a[0] = 1;
  rsMake1stOrderHighShelf(dampOmega, dampGain, &b[0], &b[1], &a[1]);
  TPar A[2], B[2]; 
  A[0] = TPar(a[0]);
  A[1] = TPar(a[1]);
  B[0] = TPar(b[0]);
  B[1] = TPar(b[1]);

  flt.setup(delay, feedback, 1, B, A, predelay);
}
// Clean this up to make it nice!
// I think dampOmega, dampGain should be TDly...but then the name TDly is not a good fit anymore.
// Maybe TTime - it applies to time variables and frequencies (reciprocal of time, kind of) and 
// also to the gain...which is an absolute amplitude/magnitude value. These are all values that
// can only take on non-negative real values.


void dampedCombAllpassComplex()
{
  // The code has been commented out because it doesn't compile anymore since we switched to using
  // rsDelayRounding instead of rsDelay. The problem is that instantiating rsDampedCombAllpass with
  // TPar = complex to get complex feedback interferes with using TPar for the non-integer delay 
  // parameter (and probably also with getTransferFunctionAt() which we did not yet try to use 
  // here). Maybe it can be solved by having a 3rd template parameter TFdbk for the feedback? Using
  // k = TSig as we did before may solve the issue with the delay but probably not the one with
  // getTransferFunctionAt.
  // ...OK...this has been "fixed" by reverting to using rsDelay again...but I think, this is not 
  // yet the final solution. I think, getTransferFunctionAt won't work.


  // We instantiate rsDampedCombAllpass with a complex datatype for the signals. The feedback gain
  // is also complex. Using a complex feedback gain further increases the space of the things that 
  // we can do with this filter. 


  // Define types to be used:
  using Real    = double;
  using Complex = rsComplex<Real>;
  //using Complex = std::complex<Real>;
  using VecR    = std::vector<Real>;
  using VecC    = std::vector<Complex>;
  using Allpass = rsDampedCombAllpass<Complex, Complex, Real, Real>;

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

  ap.setMaxIntDelayInSamples(delay);
  //ap.setupHighDamp(delay, fb, w, dampGain, false);
  rsSetupHighDampComplex(ap, delay, fb, w, dampGain, false);


  VecC h = impulseResponse(ap, N, Complex(1));
  plotComplexVectorReIm(h);                      // Doesn't accept rsComplex - fix that!

  // Extract real and imaginary parts:
  VecR hr(N), hi(N);
  for(int n = 0; n < N; n++)
  {
    hr[n] = rsReal(h[n]);
    hi[n] = rsImag(h[n]);
  }

  // Plot spectra:
  using SpecPlot = SpectrumPlotter<Real>;
  SpecPlot plt;
  plt.setFftSize(N);
  plt.plotSpectra(N, &hr[0], &hi[0]);

  // Test transfer function computation:
  Complex z(0.6, 0.8);
  //Complex H = ap.getCombTransferFunctionAt(z);
  // Doesn't compile because argument and return value are of type rsComplex<TPar> and TPar is
  // complex, so we'd get a nested complex type.
  int dummy = 0;


  // Observations:
  //
  // - Using complex numbers imprints and undulation onto the impulse reponses.
  //
  // - The spectra of real and imaginary part by themselves are not allpass like. They are more
  //   of a strange comb like structure. With fbPhase = pi/4, one has zeros between the peaks, the
  //   other doesn't. 
  //
  // - Feedback phases of 0 and pi give purely real outputs. With 0, it's unipolar, with pi its
  //   bipolar. That's how it has to be - we expect to get back to behavior of positive and
  //   negative signs. pi/2 gives a bipolar spike train both real and imaginary part. Same 
  //   for 3*pi/2.
  //
  //
  // ToDo:
  //
  // - Try to make it compile again. Maybe we need to introduce a 3rd template parameter for the
  //   delay in class rsDampedCombAllpass and assign that to Real. But that alone won't solve all
  //   problems. Maybe the getTransferFunctionAt() functions need to get their own template 
  //   parameter for the argument z and return type. This would have the additional advantage that
  //   we wouldn't need to commit to use RAPT::rsComplex or std::complex in the declaration of 
  //   getTransferFunctionAt(). This choice could be later made by the client code.
  //
  // - Check if, when using complex feedback and filter coeffs, we need to additionally conjugate
  //   the coeffs after reversing them to get an allpass. Look up the paper by Sebastain Schlecht:
  //   https://www.mdpi.com/2076-3417/10/1/187 ...hmm...it doesn't say anything about that. Maybe 
  //   make some derivations for simple low order filters (1st order, biquad). A reflection of a 
  //   root r in the unit circle is given by ref(r) = 1 / conj(r). Use that formula for a biquad
  //   defined in terms of poles and zeros.
  //
  // - Make sure that everything works with rsComplex and std::complex for Complex. We need to
  //   implement real/imag for rsComplex or rsReal/rsImag for both. I think, the latter way is the
  //   right way. It's consistent with how the rest of the library works.
  //
  // - I think, the feedback phase should scale with the delay to achieve a uniform undulation
  //   frequency. Figure this out!
}

void dampedCombAllpassNonLin()
{
  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Allpass = rsDampedCombAllpassNonLin<Real, Real, Real, Real>;

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
  // An experiment that lets us look at the content of the delaylines of rsDampedCombAllpass. To 
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
  using Allpass = rsDampedCombAllpass<Real, Real, Real, Real>;

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

void dampedAllpassBiComb_1p()
{
  // This is now obsolete because the (relevant) code has been moved into 
  // dampedAllpassBiCombUnitTest()

  // We want to create an allpass by using a correction filter after a parallel connection of two
  // Karplus-Strong like comb filters. The first step would be to derive the direct form transfer
  // function of the two parrallel combs. I did this with sageMath in AllpassStuff.txt in the 
  // private repo but trying to implement the dual comb in direct form doesn't yet work. I'm not 
  // sure if I have a mistake in the derivation or implementation. But before goimg any further, 
  // it's really necessarry to make the direct form work so we can be sure to build on formulas
  // that are actually correct ...TBC...


  // Define types to be used:
  using Real         = double;
  using Complex      = rsComplex<Real>;
  using Vec          = std::vector<Real>;
  using Allpass      = rsDampedAllpassBiComb_1p<Real, Real, Real>;
  using SparseFilter = rsSparseFilter<Real, Real, Real>;


  int  N      = 2048;    // Number of samples to render
  int  delay1 = 23;
  int  delay2 = 29;
  Real g1     = 0.6;
  Real g2     = 0.7;
  Real k1     = 0.9;
  Real k2     = 0.8;

  //k1 = k2 = 0.0;  // For test

  // Compute the feedback filter coeffs:
  Real b10, b11, a11; rsMake1stOrderHighShelf(0.5, 0.9, &b10, &b11, &a11);
  Real b20, b21, a21; rsMake1stOrderHighShelf(0.7, 0.8, &b20, &b21, &a21);
  //b10 = 1; b11 = a11 = 0; b20 = 1; b21 = a21 = 0;   //  For test with neutral filter

  // Create and set up the allpass:
  Allpass ap;
  ap.setMaxDelayInSamples(delay2);
  ap.setup(delay1, g1, k1, b10, b11, a11,
           delay2, g2, k2, b20, b21, a21);

  // Produce the impulse response of the two parallel comb filters:
  Vec hc(N);
  hc[0] = ap.getSampleCombs(1.0);
  for(int n = 1; n < N; n++)
    hc[n] = ap.getSampleCombs(0.0);

  // Let the ap convert itself into a sparse direct form filter and check that this converted 
  // filter has the same impulse response:
  SparseFilter sf;
  ap.convertCombSumToDirectForm(&sf);
  Vec hc2 = impulseResponse(sf, N, 1.0);

  bool ok = true;
  ok &= rsIsCloseTo(hc, hc2, 1.e-14);
  rsPlotVectors(hc, hc2);

  // OK - so far, so good. We can produce the output of the weighted sum of the two comb filters by
  // two different algorithms: 
  //
  //   (1) Literally implementing two combs and mixing their outputs. 
  //   (2) Implementing the whole thing in direct form
  //
  // From a practical point of view, it makes no sense to use algorithm 2 but it is a stepping 
  // stone to turn the filter into an allpass. We now need to invert the filter of algorithm 2 by
  // swapping numerator and denominator and then reflect the zeros in the unit circle. Algo 2 leads
  // to straightforward inversion.


  // Try inverting the direct form version of the bi-comb. Applying the inverted filter to the 
  // impulse response of the original filter should give us back a unit impulse:
  sf.invert();
  Vec y = filterResponse(sf, N, hc);
  ok &= rsIsUnitImpulse(y, 1.e-14);
  //rsPlotVectors(y);

  // Now we want to reflect the zeros in the unit circle. That is eventually the corrector filter 
  // that we want to use:
  sf.reflectZeros();
  y = filterResponse(sf, N, hc);
  ok &= isAllpass(y, 1.e-4);
  //rsPlotVectors(y);

  //// Test the whole filter, i.e. the comb-sum with corrector applied:
  Vec h = impulseResponse(ap, N, 1.0);
  ok &= isAllpass(h, 1.e-4);


  // Create a filter that also reflects the zeros of the original comb-filter sum. We do this by 
  // implementing the bicomb itself in (sparse) direct form and reflect the zeros of that, too:
  SparseFilter sf1, sf2;
  ap.convertCombSumToDirectForm(&sf1);
  sf2.copySettingsFrom(sf1);
  sf1.reflectZeros();
  sf2.invert();
  sf2.reflectZeros();
  Vec z(N);
  z[0] = sf2.getSample(sf1.getSample(1.0));
  for(int n = 1; n < N; n++)
    z[n] = sf2.getSample(sf1.getSample(0.0));
  ok &= isAllpass(z, 1.e-3);
  rsPlotVectors(y, z);




  int dummy = 0;


  // ToDo:
  //
  // - Implement rsDampedAllpassBiComb_1p::getTransferFunctionAt() and test it.
  //
  // - Figure out the poles and zeros of the bi-comb - at least roughly. Are they inside or outside
  //   the unit circle? Maybe try making a plot of the magnitude (or dB values) of the z-domain
  //   transfer function. But I'm not sure how numerically stable the evaluation of the transfer 
  //   function is for z with small radius. In this case, the z^-n sequence diverges. Maybe we coul 
  //   try to find them with a polynomial root finder. But for this, we should perhaps use smaller
  //   numbers for the delays because otherwise, we have to deal with very high degree polynomials
  //   and I'm not sure about the numerical performance of the root finder for very high degree
  //   polynomials. If the zeros are all inside the unit circle, then turning the numerator into
  //   a denominator (i.e. turning them into poles) shouldn't be a problem. When they are all 
  //   outside the unit circle, then we could do a reflection in the unit circle.
  //
  // - Maybe the current filter has no causal inverse. But maybe if we sap the positions of F1(z)
  //   and z^-M1 and likewise for F2(z) and z^-M2, then it could work out?
  //
  // - Apply the inverse comb to a unit impulse just for curiosity
  //
  // - Maybe we could also take the converted bi-comb, then invert-and-reflect it and combine it
  //   with the reflected one rather than the original one. Maybe that would add even more
  //   dispersion?
}

void dampedCombAllpasses()
{
  // Various experiments with ideas for turning comb filters into allpasses.

  // Currently active experiment:
  //dampedCombFilter();
  //dampedCombAllpassResponses();
  dampedMultiCombAllpassClass();


  // All experiments:
  dampedCombFilter();                   // Stub: we want to factors out the comb
  dampedCombAllpassIdea();              // Initial explorations of the idea.
  dampedCombAllpassClass();             // ..the idea has been wrapped into a class.
  dampedCombAllpassResponses();         // Plot responses for phase, group delay, etc.
  dampedCombAllpassChainOf4();          // A chain of 4 damped allpass combs.
  dampedCombAllpassFeedbackBiquad();    // With biquad in feedback path
  dampedCombAllpassFreqDependentRT60(); // Frequency dependent RT60 decay time.
  dampedCombAllpassFractional1();       // Fractional delay by lerp in feedback path.
  dampedCombAllpassFractional2();       // Fractional delay by interpolating the delay line.
  dampedMultiCombAllpassIdea();         // Idea for allpass from multiple parallel combs. 
  dampedMultiCombAllpassClass();        // ..the idea has been wrapped into a class
  dampedCombAllpassComplex();           // Experiment with complex feedback factor.
  dampedCombAllpassNonLin();            // Stub: we want to experiment with nonlinearities.
  dampedAllpassDelayContent();          // Investigate contents of the delaylines.
  dampedSchroederAllpass();             // Schroder-like allpass from damped comb (verify!)
  dampedAllpassBiComb_1p();             // Multicomb with 2 combs and 1-pole damper.


  // ToDo:
  //
  // - Instantiate rsDampedCombAllpass with rsFloat64x2 for TSig and use a feedback of [+k, -k], 
  //   i.e. different signs for the feedback for left and right channel. Thta should give 
  //   complementary combs for left and right channel.
  //
  // - Instantiate the damped comb filter for rsFloat64x2. Interpret that data type as stereo 
  //   signal. Use a feedback gain that has a positive sign for one channel and a negative sign for
  //   the other. That should give one channel a harmonic series of e.g. 100,200,300,400,500,... 
  //   and the other 50,150,250,350,450,... That may make for a nice stereoization effect. Verify
  //   the mono compatibility! Let the user dial in the strength of the effect via the feedback
  //   gain.
  //
  // - In the multicomb allpass, have an option to not use a parallel connection of combs but a 
  //   serial one. To get a serial comb chain, we just need to adapt the code that accumulates the
  //   transfer functions (from init to zero and adding to init to one and multiplying). Try a 
  //   series connection of two combs with same settings but one with different signs for the 
  //   feedback. 
  //
  // - Try to replace the delay with other allpass structures such as a Schroeder allpass, a
  //   disperser, etc. Maybe we can even use a non-allpass if we adapt the algorithm to it? We 
  //   would probably have to use the inverse filter of A(z) in the corrector, so A(z) would need 
  //   to have a stable inverse. With an allpass, this was not an issue because our end goal is an
  //   allpass anyway, so we didn't need to invert it. Maybe with a feedforward comb, we would 
  //   venture into "notchpass" territory? I'm not sure -> figure that out. I think, to turn the 
  //   delay into a feedforward comb, all we need is to add a direct path (with adjustable gain).
  //   But  I think, the feedforward comb would have one sample less delay than the whole feedback
  //   loop. That may cause a misalignment of notch- vs resonance frequencies. Cand this be fixed 
  //   somehow?
  //
  // - Instead of one filter with high feedback, try to use multiple ones with less feedback. I 
  //   guess, the ringing response will be less spikey when doing this which may be desirable.
  //   maybe first try it with two filters and the same delay for both. That's not ideal but we 
  //   may learn something from that. This situation may show more clearly how much we need to 
  //   reduce the feedback when putting two filters in series in order to achieve the same RT60.
  //   Maybe plot the decay envelopes like we do in dampedCombAllpassFreqDependentRT60(). Maybe 
  //   factor out a function for creating such plots and then apply it to this new setting.
}

//=================================================================================================

std::vector<double> extendedFDN_3x3()
{
  // We implement the idea that is outlined in Notes/DSP/FeedbackDelayNetworkIdeas.md in the 
  // research repo

  using Real  = double;
  using Vec   = std::vector<Real>;
  using Delay = RAPT::rsDelay<Real>;


  int numSamples = 1000;

  // Delay values:
  int M1 = 17;
  int M2 = 23;
  int M3 = 29;

  int N1 = 20;
  int N2 = 26;
  int N3 = 36;


  // Rotation angles for feedback matrix (in degrees):
  Real rx, ry, rz;
  rx = ry = rz = 45;    // Best theoretically? Most diffusive?
  //rx = ry = rz = 60;  
  //rx = ry = rz = 90;  // Very bad!

  // Decay time (RT60) in samples:
  Real decay = 1000;
  //Real decay = RS_INF(Real);  // Infinite decay time means no deacy at all

  // Output vectors:
  Real c1 = +1, c2 = -1, c3 = +1;
  Real d1 = -1, d2 = +1, d3 = -1;


  // Some stuff to uncomment to test with other settings:

  // Test: Using all zeros here should reduce it to the regular 3x3 FDN. Uncomment to test without
  // the second round fo delay:
  //N1 = N2 = N3 = 0; 


  // Compute the 9 = 3*3 total roundtrip delays for inspection in debugger:
  int M1_N1 = M1 + N1;
  int M1_N2 = M1 + N2;
  int M1_N3 = M1 + N3;
  int M2_N1 = M2 + N1;
  int M2_N2 = M2 + N2;
  int M2_N3 = M2 + N3;
  int M3_N1 = M3 + N1;
  int M3_N2 = M3 + N2;
  int M3_N3 = M3 + N3;


  // Compute the feedback damping gain factors:
  Real amp = Real(0.001);
  Real a1  = rsDecayTimeToFeedbackGain(decay, Real(M1), amp);
  Real a2  = rsDecayTimeToFeedbackGain(decay, Real(M2), amp);
  Real a3  = rsDecayTimeToFeedbackGain(decay, Real(M3), amp);
  Real b1  = rsDecayTimeToFeedbackGain(decay, Real(N1), amp);
  Real b2  = rsDecayTimeToFeedbackGain(decay, Real(N2), amp);
  Real b3  = rsDecayTimeToFeedbackGain(decay, Real(N3), amp);


  // Helper function:
  auto setDelay = [](Delay& delay, int amount)
  {
    delay.setMaxDelayInSamples(amount);
    delay.setDelayInSamples(amount);
  };

  // Create and set up the delay lines:
  Delay A1, A2, A3;  // Delaylines left to the feedback matrix
  Delay B1, B2, B3;  // Delaylines right to the feedback matrix
  setDelay(A1, M1);
  setDelay(A2, M2);
  setDelay(A3, M3);
  setDelay(B1, N1);
  setDelay(B2, N2);
  setDelay(B3, N3);

  // Create and set up the feedback matrix:
  rsRotationXYZ<Real> F;
  Real toRad = PI/180;
  F.setAngles(toRad*rx, toRad*ry, toRad*rz);

  //rsMatrix3x3<Real> F;
  ////Real s = Real(1) / rsSqrt(3.0);
  //Real s = 0.5015;
  //F.setValues(+s,+s,-s, +s,-s,+s, -s,+s,+s);
  // This makes the signal blow up when feedback is set to 1! Isn't it supposed to be unitary 
  // when we divide by sqrt(3)? Oh! no! The rows are actually not orthogonal! Figure out what 
  // scaler we need to use! 0.5 seems to be stable but may be slightly decaying. 0.51 blows up.
  // 0.505 also. 0.502 may be very slowly growing. 0.5015 seems to be (close to) the stability 
  // limit. ToDo: figure out theoretically what the value should be! In think, we need to compute
  // the eigenvalues of the matrix [+1,+1,-1; +1,-1,+1; -1,+1,+1] and then take the reciprocal of
  // the magnitude of the one with maximum magnitude. They are probably all the same in magnitude
  // anyway.
  //
  // The results with the rotation matrix look actually better - more complex and less regular.
  // Figure out why this is the case. It seems counterintutive.
    


  int N = numSamples;
  Vec y1(N), y2(N), y3(N);  // Outputs after the the A delaylines
  Vec z1(N), z2(N), z3(N);  // Outputs after the the B delaylines


  // Helper function:
  auto produceInternalSamples = [&](Real x, int n)
  {
    // Establish inputs to the A filters:
    Real u1 = x + b1 * B1.readOutput();
    Real u2 = x + b2 * B2.readOutput();
    Real u3 = x + b3 * B3.readOutput();
    // In a realtime implementation, we would probably just use "x + b1 * z[n-1]". But that 
    // wouldn't work for n=0 here. A realtime implementation would just keep a vector of z-values 
    // as state which could be read from in a getSample() call.

    // Apply A filters:
    y1[n] = a1 * A1.getSample(u1);
    y2[n] = a2 * A2.getSample(u2);
    y3[n] = a3 * A3.getSample(u3);

    // Apply feedback matrix to outputs of A filters:
    Real v1 = y1[n];
    Real v2 = y2[n];
    Real v3 = y3[n];
    F.apply(&v1, &v2, &v3);

    // Feed those into the B filters:
    z1[n] = b1 * B1.getSample(v1);
    z2[n] = b2 * B2.getSample(v2);
    z3[n] = b3 * B3.getSample(v3);
  };


  // Create all the internal signals and their sum:
  produceInternalSamples(1.0, 0);
  for(int n = 1; n < N; n++)
    produceInternalSamples(0.0, n);

  Vec ySum = c1*y1 + c2*y2 + c3*y3;
  Vec zSum = d1*z1 + d2*z2 + d3*z3;
  Vec sum  = ySum + zSum;

  // Plot the signals:
  //rosic::writeToMonoWaveFile("BiFDN_3x3.wav", &sum[0], N, 44100);
  //rsPlotVectors(sum, ySum, zSum);
  rsPlotVectors(sum);
  //rsPlotVectors(sum, y1, y2, y3, z1, z2, z3);

  return sum;


  // Observations:
  //
  // - It looks like we can create complexity pretty fast with this approach. It still sounds bad,
  //   though - but that is probably due to the unrealistically short delays and the fact that we 
  //   have just a 3x3 matrix in this toy example.
  //
  // - When setting the N values to all zeros, the spikes in zSum and ySum coincide. That means,
  //   using nonzero values (i.e. actually using the second delay layer rather than bypassing it)
  //   does indeed add complexity.
  //
  // - More experiments are needed to optimize the matrix and the delays. 
  //
  // - I think, it may make sense to have all M values odd and all n values even. Then we avoid any
  //   of the sums Mi+Nj to be even (because odd + even = odd). We want them all odd because 
  //   otherwise, two even sums have a gcd of 2. Although maybe a gcd of 2 is not as bad as a gcd 
  //   of 5, say. Ideally we want the gcd to be 1. If that's not possible, we may want it to be as
  //   small as possible?
  //
  // - Maybe when using real delay values (i.e. interpolated delaylines), it could make sense to 
  //   use power of the golden ratio phi (wrapped back into the interval 1..2 by dividing by an 
  //   appropriate power of 2)? Maybe look at the continued fraction expansions of the powers of
  //   phi.
  //
  // - Try it with more realistic values for the delays. Maybe write a 2nd function that does the
  //   same thing with 16 delaylines using a fast Kronecker trafo (FKT) instead of the explicit 
  //   feedback matrix.
  //
  // - But maybe when we use 2*16 = 32 delaylines, it may actually give us more bang for the buck
  //   to just use the 32 delaylines with a 32x32 matrix? The FKT would have to compute one level
  //   more. Maybe compare the results for a smaller case: 4x4 with 2nd round of delays vs 8x8 with
  //   a single round of delays. Maybe try to re-express the smaller case with a particular 
  //   feedback matrix of the bigger case. Or take an 8x8 network with 2 rounds vs a 16x16 network
  //   with 1 round of delays: Let's assume that evaluating 1 delayline takes a computational cost 
  //   of c and one elementary operation inside matrix computation takes a computational amount of 
  //   1. With the FKT, The 8x8 case with two delays would have a total cost of 
  //   8*c + 8*log2(8) + 8*c = 16*c + 8*3 = 16*c + 24. The 16x16 case with one round of delays 
  //   would take 16*c + 16*log2(16) = 16*c + 16*4 = 16*c + 64. It's probably realistic to assume 
  //   that c > 1, i.e. the evaluation of one delayline (including damping, possibly interpolation, 
  //   maybe additional allpass or notchpass, etc.) is more costly than the inner matrix operation 
  //   which we have normalized to have a cost of 1. OK - so the 16*c term is equal in both costs.
  //   They differ only the N*log2(N) term which is 24 in the 8x8 bi-delay case and 64 in the 16x16
  //   single delay case. So, indeed, the cost of a 8x8 bi-delay network is less. But maybe the 
  //   16x16 network will give better result for the higher cost? At the moment, we considered only
  //   cost, i.e. the buck. What about the bang? Maybe consider the total number of modes and/or 
  //   the initial echo density as a figure of merit. The first echoes in the 8x8 case arrive at 
  //   M1,..,M8 and then at M1+N1,...,M1+N8,.....,M8+N1,...,M8+N8. In the 16x16 case, they arrive 
  //   at M1,...,M16. If we assume that N1,..,N8 ins the 1st case are identical to M8,...,M16 in 
  //   the 2nd case, we would have more early echoes in the 2nd case. Maybe that means the 16x16 
  //   network is better. For the total number of modes, I need to write down the transfer 
  //   functions, I guess. In the 16x16 case, it's also easier to tune the delays to be all 
  //   mutually prime because we don't get these nasty sums
  //
  // - Maybe apply different types of waveshaping in the feedback path: 
  //   -Saturation: should make early decay faster
  //   -Antisaturation: should make late decay faster. But be careful: this may break stability
  //    when we do not limit the slope to unity.
  //   -Bitcrush: should noisify the tail, relative noisification is stronger for quiet signals
  //   -Floatcrush: should noisify the tail, relative nosification is independent from level
  //   -Soft-Bitcrush: should be a bitcrush with adjustable smoothing of the stairsteps into softer
  //    sigmoid shapes
  //
  // - Maybe use a function x + a*x^3 with adjustable a > 0 and its inverse (used with |a| when 
  //   a < 0). This waveshaping in the feedback path could be used to dial in a balance between
  //   early and late decay, I guess. But maybe for the a > 0 case, we need to limit the slope of
  //   the function somehow to keep the feedback loop gain stable.
}

// Move to somwhere else - maybe TestUtilities.h/cpp or something
template<class T>
void rsRotationMatrixFromEulerAngles(T rx, T ry, T rz, rsMatrix<T>* R)
{
  rsAssert(R->hasShape(3,3), "Matrix needs to be 3x3 in rsRotationMatrixFromEulerAngles");

  // Sines/cosines:
  T sx = sin(rx); T cx = cos(rx);
  T sy = sin(ry); T cy = cos(ry);
  T sz = sin(rz); T cz = cos(rz);

  // Rotation matrix coeffs:
  (*R)(0,0) =  cz*cy;
  (*R)(0,1) = -sz*cx + cz*sy*sx;
  (*R)(0,2) =  sz*sx + cz*sy*cx;
  (*R)(1,0) =  sz*cy;
  (*R)(1,1) =  cz*cx + sz*sy*sx;
  (*R)(1,2) = -cz*sx + sz*sy*cx;
  (*R)(2,0) = -sy;
  (*R)(2,1) =  cy*sx;
  (*R)(2,2) =  cy*cx;

  // See also: rsRotationXYZ<T>::updateCoeffs(). We actually duplicate the code from there here.
  // Maybe try to refactor to get rid of that duplication.
}

void protoFDN1()
{
  // We reproduce the results from the function prePostDelayFDN_3x3 above but this time using the 
  // rsProtoFDN class. 

  bool ok = true;

  using Real    = double;
  using Complex = rsComplex<Real>;
  using VecI    = std::vector<int>;
  using VecR    = std::vector<Real>;
  using MatC    = rsMatrix<Complex>;
  using FDN     = rsProtoFDN<Real, Real>;
  using TF      = rsSparseTransferFunction<Real, Real>;

  int numSamples = 1000;
  int numChans   = 3;

  // Decay time (RT60) in samples:
  Real decay = 1000;

  // Create the vector of delay values:
  VecI delays( { 17, 23, 29 });

  // Create the feedback matrix:
  Real rx, ry, rz;
  rx = ry = rz = 45;
  //rx = ry = rz = 0;                     // Leads to identity matrix for feedback.
  Real toRad = PI/180;
  rsMatrix<Real> fbMatrix(numChans, numChans);
  rsRotationMatrixFromEulerAngles(toRad*rx, toRad*ry, toRad*rz, &fbMatrix);

  // Create the input and output matrices:
  rsMatrix<Real> inMatrix( numChans, 1, {+1, +1, +1});
  rsMatrix<Real> outMatrix(1, numChans, {+1, -1, +1});
  //rsMatrix<Real> inMatrix( numChans, 1, {+1.1, +1.3, +1.2});
  //rsMatrix<Real> outMatrix(1, numChans, {+1.4, -1.2, +1.1});

  // Compute the damping factors from the desired decay time and delay lengths:
  VecR dampFactors(numChans);
  Real amp = Real(0.001);
  for(int i = 0; i < numChans; i++)
    dampFactors[i] = rsDecayTimeToFeedbackGain(decay, Real(delays[i]), amp);

  // Create and set up the FDN:
  FDN fdn;
  fdn.setFeedbackMatrix(fbMatrix);
  fdn.setDelays(        delays);
  fdn.setInputMatrix(   inMatrix);
  fdn.setOutputMatrix(  outMatrix);
  fdn.setDampFactors(   dampFactors);
  rsAssert(fdn.areSettingsConsistent());  // Sanity check
  // We can set the settings in any order but it's important that after calling all the setters, 
  // all the vectors and matrices in the FDN are consistent with respect to their sizes and shapes.
  // They may be inconsistent at an intermediate stage, i.e. in between the calls.

  // Create input impulse signal:
  int N = numSamples;
  VecR x(N);
  x[0] = 1.0;

  // Helper function to produce one output sample at a time. We need it because the API of the FDN 
  // class uses std::vector for supporting multichannel I/O:
  auto getSample = [&](Real in)
  {
    VecR tmpIn(1), tmpOut(1);
    tmpIn[0] = in;
    fdn.processFrame(tmpIn, tmpOut);
    return tmpOut[0];
  };

  // Produce FDN output signal, i.e. the impulse response:
  VecR y(N);
  y[0] = getSample(1.0);
  for(int n = 1; n < N; n++)
    y[n] = getSample(0.0);

  // Produce the impulse response using the convenience function getSample():
  VecR y2 = impulseResponse(fdn, N, 1.0);
  VecR err = y2 - y;
  ok &= y2 == y;

  // Compute transfer function numerically:
  Complex z(0.9, 0.8);
  Complex Hn = rsEvaluateTransferFunctionNumerically(fdn, z, 5000);

  // Evaluate transfer function via getTransferFunctionAt:
  MatC    H  = fdn.getTransferFunctionAt(z);
  Complex errH = H(0,0) - Hn;
  ok &= rsIsNegligible(errH, 1.e-13);        // new
  //ok &= rsIsCloseTo(H(0,0), Hn, 1.e-13);   // old

  // Create symbolic transfer function and evaluate that at z:
  //rsMatrix<TF> tf = fdn.getTransferFunction(1.e-13);
  rsMatrix<TF> tf = fdn.getTransferFunction(0.0);
  Complex H2 = tf(0,0)(z);
  errH = H2 - Hn;
  //ok &= rsIsNegligible(errH, 1.e-13);  // FAILS!!!
  // Apparently, the Gaussian elimination algo with matrices of transfer functions does not yet 
  // work right. 
  // ToDo: Implement unit tests for linear algebra with matrices of rational functions of type
  // rsSparseRationalFunction and rsSparseTransferFunctions - maybe also rsSparsePolynomial, if
  // possible.



  // Plot the generated signal together with the reference signal:
  //rsPlotVectors(y, y2);
  rsAssert(ok);
}

void protoFDNvsSSF()
{
  // We compare the simplemost case of the rsProtoFDN with just one single channel and with unit 
  // delay against the state space filter (SSF) implementation. FDNs are generalizations of SSFs in
  // the sense that if you just set all the delays to unit delays in an FDN, you get an SSF. We 
  // want to compare the computation of the transfer function of both classes because the one in 
  // the FDN seems to be buggy and the one in the SSF is known to be working. I hope to find the 
  // bug in the FDN that way. ...TBC...

  bool ok = true;

  using Real    = double;
  using Complex = rsComplex<Real>;
  using VecI    = std::vector<int>;
  using VecR    = std::vector<Real>;
  using MatR    = rsMatrix<Real>;
  using MatC    = rsMatrix<Complex>;
  using FDN     = rsProtoFDN<Real, Real>;
  using SSF     = rsStateSpaceFilter<Real>;

  // Setup:
  int  N        = 100;
  Real feedback = 0.8;
  Real inGain   = 1.0;
  Real outGain  = 1.0;
  Real thruGain = 0.0;

  // Create the 1-vectors and 1x1-matrices needed to set up the FDN and SSF:
  VecI delays(     1,    { 1        });
  VecR dampFactors(1,    { 1.0      });   // Unity because damping is done via the matrix here.
  MatR fbMatrix(   1, 1, { feedback });
  MatR inMatrix(   1, 1, { inGain   });
  MatR outMatrix(  1, 1, { outGain  });
  MatR thruMatrix (1, 1, { thruGain });

  // Create and set up the FDN:
  FDN fdn;
  fdn.setFeedbackMatrix(fbMatrix);
  fdn.setDelays(        delays);
  fdn.setInputMatrix(   inMatrix);
  fdn.setOutputMatrix(  outMatrix);
  fdn.setDampFactors(   dampFactors);
  rsAssert(fdn.areSettingsConsistent());  // Sanity check

  // Create and set up the SSF:
  SSF ssf;
  ssf.setup(fbMatrix, inMatrix, outMatrix, thruMatrix);

  // Produce and plot impulse responses of both filters:
  VecR h_FDN = impulseResponse(fdn, N, 1.0);
  VecR h_SSF = impulseResponse(ssf, N, 1.0); 
  //rsPlotVectors(h_FDN, h_SSF);

  // Compute transfer functions at a given z of SSF and FDN:
  Complex z(0.9, 0.7);
  MatC    H_SSF   = ssf.getTransferFunctionAt(z);
  MatC    H_FDN   = fdn.getTransferFunctionAt(z);
  Complex H_FDN_n = rsEvaluateTransferFunctionNumerically(fdn, z, N);
  //ok &= rsIsCloseTo(H_SSF,      H_FDN,   1.e-13); // Doesn't compile
  ok &= rsIsCloseTo(H_FDN(0,0), H_FDN_n, 1.e-13);

  // Now let thr FDN do the damping via the actual damping coeff:
  dampFactors[0] = feedback;
  fbMatrix(0, 0) = 1.0;
  fdn.setFeedbackMatrix(fbMatrix);
  fdn.setDampFactors(   dampFactors);
  VecR h_FDN2 = impulseResponse(fdn, N, 1.0);
  H_FDN   = fdn.getTransferFunctionAt(z);
  H_FDN_n = rsEvaluateTransferFunctionNumerically(fdn, z, N);
  ok &= rsIsCloseTo(H_FDN(0,0), H_FDN_n, 1.e-13);
  //rsPlotVectors(h_FDN, h_SSF, h_FDN2);
  // The impulse response is now scaled by the damping factor compared to what it was before. I 
  // think, that is the correct behavior. Or is it? Verify!

  rsAssert(ok);
}

void extendedProtoFDN1()
{
  // We reproduce the results from the function prePostDelayFDN_3x3 above but this time using the 
  // rsProtoFDN class. 

  bool ok = true;

  using Real = double;
  using VecI = std::vector<int>;
  using VecR = std::vector<Real>;
  using EFDN = rsExtendedProtoFDN<Real, Real>;

  int numSamples = 1000;
  int numChans   = 3;

  // Decay time (RT60) in samples:
  Real decay = 1000;

  // Create the vectors of delay values:
  VecI delaysPre( { 17, 23, 29 });
  VecI delaysPost({ 20, 26, 36 });

  // Create the feedback matrix:
  Real rx, ry, rz;
  rx = ry = rz = 45;
  Real toRad = PI/180;
  rsMatrix<Real> fbMatrix(numChans, numChans);
  rsRotationMatrixFromEulerAngles(toRad*rx, toRad*ry, toRad*rz, &fbMatrix);

  // Create the input and output matrices:
  rsMatrix<Real> inMatrixPre(  numChans, 1, {+1, +1, +1});
  rsMatrix<Real> inMatrixPost( numChans, 1, { 0,  0,  0});
  rsMatrix<Real> outMatrixPre( 1, numChans, {+1, -1, +1});
  rsMatrix<Real> outMatrixPost(1, numChans, {-1, +1, -1}); 

  // Compute the damping factors from the desired decay time and delay lengths:
  VecR dampFactorsPre( numChans);
  VecR dampFactorsPost(numChans);
  Real amp = Real(0.001);
  for(int i = 0; i < numChans; i++)
  {
    dampFactorsPre[i]  = rsDecayTimeToFeedbackGain(decay, Real(delaysPre[i]),  amp);
    dampFactorsPost[i] = rsDecayTimeToFeedbackGain(decay, Real(delaysPost[i]), amp);
  }

  // Create and set up the FDN:
  EFDN efdn;
  efdn.setFeedbackMatrix(fbMatrix);
  efdn.setDelays(        delaysPre,      delaysPost);
  efdn.setInputMatrices( inMatrixPre,    inMatrixPost);
  efdn.setOutputMatrices(outMatrixPre,   outMatrixPost);
  efdn.setDampFactors(   dampFactorsPre, dampFactorsPost);
  rsAssert(efdn.areSettingsConsistent());                     // Sanity check
  // We can set the settings in any order but it's important that after calling all the setters, 
  // all the vectors and matrices in the FDN are consistent with respect to their sizes and shapes.
  // They may be inconsistent at an intermediate stage, i.e. in between the calls.

  // Create input impulse signal:
  int N = numSamples;
  VecR x(N);
  x[0] = 1.0;

  // Helper function to produce one output sample at a time. We need it because the API of the FDN 
  // class uses std::vector for supporting multichannel I/O:
  auto getSample = [&](Real in)
  {
    VecR tmpIn(1), tmpOut(1);
    tmpIn[0] = in;
    efdn.processFrame(tmpIn, tmpOut);
    return tmpOut[0];
  };

  // Produce FDN output signal, i.e. the impulse response:
  VecR y(N);
  y[0] = getSample(1.0);
  for(int n = 1; n < N; n++)
    y[n] = getSample(0.0);

  // Produce reference signal with code from the previous experiment:
  VecR yr  = extendedFDN_3x3();  // This call also fires up a plot which is a bit inconvenient here.
  VecR err = yr - y;
  ok &= rsIsCloseTo(y, yr, 1.e-16);  // There's a tiny roundoff error! Why?

  // Plot the generated signal together with the reference signal:
  //rsPlotVectors(y, yr, err);
  rsAssert(ok);


  // ToDo:
  //
  // - Add damping factors to the rsProtoFDN class. The computation of them should be left to 
  //   client code for flexibility.
  //
  // - Include a self-test function areSettingsConsistent() and call it at the start of 
  //   processFrame(). It should check, if the input and output matrices, damping factor arrays,
  //   feedback matrix, etc. all have lengths and shapes that fit together.
  //
  // - Figure out, how a 3x3 FDN with post-matrix delaylines can be re-expressed as 6x6 FDN without
  //   it. I think, if we arrange the vectors u,v of pre- and post- delayline outputs as single 
  //   vector by concatenation, the 6x6 feedback matrix would have to look like in this update 
  //   equation:
  //
  //     u1     0   0   0   a11 a12 a13     u1
  //     u2     0   0   0   a21 a22 a23     u2
  //     u3  =  0   0   0   a31 a32 a33  *  u3
  //     v1     1   0   0    0   0   0      v1d
  //     v2     0   1   0    0   0   0      v2d
  //     v3     0   0   1    0   0   0      v3d
  //
  //   where the top-right 3x3 submatrix would be the feedback matrix used in the 3x3 FDN. Here,
  //   the d in v1d, etc means (unit) delayed. Or maybe this submatrix should sit in the bottom 
  //   left (swapped with the current bottom left)? Hmm - yes - I think so. The current u-vector 
  //   just gets the old v-vector (from the previous sample) fed back whereas the current v-vector 
  //   is obtained from the current u-vector via the feedback matrix. So, yeah - I think, that 
  //   means we need to swap the top-right and bottom-left sub-matrices. Try it in practice!
  //   Implement a class rsProtoFDN and let rsExtendedProtoFDN have a method that converts the
  //   extended structure into an equivalent basic structure (of twice the order). The API should
  //   be something like:  rsProtoFDN rsExtendedProtoFDN::getAsBasicFDN()  .or maybe it should take
  //   the basic FDN as pointer parameter and be called convertToBasicFDN() or something.
}

void feedbackDelayNetworks()
{
  protoFDN1();
  protoFDNvsSSF();


  extendedFDN_3x3();
  protoFDN1();
  protoFDNvsSSF();
  extendedProtoFDN1();

  // Notes:
  //
  // - There are some more, older FDN experiments in rosic_EffectsTest.cpp. Maybe drag them over to
  //   here. there is also allpassFDN() in FilterExperiments.cpp.
  //
  //
  // Resources:
  //
  // - Modal Decomposition of Feedback Delay Networks  https://arxiv.org/pdf/1901.08865
  //
  // - Choosing Optimal Delays for Feedback Delay Networks 
  //   https://pub.dega-akustik.de/DAGA_2014/data/articles/000025.pdf
  //
  // - Optimization Approaches for Feedback Delay Networks
  //   https://www.diva-portal.org/smash/get/diva2:1668199/FULLTEXT01.pdf
  //
  // - DIFFERENTIABLE FEEDBACK DELAY NETWORK FOR COLORLESS REVERBERATION
  //   https://www.dafx.de/paper-archive/2023/DAFx23_paper_32.pdf
  //
  // - FDNTB: THE FEEDBACK DELAY NETWORK TOOLBOX
  //   https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_53.pdf
  //
  // - Feedback Delay Network Optimization  https://arxiv.org/pdf/2402.11216v1
}
