// todo: SinusoidalSynthesizer to RAPT to make it available here

using namespace RAPT;


// convenience function:
std::vector<double> synthesizeSinusoidal(
  const RAPT::rsSinusoidalModel<double>& model, double sampleRate, double fadeTime)
{
  std::cout << "Synthesizing sinusoidal model with " << model.getNumPartials() << " partials...";
  typedef RAPT::rsSinusoidalSynthesizer<double> SS;
  typedef SS::PhaseInterpolationMethod PIM;
  SS synth;
  synth.setSampleRate(sampleRate);

  //synth.setCubicAmplitudeInterpolation(true);

  //synth.setPhaseInterpolation(PIM::tweakedFreqIntegral);
  //synth.setPhaseInterpolation(PIM::cubicHermite);
  synth.setPhaseInterpolation(PIM::quinticHermite);

  std::vector<double> x = synth.synthesize(model);
  if(fadeTime > 0.0)
    applyFadeInAndOut( &x[0], (int) x.size(), int (fadeTime*sampleRate));
  std::cout << "Done\n";
  return x;
}

void setupHarmonicAnalyzerFor(RAPT::rsHarmonicAnalyzer<double>& analyzer, 
  const std::string& sampleName, double fs, double f0)
{
  // First, set up some standard settings:

  typedef rsWindowFunction::WindowType WT;
  analyzer.setSampleRate(fs);
  analyzer.setSincInterpolationLength(64);
  analyzer.setNumCyclesPerBlock(4);
  //analyzer.setWindowType(WT::hamming);
  analyzer.setWindowType(WT::blackman);
  analyzer.setSpectralOversampling(8);  // zero padding - maybe use 4 - 8 seems a bit excessive
  analyzer.setAllowInharmonics(true);
  analyzer.setSpectralPeakSearchWidth(0.5);       // default: 1 - blackman needs a value less than 1
  analyzer.setMinPeakToMainlobeWidthRatio(0.75);  // default: 0.75


  //analyzer.setMinPeakToMainlobeWidthRatio(2.5);
  //analyzer.setMinPeakToMainlobeWidthRatio(0.125);
  //analyzer.setMinPeakToMainlobeWidthRatio(0.0);
  // test - trying to reduce gap-artifacts - smaller values for this should reduce the gappiness - 
  // todo:
  // figure out, how low we can go without introducing other artifacts due to picking up on 
  // sidelobes - wtf - this doesn't seem to have much influence
  // -even with a zero value, we still get gaps with the rhodes - is there some other criterion 
  //  that discards the partial?
  // -even if we immediately retrun true in isPeakPartial, we get a gap - but it's narrower
  //  -> check all branches that return -1 in findPeakBinNear
  // -test this as follows: create a sine with strong beating - the amp-envelope of the partial 
  //  that shows this artifact really drops to almost zero - maybe use sines at 100, 199, 201, 300

  //analyzer.setFreqsByPhaseDerivative(true);
  //analyzer.setFreqPhaseConsistency(true);
  // todo: maybe provide different freq-refinement methods (not necessarily mutually exclusive)

  // set up settings of the embedded cycle-mark finder:
  rsCycleMarkFinder<double>& cmf = analyzer.getCycleFinder();
                                              // defaults:
  cmf.setRelativeBandpassWidth(0.5);          // 1.0
  cmf.setBandpassSteepness(5);                // 3
  cmf.setSubSampleApproximationPrecision(2);  // 1, 0: linear, 1: cubic, 2: quintic, ...
  cmf.setFundamentalRange(50., 1000.);        // 20, 5000 
  cmf.setFundamental(f0);                     // 0 -> auto-detect
  //cmf.setAlgorithm(cmf.F0_ZERO_CROSSINGS);
  cmf.setAlgorithm(cmf.CYCLE_CORRELATION);


  // Second: override some settings for specific samples:

  if(sampleName == "flute-C-octave2")
  {

  }
  else if(sampleName == "piano_E2")
  {
    analyzer.setNumCyclesPerBlock(8);
    analyzer.setSpectralOversampling(4);
    analyzer.setSpectralPeakSearchWidth(1.0);
    analyzer.setMinPeakToMainlobeWidthRatio(0.5);
  }
  else if(sampleName == "Rhodes_F3")
  {
    //cmf.setAlgorithm(cmf.F0_ZERO_CROSSINGS);
    // because CYCLE_CORRELATION hangs - but why?
  }
  // etc...

}

// maybe it will soon make sense to wrap this into a class
void testHarmonicResynthesis(const std::string& name, std::vector<double>& input, 
  double fs, double f0, bool writeWaveFiles, bool plotResults)
{
  // analyze, resynthesize and create resiudal signal:


  double* x = &input[0];   // pointer to first sample (for convenience)
  int Nx = (int) input.size();


  // Analysis:

  RAPT::rsHarmonicAnalyzer<double> analyzer;
  setupHarmonicAnalyzerFor(analyzer, name, fs, f0);
  RAPT::rsSinusoidalModel<double> mdl = analyzer.analyze(x, Nx);

  rsAssert(mdl.isDataValid());


  // Manipulations:

  //mdl.removePartial(0);                      // remove DC
  mdl.removePartialsWithMeanFreqAbove(fs/2); // anti-alias

  //mdl.keepOnly({0, 9});  // for test with TwoSines_Freq1=200_Freq2=2025
  //mdl.removePartial(0);    // test: resynthesize without fundamental

  //mdl.keepOnly({0, 99});  // fom TwoSines_Freq1=100_Freq2=10020

  //mdl.removePartialsAbove(9); 
  // for testing artifacts (spurious high partials) in flute-c1 - are they analysis or synthesis 
  // artifacts? if they occur evene though we remove the partials above 9, they occur in synthesis
  // ->they seem to be analysis artifacts

  //double splitFreq = 1000;
  //RAPT::rsSinusoidalModel<double> lp = rsSinusoidalProcessor<double>::extractLowpassPart(mdl, splitFreq);
  //RAPT::rsSinusoidalModel<double> hp = rsSinusoidalProcessor<double>::extractHighpassPart(mdl, splitFreq);


  // Resynthesis:

  typedef RAPT::rsSinusoidalSynthesizer<double> SS;
  typedef SS::PhaseInterpolationMethod PIM;
  SS synth;
  synth.setSampleRate(fs);
  //synth.setCubicAmplitudeInterpolation(true);
  synth.setPhaseInterpolation(PIM::tweakedFreqIntegral);
  //synth.setPhaseInterpolation(PIM::cubicHermite);
  //synth.setPhaseInterpolation(PIM::quinticHermite);
  //synth.setPhaseInterpolation(PIM::linear);  // does not yet exist
  std::vector<double> output = synth.synthesize(mdl);

  std::vector<double> error = output-input;
  double* y = &output[0]; int Ny = (int) output.size(); // again, for convenience
  double* e = &error[0];  int Ne = (int) error.size();  // dito




  // Plotting and file output:

  // write original, resynthesized and error signals to files, if desired:
  if(writeWaveFiles == true) {
    rosic::writeToMonoWaveFile((name + "Original.wav").c_str(),      x, Nx, (int)fs);
    rosic::writeToMonoWaveFile((name + "Resynthesized.wav").c_str(), y, Ny, (int)fs);
    rosic::writeToMonoWaveFile((name + "Residual.wav").c_str(),      e, Ne, (int)fs);
  }


  // move the required functions to rs_testing module
  // plot model data, if desired....
  bool plotModel = false;      // make user parameter
  if(plotModel == true)
    plotSineModel(mdl, fs);


  // plot original, resynthesized and error signals, if desired:
  if(plotResults == true) {
    GNUPlotter plt;
    plt.addDataArrays(Nx, x);
    plt.addDataArrays(Ny, y);
    plt.addDataArrays(Ne, e);

    std::vector<double> marks = analyzer.getOriginalTimeStamps();
    std::vector<double> zeros(marks.size());    // y values for plotting (all zero)
    RAPT::rsArrayTools::fillWithZeros(&zeros[0], (int) marks.size());
    plt.addDataArrays((int) marks.size(), &marks[0], &zeros[0]);

    plt.setGraphStyles("lines", "lines", "lines", "points");
    plt.setPixelSize(1000, 300);
    plt.plot();
  }
}

void testMakeHarmonic(const std::string& name, std::vector<double>& input,
  double fs, double f0Out, double inharmonicity,  double f0In)
{
  //double* x = &input[0];   // pointer to first sample (for convenience)
  //int Nx = (int) input.size();

  // analyze:

  RAPT::rsHarmonicAnalyzer<double> analyzer;
  setupHarmonicAnalyzerFor(analyzer, name, fs, f0In);
  RAPT::rsSinusoidalModel<double> mdl = analyzer.analyze(&input[0], (int) input.size());

  // plotSineModel(mdl, fs); // move to rapt

  // process model data:

  mdl.removePartial(0); 
  //mdl.keepOnly({ 0 });  // for test
  rsSinusoidalProcessor<double>::makeStrictlyHarmonic(mdl, f0Out, inharmonicity, 0.5);
  mdl.removePartialsWithMeanFreqAbove(fs/2); // anti-alias - maybe use min- insetad of mean-freq to be safe...


  //plotSineModel(mdl, fs);

  // (re)synthesize:
  typedef RAPT::rsSinusoidalSynthesizer<double> SS;
  typedef SS::PhaseInterpolationMethod PIM;
  SS synth;
  synth.setSampleRate(fs);
  synth.setCubicAmplitudeInterpolation(true);
  //synth.setPhaseInterpolation(PIM::tweakedFreqIntegral);
  synth.setPhaseInterpolation(PIM::cubicHermite);
  std::vector<double> output = synth.synthesize(mdl);

  // due to new phase-relationships, the maximum output amplitude may be different from the 
  // original sound - renormalize:
  rsArrayTools::normalize(&output[0], (int) output.size(), 1.0);


  std::string name2 = name + "_F=" + std::to_string((int)f0Out) + "Hz";
  if(inharmonicity != 0)
    name2 += "_B=" + std::to_string(inharmonicity);
  name2 += ".wav";
  // the formatting of floating point numbers std::to_string sucks -> write a better rsToString 
  // function
  rosic::writeToMonoWaveFile(name2.c_str(), &output[0], (int) output.size(), (int)fs);
}

void getPaddedSignals(double* xIn, int Nx, 
  const RAPT::rsSinusoidalModel<double>& model,
  const RAPT::rsSinusoidalSynthesizer<double>& synth,
  std::vector<double>& x, std::vector<double>& y)
{
  typedef std::vector<double> Vec;

  // synthesize y from the model:
  y = synth.synthesize(model);
  int nf = model.getStartSampleIndex(synth.getSampleRate()); // # fade-in samples

  int n;
  if(nf < 0) {
    // obtain a version of x with an appropriate number of zeros prepended
    nf = -nf;
    x.resize(Nx+nf);
    for(n = 0; n < nf; n++)        x[n] = 0;
    for(n = nf; n < (int) x.size(); n++) x[n] = xIn[n-nf];
  }
  else {
    RAPT::rsPadLeft(y, nf, 0.0);       // prepend zeros to y, if nf > 0
    x = RAPT::toVector(xIn, Nx);
  }

  // extend the shorter of both signals with zeros:
  size_t nx = x.size();
  size_t ny = y.size();
  if(nx > ny)
    RAPT::rsResizeWithInit(y, nx, 0.0);
  else if(ny > nx)
    RAPT::rsResizeWithInit(x, ny, 0.0);
}

void testModalResynthesis(const std::string& name, std::vector<double>& x,
  double fs, double f0)
{
  int N = (int) x.size();

  RAPT::rsHarmonicAnalyzer<double> sineAnalyzer;
  setupHarmonicAnalyzerFor(sineAnalyzer, name, fs, f0);
  RAPT::rsSinusoidalModel<double> sineModel = sineAnalyzer.analyze(&x[0], N);
  sineModel.removePartial(0);  // we don't want DC
  //plotSineModel(sineModel, fs);
  //std::vector<double> ys = synthesizeSinusoidal(sineModel, fs);
  //rosic::writeToMonoWaveFile("ModalSineOutput.wav", &ys[0], (int) ys.size(), (int)fs);
  // ...maybe also remove all partials whose (mean) frequency is above the nyquist freq

  rsModalAnalyzer<double> modeAnalyzer;
  std::vector<rsModalFilterParameters<double>> modeModel
    = modeAnalyzer.getModalModel(sineModel);
  plotModalAmplitudes(modeModel);
  std::vector<double> y = synthesizeModal(modeModel, fs, N);

  rosic::writeToMonoWaveFile("ModalInput.wav",  &x[0],  N, (int)fs);
  rosic::writeToMonoWaveFile("ModalOutput.wav", &y[0],  N, (int)fs);
}

void testDeBeating(const std::string& name, std::vector<double>& x, double fs, double f0)
{
  // for reference, write input into a wavefile:
  rosic::writeToMonoWaveFile(name + "DeBeatInput.wav",  &x[0], (int)x.size(), (int)fs);

  // create and set up analyzer and obtain sinusoidal model and write a resynthesized unmodified
  // signal into a wavefile:
  int N = (int) x.size(); // x is the input signal

  RAPT::rsHarmonicAnalyzer<double> analyzer;

  // temporary - to make experimentation faster:
  //analyzer.setMinPartialIndex(0);  // not yet working
  analyzer.setMaxPartialIndex(6);  // for quicker tests

  //analyzer.getCycleFinder().setAlgorithm(rsCycleMarkFinder<double>::F0_ZERO_CROSSINGS);
  // for test with Rhodes Tuned F3 V12TX -16.4 10-17-16 shorter

  setupHarmonicAnalyzerFor(analyzer, name, fs, f0);
  std::cout << "Analyzing...\n";  // maybe write function that prints the time it took to analyze
                                  // takes a string and a function to call - callAndEcho
  RAPT::rsSinusoidalModel<double> mdl = analyzer.analyze(&x[0], N);
  rsAssert(mdl.isDataValid());
  //plotSineModel(mdl, fs);
  //plotSineModelAmplitudes(mdl, {1,2,3,4,5,6});
  //std::cout << "Resynthesizing...\n";
  std::vector<double> y = synthesizeSinusoidal(mdl, fs);
  rosic::writeToMonoWaveFile(name + "DeBeatOutputUnmodified.wav", &y[0], (int)y.size(), (int)fs);

  // mdl now contains the sinusoidal model for the sound. Now, we set up a rsPartialBeatingRemover 
  // and apply the de-beating to model data, then synthesize the de-beated signal and write into
  // wavefile:
  rsPartialBeatingRemover<double> deBeater;
  deBeater.setPhaseSmoothingParameters(5.0, 1, 4); // cutoff = 10 leaves a vibrato

  // testing the new peak-picker:
  deBeater.envExtractor.peakPicker.setShadowWidths(0.0, 0.5);
  //deBeater.envExtractor.peakPicker.setWorkOnPeaksOnly(true); // default is false
  //deBeater.envExtractor.peakPicker.setShadowWidths(0.0, 2.7);
  deBeater.envExtractor.peakPicker.setNumNeighbors(3);


  deBeater.envExtractor.setMaxSampleSpacing(0.5);
  //deBeater.envExtractor.setInterpolationMode(
  //  rsInterpolatingFunction<double, double>::CUBIC_HERMITE); // cubic interpolation produces overshoot

  //deBeater.envExtractor.maxSpacing = 0.5;  // temporary, during development
  // -the optimal settings need to be figured out
  // -i think, we may have to set it up per partial and not just use one set of settings for all 
  //  partial
  // -i think, the decay-rate of the shadowing algorithm needs to be less (but not much less) than
  //  the overall decay-rate of the envelope - maybe we should estimate it from the area or energy
  //  in the signal after the peak - but what about non-decaying sounds?
  // -rhodes, 1st partial: 
  //  -decay should be in the range 2.6 < width < 2.7 (note that this is the half-height width (in 
  //   seconds) - may need to be converted to the 1/e width), with decay-width = 3.0, the 
  //   sluggishness can be clearly seen
  //  -the release shoud be less than that - could even be zero
  // -rhodes, 2nd - 4th partial: 0.7 < width < 0.8
  // -rhodes, 5th and 6th partial - something around 0.3
  // -i think, the 6th partial is a good example to try the pre-filtering by DC-blocking idea
  // -the 0th partial (DC) is again special - reasonable settings are hard to find here because it
  //  doesn't really behave like an envelope - more like noise with a big spike at the beginning, 
  //  the noise even going negative - but we may not have to waste too much time with it - in the 
  //  end, we may not want to resynthesize it anyway
  // -in order to not be too sluggish, i think, the decay-time constant for rightward shadows 
  //  should be less that the decay time of the envelope - we may have to estimate that - maybe the
  //  user parameter should be something like a scaler for the estimated decay time - or maybe let 
  //  the user decide between and absolute, fixed setting or a setting based on estimated decay 
  //  time - the fixed setting may include some factor by which the shadow widths decreases with
  //  increasing partial frequency - maybe a power rule could be suitable (i.e. 
  //  width = baseWidth * (freq/fundamental)^p where p is the adjustable power/exponent
  // -maybe to estimate the decay of the envelope, we could use the code to fit the exponential 
  //  envelope (where we still have to select the second estimation-point - maybe that could be 
  //  automated as well)
  // -the leftward shadow-widths may be set to some scale factor times the rightward shadow - in 
  //  general, we want them to be smaller - even a factor of 0 could make sense
  // -maybe we should indeed implement the idea of removing "undulating DC", i.e. apply a highpass
  //  at sub-tremolo frequencies (maybe 0.5 Hz or something)
  // -maybe try prominence thresholding, too (i've almost forgotten that we can do that, too)

  // -why are the final estimated envelopes still bad? the densification does not seem to work yet
  //  ...is it still commented out? when we use a too large setting for the widths, the beating 
  //  re-appears - maybe it's some post-processing step (no-stickout something) that adds the peaks
  //  from the original envelope - nope - it's the densification - apparently, it uses a very small
  //  maxSpacing in these cases
  // -maxSpacing = 4.43 in rsEnvelopeExtractor<T>::fillSparseAreasNew with the 1st partial - that
  //  distance is far too large, if it's in seconds
  // -so, sometimes the densification uses a too large value and sometimes a too small one
  // -estimating the densification parameter maxSpacing does not yet work properly
  // -maybe let the user choose a fixed maxSpacing (in seconds) for the moment, while trying to 
  //  make it adaptive
  // -maybe we somehow need to detect, if we are in a decaying section and densify only there
  // -this could be done by checking, if the env is monotonic between two peaks



  //---------------------------------------------------------------------------------------
  // experimental - set up the minimum distance between de-beated amplitude envelope peaks
  // 10 is ad hoc - at least one sample per 10 cycles:
  if(f0 == 0) 
    f0 = mdl.getPartial(1).getMeanFreq();

  // only during development - look only at a single partial - comments relate to the sample
  // Rhodes_F3_Medium.wav in the test repo

  // partial 1 exposes the problem of not de-beating at all
  //mdl.keepOnly({1});  // the fundamental exposes the not-debeating-at-all problem
  //deBeater.envExtractor.setMaxSpacingMultiplier(1.0);
  // maxSpacing = 0.32 when the multiplier is 1 - maybe we also need a min-spacing 
  // -the problem is that the valleys in between the relevant maxima show tiny local maxima so the
  //  envelope gets sampled at these irrelevant local maxima, too 
  //  -to fix it, maybe we need a stricter criterion for what is considered a relevant maximum 
  //   and/or we look for maxima not in the envelope itself but in a smoothed version of it

  // partial 2 exposes the "straight-line" problem
  //mdl.keepOnly({2});
  // ...maybe the max-spacing is too large? figure out what fillSparseAreas does in this case
  // -yep: maxSapcing is 4.74 seconds - faaar too large 
  // -there is some hardly visible shallow maximum toward the end of the envelope which is picked 
  //  up by the peak-finder
  // -in both cases, it is the presence of spurious maxima, that peak finder should ignore but 
  //  doesn't
  // -how can we prevent the peak-picker from picking up these peaks undesired peaks?
  // -we need additional criteria - but which ones? 
  //  -maybe, we should use smoothing before the peak-picker?
  //  -maybe, we should accept a peak only, if it's greater than two neighbours left and right?
  //   -maybe that should be a user-option: accept maxima only, if the are larger than N_l 
  //    neighbours to the left and N_r neighbours to the right with default setting: 2,2
  //    or maybe 3,2 - i think, to the left, we should include more samples
  //  -maybe make a class rsPeakPicker that let's the user set up various criteria for considering
  //   a peak relevant ...maybe there should also be some threshold by which it must be above a 
  //   local average - and the size of that average can be set by the user


  // maybe, we should give the user the option to manually set a maximum spacing that is used in
  // addition to the computed one (the algo should use the minimum of both values)

  //mdl.keepOnly({2});


  //deBeater.setMaxEnvelopeSampleSpacing(16.0/f0);     // works reasonably
  //deBeater.setMaxEnvelopeSampleSpacing(1.0/f0);    // hangs
  //deBeater.setMaxEnvelopeSampleSpacing(2.0/f0);    // works but doesn't actually de-beat
  // ...maybe use the maximum distance between the found peaks as the minimum distance between 
  // datapoints in the de-beated envelope

  // that is wrong - we need to set it to a value a bit above the beat-period expressed in the 
  // frame-rate ...or - wait - is that actually true

  //deBeater.removeAmplitudeBeating(mdl.getModifiablePartialRef(5)); //for debugging

  // Ideas:
  // -maybe we should have a minimum allowed disatance (to avoid the straight-lines problem) and
  //  a maximum distance (to ensure that it does some de-beating at all)
  // -the minimum allowed distance should

  // end of exp experimental code to set up minimum distance of amp-env samples
  //-------------------------------------------------------------------------------------



  //mdl.removePartial(0);  // test - remove DC - the DC component crashes with Rhodes Tuned F3 V12TX -16.4 10-17-16 short
  std::cout << "De-Beating...\n";
  deBeater.processModel(mdl);
  rsAssert(mdl.isDataValid());
  std::vector<rsInstantaneousSineParams<double>> invalidDataPoints = mdl.getInvalidDataPoints();
  //plotSineModel(mdl, fs);



  // these partials with the Rhodes_F3_Medium sample...
  plotSineModelAmplitudes(mdl, {1,2,3,4,5,6});

  plotSineModelAmplitudes(mdl, {1}); // ...shows problem with detecting min-peaks
  plotSineModelAmplitudes(mdl, {2}); // ...shows problem with straight lines (to small peak-densisty)
  //plotSineModelAmplitudes(mdl, {3}); // shows problem
  plotSineModelAmplitudes(mdl, {4}); // ...shows problem with drop in/out



  y = synthesizeSinusoidal(mdl, fs);
  rosic::writeToMonoWaveFile(name + "DeBeatOutput.wav", &y[0], (int)y.size(), (int)fs);

  // todo: the peak-picking algo in the amplitude de-beating is still not working as it should
  // -sometimes, there is a too great distance between two peaks (we need more density)
  //  -example: 2nd partial of Rhodes_F3_Medium
  // -sometimes, it also detects small mini-peaks within valleys which should be ignored
  //  -example: fundamental of Rhodes_F3_Medium
}

void testEnvelopeMatching(std::vector<double>& x1, std::vector<double>& x2)
{
  // todo: 
  // (1) extract envelopes of both signals
  // (2) find best match offset

  //rsEnvelopeExtractor<double> ee;

  //rsEnvelopeFollower2<double> ef;
  //ef.setSampleRate(44100);  // make this a function parameter
  //ef.setInputFrequency(85); // this too

  rsEnvelopeFollower<double, double> ef;
  ef.setSampleRate(44100);  // make this a function parameter
  ef.setAttackTime(0.0);    // in ms?
  ef.setReleaseTime(200.0);
  //ef.setInputFrequency(85); // this too



  // todo: use simpler envelope follower with faster attack - and then use a smaller initial ignore
  // section for the shiftee (the need for an initial ignore section is an artifact of the slow 
  // attack of the env follower)

  // exctract envelopes:
  std::vector<double> e1(x1.size()), e2(x2.size());
  int n;
  for(n = 0; n < x1.size(); n++) e1[n] = ef.getSample(x1[n]);
  ef.reset();
  for(n = 0; n < x2.size(); n++) e2[n] = ef.getSample(x2[n]);

  //rsPlotVectors(x2, e2);

  double thresh = -65;

  RAPT::rsExponentialEnvelopeMatcher<double> em;
  em.setMatchLevel(-55);               // make function parameter

  //em.setInitialIgnoreSection1(16000);  // reference signal has 2-stage decay
  em.setInitialIgnoreSection1(60000);  // ..or actually mor like a 3-stage decay
  em.setInitialIgnoreSection2( 1000);

  // good for when the tail is actually cut out from the full signal:
  //em.setInitialIgnoreSection1(60000); 
  //em.setInitialIgnoreSection2(    0);

  em.setIgnoreThreshold1(thresh);
  em.setIgnoreThreshold2(thresh);

  int dt = (int) em.getMatchOffset(&e1[0], (int) e1.size(), &e2[0], (int) e2.size());
  // hmm...maybe we should pass references to the enve-follower and env-matcher, so the caller can
  // set them up - it would be too many function parameters otherwise




  // decimate enevlopes for plotting (GNUPlot doesn't like big datasets):
  int decimation = 16;
  std::vector<double> e1d = rsDecimate(e1, decimation);
  std::vector<double> e2d = rsDecimate(e2, decimation);

  // convert to decibels:
  std::vector<double> db1d(e1d.size()), db2d(e2d.size());
  for(n = 0; n < e1d.size(); n++)  db1d[n] = rsMax(rsAmp2dB(e1d[n]), thresh);
  for(n = 0; n < e2d.size(); n++)  db2d[n] = rsMax(rsAmp2dB(e2d[n]), thresh);

  // create the two time axes:
  std::vector<double> t1d(e1d.size()), t2d(e2d.size());
  for(n = 0; n < t1d.size(); n++)  t1d[n] = n * decimation;
  for(n = 0; n < t2d.size(); n++)  t2d[n] = n * decimation + dt;

  // plot:
  GNUPlotter plt;
  //plt.addDataArrays((int) t1d.size(), &t1d[0], &e1d[0]);
  //plt.addDataArrays((int) t2d.size(), &t2d[0], &e2d[0]);
  plt.addDataArrays((int) t1d.size(), &t1d[0], &db1d[0]);
  plt.addDataArrays((int) t2d.size(), &t2d[0], &db2d[0]);
  plt.plot();

  // todo: maybe plot the two regression lines as well
  // try to cut an actual section of the tail from the 1st signal and match that to the full
  // length signal
}

void testEnvelopeMatching2(std::vector<double>& x1, std::vector<double>& x2)
{
  typedef std::vector<double> Vec;

  // exctract envelopes:
  rsEnvelopeFollower<double, double> ef;
  ef.setSampleRate(44100);  // make this a function parameter
  ef.setAttackTime(0.0);    // in ms?
  ef.setReleaseTime(200.0);
  Vec e1(x1.size()), e2(x2.size());
  size_t n;
  for(n = 0; n < (int)x1.size(); n++) e1[n] = ef.getSample(x1[n]);
  ef.reset();
  for(n = 0; n < (int)x2.size(); n++) e2[n] = ef.getSample(x2[n]);

  // find match offset:
  double dt = 0;
  int D = 100;  // decimation factor
  dt = rsEnvelopeMatchOffset(&e1[0], (int) e1.size(), &e2[0], (int) e2.size(), D);

  // create the two time axes and decimated enveloeps for plotting (using the same decimation
  // factor as for matching):
  Vec e1d = rsDecimateViaMean(e1, D);
  Vec e2d = rsDecimateViaMean(e2, D);  
  Vec t1d(e1d.size()), t2d(e2d.size());
  for(n = 0; n < (int)t1d.size(); n++)  t1d[n] = double(n*D);
  for(n = 0; n < (int)t2d.size(); n++)  t2d[n] = double(n*D) + dt;

  // plot:
  GNUPlotter plt;
  plt.addDataArrays((int) t1d.size(), &t1d[0], &e1d[0]);
  plt.addDataArrays((int) t2d.size(), &t2d[0], &e2d[0]);
  plt.plot();
}

void testTimeWarping1(const std::vector<double>& x, double fs, double f0)
{
  //rsTimeWarper<double, double> tw;
  //tw.timeWarpSinc(

  //rsPitchFlattener<double, double> pf;


  int dummy = 0;
}



