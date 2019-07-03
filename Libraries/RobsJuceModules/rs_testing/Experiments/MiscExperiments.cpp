// todo: SinusoidalSynthesizer to RAPT to make it available here

using namespace RAPT;


// convenience function:
std::vector<double> synthesizeSinusoidal(
  const RAPT::rsSinusoidalModel<double>& model, double sampleRate, double fadeTime)
{
  typedef RAPT::rsSinusoidalSynthesizer<double> SS;
  typedef SS::PhaseInterpolationMethod PIM;
  SS synth;
  synth.setSampleRate(sampleRate);
  //synth.setCubicAmplitudeInterpolation(true);
  synth.setPhaseInterpolation(PIM::tweakedFreqIntegral);
  std::vector<double> x = synth.synthesize(model);
  if(fadeTime > 0.0)
    applyFadeInAndOut( &x[0], (int) x.size(), int (fadeTime*sampleRate));
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
  analyzer.setSpectralOversampling(8);  // zero padding
  analyzer.setAllowInharmonics(true);
  analyzer.setSpectralPeakSearchWidth(0.5);       // default: 1 - blackman needs a value less than 1
  analyzer.setMinPeakToMainlobeWidthRatio(0.75);  // default: 0.75

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


  // Manipulations:

  mdl.removePartial(0);                      // remove DC
  mdl.removePartialsWithMeanFreqAbove(fs/2); // anti-alias

  //mdl.keepOnly({0, 9});  // for test with TwoSines_Freq1=200_Freq2=2025
  //mdl.removePartial(0);    // test: resynthesize without fundamental

  //mdl.keepOnly({0, 99});  // fom TwoSines_Freq1=100_Freq2=10020

  //mdl.removePartialsAbove(9); 
  // for testing artifacts (spurious high partials) in flute-c1 - are they analysis or synthesis 
  // artifacts? if they occur evene though we remove the partials above 9, they occur in synthesis
  // ->they seem to be analysis artifacts


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
    RAPT::rsArray::fillWithZeros(&zeros[0], (int) marks.size());
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
  rsArray::normalize(&output[0], (int) output.size(), 1.0);


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

  size_t n;
  if(nf < 0) {
    // obtain a version of x with an appropriate number of zeros prepended
    nf = -nf;
    x.resize(Nx+nf);
    for(n = 0; n < nf; n++)        x[n] = 0;
    for(n = nf; n < x.size(); n++) x[n] = xIn[n-nf];
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

template<class T>
rsModalFilterParameters<T> getModalModel(const RAPT::rsSinusoidalPartial<T>& partial)
{
  // todo: estimate the actual parameters:
  // frequency: take average freq of partial
  // amplitude: take peak/max amplitude fo partial
  // attack:    time-instant of peak
  // decay:     compare average amplitudes of 1st and 2nd half after the peak
  // phase:     start phase in model

  rsModalFilterParameters<T> params;
  params.freq  = partial.getMeanFreq();  
  // maybe we should cut off the transient before taking the mean?

  params.phase = partial.getFirstDataPoint().getWrappedPhase();
  // maybe this is not the best strategy either - maybe take the phase at the amplitude peak and
  // compute the start-phase from that and the mean frequency ...it seems, it would be a good idea
  // to implement a class where the

  int peakIndex = partial.getMaxAmpIndex();
  params.amp = partial.getDataPoint(peakIndex).getAmplitude();
  params.att = partial.getDataPoint(peakIndex).getTime();

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
  for(int i = peakIndex+1; i < (int) partial.getNumDataPoints(); i++)
  {
    if(partial.getDataPoint(i).gain < targetAmp)
    {
      params.dec = partial.getDataPoint(i).time - params.att;
      // this is very coarse - todo: interpolate (linearly on the dB-scale)

      break;
    }
  }
  // maybe sometimes, the sample isn't long enough to contain the data, where it has decayed to
  // peak/e (because of an early fade-out or palm-muting, whatever) - then we should use 
  // c*peak/e for some c < 1 and multiply the decay-time by that same c - maybe have a loop
  // with exponentially decreasing c, i.e. c = 1,0.5,0.25,0.125,...
  // what about samples with no decay at all, i.e. sustained sounds?





  // maybe instead of averaging, take the maximum value of the section that starts halfway
  // after the peak - maybe the function getMaxAmpIndex can take a start-search index as parameter




  return params;
  //return rsModalFilterParameters<T>(); // preliminary
}

template<class T>
std::vector<rsModalFilterParameters<T>> getModalModel(const RAPT::rsSinusoidalModel<T>& model)
{
  std::vector<rsModalFilterParameters<T>> p(model.getNumPartials());
  for(int i = 0; i < model.getNumPartials(); i++)
    p[i] = getModalModel(model.getPartial(i));
  return p;
}
template std::vector<rsModalFilterParameters<double>> 
  getModalModel(const RAPT::rsSinusoidalModel<double>& model);

std::vector<double> synthesizeModal(
  const std::vector<rsModalFilterParameters<double>>& p, double fs, int N)
{
  std::vector<double> x(N);
  rsArray::fillWithZeros(&x[0], N);
  rosic::rsModalFilterWithAttackDD flt;
  for(size_t i = 0; i < p.size(); i++) {
    flt.setModalParameters(p[i].freq, p[i].amp, p[i].att, p[i].dec, p[i].phase, fs);
    flt.reset();
    x[0] += flt.getSample(1);
    for(int n = 1; n < N; n++)
      x[n] += flt.getSample(0); }
  return x;
}