// todo: SinusoidalSynthesizer to RAPT to make it available here

/*
template<class T>
void applyFadeIn(T* x, int N, int numFadeSamples)
{
  int nf = rsMin(numFadeSamples, N);
  for(int n = 0; n < nf; n++) {
    T t = T(n) / T(nf);
    x[n] *= t;
  }
}
template<class T>
void applyFadeOut(T* x, int N, int numFadeSamples)
{
  int nf = rsMin(numFadeSamples, N);
  for(int n = 0; n < nf; n++) {
    T t = T(n) / T(nf);
    x[N-n-1] *= t;
  }
}
template<class T>
void applyFadeInAndOut(T* x, int N, int numFadeSamples)
{
  applyFadeIn( &x[0], N, numFadeSamples);
  applyFadeOut(&x[0], N, numFadeSamples);
}

// convenience function:
std::vector<double> synthesizeSinusoidal(
  const RAPT::rsSinusoidalModel<double>& model, double sampleRate, double fadeTime = 0.0)
{
  RAPT::SinusoidalSynthesizer<double> synth;
  synth.setSampleRate(sampleRate);
  //synth.setCubicAmplitudeInterpolation(true);
  std::vector<double> x = synth.synthesize(model);
  if(fadeTime > 0.0)
    applyFadeInAndOut( &x[0], (int) x.size(), int (fadeTime*sampleRate));
  return x;
}

void testHarmonicResynthesis(const std::string& name, std::vector<double>& input, 
  double fs, bool writeWaveFiles, bool plotResults)
{
  // analyze, resynthesize and create error signal:
  double* x = &input[0];   // pointer to first sample (for convenience)
  int Nx = (int) input.size();
  RAPT::rsHarmonicAnalyzer<double> analyzer;
  analyzer.setSampleRate(fs);
  analyzer.setSincInterpolationLength(512);
  RAPT::rsSinusoidalModel<double> mdl = analyzer.analyze(x, Nx);
  //plotSineModel(mdl, fs);
  std::vector<double> output = synthesizeSinusoidal(mdl, fs); 
  std::vector<double> error = output-input;
  double* y = &output[0]; int Ny = (int) output.size(); // again, for convenience
  double* e = &error[0];  int Ne = (int) error.size();  // dito

  // write original, resynthesized and error signals to files, if desired:
  if(writeWaveFiles == true) {
    rosic::writeToMonoWaveFile((name + "Original.wav").c_str(),      x, Nx, (int)fs);
    rosic::writeToMonoWaveFile((name + "Resynthesized.wav").c_str(), y, Ny, (int)fs);
    rosic::writeToMonoWaveFile((name + "Error.wav").c_str(),         e, Ne, (int)fs);
  }

  // plot original, resynthesized and error signals, if desired:
  if(plotResults == true) {
    GNUPlotter plt;
    plt.addDataArrays(Nx, x);
    plt.addDataArrays(Ny, y);
    plt.addDataArrays(Ne, e);
    //plt.addDataArrays(Ne-2000, &e[1000]);  // middle part of error
    plt.plot();
  }
}
*/