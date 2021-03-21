#pragma once

// convenience function:
std::vector<double> synthesizeSinusoidal(
  const RAPT::rsSinusoidalModel<double>& model, double sampleRate, double fadeTime = 0.0);

// maybe make a getSineModel convenience function

/** ...
If plotResults is true, the function will plot the original, resynthesized and residual signals 
along with markers at the cycle-marks. The analysis datapoints are always in the middle between two
cycle-marks
todo: maybe also optionally plot the model..maybe use thickness and/or color of the lines to 
indicate the amplitude of the partials */
void testHarmonicResynthesis(const std::string& name, std::vector<double>& input, 
  double sampleRate, double fundamental = 0, bool writeWaveFiles = true, bool plotResults = false);

void testMakeHarmonic(const std::string& name, std::vector<double>& input, 
  double sampleRate, double targetFundamental, double inharmonicity = 0, 
  double originalFundamental = 0);

// maybe make a class rsWaveFileBuffer that has fields "name", "sampleData" "sampleRate" and maybe 
// some metadata - to avoid having to pass so many saparate function parameters.

// maybe make a class rsSinusoidalModelEffect that wraps the analyzer, a chain of processors and 
// synthesizer into a single, convenient object - this class could then accept rsWaveFileBuffer
// objects as input for processing
// maybe we should have separate (base)classes for processing the model data and for processing
// the interpolated instantaneous phase- and amplitude trajectories

/** Convenience function to produce a vector x from the array xIn of length N and a vector y from 
the model in such a way that they are time aligned even in cases when the lengths dont match by 
padding an appropriate number of zero samples where necessarry
x: original input of length N, model: a sinusoidal model 
This is useful for creating time aligned signal vectors from a reference signal x and a model that
is supposed to resynthesize x for comparing original to resynthesized sound. */
void getPaddedSignals(double* xIn, int Nx,  const RAPT::rsSinusoidalModel<double>& model,
  const RAPT::rsSinusoidalSynthesizer<double>& synth, 
  std::vector<double>& x, std::vector<double>& y);
// actually, we may need this later also to obtain a residual...function should be renamed


void testModalResynthesis(const std::string& name, std::vector<double>& input, 
  double sampleRate, double fundamental = 0);

void testDeBeating(const std::string& name, std::vector<double>& input, 
  double sampleRate, double fundamental = 0);


void testEnvelopeMatching(std::vector<double>& input1, std::vector<double>& input2);
void testEnvelopeMatching2(std::vector<double>& input1, std::vector<double>& input2);


void testTimeWarping1(const std::vector<double>& x, double fs, double f0);

// move to RAPT::rsArrayTools
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
