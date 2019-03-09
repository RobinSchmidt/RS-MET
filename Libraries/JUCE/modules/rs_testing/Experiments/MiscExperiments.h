#pragma once

// convenience function:
std::vector<double> synthesizeSinusoidal(
  const RAPT::rsSinusoidalModel<double>& model, double sampleRate, double fadeTime = 0.0);


void testHarmonicResynthesis(const std::string& name, std::vector<double>& input, 
  double fs, bool writeWaveFiles, bool plotResults);

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


// move to RAPT::rsArray
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
