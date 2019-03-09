#pragma once

//void testHarmonicResynthesis(const std::string& name, std::vector<double>& input, 
//  double fs, bool writeWaveFiles, bool plotResults);

/** Convenience function to produce a vector x from the array xIn of length N and a vector y from 
the model in such a way that they are time aligned even in cases when the lengths dont match by 
padding an appropriate number of zero samples where necessarry
x: original input of length N, model: a sinusoidal model 
This is useful for creating time aligned signal vectors from a reference signal x and a model that
is supposed to resynthesize x for comparing original to resynthesized sound. */
void getPaddedSignals(double* xIn, int Nx,  const RAPT::rsSinusoidalModel<double>& model,
  const RAPT::SinusoidalSynthesizer<double>& synth, std::vector<double>& x, std::vector<double>& y);
// actually, we may need this later also to obtain a residual...function should be renamed