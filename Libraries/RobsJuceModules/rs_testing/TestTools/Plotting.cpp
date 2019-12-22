#ifndef RAPT_PLOTTING_CPP // why? i think, it may be deleted
#define RAPT_PLOTTING_CPP

#include "Plotting.h"

//#include "rapt/rapt.h"

using namespace RAPT;

void createTimeAxis(int numSamples, float *timeAxis, float sampleRate)
{
  for(int n = 0; n < numSamples; n++)
    timeAxis[n] = n / sampleRate;
}

void createTimeAxis(int numSamples, double *timeAxis, double sampleRate)
{
  for(int n = 0; n < numSamples; n++)
    timeAxis[n] = n / sampleRate;
}

//void plotData(int N, float *x, float *y1, float *y2, float *y3, float *y4, float *y5)
//{
//  GNUPlotter plt;
//  plt.addDataArrays(N, x, y1, y2, y3, y4, y5);
//  plt.plot();
//}

/* // commented because inlined in .h now
template<class TSig, class TFlt>
vector<TSig> impulseResponse(TFlt &filter, int N, TSig scale)
{
  vector<TSig> y(N);
  filter.reset();
  y[0] = filter.getSample(scale);
  for(int n = 1; n < N; n++)
    y[n] = filter.getSample(0.0);
  return y;
}

template<class TSig, class TFlt>
void plotImpulseResponse(TFlt &filter, int N, TSig scale)
{
  vector<TSig> y = impulseResponse(filter, N, scale);
  GNUPlotter plt;
  plt.addDataArrays(N, &y[0]);
  plt.plot();
}
*/

// move this code to FilterPlotter (or maybe there already is something similar?) - allows user
// to create an object and access further settings
void plotFrequencyResponse(std::vector<double>& f, std::vector<double>& dB,
  std::vector<double>& degrees, bool logFreq)
{
  //rsAssert(size(f) == size(dB) && size(f) == size(degrees));
  int N = (int) f.size();
  GNUPlotter p;
  p.addDataArrays(N, &f[0], &dB[0]);
  p.addDataArrays(N, &f[0], &degrees[0]);
  p.setPixelSize(1200, 400);
  p.setTitle("Filter Frequency Response");
  //p.setGraphColors("A00000", "909000", "008000", "0000A0", "800080",
  //  "A00000", "909000", "008000", "0000A0", "800080" );
  if(logFreq)
    p.addCommand("set logscale x");
  //p.addCommand("set xrange  [0.0625:16]");
  //p.addCommand("set yrange  [-100:0]");
  //p.addCommand("set y2range [-450:0]");
  p.addCommand("set xlabel \"Frequency in Hz\"");
  p.addCommand("set ylabel \"Magnitude in dB\"");
  p.addCommand("set y2label \"Phase in Degrees\"");
  //p.addCommand("set xtics 2");    // factor 2 between (major) frequency axis tics
  //p.addCommand("unset mxtics");   // no minor tics for frequency axis
  p.addCommand("set ytics 10");   // 10 dB steps for magnitude axis
  p.addCommand("set y2tics 45");  // 45° steps for phase axis

  // add magnitude and phase graphs:
  p.addGraph("i 0 u 1:2 w lines lw 1.5 axes x1y1 notitle");
  p.addGraph("i 1 u 1:2 w lines lw 1.5 axes x1y2 notitle");
  p.plot();
}



template<class T>
void plotArrays(int N, T *y1, T *y2, T *y3, T *y4, T *y5, T *y6, T *y7, T *y8, T *y9)
{
  GNUPlotter plt;
  plt.plotArrays(N, y1, y2, y3, y4, y5, y6, y7, y8, y9);
}
template void  plotArrays(int N, double *y1, double *y2, double *y3, double *y4, double *y5, 
  double *y6, double *y7, double *y8, double *y9);

// new from RSLib test - needs testing:

void plotData(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  GNUPlotter plt;
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
}

void plotData(int N, double x0, double dx, double *y1, double *y2, double *y3, double *y4,
  double *y5)
{
  double *x = new double[N];
  for(int n = 0; n < N; n++)
    x[n] = x0 + n*dx;
  GNUPlotter plt;
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
  delete[] x;
}

void plotDataLogX(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  GNUPlotter plt;
  plt.setLogScale("x");
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
}

// replace by RAPT::rsPlotVector
void plotVector(std::vector<double> v)
{
  int N = (int) v.size();
  double *y = new double[N];
  for(int n = 0; n < N; n++)
    y[n] = v[n];
  GNUPlotter plt;
  plt.plotArrays(N, y);
  delete[] y;
}

void plotComplexVectorReIm(std::vector<std::complex<double>> v)
{
  GNUPlotter plt;
  int N = (int) v.size();
  double *y = new double[N];
  for(int n = 0; n < N; n++) 
    y[n] = real(v[n]);
  plt.addDataArrays(N, y);
  for(int n = 0; n < N; n++) 
    y[n] = imag(v[n]);
  plt.addDataArrays(N, y);
  delete[] y;
  plt.plot();
}


// use shorter variable names here...
void plotSpectrogram(int numFrames, int numBins, double **s, double fs, int H, 
  double dbMin, double dbMax)
{
  // todo:
  // code copied to class SpectrogramPlotter - call it here and delete code below

  // fs: sample rate, H: hop size

  GNUPlotter p;

  // create time- and frequency axis and add the data:
  double tMax = H*(numFrames-1)/fs;
  double fMax =  0.5*fs*(numBins-1)/numBins;
  double *t = new double[numFrames];
  double *f = new double[numBins];
  p.rangeLinear(t, numFrames, 0.0, tMax);
  p.rangeLinear(f, numBins,   0.0, fMax);


  p.setDataPrecision(4);                                  // make temp files smaller
  p.addDataMatrix(numFrames, numBins, t, f, s);

  // set up style settings:
  //p.setAxisLabels("Time in s", "Frequency in Hz", "Level in dB");
  p.setPixelSize(800, 400);
  p.addGraph("i 0 nonuniform matrix w image notitle");

  //p.addCommand("set palette color");                  // this is used by default
  //p.addCommand("set palette color negative");         // reversed colors
  //p.addCommand("set palette gray negative");          // maximum is black
  //p.addCommand("set palette gray");                   // maximum is white
  p.addCommand("set palette rgbformulae 30,31,32");     // colors printable as grayscale


  // maybe use signal-dependent values later:
  //double dbMin = -100.0;
  //double dbMax =   10.0;
  p.setRange(0.0, tMax, 0.0, fMax, dbMin, dbMax);
  // doesn't seem to have any effect on the color axis (z). it seems like GNUPlot uses
  // autoscaling for it, no matter what.

  p.plot();
  delete[] t;
  delete[] f;
}
// maybe factor out a matrix plotting function

void plotSpectrogram(int numFrames, int numBins, const rsMatrix<std::complex<double>>& spec, 
  double sampleRate, int hopSize, double dbMin, double dbMax)
{
  double **dB;
  MatrixTools::rsAllocateMatrix(dB, numFrames, numBins);
  for(int i = 0; i < numFrames; i++)
    for(int  j = 0; j < numBins; j++)
      dB[i][j]  = rsMax( rsAmp2dB(abs(spec(i,j))) , dbMin);
  plotSpectrogram(numFrames, numBins, dB, sampleRate, hopSize, dbMin, dbMax);
  MatrixTools::rsDeAllocateMatrix(dB, numFrames, numBins);
}


void plotPhasogram(int numFrames, int numBins, double **phases, double sampleRate, int hopSize)
{
  plotSpectrogram(numFrames, numBins, phases, sampleRate, hopSize, -PI, PI);
  // a bit dirty to abuse the spectrogram...clean up..
}

void plotMagnitudeResponse(const RAPT::rsFilterSpecificationBA<double>& specBA)
{
  FilterPlotter<double> plt;
  plt.addFilterSpecificationBA(specBA);
  plt.addCommand("set xtics ('0' 0, 'sr/4' pi/2, 'sr/2' pi)");
  plt.plotMagnitude(1000, 0.0, PI, false, false);
}

void plotPolesAndZeros(const RAPT::rsFilterSpecificationBA<double>& specBA)
{
  FilterPlotter<double> plt;
  plt.addFilterSpecificationBA(specBA);
  plt.plotPolesAndZeros();
}

void showFilterPlots(const RAPT::rsFilterSpecificationBA<double>& specBA)
{
  plotMagnitudeResponse(specBA);
  plotPolesAndZeros(specBA);
}

void stemPlot(int N, double *x, double *y)
{
  GNUPlotter plt;
  plt.addDataArrays(N, x, y);
  plt.addDataArrays(N, x, y); // can probably be done without adding the data twice
  plt.setGraphStyles("impulses", "points pt 7 ps 1.2");
  plt.plot();
}

void plotSinusoidalAnalysisResult(RAPT::rsSinusoidalAnalyzer<double>& sa, double* x, int N, double fs)
{
  SinusoidalModelPlotter<double> plt;
  plt.plotAnalysisResult(sa, x, N, fs);
}

void plotSineModel(const RAPT::rsSinusoidalModel<double>& model, double fs)
{
  SinusoidalModelPlotter<double> plt;
  plt.plotModel(model, fs);
}

void plotTwoSineModels(
  const RAPT::rsSinusoidalModel<double>& model1, 
  const RAPT::rsSinusoidalModel<double>& model2, 
  double fs)
{
  SinusoidalModelPlotter<double> plt;
  plt.plotTwoModels(model1, model2, fs);
}

std::vector<int> rsIndexFilledVector(int N)
{
  std::vector<int> v(N);
  for(int i = 0; i < N; i++)
    v[i] = i;
  return v;
}


void plotSineModelAmplitudes(
  const RAPT::rsSinusoidalModel<double>& model, std::vector<int> partialIndices)
{
  GNUPlotter plt; 
  std::vector<double> t, a;  // time and amplitude arrays (x- and y-values for the plot)

  if(partialIndices.empty()) 
    partialIndices = rsIndexFilledVector((int)model.getNumPartials());

  for(int i = 0; i < (int)partialIndices.size(); i++) {
    int index = partialIndices[i];
    t = model.getPartial(index).getTimeArray();
    a = model.getPartial(index).getAmplitudeArray();
    plt.addDataArrays((int)t.size(), &t[0], &a[0]);
  }

  plt.setPixelSize(1200, 300);
  plt.plot();
}

void plotSineModelPhases(
  const RAPT::rsSinusoidalModel<double>& model,
  const std::vector<int>& partialIndices, bool derivative)
{
  GNUPlotter plt;
  std::vector<double> t, p, f, pd;

  for(int i = 0; i < (int)partialIndices.size(); i++) {
    int index = partialIndices[i];
    t = model.getPartial(index).getTimeArray();
    p = model.getPartial(index).getPhaseArray();
    f = model.getPartial(index).getFrequencyArray();
    p = rsSinusoidalProcessor<double>::unwrapPhase(t, f, p);          // unwrap
    if(derivative) {
      pd.resize(p.size());
      rsNumericDerivative(&t[0], &p[0], &pd[0], (int)p.size(), true); // take derivative
      plt.addDataArrays((int)t.size(), &t[0], &pd[0]); }
    else {
      rsDeTrender<double> dtr;
      dtr.removeTrendAndOffset((int)p.size(), &t[0], &p[0], &p[0]);
      plt.addDataArrays((int)t.size(), &t[0], &p[0]); 
    }
  }
  plt.setPixelSize(1600, 400);
  plt.plot();
}

void plotSineResynthesisResult(const RAPT::rsSinusoidalModel<double>& model, 
  const RAPT::rsSinusoidalSynthesizer<double>& synth, double* x, int Nx)
{
  // Synthesize model output and create residual:
  typedef std::vector<double> Vec;
  Vec xp, yp;        // xp: padded x, yp: padded model output
  getPaddedSignals(x, Nx, model, synth, xp, yp);
  Vec r = xp - yp;   // residual

  // write results to wavefile, if desired:
  int N = (int)yp.size();
  bool writeToWaveFile = false;  // make user parameter
  if(writeToWaveFile == true) {
    int fs = (int) synth.getSampleRate();
    rosic::writeToMonoWaveFile("SineModeling_Original.wav",      &xp[0], N, fs, 16);
    rosic::writeToMonoWaveFile("SineModeling_Resynthesized.wav", &yp[0], N, fs, 16);
    rosic::writeToMonoWaveFile("SineModeling_Residual.wav",      &r[0],  N, fs, 16);
  }

  // plot original, resynthesized and residual signals:
  Vec t(N);
  createTimeAxis(N, &t[0], synth.getSampleRate());
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0],  &xp[0]);
  plt.addDataArrays(N, &t[0], &yp[0]);
  plt.addDataArrays(N, &t[0], &r[0]);
  plt.setPixelSize(1200, 400);
  plt.plot();
}

void plotModelOutputComparison(
  const RAPT::rsSinusoidalModel<double>& model1,
  const RAPT::rsSinusoidalModel<double>& model2,
  const RAPT::rsSinusoidalSynthesizer<double>& synth)
{
  std::vector<double> x = synth.synthesize(model1); // reference signal
  plotSineResynthesisResult(model2, synth, &x[0], (int)x.size());
}

void plotModalAmplitudes(const std::vector<rsModalFilterParameters<double>>& modeModel)
{
  int N = (int)modeModel.size();
  std::vector<double> f(N), a(N), d(N);
  for(int n = 0; n < N; n++) {
    f[n] = modeModel[n].freq;
    a[n] = modeModel[n].amp;
    d[n] = modeModel[n].dec;
  }
  GNUPlotter plt;
  plt.addDataArrays(N, &f[0], &a[0]);
  //plt.addDataArrays(N, &f[0], &a[0], &d[0]);
  plt.plot();
  // todo: plot with impulses, maybe plot attack and decay times
}


void plotModeVsSineAmpEnv(
  rsModalFilterParameters<double>& m, RAPT::rsSinusoidalPartial<double>& s)
{
  std::vector<double> ts, as;
  ts = s.getTimeArray();
  as = s.getAmplitudeArray();
  int N = (int) ts.size();


  double fs = N / (ts[N-1] - ts[0]); // average sample rate
  std::vector<double> tm(N), am(N);
  rsArrayTools::fillWithRangeLinear(&tm[0], N, 0.0, ts[N-1]);

  // factor out (synthesizeModalPartial or something)
  RAPT::rsModalFilterWithAttack<double, double> flt;
  flt.setModalParameters(0.0, m.amp, m.att, m.dec, 90, fs);
  am[0] = flt.getSample(1.0);
  for(int n = 1; n < N; n++)
    am[n] = flt.getSample(0.0);

  GNUPlotter plt;
  plt.addDataArrays(N, &ts[0], &as[0]);
  plt.addDataArrays(N, &tm[0], &am[0]);
  plt.plot();
}


#endif
