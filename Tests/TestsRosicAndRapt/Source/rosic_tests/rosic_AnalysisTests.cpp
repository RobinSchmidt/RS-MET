//#include "rosic_AnalysisTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
//#include "../Shared/Plotting/rosic_Plotter.h"
using namespace rosic;

void rotes::testOscilloscopeBuffer()
{
  // test cases: W = 500
  // p = 450:  N=400:y, N=500:y, N=600:y
  // p = 500:  N=400:y, N=500:y, N=600:y
  // p = 550:  N=400:y, N=500:y, N=600:y


  static const int    L   = 10000;   // total signal length in samples
  static const int    W   = 882;     // display width in pixels
  static const int    n0  = 882;    // sample-instant at which we obtain the display buffer

  static const double p   = 441.0;   // period length in samples
  static const double N   = 882.0;  // number of samples inside time window

  static const double fs  = 44100.0; // sample rate
  static const double f   = fs/p;   // frequency of the test sinusoid
  static const double tau = 0.05;    // damping constant for the test sinusoid

  // create the absolute time-axis and the input signal: 
  int n;
  double t[L];
  RAPT::rsArrayTools::fillWithRangeLinear(t, L, 0.0, (double) (L-1) / fs);
  //fillWithRangeLinear(t, L, 0.0, (double) L);  // test - show time-axis in samples
  double xL[L], xR[L], xM[L], xS[L];
  double w = 2.0*PI*f/fs;
  double d = 1.0/(tau*fs);
  for(n = 0; n < L; n++)
  {
    xL[n] = sin(w*n) * exp(-d*n); 
    xR[n] = cos(w*n) * exp(-d*n); 
    xM[n] = (xL[n]+xR[n]) / sqrt(2.0);
    xS[n] = (xL[n]-xR[n]) / sqrt(2.0);
  }
  //Plotter::plotData(L, t, xL);
  //Plotter::plotData(L, t, xL, xR, xM, xS);


  // create and set up the OscilloscopeBuffer object:
  rosic::SyncedWaveformDisplayBuffer oscBuf; 
  oscBuf.setSampleRate(fs);
  oscBuf.setTimeWindowLength((double) N / fs);
  oscBuf.setDisplayWidth(W);
  oscBuf.setSyncMode(SyncedWaveformDisplayBuffer::ZEROS);
  //oscBuf.setSyncBackward(true);

  // create oscilloscope buffer's relative time axis:
  //double tb[2*W];
  //fillWithRangeLinear(tb, 2*W, 0.0, (double) (2*W-1) / fs);
  //fillWithRangeLinear(tb, 2*W, 0.0, (double) (N-1)); // test - in samples


  // create the display buffer:
  double xb[2*W];
  for(n = 0; n < L; n++)
  {
    oscBuf.feedInputBuffer(&xL[n], 1);
    if( n == n0 )
    {
      // obtain the buffer to be shown:
      //oscBuf.updateDisplayBuffer();
      double *buf = oscBuf.getDisplayBuffer();
      RAPT::rsArrayTools::copy(buf, xb, 2*W);
      int dummy = 0;
    }
  }


  plotData(2*W, oscBuf.getTimeAxis(), xb);
  //plotData(2*W, tb, xb);
}
