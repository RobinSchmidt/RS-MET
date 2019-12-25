#include "romos_GenerateDesiredOutput.h"
using namespace romos;

void GenerateDesiredOutput::forImpulse(int N, double *x)
{
  x[0] = 1.0;
  for(int n = 1; n < N; n++)
    x[n] = 0.0;
}

void GenerateDesiredOutput::forWhiteNoiseUniform(int N, double  *d, unsigned long seed)
{ 
  unsigned long state = seed;
  for(int n = 0; n < N; n++)
  {
    state  = 1664525 * state + 1013904223;   // mod implicitely by integer overflow
    d[n] = -1.0 + ((2.0/4294967296.0) * state);
  }
}

void GenerateDesiredOutput::forUnitDelay(int N, double *x, double *d)
{
  d[0] = 0.0;
  for(int n = 1; n < N; n++)
    d[n] = x[n-1];
}

void GenerateDesiredOutput::forSummedDiffs(int N, double **x, double **d)
{
  double d1, d2, d3, d4;
  for(int n = 0; n < N; n++)
  {
    d1 = x[0][n] - x[1][n];
    d2 = x[1][n] - x[0][n];
    d3 = x[2][n] - x[1][n];
    d4 = x[1][n] - x[2][n];
    d[0][n] = d1 + d2 + d4;
    d[1][n] = d2 + d4;
    d[2][n] = d1 + d3;
    d[3][n] = d1 + d3 + d4;
  }
}

void GenerateDesiredOutput::forMovingAverage(int N, double *x, double *b0, double *b1, double *d)
{
  d[0] = b0[0]*x[0] + b1[0]*0.0;  // assumes 0.0 as initial state of x[n-1] buffer
  for(int n = 1; n < N; n++)
    d[n] = b0[n]*x[n] + b1[n]*x[n-1]; 
}

void GenerateDesiredOutput::forLeakyIntegrator(int N, double *x, double *c, double *d)
{
  d[0] = c[0]*x[0] + (1.0-c[0])*0.0;  // assumes 0.0 as initial state of y[n-1] buffer
  for(int n = 1; n < N; n++)
    d[n] = c[n]*x[n] + (1.0-c[n])*d[n-1]; 
}

void GenerateDesiredOutput::forLeakyIntegratorDoubleDelay(int N, double *x, double *c, double *d)
{
  d[0] = c[0] * x[0] + (1.0 - c[0]) * 0.0;  // assumes 0.0 as initial state of y[n-1] buffer
  d[1] = c[1] * x[1] + (1.0 - c[1]) * 0.0;  // assumes 0.0 as initial state of y[n-2] buffer
  for(int n = 2; n < N; n++)
    d[n] = c[n]*x[n] + (1.0-c[n])*d[n-2]; 
}

void GenerateDesiredOutput::forTestFilter1(int N, double *x, double *b0, double *b1, double *c, 
                                            double *dSum, double *dDiff, double *dProd)
{  
  double *dMovAv   = new double[N];
  double *dLeakInt = new double[N];
  GenerateDesiredOutput::forMovingAverage(  N, x, b0, b1, dMovAv);
  GenerateDesiredOutput::forLeakyIntegrator(N, x, c,      dLeakInt);
  RAPT::rsArrayTools::add(     dMovAv, dLeakInt, dSum,  N);
  RAPT::rsArrayTools::subtract(dMovAv, dLeakInt, dDiff, N);
  RAPT::rsArrayTools::multiply(dMovAv, dLeakInt, dProd, N);
  delete[] dLeakInt;
  delete[] dMovAv;
}

void GenerateDesiredOutput::forBiquad(int N, double *x, double *b0, double *b1, double *b2, double *a1, double *a2, double *y)
{
  /*
  // assumes 0.0 as initial values of all buffers:
  y[0] = b0[0]*x[0];  
  y[1] = b0[1]*x[1] + b1[1]*x[0] - a1[1]*y[0];
  for(int n = 2; n < N; n++)
    y[n] = b0[n]*x[n] + b1[n]*x[n-1] + b2[n]*x[n-2] - a1[n]*y[n-1] - a2[n]*y[n-2];
    */

  double x1   = 0.0;
  double x2   = 0.0;
  double y0   = 0.0;
  double y1   = 0.0;
  double y2   = 0.0;
  for(int n = 0; n < N; n++) {
    y0   = b0[n] * x[n] + b1[n] * x1 + b2[n] * x2 - a1[n] * y1 - a2[n] * y2;
    x2   = x1;
    x1   = x[n];
    y2   = y1;
    y1   = y0;
    y[n] = y0;
    //int dummy = 0;
  }
}


void GenerateDesiredOutput::forBiquadWithFixedCoeffs(int N, double  *x, double  b0, double  b1, 
  double  b2, double  a1, double  a2, double *d)
{
  double *b0a = new double[N];  RAPT::rsArrayTools::fillWithValue(b0a, N, b0);
  double *b1a = new double[N];  RAPT::rsArrayTools::fillWithValue(b1a, N, b1);
  double *b2a = new double[N];  RAPT::rsArrayTools::fillWithValue(b2a, N, b2);
  double *a1a = new double[N];  RAPT::rsArrayTools::fillWithValue(a1a, N, a1);
  double *a2a = new double[N];  RAPT::rsArrayTools::fillWithValue(a2a, N, a2);

  GenerateDesiredOutput::forBiquad(N, x, b0a, b1a, b2a, a1a, a2a, d); 

  delete[] b0a;
  delete[] b1a;
  delete[] b2a;
  delete[] a1a;
  delete[] a2a;
}

void GenerateDesiredOutput::forFormula1In1Out(int N, double *x, double *d)
{
  for(int n = 0; n < N; n++) 
    d[n] = tanh(2 * x[n]*x[n]); // tanh(2*x^2) is our example formula
}


void GenerateDesiredOutput::forFilterBlip(int N, double frequency, double q, double *desiredOutput)
{
  double *x = new double[N];  
  double coeffs[5];  
  romos::biquadBandpassConstSkirtCoeffs(coeffs, frequency, q);  
  GenerateDesiredOutput::forImpulse(              N, x);
  GenerateDesiredOutput::forBiquadWithFixedCoeffs(N, x, coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], desiredOutput);
  delete[] x;
}

/*
void GenerateDesiredOutput::forGatedNoteFrequencies(int N, std::vector<NoteEvent> *events, double ***desiredOutputs,                                           
                                                     bool containerIsPolyphonic, bool noteFreqModuleIsPolyphonic)
{
  std::vector<NoteEvent> noteOns;
  std::vector<int>       durations;
  convertNoteEventsToStartsAndDurations(*events, noteOns, durations);


  zeroDesiredOutputs();

  // in the generation of the desired output signals, we assume here that the note-ons do not exceed the number of available voices, i.e. 
  // we don't consider voice-stealing and we also use a new voice for each of the notes, regardless whether the notes are actually 
  // simultaneous or not
  // \todo use the romos::ProcessingStatus object to do actual voice-allocation as the synth would
  // write separate tests for the ProcessingStatus object


  if( containerIsPolyphonic && noteFreqModuleIsPolyphonic )
  {
    for(unsigned int eventIndex = 0; eventIndex < noteOns.size(); eventIndex++)
    {
      NoteEvent noteOn        = noteOns.at(eventIndex);
      double    noteFrequency = pitchToFreq(noteOn.key);
      int       noteStart     = noteOn.deltaFrames;
      int       noteEnd       = rmin(noteStart + durations[eventIndex], N);
      int       voiceToUse    = eventIndex;
      for(int frameIndex = noteStart; frameIndex < noteEnd; frameIndex++)
        desiredOutputs[voiceToUse][0][frameIndex] += noteFrequency;
    }
  }
  else if( containerIsPolyphonic && !noteFreqModuleIsPolyphonic )
  {
    for(unsigned int eventIndex = 0; eventIndex < noteOns.size(); eventIndex++)
    {
      NoteEvent noteOn        = noteOns.at(0);
      double    noteFrequency = pitchToFreq(noteOn.key);
      int       noteStart     = noteOn.deltaFrames;
      int       noteEnd       = rmin(noteStart + durations[0], N);
      int       voiceToUse    = eventIndex;
      for(int frameIndex = noteStart; frameIndex < noteEnd; frameIndex++)
        desiredOutputs[voiceToUse][0][frameIndex] += noteFrequency;
    }
  }
  else if( !containerIsPolyphonic && noteFreqModuleIsPolyphonic )
  {
    for(unsigned int eventIndex = 0; eventIndex < noteOns.size(); eventIndex++)
    {
      NoteEvent noteOn        = noteOns.at(eventIndex);
      double    noteFrequency = pitchToFreq(noteOn.key);
      int       noteStart     = noteOn.deltaFrames;
      int       noteEnd       = rmin(noteStart + durations[eventIndex], N);
      for(int frameIndex = noteStart; frameIndex < noteEnd; frameIndex++)
        desiredOutputs[0][0][frameIndex] += noteFrequency;
    }
  }
  else if( !containerIsPolyphonic && !noteFreqModuleIsPolyphonic )
  {
    for(unsigned int eventIndex = 0; eventIndex < 1; eventIndex++)
    {
      NoteEvent noteOn        = noteOns.at(0);
      double    noteFrequency = pitchToFreq(noteOn.key);
      int       noteStart     = noteOn.deltaFrames;
      int       noteEnd       = rmin(noteStart + durations[0], N);
      for(int frameIndex = noteStart; frameIndex < noteEnd; frameIndex++)
        desiredOutputs[0][0][frameIndex] += noteFrequency;
    }
  }
}
*/



