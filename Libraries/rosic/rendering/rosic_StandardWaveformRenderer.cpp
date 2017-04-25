#include "rosic_StandardWaveformRenderer.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

StandardWaveformRenderer::StandardWaveformRenderer()
{
  waveform = SINE;
}

//-------------------------------------------------------------------------------------------------
// waveform rendering:

void StandardWaveformRenderer::renderWaveform(double *targetBuffer, int length)
{
  // dispatch the rendering request to the appropriate renderer:
  switch( waveform )
  {
  case SINE:     renderSineWaveform(    targetBuffer, length);  break; 
  case SAW:      renderSawWaveform(     targetBuffer, length);  break; 
  case SQUARE:   renderSquareWaveform(  targetBuffer, length);  break; 
  case TRIANGLE: renderTriangleWaveform(targetBuffer, length);  break; 
  default:       fillWithZeros(         targetBuffer, length);
  }
}

void StandardWaveformRenderer::renderSineWaveform(double *buffer, int length)
{
  for(int i=0; i<length; i++)
    buffer[i] = sin( (2.0*PI*i) / (double) (length) );
}

void StandardWaveformRenderer::renderSawWaveform(double *buffer, int length)
{
  int    N  = length;
  double k  = 0.5;  // more general: k = symmetry
  int    N1 = clip(roundToInt(k*(N-1)), 1, N-1);
  int    N2 = N-N1;
  double s1 = 1.0 / (N1-1);
  double s2 = 1.0 / N2;
  for(int n=0; n<N1; n++)
    buffer[n] = s1*n;
  for(int n=N1; n<N; n++)
    buffer[n] = -1.0 + s2*(n-N1);
}

void StandardWaveformRenderer::renderSquareWaveform(  double *buffer, int length)
{
  int    N  = length;
  double k  = 0.5;  // more general: k = symmetry
  int    N1 = clip(roundToInt(k*(N-1)), 1, N-1);
  int    N2 = N-N1;
  for(int n=0; n<N1; n++)
    buffer[n] = +1.0;
  for(int n=N1; n<N; n++)
    buffer[n] = -1.0;
}

void StandardWaveformRenderer::renderTriangleWaveform(double *buffer, int length)
{
  int i;
  for (i=0; i<(length/4); i++)
    buffer[i] = (double)(4*i) / (double)(length);
  for (i=(length/4); i<(3*length/4); i++)
    buffer[i] = 2.0 - ((double)(4*i) / (double)(length));
  for (i=(3*length/4); i<(length); i++)
    buffer[i] = -4.0+ ((double)(4*i) / (double)(length));
}





