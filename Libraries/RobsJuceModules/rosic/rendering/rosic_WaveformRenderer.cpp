//#include "rosic_WaveformRenderer.h"
//using namespace rosic;

WaveformRenderer::WaveformRenderer()
{
  mode           = STANDARD_WAVEFORM;
  //mode           = ALGORITHM;  // test
  dcRemove       = false;
  normalize      = false; 
  fitToUnitRange = false;
}

WaveformRenderer::~WaveformRenderer()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void WaveformRenderer::renderWaveForm(double *targetBuffer, int length)
{
  // dispatch the rendering request to the appropriate renderer:
  switch( mode )
  {
  case STANDARD_WAVEFORM: standardRenderer.renderWaveform(    targetBuffer, length); break;
  case AUDIO_FILE:        waveBuffer.getWaveform(             targetBuffer, length); break;
  case ALGORITHM:         algorithmicRenderer.renderWaveform( targetBuffer, length); break;
  case MULTI_SEGMENT:     multiSegmentRenderer.renderWaveform(targetBuffer, length); break;
  default:                RAPT::rsArray::fillWithZeros(                      targetBuffer, length);
  }

  // some optional post processing steps:
  if( fitToUnitRange == true )
    RAPT::rsArray::fitIntoRange(targetBuffer, length, -1.0, 1.0);
  else
  {
    if( dcRemove == true )
      RAPT::rsArray::removeMean(targetBuffer, length);
    if( normalize == true )
      RAPT::rsArray::normalize(targetBuffer, length, 1.0);
  }
}
    
