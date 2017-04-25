#include "rosic_MultiSegmentWaveformRenderer.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MultiSegmentWaveformRenderer::MultiSegmentWaveformRenderer()
{

}

//-------------------------------------------------------------------------------------------------
// waveform rendering:

void MultiSegmentWaveformRenderer::renderWaveform(double *targetBuffer, int length)
{
  breakpointModulator.fillBufferWithEnvelope(targetBuffer, length, 
    breakpointModulator.getLoopMode() != BreakpointModulator::NO_LOOP);
}






