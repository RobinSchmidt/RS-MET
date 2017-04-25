#include "rosic_AlgorithmicWaveformRenderer.h"
using namespace rosic;

UnaryFunctionPointer rosic::getWaveformFunction(int waveformIndex)
{
  switch( waveformIndex )
  {
  case SINE:     return &sin;          break;
  case SAW:      return &sawWave;      break;
  case SQUARE:   return &sqrWave;      break;
  case TRIANGLE: return &triWave;      break;
  default:       return &zeroFunction; break;
  }
}


//=================================================================================================
// class ModulationWaveformRenderer:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulationWaveformRenderer::ModulationWaveformRenderer()
{
  cWave    = SINE;
  mWave    = SINE;
  cFreq    = 1;
  mFreq    = 1;       
  cPhase   = 0.0;
  mPhase   = 0.0;
  modIndex = 1.0;
}

//-------------------------------------------------------------------------------------------------
// others:

void ModulationWaveformRenderer::assignRenderingParameters(UnaryFunctionPointer &cFun, 
                                                           UnaryFunctionPointer &mFun, 
                                                           double &wc, double &wm, 
                                                           double &pc, double &pm,
                                                           int tableLength)
{
  double tmp = 2.0*PI/tableLength;         // fundamental frequency w0
  wc         = tmp*cFreq;                  // carrier frequency
  wm         = tmp*mFreq;                  // modulator frequency
  tmp        = 2.0*PI/360.0;               // phase conversion factor
  pc         = tmp*cPhase;                 // carrier phase
  pm         = tmp*mPhase;                 // modulator phase
  cFun       = getWaveformFunction(cWave); // pointer to carrier waveform function
  mFun       = getWaveformFunction(mWave); // pointer to modulator waveform fucntion
}


//=================================================================================================
// class AlgorithmicWaveformRenderer:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AlgorithmicWaveformRenderer::AlgorithmicWaveformRenderer()
{
  algorithm = PHASE_MODULATION;
}

AlgorithmicWaveformRenderer::~AlgorithmicWaveformRenderer()
{

}

//-------------------------------------------------------------------------------------------------
// waveform rendering:

void AlgorithmicWaveformRenderer::renderWaveform(double *targetBuffer, int length)
{
  // dispatch the rendering request to the appropriate renderer:
  switch( algorithm )
  {
  case AMPLITUDE_MODULATION: 
    amplitudeModulationRenderer.renderWaveform(targetBuffer, length); break;
  case PHASE_MODULATION:     
    phaseModulationRenderer.renderWaveform(    targetBuffer, length); break;
  case RING_MODULATION:     
    ringModulationRenderer.renderWaveform(     targetBuffer, length); break;
  default: 
    fillWithZeros(targetBuffer, length);
  }
}