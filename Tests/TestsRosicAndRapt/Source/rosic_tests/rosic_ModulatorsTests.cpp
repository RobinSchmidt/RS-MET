//#include "rosic_ModulatorsTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
//#include "../Shared/Plotting/rosic_Plotter.h"
using namespace rosic;

void rotes::testConsecutiveExponentialDecay()
{
  static const int numSamples = 3000;

  double indices[numSamples];
  double impulseResponse[numSamples];
  double stepResponse[numSamples];
  RAPT::rsArrayTools::fillWithIndex(indices, numSamples);



  AttackDecayEnvelope envGen;
  envGen.setPeakTime(1000 * 500.0 / 44100.0);
  envGen.setDecayTimeConstant(10.0);
  envGen.trigger(false);

  int n;
  for(n=0; n<numSamples; n++)
    impulseResponse[n] = envGen.getSample();



  //RAPT::rsArrayTools::copy(impulseResponse, stepResponse, numSamples);
  //RAPT::rsArrayTools::cumulativeSum(stepResponse, numSamples, 1);

  RAPT::rsArrayTools::cumulativeSum(impulseResponse, stepResponse, numSamples, 1);


  double scaler = 1.0 / RAPT::rsArrayTools::maxValue(stepResponse, numSamples);
  RAPT::rsArrayTools::scale(stepResponse, stepResponse, numSamples, scaler);
  //ste


  plotData(numSamples, indices, impulseResponse, stepResponse);








  int dummy = 0;



  /*
  // set up the WaveTable:
  static const int prototypeLength = 2048;
  double prototypeWaveFlat[prototypeLength];
  double *prototypeWave[2];
  prototypeWave[0] = &prototypeWaveFlat[0];
  prototypeWave[1] = &prototypeWaveFlat[0];
  int n;
  for(n=0; n<prototypeLength; n++)
    prototypeWaveFlat[n] = sawWave( 2*PI*n / prototypeLength);
  MipMappedWaveTableStereo waveTable;
  waveTable.setWaveform(prototypeWave, prototypeLength);

  // set up the oscillator:
  OscillatorStereo osc;
  osc.setWaveTableToUse(&waveTable);
  osc.setStereoPhaseShift(90.0);
  osc.setStartPhase(95.0);

  // obtain the data for the plot:
  static const int plotLength = 160;
  double plotIndices[plotLength];
  fillWithIndex(plotIndices, plotLength);

  double plotDataFlat1[2*plotLength];
  double *plotData1[2];
  plotData1[0] = &plotDataFlat1[0];
  plotData1[1] = &plotDataFlat1[plotLength];
  osc.getWaveformForDisplay(plotData1, plotLength);

  osc.setStartPhase(96.0);
  double plotDataFlat2[2*plotLength];
  double *plotData2[2];
  plotData2[0] = &plotDataFlat2[0];
  plotData2[1] = &plotDataFlat2[plotLength];
  osc.getWaveformForDisplay(plotData2, plotLength);


  plotData(plotLength, plotIndices, plotData1[0], plotData2[0]);
  */
}