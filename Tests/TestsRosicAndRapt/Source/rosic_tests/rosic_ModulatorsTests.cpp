using namespace rotes;
using namespace rosic;

void rotes::testConsecutiveExponentialDecay()
{
  // We test rosic::AttackDecayEnvelope by recording and plotting its impulse response. We also
  // plot the cumulative sum of the impulse response which is (supposed to be) the step response.
  // (ToDo: verify that this is actually the case)

  // todo: rename to testAttackDecayEnv...we have some other such experiments elsewhere - maybe 
  // move the code together..but maybe not - this here is from rosic, the other is from rapt

  static const int numSamples = 3000;

  using AT = RAPT::rsArrayTools;

  double indices[numSamples];
  double impulseResponse[numSamples];
  double stepResponse[numSamples];
  AT::fillWithIndex(indices, numSamples);

  AttackDecayEnvelope envGen;
  envGen.setPeakTime(1000 * 500.0 / 44100.0);
  envGen.setDecayTimeConstant(10.0);
  envGen.trigger(false);
  for(int n = 0; n < numSamples; n++)
    impulseResponse[n] = envGen.getSample();

  AT::cumulativeSum(impulseResponse, stepResponse, numSamples, 1);
  AT::normalize(stepResponse, numSamples);
  // We need to scale it down because otherwise the values go into the thousands and totally dwarf
  // the plot of the impulse response.

  plotData(numSamples, indices, impulseResponse, stepResponse);

  // ToDo: maybe produce the step response not by using a cumulative sum but by actually feeding
  // a step...maybe verify if the results are the same - i think, they should...right?


  // what's this? it looks like the code in testOscillatorStereo. is this a copy/paste remnant?
  // if so, delete it:
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