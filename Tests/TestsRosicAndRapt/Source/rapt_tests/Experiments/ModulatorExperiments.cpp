#include "ModulatorExperiments.h"
//#include <JuceHeader.h>
using namespace RAPT;

void attackDecayEnvelope()
{
  int N = 200;

  rsAttackDecayEnvelope<double> env;

  std::vector<double> y(N);
  env.noteOn(64, 64);
  for(int n = 0; n < N; n++)
    y[n] = env.getSample();

  // maybe change the code such that the first sample is already nonzero - we need to make the 
  // attack one sample longe then...


  rsPlotVector(y);
}
