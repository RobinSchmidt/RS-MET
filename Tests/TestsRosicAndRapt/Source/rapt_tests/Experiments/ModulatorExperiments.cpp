//#include <JuceHeader.h>
using namespace RAPT;

void attackDecayEnvelope()
{
  int N    = 1000;
  int nOff = 800;     // note-off sample instant
  int key  = 64;
  int vel  = 64;


  std::vector<double> y(N);


  rsAttackDecayFilter<double> flt;
  flt.setAttackSamples(5);
  flt.setDecaySamples(15);

  // plot DC response:
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(1.0);
  //rsPlotVector(y);
  double dcGain = flt.getGainAtDC();  // should be 22.335... - yep, works




  rsAttackDecayEnvelope<double> env;
  flt.setAttackSamples(20);
  flt.setDecaySamples(100);
  env.setSustain(0.5);
  env.noteOn(key, vel);
  for(int n = 0; n < N; n++)
    y[n] = env.getSample();
  //env.noteOff(key, vel);  // why does note-off need key and vel?

  rsPlotVector(y);

  // Observations:
  // -the sustain level works but using nonzero sustain slightly changes the attack-time and 
  //  maximum peak level: the peak gets higher and occurs later with increasing sustain, for 
  //  example, with attack = 20, decay = 100, sustain = 0.5, the actual peak occurs at 24 samples
  //  and has a height of around 1.07 (as per the settings, it should occur at sample 20 witha 
  //  height of 1), with sustain = 1, we get apeak at sample 33 with height 1.16
  
}
