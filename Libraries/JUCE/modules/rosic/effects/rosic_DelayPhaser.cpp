//#include "rosic_DelayPhaser.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DelayPhaser::DelayPhaser()
{
  fb1 = 0.0;
  fb2 = 0.0;
  fb3 = 0.0;
  reset();
}

DelayPhaser::~DelayPhaser()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void DelayPhaser::setSampleRate(double newSampleRate)
{
  phaser1.setSampleRate(newSampleRate);
  phaser2.setSampleRate(newSampleRate);
  delay.setSampleRate(newSampleRate);
}

void DelayPhaser::setTempoInBPM(double newTempoInBPM)
{
  phaser1.setTempoInBPM(newTempoInBPM);
  phaser2.setTempoInBPM(newTempoInBPM);
  delay.setTempoInBPM(newTempoInBPM);
}

//-------------------------------------------------------------------------------------------------
// others:

void DelayPhaser::reset()
{
  dL = dR = p2L = p2R = 0.0;
  phaser1.reset();
  phaser2.reset();
  delay.reset();
}  

