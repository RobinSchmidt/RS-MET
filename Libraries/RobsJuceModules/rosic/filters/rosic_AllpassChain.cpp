//#include "rosic_AllpassChain.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AllpassChain::AllpassChain()
{
  sampleRate        = 44100.0;
  oneOverSampleRate = 1.0 / sampleRate;
  mode              = FIRST_ORDER_ALLPASS;
  frequency         = 1000.0;
  q                 = 2.0;
  numStages         = 4;
  updateCoeffs();
  reset();
}

AllpassChain::~AllpassChain()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AllpassChain::setSampleRate(double newSampleRate, bool updateCoefficients)
{
  if( newSampleRate <= 0.0 )
  {
    DEBUG_BREAK;
    return;
  }
  sampleRate        = newSampleRate;
  oneOverSampleRate = 1.0 / sampleRate;
  if( updateCoefficients == true )
    updateCoeffs();
}

void AllpassChain::setMode(int newMode, bool updateCoefficients)
{
  if( newMode >= FIRST_ORDER_ALLPASS && newMode < NUM_FILTER_MODES )
    mode = newMode;
  if( updateCoefficients == true )
    updateCoeffs();
}

void AllpassChain::setFrequency(double newFrequency, bool updateCoefficients)
{
  frequency = RAPT::rsClip(newFrequency, 2.0, 20000.0);
  if( updateCoefficients == true )
    updateCoeffs();
}

void AllpassChain::setQ(double newQ, bool updateCoefficients)
{
  q = RAPT::rsMax(newQ, 0.0);
  if( updateCoefficients == true )
    updateCoeffs();
}

void AllpassChain::setNumStages(int newNumStages)
{
  RAPT::rsAssert(newNumStages <= maxNumStages, "We currently have a limit of 24 stages");
  if( newNumStages >= 0 && newNumStages <= maxNumStages )
    numStages = newNumStages;
}

//-------------------------------------------------------------------------------------------------
// others:

void AllpassChain::reset()
{
  for(int i=0; i<maxNumStages; i++)
  {
    x1[i] = 0.0;
    x2[i] = 0.0;
    y1[i] = 0.0;
    y2[i] = 0.0;
  }
}

void AllpassChain::updateCoeffs()
{
  switch(mode)
  {
  case FIRST_ORDER_ALLPASS:
    BiquadDesigner::calculateFirstOrderAllpassCoeffs(b0, b1, b2, a1, a2, oneOverSampleRate, 
      frequency);  
    break;
  case SECOND_ORDER_ALLPASS:
    BiquadDesigner::calculateCookbookAllpassCoeffs(b0, b1, b2, a1, a2, oneOverSampleRate, 
      frequency, q);  
    break;
  default:
    BiquadDesigner::makeBypassBiquad(b0, b1, b2, a1, a2);  
    break;
  }
}
