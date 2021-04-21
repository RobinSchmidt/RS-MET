//#include "rosic_PitchShifter.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

PitchShifter::PitchShifter()
{
  // initialize parameters:
  detuneCoarse      = 7.0;
  detuneFine        = 0.0;
  sampleRate        = 44100.0;
  //grainLength       = 20.0;

  tapIn             = 1;
  distanceIncrement = 0.0;
  //incrementSign     = 1.0;
  feedbackFactor    = 0.0;
  dry               = 0.0;
  wet               = 1.0;
  wetPolarity       = 1.0;
  reversePlayback   = false;
  antiAliasingIsOn  = false;
  //mono              = false;

  //delayLineLength   = 32768;
  delayLineLength   = 262144;
  delayLineL        = NULL;
  delayLineR        = NULL;
  delayLineL        = new double[delayLineLength];
  delayLineR        = new double[delayLineLength];


  setGrainLength(20.0);
  updateIncrement();
  initTapDistances();
  reset();

  antiAliasFilterL.setSubDivision(1.0);
  antiAliasFilterR.setSubDivision(1.0);
}

PitchShifter::~PitchShifter()
{
  if( delayLineL != NULL )
    delete[] delayLineL;
  if( delayLineR != NULL )
    delete[] delayLineR;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void PitchShifter::setSampleRate(double /*newSampleRate*/)
{
  //antiAliasFilterL.setSampleRate(newSampleRate);
  //antiAliasFilterR.setSampleRate(newSampleRate);
  updateDistanceVariables();
  updateIncrement();
  //initTapDistances();
}

void PitchShifter::setDetuneCoarse(double newDetuneCoarse)
{
  detuneCoarse = newDetuneCoarse;
  updateIncrement();
  //initTapDistances();
}

void PitchShifter::setDetuneFine(double newDetuneFine)
{
  detuneFine = newDetuneFine;
  updateIncrement();
  //initTapDistances();
}

void PitchShifter::setGrainLength(double newGrainLength)
{
  if( newGrainLength >= 1.0 && newGrainLength != grainLength )
  {
    if( 0.001*newGrainLength*sampleRate < delayLineLength-10 ) // some arbitrary safety-margin
    {
      grainLength = newGrainLength;
      updateDistanceVariables();
    }
  }
}

void PitchShifter::setFeedback(double newFeedback)
{
  feedbackFactor = 0.01 * newFeedback;
}

void PitchShifter::setDryWet(double newDryWet)
{
  wet = 0.01 * newDryWet;
  dry = 1.0 - wet;
}

void PitchShifter::setReversePlayback(bool shouldBeReverse)
{
  reversePlayback = shouldBeReverse;
}

void PitchShifter::setNegativePolarity(bool shouldBeNegative)
{
  if( shouldBeNegative == true )
    wetPolarity = -1.0;
  else
    wetPolarity = +1.0;
}

void PitchShifter::setAntiAliasing(bool shouldAntiAlias)
{
  antiAliasingIsOn = shouldAntiAlias;
  antiAliasFilterL.reset();
  antiAliasFilterR.reset();
}

//-------------------------------------------------------------------------------------------------
// inquiry:
/*
double PitchShifter::getDetuneCoarse()
{
  return detuneCoarse;
}

double PitchShifter::getDetuneFine()
{
  return detuneFine;
}

double PitchShifter::getGrainLength()
{
  return grainLength;
}

double PitchShifter::getFeedback()
{
  return 100.0 * feedbackFactor;
}

double PitchShifter::getDryWet()
{
  return 100.0 * wet;
}

bool PitchShifter::isPlaybackReverse()
{
  return reversePlayback;
}

bool PitchShifter::isPolarityNegative()
{
  if( wetPolarity == -1.0 )
    return true;
  else
    return false;
}

bool PitchShifter::getAntiAliasing()
{
  return antiAliasingIsOn;
}
*/
//-------------------------------------------------------------------------------------------------
// others:

void PitchShifter::reset()
{
  for(int i=0; i<delayLineLength; i++)
  {
    delayLineL[i] = 0.0;
    delayLineR[i] = 0.0;
  }
  antiAliasFilterL.reset();
  antiAliasFilterR.reset();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void PitchShifter::updateIncrement()
{
  // calculate the total detune in semitones
  double detune    = detuneCoarse + 0.01*detuneFine;

  // calculate, to which tap-increment this corresponds:
  double increment = RAPT::rsPitchOffsetToFreqFactor(detune);

  // the difference between the hypotheticial increment for the read pointers and a unit increment
  // (which the write pointer has) gives the increment for the distance between the read and write
  // pointers:
  distanceIncrement = -(increment-1.0);

  // adjust the cutoff frequency of the anti-aliasing filter:
  if( increment > 1.0 )
  {
    antiAliasFilterL.setSubDivision(increment);
    antiAliasFilterR.setSubDivision(increment);
  }
  else
  {
    antiAliasFilterL.setSubDivision(1.0);
    antiAliasFilterR.setSubDivision(1.0);
  }
}

void PitchShifter::initTapDistances()
{
  distance1 = maxDistance;
  distance2 = 0.5*maxDistance;
}

void PitchShifter::updateDistanceVariables()
{
  maxDistance     = 0.001*grainLength*sampleRate;
  maxDistanceHalf = 0.5*maxDistance;
  maxDistanceRec  = 1.0/maxDistance;
}




/*
Ideas:
-instead using just two tap-out pointers, let the user select, how many
-let the user select the window function, too
-the overlap-to-1 condition of the windows can be enforced in realtime for any window and any 
 overlap factor by just keeping track of the sum of the weights during output accumulation and 
 diving the output by that number afterwards. of course, that may change the effective shape of the 
 window (unless it already obeys the sum-to-1 condition for given overlap factor)....but in a 
 desirable way for this application
-currently, we use a cosine-squared window with an overlap of 1/2 because that naturally gives the
 sum-to-one condition - but it may be advatageous to take larger overlaps: 3/2, 3/4, 4/5 (the 
 distance between the taps is 1/2, 1/3, 1/4, 1/5, ...)
-do tests: with a stable sinusoidal input, they should be all perfect - but differences may become 
 obvious when taking an input that features vibrato and/or an envelope. Test also with DC input. 
 That should also be perfect, i.e. change nothing (except transient artifacts, i.e. some sort 
 of fade-in at the very beginning). Test it also with material that features strong transients. 
 Maybe try an impulse (train), too.

*/