#include "rosic_VectorSamplePlayer.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

VectorSamplePlayer::VectorSamplePlayer()
{
  topLeftGainFactor = topRightGainFactor = bottomLeftGainFactor = bottomRightGainFactor = 0.5;
}

VectorSamplePlayer::~VectorSamplePlayer()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings: 

void VectorSamplePlayer::setSampleRate(double newSampleRate)
{
  samplePlayerTopLeft.setSampleRate(    newSampleRate);
  samplePlayerTopRight.setSampleRate(   newSampleRate);
  samplePlayerBottomLeft.setSampleRate( newSampleRate);
  samplePlayerBottomRight.setSampleRate(newSampleRate);
  xLfo.setSampleRate(newSampleRate);
  yLfo.setSampleRate(newSampleRate);
}

void VectorSamplePlayer::setBeatsPerMinute(double newBpm)
{
  xLfo.setBeatsPerMinute(newBpm);
  yLfo.setBeatsPerMinute(newBpm);
}

void VectorSamplePlayer::setKeyAndVel(int newKey, int newVel)
{ 
  samplePlayerTopLeft.setKeyAndVel(newKey, newVel);
  samplePlayerTopRight.setKeyAndVel(newKey, newVel);
  samplePlayerBottomLeft.setKeyAndVel(newKey, newVel);
  samplePlayerBottomRight.setKeyAndVel(newKey, newVel);
}

void VectorSamplePlayer::setPlaybackFrequencyNominal(double newFrequency)
{
  samplePlayerTopLeft.setPlaybackFrequencyNominal(    newFrequency);
  samplePlayerTopRight.setPlaybackFrequencyNominal(   newFrequency);
  samplePlayerBottomLeft.setPlaybackFrequencyNominal( newFrequency);
  samplePlayerBottomRight.setPlaybackFrequencyNominal(newFrequency);
}



void VectorSamplePlayer::reset()
{
  samplePlayerTopLeft.reset();
  samplePlayerTopRight.reset();
  samplePlayerBottomLeft.reset();
  samplePlayerBottomRight.reset();
  xLfo.reset();
  yLfo.reset();
}