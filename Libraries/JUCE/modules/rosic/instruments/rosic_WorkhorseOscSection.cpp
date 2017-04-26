#include "rosic_WorkhorseOscSection.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WorkhorseOscSection::WorkhorseOscSection()
{

}

WorkhorseOscSection::~WorkhorseOscSection()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings: 

void WorkhorseOscSection::setPlaybackFrequencyNominal(double newFrequency)
{
  samplePlayerTopLeft.setPlaybackFrequencyNominal(    newFrequency);
  samplePlayerTopRight.setPlaybackFrequencyNominal(   newFrequency);
  samplePlayerBottomLeft.setPlaybackFrequencyNominal( newFrequency);
  samplePlayerBottomRight.setPlaybackFrequencyNominal(newFrequency);
}

void WorkhorseOscSection::setSampleRate(double newSampleRate)
{
  samplePlayerTopLeft.setSampleRate(    newSampleRate);
  samplePlayerTopRight.setSampleRate(   newSampleRate);
  samplePlayerBottomLeft.setSampleRate( newSampleRate);
  samplePlayerBottomRight.setSampleRate(newSampleRate);
}

void WorkhorseOscSection::reset()
{
  samplePlayerTopLeft.reset();
  samplePlayerTopRight.reset();
  samplePlayerBottomLeft.reset();
  samplePlayerBottomRight.reset();
}