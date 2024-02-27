//#include "rosic_Equalizer.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

Equalizer::Equalizer()
{
  sampleRate       = 44100.0;
  globalGainFactor = 1.0;
  mono             = false;
  bands.reserve(8);              // Added 2024/02/27 - might be a good idea.
}

Equalizer::~Equalizer()
{
  bands.clear();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void Equalizer::setSampleRate(double newSampleRate)
{
  if( newSampleRate <= 0.0 )
  {
    DEBUG_BREAK;
    return;
  }
  sampleRate = newSampleRate;
  for(unsigned int i=0; i<bands.size(); i++)
    bands[i].setSampleRate(sampleRate);
}

int Equalizer::addBand(int newMode, double newFrequency, double newGain, double newBandwidth)
{
  int result = -1;

  TwoPoleFilter newBand;

  newBand.setParameters(newMode, newFrequency, newGain, newBandwidth, false);
  newBand.setSampleRate(sampleRate);

  bands.push_back(newBand);
  result = (int) bands.size()-1;

  return result;
}

bool Equalizer::removeBand(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;   // index out of range
    return false;
  }
  bands.erase(bands.begin()+index);
  return true;
}

bool Equalizer::removeAllBands()
{
  bool result = (bands.size() != 0);
  bands.clear();
  return result;
}

bool Equalizer::modifyBand(unsigned int index, int newMode, double newFrequency, double newGain, double newBandwidth)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return false;
  }
  bands[index].setParameters(newMode, newFrequency, newGain, newBandwidth, true);
  return true;
}

bool Equalizer::setBandMode(int index, int newMode)
{
  return modifyBand(index, newMode, getBandFrequency(index), getBandGain(index), getBandBandwidth(index));
}

bool Equalizer::setBandFrequency(int index, double newFrequency)
{
  return modifyBand(index, getBandMode(index), newFrequency, getBandGain(index), getBandBandwidth(index));
}

bool Equalizer::setBandGain(int index, double newGain)
{
  return modifyBand(index, getBandMode(index), getBandFrequency(index), newGain, getBandBandwidth(index));
}

bool Equalizer::setBandBandwidth(int index, double newBandwidth)
{
  return modifyBand(index, getBandMode(index), getBandFrequency(index), getBandGain(index), newBandwidth);
}

bool Equalizer::setLowerBandedgeFrequency(unsigned int index, double newLowerBandedgeFrequency)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return false;
  }
  bands[index].setLowerBandedgeFrequency(newLowerBandedgeFrequency);
  return true;
}

bool Equalizer::setUpperBandedgeFrequency(unsigned int index, double newUpperBandedgeFrequency)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return false;
  }
  bands[index].setUpperBandedgeFrequency(newUpperBandedgeFrequency);
  return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

int Equalizer::getNumBands()
{
  return (int) bands.size();
}

int Equalizer::getBandMode(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0;
  }
  return bands[index].getMode();
}

double Equalizer::getBandFrequency(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  return bands[index].getFrequency();
}

double Equalizer::getBandGain(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  return bands[index].getGain();
}

double Equalizer::getBandBandwidth(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  return bands[index].getBandwidth();
}

double Equalizer::getLowerBandedgeFrequency(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  return bands[index].getLowerBandedgeFrequency();
}

double Equalizer::getUpperBandedgeFrequency(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  return bands[index].getUpperBandedgeFrequency();
}

bool Equalizer::doesModeSupportBandwidth(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  return bands[index].doesModeSupportBandwidth();
}

bool Equalizer::doesModeSupportGain(unsigned int index)
{
  if( /*index < 0 ||*/ index >= bands.size() )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  return bands[index].doesModeSupportGain();
}

void Equalizer::getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins)
{
  for(int k = 0; k < numBins; k++)
  {
    magnitudes[k] = globalGainFactor;
    if( frequencies[k] < 0.5*sampleRate )
    {
      for(unsigned int s=0; s<bands.size(); s++)
        magnitudes[k] *= bands[s].getMagnitudeAt(frequencies[k]);
      magnitudes[k] = RAPT::rsAmpToDbWithCheck(magnitudes[k], 0.0001);
    }
    else
      magnitudes[k] = -120.0;
  }

//  // debug:
//  double m[1000];
//  copy(magnitudes, m, rmin(numBins, 1000));
//  int dummy = 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void Equalizer::reset()
{
  for(unsigned int s=0; s<bands.size(); s++)
    bands[s].reset();
}
