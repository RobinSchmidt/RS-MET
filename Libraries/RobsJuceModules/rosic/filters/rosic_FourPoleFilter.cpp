//#include "rosic_FourPoleFilter.h"
//using namespace rosic;

// Construction/Destruction:

FourPoleFilter::FourPoleFilter()
{
  parameters = new FourPoleFilterParameters;
  isMaster   = true;
  freq       = 1000.0;
  preGain = 1.0;
  reset();
  updateFilterCoefficients();
}

FourPoleFilter::~FourPoleFilter()
{
  if( isMaster == true )
    delete parameters;
}

// Setup:

void FourPoleFilter::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    parameters->sampleRate    = newSampleRate;
    parameters->sampleRateRec = 1.0 / newSampleRate;
  }
  updateFilterCoefficients();
  for(unsigned int s=0; s<slaves.size(); s++)
    slaves[s]->updateFilterCoefficients();
}

void FourPoleFilter::setMode(int newMode)
{
  parameters->mode = newMode;
  //secondStageIsActive = parameters->twoStageSwitch;

  unsigned int s;
  for(s = 0; s < slaves.size(); s++)
    slaves[s]->secondStageIsActive = this->secondStageIsActive;
  updateFilterCoefficients();
  for(s = 0; s < slaves.size(); s++)
    slaves[s]->updateFilterCoefficients();
}

void FourPoleFilter::useTwoStages(bool shouldUseTwoStages)
{
  parameters->twoStageSwitch = shouldUseTwoStages;
  secondStageIsActive        = parameters->twoStageSwitch;

  unsigned int s;
  reset();
  updateFilterCoefficients();
  for(s=0; s<slaves.size(); s++)
  {
    slaves[s]->secondStageIsActive = this->secondStageIsActive;
    slaves[s]->reset();
    slaves[s]->updateFilterCoefficients();
  }
}

// Inquiry:

/*
int FourPoleFilter::getMode()
{
  return parameters->mode;
}

bool FourPoleFilter::usesTwoStages()
{
  return parameters->twoStageSwitch;
}

double FourPoleFilter::getFrequency()
{
  return freq;
}

double FourPoleFilter::getQ()
{
  return parameters->q;
}

double FourPoleFilter::getGain()
{
  return parameters->gainDb;
}

double FourPoleFilter::getMorph()
{
  return parameters->morph;
}
*/

double FourPoleFilter::getMagnitudeAt(double frequency)
{
  double omega = 2.0 * PI * frequency * parameters->sampleRateRec; // optimize to one mul
  double c1  = cos(omega);
  double c2  = cos(2.0*omega);

  double b0  =  b0_s1;
  double b1  =  b1_s1;
  double b2  =  b2_s1;
  double a1  = -a1_s1; // we use the other sign-convention in the getSample-function
  double a2  = -a2_s1;

  double num = b0*b0 + b1*b1 + b2*b2 + 2.0*(b0*b1 + b1*b2)*c1 + 2.0*b0*b2*c2;
  double den = 1.0   + a1*a1 + a2*a2 + 2.0*(   a1 + a1*a2)*c1 + 2.0*   a2*c2;
  double mag = sqrt(num/den);

  if( secondStageIsActive )
  {
    b0   =  b0_s2;
    b1   =  b1_s2;
    b2   =  b2_s2;
    a1   = -a1_s2;
    a2   = -a2_s2;

    num  = b0*b0 + b1*b1 + b2*b2 + 2.0*(b0*b1 + b1*b2)*c1 + 2.0*b0*b2*c2;
    den  = 1.0   + a1*a1 + a2*a2 + 2.0*(   a1 + a1*a2)*c1 + 2.0*   a2*c2;
    mag *= sqrt(num/den);
  }

  return preGain*mag;
}

void FourPoleFilter::getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, 
                                          bool inDecibels, bool accumulate)
{
  int k;
  if( inDecibels == false )
  {
    if( accumulate == false )
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] = getMagnitudeAt(frequencies[k]);
    }
    else
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] *= getMagnitudeAt(frequencies[k]);
    }
  }
  else
  {
    if( accumulate == false )
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] = RAPT::rsAmpToDbWithCheck(getMagnitudeAt(frequencies[k]), 0.000001);
    }
    else
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] += RAPT::rsAmpToDbWithCheck(getMagnitudeAt(frequencies[k]), 0.000001);
    }
  }
}

// Polyphony handling:

void FourPoleFilter::addSlave(FourPoleFilter* newSlave)
{
  // add the new slave to the vector of slaves:
  slaves.push_back(newSlave);

  // delete the original parameter-set of the new slave and redirect it to ours (with some safety 
  // checks):
  if( newSlave->parameters != NULL && newSlave->parameters != this->parameters )
  {
    delete newSlave->parameters;
    newSlave->parameters = this->parameters;
  }
  else
  {
    DEBUG_BREAK; 
    // the object to be added as slave did not contain a valid parameter-pointer - maybe it has 
    // been already added as slave to another master?
  }

  // set the isMaster-flag of the new slave to false: 
  newSlave->isMaster = false;

  // this flag will prevent the destructor of the slave from trying to delete the parameter-set 
  // which is now shared - only masters delete their parameter-set on destruction
}

// Miscellaneous:

void FourPoleFilter::reset()
{
  xL_s1_d1 = 0.0;
  xR_s1_d1 = 0.0;
  xL_s1_d2 = 0.0;  
  xR_s1_d2 = 0.0;   
  yL_s1_d1 = 0.0;  
  yR_s1_d1 = 0.0;   
  yL_s1_d2 = 0.0;    
  yR_s1_d2 = 0.0;

  xL_s2_d1 = 0.0;
  xR_s2_d1 = 0.0;
  xL_s2_d2 = 0.0;  
  xR_s2_d2 = 0.0;   
  yL_s2_d1 = 0.0;  
  yR_s2_d1 = 0.0;   
  yL_s2_d2 = 0.0;    
  yR_s2_d2 = 0.0;
}

void FourPoleFilter::resetBuffersStage2()
{
  xL_s2_d1 = 0.0;
  xR_s2_d1 = 0.0;
  xL_s2_d2 = 0.0;  
  xR_s2_d2 = 0.0;   
  yL_s2_d1 = 0.0;  
  yR_s2_d1 = 0.0;   
  yL_s2_d2 = 0.0;    
  yR_s2_d2 = 0.0;
}
