//#include "rosic_OscillatorBank.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

OscillatorBank::OscillatorBank(int newMaxNumSines)
{
  sampleRate           = 44100.0; 
  fundamentalFrequency = 50.0;
  masterAmplitude      = 1.0;
  if( newMaxNumSines >= 1 )
    maxNumSines = newMaxNumSines;
  else
    maxNumSines = 512;
  numSines = 400;

  // create and initialize all the required arrays:
  a1 = new double[maxNumSines];
  s1 = new double[maxNumSines];
  s2 = new double[maxNumSines];
  c1 = new double[maxNumSines];
  c2 = new double[maxNumSines];

  relativeFrequencies = new double[maxNumSines];
  relativeAmplitudes  = new double[maxNumSines];
  phases              = new double[maxNumSines];

  initToSawWave();
  calculateAllCoefficients();
  triggerAll();
}

OscillatorBank::~OscillatorBank()
{
  delete[] a1;
  delete[] s1;
  delete[] s2;
  delete[] c1;
  delete[] c2;
  delete[] relativeFrequencies;     
  delete[] relativeAmplitudes;  
  delete[] phases;      
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void OscillatorBank::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.0)
    sampleRate = newSampleRate;

  double omega;
  for(int i=0; i<maxNumSines; i++)
  {
    omega = 2.0*PI*fundamentalFrequency*relativeFrequencies[i] / sampleRate;
    a1[i] = 2.0*cos(omega);
  }
  triggerAll();
}

void OscillatorBank::setNumSines(int newNumSines)
{
  if( newNumSines >= 0 && newNumSines <= maxNumSines )
    numSines = newNumSines;
  else 
    numSines = 1;
}

void OscillatorBank::setFundamentalFrequency(double newFundamentalFrequency)
{
  fundamentalFrequency = newFundamentalFrequency;
  calculateAllCoefficients();
}

void OscillatorBank::setMasterAmplitude(double newMasterAmplitude)
{
  masterAmplitude = newMasterAmplitude;
}

void OscillatorBank::setRelativeFrequencies(double *newRelativeFrequencies)
{
  for(int i=0; i<maxNumSines; i++)
    relativeFrequencies[i] = newRelativeFrequencies[i];
  calculateAllCoefficients();
}

void OscillatorBank::setRelativeAmplitudes(double *newRelativeAmplitudes)
{
  for(int i=0; i<maxNumSines; i++)
    relativeAmplitudes[i] = newRelativeAmplitudes[i];
}

//-------------------------------------------------------------------------------------------------
// event handling:

void OscillatorBank::triggerAll()
{
  double omega;
  for(int i=0; i<maxNumSines; i++)
  {
    omega = 2.0*PI*fundamentalFrequency*relativeFrequencies[i] / sampleRate;
    s1[i] = sin(phases[i] -     omega);
    s2[i] = sin(phases[i] - 2.0*omega);
    c1[i] = cos(phases[i] -     omega);
    c2[i] = cos(phases[i] - 2.0*omega);
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void OscillatorBank::initToSawWave()
{
  for(int i=0; i<maxNumSines; i++)
  {
    relativeFrequencies[i] = (double) (i+1);
    relativeAmplitudes[i]  = 1.0 / ((double) (i+1));
    phases[i]              = 0.0;
  }
}

void OscillatorBank::calculateAllCoefficients()
{
  double omega;
  for(int i=0; i<maxNumSines; i++)
  {
    omega = 2.0*PI*fundamentalFrequency*relativeFrequencies[i] / sampleRate;
    a1[i] = 2.0*cos(omega);
  }
}


