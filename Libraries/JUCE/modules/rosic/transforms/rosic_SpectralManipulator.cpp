//#include "rosic_SpectralManipulator.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SpectralManipulator::SpectralManipulator()
{
  sampleRate = 44100.0;
}

SpectralManipulator::~SpectralManipulator()
{

}

//-------------------------------------------------------------------------------------------------
// magnitude spectrum manipulations (static):

void SpectralManipulator::applyContrast(double *magnitudes, int length, double power)
{
  int k;

  double maxMagnitude = 0.0;
  for(k=0; k<length; k++)
  {
    if( magnitudes[k] > maxMagnitude )
      maxMagnitude = magnitudes[k];
  }
  double normalizer = maxMagnitude / pow(maxMagnitude, power);

  for(k=0; k<length; k++)
    magnitudes[k] = normalizer * pow(magnitudes[k], power);
}

void SpectralManipulator::applySlope(double *magnitudes, int length, double slope /*,int cutoffBin*/)
{
  // this normalizer is good for waveforms with 1/n falloff of the harmonics:
  double normalizer = 1.0;
  if( slope > 0.0 )
    normalizer = 1.0 / RAPT::rsDbToAmp( 0.5*slope*log2((double)length) )   ;

  for(int k=0; k<length; k++)
    magnitudes[k] *= normalizer * RAPT::rsDbToAmp(slope*log2(k));
}

void SpectralManipulator::applyBrickwallLowpass(double *magnitudes, int length,
                                                int highestBinToKeep)
{
  for(int k=highestBinToKeep+1; k<length; k++)
    magnitudes[k] = 0.0;
}

void SpectralManipulator::applyBrickwallHighpass(double *magnitudes, int length,
                                                 int lowestBinToKeep)
{
  if(lowestBinToKeep > length)
    lowestBinToKeep = length;
  for(int k=0; k<lowestBinToKeep; k++)
    magnitudes[k] = 0.0;
}

void SpectralManipulator::applyEvenOddBalance(double *magnitudes, int length, double balance)
{
  double evenAmp, oddAmp;
  if( balance > 0.5 )
  {
    oddAmp  = 1.0;
    evenAmp = 2.0 - 2.0*balance;
  }
  else
  {
    evenAmp = 1.0;
    oddAmp  = 2.0 * balance;
  }

  int k;
  for(k=1; k<length; k+=2)
    magnitudes[k] *= oddAmp;
  for(k=2; k<length; k+=2)
    magnitudes[k] *= evenAmp;

  // starting at k=2 for the even harmonics leaves DC alone - this is particularly desirable for
  // very short blocklengths where the (moving) DC represents low frequency content
}

void SpectralManipulator::applyMagnitudeRandomization(double* /*magnitudes*/, int /*length*/,
                                                      double /*amount*/, int seed)
{
  randomUniform(0.0, 1.0, seed);
  // \todo: write a class for this PNRG and create an instance on the stack here to make it
  // thread-safe (i.e. avoid that other threads spoil the state of PNRG during we are using it. */

  //....
}

//-------------------------------------------------------------------------------------------------
// phase spectrum manipulations (static):

void SpectralManipulator::applyEvenOddPhaseShift(double *phases, int length, double shiftInDegrees)
{
  double phi = RAPT::rsDegreeToRadiant(0.5*shiftInDegrees);
  int k;
  for(k=1; k<length; k+=2)
    phases[k] += phi;
  for(k=2; k<length; k+=2)
    phases[k] -= phi;
}

void SpectralManipulator::applyStereoPhaseShift(double *phasesL, double *phasesR, int length,
                                                double shiftInDegrees)
{
  double phi = RAPT::rsDegreeToRadiant(0.5*shiftInDegrees);
  int k;
  for(k=1; k<length; k+=2)
  {
    phasesL[k] -= phi;
    phasesR[k] += phi;
  }
}

void SpectralManipulator::applyEvenOddStereoPhaseShift(double *phasesL, double *phasesR,
                                                       int length, double shiftInDegrees)
{
  double phi = RAPT::rsDegreeToRadiant(0.5*shiftInDegrees);
  int k;
  for(k=1; k<length; k+=2)
  {
    phasesL[k] += phi;
    phasesR[k] -= phi;
  }
  for(k=2; k<length; k+=2)
  {
    phasesL[k] -= phi;
    phasesR[k] += phi;
  }
}
