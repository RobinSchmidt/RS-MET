//#include "rosic_WaveformGenerator.h"
//using namespace rosic;

WaveformGenerator::WaveformGenerator()
{
  algorithmIndex = 0;

  int i;
  for(i=0; i<tableLength; i++)
    wave[i] = 0.0;
  for(i=0; i<numParameters; i++)
    p[i] = 0.0;

  wave1 = sin;
  wave2 = sin;

  fourierTransformer.setBlockSize(tableLength);
  fourierTransformer.setDirection(FourierTransformerRadix2::INVERSE);
}

WaveformGenerator::~WaveformGenerator()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void WaveformGenerator::setAlgorithm(int newAlgorithmIndex)
{
  if( newAlgorithmIndex >= 0 && newAlgorithmIndex < NUM_ALGORITHMS )
  {
    algorithmIndex = newAlgorithmIndex;
    renderWaveform();
  }
}
    
void WaveformGenerator::setBasicWaveshape(int index, int newShape)
{
  double (*waveTmp) (double); 
  switch( newShape )
  {
  case SINE:     waveTmp = sin;     break;
  case SAW:      waveTmp = sawWave; break;
  case SQUARE:   waveTmp = sqrWave; break;
  case TRIANGLE: waveTmp = triWave; break;
  default:       waveTmp = sin;  
  }
  switch( index )
  {
  case 0: wave1 = waveTmp; break;
  case 1: wave2 = waveTmp; break;
  }
  renderWaveform();
}

void WaveformGenerator::setParameter(int index, double value)
{
  if( index >= 0 && index < numParameters )
  {
    p[index] = value;
    renderWaveform();
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void WaveformGenerator::clearWaveform()
{
  for(int i=0; i<tableLength; i++)
    wave[i] = 0.0;
}

void WaveformGenerator::clearMagnitudes()
{
  for(int i=0; i<tableLength/2; i++)
    magnitudes[i] = 0.0;
}

void WaveformGenerator::initPhases(double phaseValue)
{
  phases[0] = 0.0;  // this slot has the special role of representing the imaginary part at the 
                    // Nyquist-frequency
  for(int i=1; i<tableLength/2; i++)
    phases[i] = phaseValue;
}

void WaveformGenerator::renderWaveform()
{
  switch( algorithmIndex )
  {
  case AMPLITUDE_MODULATION: renderWaveformAmplitudeModulation();    break;
  case OCTAVES_OF_PRIMES:    renderWaveformOctavesOfPrimes();        break;
  case PHASE_MODULATION:     renderWaveformPhaseModulation();        break;
  case WINDOW_FUNCTION:      renderWaveformWindowFunction();         break;
  }
}

void WaveformGenerator::renderWaveformAmplitudeModulation()
{
  double w0 = 2.0*PI/tableLength;   // fundamental radian frequency  
  double wc = w0*p[0];              // radian frequency of the carrier  
  double wm = w0*p[1];              // radian frequency of the modulator
  double mi = p[2];                 // modulation index
  for(int n=0; n<tableLength; n++)
    wave[n] = wave1(wc*n) * (1+mi*wave2(wm*n));
  removeMean(wave, tableLength);
  normalize( wave, tableLength, 1.0);
}

void WaveformGenerator::renderWaveformOctavesOfPrimes()
{
  int lo = PrimeNumbers::findClosestLowerPrimeIndex(rmax((int)p[0],2));  // lowest prime index
  int hi = PrimeNumbers::findClosestLowerPrimeIndex(rmax((int)p[1],2));  // highest prime index

  if( wave1 == sin )
  {
    clearMagnitudes();
    initPhases(-PI/2);
    magnitudes[1] = 1.0;
    for(int i=lo; i<=hi; i++)
    {
      int prime = PrimeNumbers::getPrime(i);
      int h     = prime;                     // harmonic number
      while( h < tableLength/2 )
      {
        magnitudes[h] = 1.0 / (double) h;
        h *= 2;
      }
    }
    fourierTransformer.getRealSignalFromMagnitudesAndPhases(magnitudes, phases, wave);
  }
  else
  {
    clearWaveform();
    double w0 = 2.0*PI/tableLength;           // radian frequency of the fundamental
    if( p[0] == 1 )
    {
      // create the fundamental waveform:
      for(int n=0; n<tableLength; n++)        
        wave[n] += wave1(w0*n);
    }
    for(int i=lo; i<=hi; i++)
    {
      int prime = PrimeNumbers::getPrime(i);
      int h     = prime;                      // harmonic number
      while( h < tableLength/2 )
      {
        double a  = 1.0 / (double) h;         // amplitude
        for(int n=0; n<tableLength; n++)
          wave[n] += a * wave1(w0*n*h);
        h *= 2;
      }
    }
  }

  removeMean(wave, tableLength);
  normalize( wave, tableLength, 1.0);
}

void WaveformGenerator::renderWaveformPhaseModulation()
{
  double w0 = 2.0*PI/tableLength;   // fundamental radian frequency  
  double wc = w0*p[0];              // radian frequency of the carrier  
  double wm = w0*p[1];              // radian frequency of the modulator
  double mi = p[2];                 // modulation index
  for(int n=0; n<tableLength; n++)
    wave[n] = wave1( wc*n + mi*wave2(wm*n) );
  removeMean(wave, tableLength);
  normalize( wave, tableLength, 1.0);
}

void WaveformGenerator::renderWaveformWindowFunction()
{
  double p0 = rmax(p[0], 1.0);        // number of cycles of the carrier that fit into the table
  double p1 = rmax(p[1], 1.0);        // number of cycles of the window that fit into the table
  double wc = p0*2.0*PI/tableLength;  // carrier frequency
  double wm = wc/p1;                  // frequency of the widowing cosine
  double w  = p0/p1;                  // bandwidth parameter
  int    N  = (int) (tableLength/w);  // number of nonzero values

  clearWaveform();
  for(int n=0; n<N; n++)
    wave[n] = wave1(wc*n) * 0.5*( 1-cos(wm*n) );

  removeMean(wave, tableLength);
  normalize( wave, tableLength, 1.0);
}