//#include "rosic_NoiseGenerator.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

NoiseGenerator::NoiseGenerator(int bufferLengthToAllocate)
{
  mutex.lock();
  if( bufferLengthToAllocate > 0 )
  {
    length = RAPT::rsNextPowerOfTwo(bufferLengthToAllocate); 
    buffer = new(std::nothrow) double[length];
    if( buffer == NULL )
    {
      pointerInvalid = true;
      memoryIsShared = true;
    }
    else
    {
      pointerInvalid = false;
      memoryIsShared = false;
    }
  }
  else
  {
    length         = 0;
    buffer         = NULL;
    pointerInvalid = true;
    memoryIsShared = true;
  }
  mutex.unlock();

  sampleRate  = 44100.0;
  slope       = 0.0;
  lowestFreq  = 20.0;
  highestFreq = 20000.0;
  seed        = 0;
  createNoiseSequence();
  trigger();
}

NoiseGenerator::~NoiseGenerator()
{
  freeMemoryIfNotShared();
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void NoiseGenerator::setSharedMemoryAreaToUse(void *startAddress, int sizeInBytes)
{
  mutex.lock();

  freeMemoryIfNotShared();

  int sampleSizeInBytes = sizeof(double);

  length = sizeInBytes/sampleSizeInBytes; 
  if( !RAPT::rsIsPowerOfTwo(length) )
    length = RAPT::rsNextPowerOfTwo(length) / 2;  // power of 2 and <= available space
  
  buffer = (double*) startAddress;

  if( buffer != NULL )
    pointerInvalid = false;
  else
    pointerInvalid = true;
  memoryIsShared  = true;

  createNoiseSequence();

  mutex.unlock();
}

void NoiseGenerator::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;
    createNoiseSequence();
  }
}

void NoiseGenerator::setSpectralSlope(double newSlope)
{
  if( newSlope != slope )
  {
    slope = newSlope;
    createNoiseSequence();
  }
}

void NoiseGenerator::setLowestFrequency(double newLowestFrequency)
{
  if( newLowestFrequency >= 0.0 && newLowestFrequency != lowestFreq )
  {
    lowestFreq = newLowestFrequency;
    createNoiseSequence();
  }
}

void NoiseGenerator::setHighestFrequency(double newHighestFrequency)
{
  if( newHighestFrequency >= 0.0 && newHighestFrequency != highestFreq )
  {
    highestFreq = newHighestFrequency;
    createNoiseSequence();
  }
}

void NoiseGenerator::setRandomPhaseSeed(int newSeed)
{
  if( seed != newSeed )
  {
    seed = newSeed;
    createNoiseSequence();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void NoiseGenerator::trigger()
{
  readIndex = 0;
}

void NoiseGenerator::createNoiseSequence()   
{
  mutex.lock();
  if( pointerInvalid )
  {
    mutex.unlock();
    return;
  }

  double *magnitudes = new double[length/2];
  double *phases     = new double[length/2];
  magnitudes[0]      = 0.0;
  phases[0]          = 0.0;

  RAPT::rsRandomUniform(0.0, 1.0, seed);  // init random number generator
  int k;
  for(k=1; k<length/2; k++)
  {
    magnitudes[k] = RAPT::rsDbToAmp(slope*log2(k));
    phases[k]     = RAPT::rsRandomUniform(0.0, 2.0*PI);
  }

  // zero out the magnitudes below the lowest frequency and above the highest frequency:
  int lowBin  = (int) ceil( lowestFreq *length/sampleRate);
  int highBin = (int) floor(highestFreq*length/sampleRate);
  for(k=0; k<=lowBin; k++)
    magnitudes[k] = 0.0;
  for(k=highBin; k<length/2; k++)
    magnitudes[k] = 0.0;

  // create the noise by means of IFFT:
  FourierTransformerRadix2 transformer;
  transformer.setBlockSize(length);
  transformer.setDirection(FourierTransformerRadix2::INVERSE);
  transformer.getRealSignalFromMagnitudesAndPhases(magnitudes, phases, buffer);

  RAPT::rsArrayTools::normalize(buffer, length, 1.0);

  mutex.unlock();
}

void NoiseGenerator::freeMemoryIfNotShared()
{
  mutex.lock();
  if( memoryIsShared == false )
  {
    if( buffer != NULL )
    {
      delete[] buffer;
      buffer = NULL;
    }
    length         = 0;
    pointerInvalid = true;
  }
  mutex.unlock();
}
