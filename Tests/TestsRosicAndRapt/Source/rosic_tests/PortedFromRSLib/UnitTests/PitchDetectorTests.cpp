//#include "PitchDetectorTests.h"

template <class T>
bool isRangeLinear(T *buffer, int length, T min, T max)
{
  double factor = (max-min) / (double) (length-1);
  for(int i=0; i<length; i++)
  {
    if( buffer[i] != (T) (factor * T(i) + min) )
      return false;
  }
  return true;
}

bool testAutoCorrelationPitchDetector()
{
  bool testResult = true;

  //static const int bufferSize = 16;
  //double buffer[bufferSize];
  double signal[1000];
  RAPT::rsArrayTools::fillWithIndex(signal, 1000);

  rsAutoCorrelationPitchDetectorD pd;
  pd.setBufferSize(16);
  pd.setUpdateInterval(5);

  int start = 1;
  int size  = 3;

  /*
  pd.processBlock(&signal[start], size);

  rsCopyBuffer(pd.circularBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0],  3, 1.0, 3.0);
  testResult &= isRangeLinear(&buffer[3], 13, 0.0, 0.0);
  // should be 1, 2, 3, 0, 0, 0, ....

  rsCopyBuffer(pd.linearBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 16, 0.0, 0.0);
  // should be all zeros

  start += size;
  size   = 4;
  pd.processBlock(&signal[start], size);

  rsCopyBuffer(pd.circularBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 7, 1.0, 7.0);
  testResult &= isRangeLinear(&buffer[7], 9, 0.0, 0.0);
  // 1, ... , 7, 0, 0, ....

  rsCopyBuffer(pd.linearBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 9, 0.0, 0.0);
  testResult &= isRangeLinear(&buffer[9], 7, 1.0, 7.0);
  // 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7

  start += size;
  size   = 5;
  pd.processBlock(&signal[start], size);

  rsCopyBuffer(pd.circularBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 12, 1.0, 12.0);
  testResult &= isRangeLinear(&buffer[12], 4, 0.0,  0.0);
  // 1, ... , 12, 0, 0, ....

  rsCopyBuffer(pd.linearBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0],  4, 0.0,  0.0);
  testResult &= isRangeLinear(&buffer[4], 12, 1.0, 12.0);
  // 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12

  start += size;
  size   = 6;
  pd.processBlock(&signal[start], size);

  rsCopyBuffer(pd.circularBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0],  2, 17.0, 18.0);
  testResult &= isRangeLinear(&buffer[2], 14,  3.0, 16.0);
  // 17, 18, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16

  rsCopyBuffer(pd.linearBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 16, 3.0, 18.0);
  // 3, 4, ... , 17, 18

  start += size;
  size   = 16;
  pd.processBlock(&signal[start], size);

  rsCopyBuffer(pd.circularBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0],  2, 33.0, 34.0);
  testResult &= isRangeLinear(&buffer[2], 14, 19.0, 32.0);
  // 33, 34, 19, 20, ... , 32

  rsCopyBuffer(pd.linearBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 16, 19.0, 34.0);
  // 19, 20, ... , 33, 34

  start += size;
  size   = 20;
  pd.processBlock(&signal[start], size);

  rsCopyBuffer(pd.circularBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 16, 39.0, 54.0);
  // 39, ... , 54

  rsCopyBuffer(pd.linearBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 16, 39.0, 54.0);
  // 39, ... , 54

  start += size;
  size   = 8;
  pd.processBlock(&signal[start], size);

  rsCopyBuffer(pd.circularBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 8, 55.0, 62.0);
  testResult &= isRangeLinear(&buffer[8], 8, 47.0, 54.0);
  // 55,...,62,47,...,54

  rsCopyBuffer(pd.linearBuffer, buffer, bufferSize);
  testResult &= isRangeLinear(&buffer[0], 16, 47.0, 62.0);
  // 47,...,62
  */



  // now test with some more realistic input:
  double sampleRate = 44100;
  static const int numSamples = 10000;
  pd.setBufferSize(2048);
  pd.setUpdateInterval(256);
  pd.setSampleRate(sampleRate);

  double x[numSamples];
  double f = 430.0;

  double w = 2*PI*f/sampleRate;
  int n;
  for(n = 0; n < numSamples; n++)
    x[n] = sin(w*n);

  double estimate;
  start = 0;
  size  = 120;
  estimate = pd.processBlock(&x[start], size);
  testResult &= (estimate == 1000.0);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  testResult &= (estimate == 1000.0);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  size   = 250;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  size   = 256;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);
  start += size;
  estimate = pd.processBlock(&x[start], size);

  return testResult;
}
