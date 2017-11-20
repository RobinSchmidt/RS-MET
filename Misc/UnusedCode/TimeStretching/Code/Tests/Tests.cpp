#include "TimeStretcher.h"
#include "../Bernsee/smbPitchShiftModified.cpp"

// x:  input signal (unprocessed)
// xN: length of x in samples
// y:  output signal
// yN: length of y in samples
void writeOriginalAndProcessedToWaveFile(float *x, int xN, float *y, int yN, char *fileName)
{
  // find out total number of frames in file:
  int N = max(xN, yN);;                       

  // allocate and initialize stereo buffer:
  float *xy = new float[2*N];  
  int n;
  for(n = 0; n < 2*N; n++)
    xy[n] = 0.f;

  // fill stereo buffer:
  for(int n = 0; n < xN; n++)
    xy[2*n] = x[n];
  for(int n = 0; n < yN; n++)
    xy[2*n+1] = y[n];

  // write to file and clean up:
  WavOutFile wavOutFile(fileName, 44100, 16, 2);
  wavOutFile.write(xy, 2*N);
  delete[] xy;
}


void testRandomBufferSizes(TimeStretcher *ts, float *y, int yN)
{
  int minBufferSize  = 1;
  int maxBufferSize  = 1500;
  int framesProduced = 0;
  int bufferSize;
  while( framesProduced < yN - maxBufferSize )
  {
    bufferSize = (int) random(minBufferSize, maxBufferSize);
    ts->getBuffer(&y[framesProduced], bufferSize);
    framesProduced += bufferSize;
  }
}

void testPositionJump(TimeStretcher *ts, float *y, int yN)
{
  ts->setPitchFactor(1.2f);
  ts->getBuffer(&y[0], 88200);        // retrieve 1st bar of output
  ts->jumpToInputFrame(176400);       // jump to 3rd bar
  ts->getBuffer(&y[88200], 88200);    // retrieve 3rd bar of output
  ts->jumpToInputFrame(88200);        // jump to 2nd bar
  ts->getBuffer(&y[176400], 88200);   // retrieve 2nd bar of output
  ts->jumpToInputFrame(0);            // jump to 1st bar
  ts->getBuffer(&y[264600], 88200);   // retrieve 1st bar of output
}

void testPitchJump(TimeStretcher *ts, float *y, int yN)
{
  // drag the pitch of the saw to a constant value by requesting an appropriate pitch factor
  // for each bar:
  ts->setStretchAndPitchFactor(1.f, (float) pitchOffsetToFreqFactor( 0.0));
  ts->getBuffer(&y[0], 88200);
  ts->setStretchAndPitchFactor(1.f, (float) pitchOffsetToFreqFactor(+2.0));
  ts->getBuffer(&y[88200], 88200);
  ts->setStretchAndPitchFactor(1.f, (float) pitchOffsetToFreqFactor(-3.0));
  ts->getBuffer(&y[176400], 88200);
  ts->setStretchAndPitchFactor(1.f, (float) pitchOffsetToFreqFactor(-5.0));
  ts->getBuffer(&y[264600], 88200);
}

void testStretchJump(TimeStretcher *ts, float *y, int yN)
{
  // Switch the stretch factor from 1.5 to 0.5 after half of the output has been produced. This 
  // should result in an output that stretches the 1st half of the input such that it spans three
  // quarters of the total length, and the second half of the input is cramped into the last 
  // quarter:
  int xN         = ts->getInputLength();
  int numSamples = yN;  
  int n          = 0;
  ts->setStretchAndPitchFactor(1.5f, 1.0f);
  while(numSamples > 0)
  {
    int amount = std::min(1000, numSamples);
    //int amount = std::min(500, numSamples);
    //int amount = std::min(512, numSamples);
    if( n > 0.75*xN )
      ts->setStretchAndPitchFactor(0.5f, 1.0f);
    ts->getBuffer(&y[n], amount);
    n += amount;
    numSamples -= amount;
  }
}

void testStretchAndPitchWobble(TimeStretcher *ts, float *y, int yN)
{
  int numSamples  = yN;
  float modDepth  = 5;
  float pitch     = -modDepth;
  float increment = 0.1f;
  while (numSamples > 0)
  {
    //int amount = std::min(2000, numSamples);
    //int amount = std::min(100, numSamples);
    int amount = std::min(500, numSamples);
    //int amount = std::min(512, numSamples);

    float pitchFactor   = (float) pitchOffsetToFreqFactor(pitch);
    float stretchFactor = pitchFactor;
    ts->setStretchAndPitchFactor(stretchFactor, pitchFactor);
    //stretcher.setStretchAndPitchFactor(stretchFactor, 1.0f);

    ts->getBuffer(y, amount);
    y += amount;

    pitch += increment;
    if (pitch > modDepth)
    {
      increment = - increment;
      pitch = modDepth + increment;
    }
    else if (pitch < -modDepth)
    {
      increment = - increment;
      pitch = -modDepth + increment;
    }
    numSamples -= amount;
  }
}

void testBernseePitchShift(float *x, float *y, int N, float pitchFactor)
{
  int fftBlockSize = 1024;
  int overlap      = 8;
  smbPitchShift(pitchFactor, N, fftBlockSize, overlap, 44100.f, x, y);

  // for overlap = 2, smbPitchShift introduces amplitude-modulation with a sinusoid corrsponding
  // to the blocksize, this disappears for overlap = 4
}


int main(int argc, char** argv)
{
  //float stretchFactor = 1.0f;
  float stretchFactor = 1.0f;
  float pitchFactor   = 1.5f;

  // read the input wave file inot a buffer:
  //WavInFile wavInFile("..\\Data\\Input\\Sine100Hz.wav");
  //WavInFile wavInFile("..\\Data\\Input\\Sine100HzImpulseAt25k.wav");
  WavInFile wavInFile("..\\Data\\Input\\TestLoop1.wav");
  //WavInFile wavInFile("..\\Data\\Input\\DirectCurrentHalfAmplitude.wav");
  int inputLength = wavInFile.getNumSamples();
  float *original  = new float[inputLength];
  wavInFile.read(original, inputLength);
  
  // create the processed (time-stretched) output buffer:
  int padding      = 10000; // should be >= the I/O latency of the stretcher
  int outputLength = max(inputLength, (int) ceil(stretchFactor*inputLength)) + padding;
  float *stretched = new float[outputLength];
  clearBuffer(stretched, outputLength);

  // create the TimeStretcher object (uncomment one of the lines at a time):
  //TimeStretcher *ts = new TimeStretcherElastique;
  TimeStretcher *ts = new TimeStretcherSoundTouch;
  ts->setInputStream(original, inputLength);
  ts->setStretchAndPitchFactor(stretchFactor, pitchFactor);


  // do one of the tests with the timestretcher object (uncomment one of the lines at a time):
  //testRandomBufferSizes(ts, stretched, outputLength);
  //testPositionJump(ts, stretched, outputLength);
  //testPitchJump(ts, stretched, outputLength);
  //testStretchJump(ts, stretched, outputLength);
  //testStretchAndPitchWobble(ts, stretched, outputLength);



  testBernseePitchShift(original, stretched, inputLength, 1.5);



  // write the original data into the left and the the processed data into the right channel
  // of a stereo output wavefile:
  writeOriginalAndProcessedToWaveFile(original, inputLength, stretched, outputLength, 
    "..\\Data\\Output\\TimeStretchTest.wav");

  delete[] original;
  delete[] stretched;
  delete ts;
  printf("%s", "Done. Press any key.\n");
  getchar();
}




/*
// currently not used:
void stretchTest7(float *x, int xN, float *y, int yN, float stretchFactor)
{
  // start processing from differnt points in the data - check from which sample onwards the 
  // outputs agree.

  float *tmp1 = new float[yN];
  float *tmp2 = new float[yN];
  clearBuffer(tmp1, yN);
  clearBuffer(tmp2, yN);
  int startOffset = 10000;
  //int startOffset = 9702;  // a zero crossing of the 100 Hz sine

  TimeStretcherElastique stretcher1(512);
  stretcher1.setInputStream(x, xN);
  stretcher1.getBuffer(tmp1, yN);

  TimeStretcherElastique stretcher2(512);
  stretcher2.setInputStream(&x[startOffset], xN-startOffset);
  stretcher2.getBuffer(&tmp2[startOffset], yN-startOffset);

  writeOriginalAndProcessedToWaveFile(tmp1, yN, tmp2, yN, 
    "..\\Data\\Output\\StartOffsetTest.wav");

  delete[] tmp1;
  delete[] tmp2;
}
*/

