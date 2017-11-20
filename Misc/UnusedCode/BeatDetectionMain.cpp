#include "../rosic/infrastructure/rosic_AudioFileHandling.h"
#include "../rosic/analysis/rosic_BeatDetector.h"
using namespace rosic;

int main(int argc, char* argv[])
{
  // read in the test-signal and convert it to single-precision, mono:
  int numFrames;
  int numChannels;
  int sampleRate;
  double** tmp = readFromWaveFile("D:/TmpAudio/ElektroBeatDry140BPM.wav", numChannels, numFrames, sampleRate);
  //double** tmp = readFromWaveFile("D:/TmpAudio/Drumloop.wav", numChannels, numFrames, sampleRate);
  //double** tmp = readFromWaveFile("D:/TmpAudio/RebelYellChunk2.wav", numChannels, numFrames, sampleRate);
  //double** tmp = readFromWaveFile("D:/TmpAudio/BinaryFinaryChunk.wav", numChannels, numFrames, sampleRate);
  //double** tmp = readFromWaveFile("D:/TmpAudio/test.wav", numChannels, numFrames, sampleRate);

  float* x = NULL;  // the test-signal in /mono format
  int n;
  if( tmp != NULL )
  {
    x = new float[numFrames];
    if( numChannels == 2 )
    {
      for(n=0; n<numFrames; n++)
        x[n] = (float) (0.5 * (tmp[0][n] + tmp[1][n]));
    }
    else
    {
      for(n=0; n<numFrames; n++)
        x[n] = (float) tmp[0][n];
    }
    delete[] tmp[0];
    delete[] tmp;
  }

  // x now contains the signal to be analyzed
  BeatDetector detector;
  detector.processSignalAtOnce(x, numFrames, sampleRate);
  detector.processOnsets();
  std::vector<Onset> onsets = detector.getOnsets();

  // for comparison, process the same signal blockwise:
  BeatDetector detector2;
  detector2.prepareForBlockProcessing(numFrames, sampleRate);
  int currentBlockStart   = 0;
  while( currentBlockStart < numFrames )
  {
    int remainingSamples = numFrames - currentBlockStart;
    int currentBlockSize = (int) randomUniform(1, 5000);
    currentBlockSize = rmin(currentBlockSize, remainingSamples);

    detector2.feedSignalBlock(&x[currentBlockStart], currentBlockSize);

    currentBlockStart += currentBlockSize;
    int dummy = 0;
  }
  detector2.finishBlockProcessing();
  detector2.processOnsets();
  std::vector<Onset> onsets2 = detector.getOnsets();

  // attenuate signal, insert clicks at the onsets and write the result out to a file::    
  for(n=0; n<numFrames; n++)
    x[n] *= 0.25f;

  for(int i=0; i < (int) onsets.size(); i++)
  {
    if( onsets[i].isBeat )
    {
      n      = onsets[i].timeInSamples;
      x[n]   =  1.f;
      x[n+1] = -1.f;
    }
  }
  writeToMonoWaveFile("D:/TmpAudio/Beats.wav", x, numFrames, sampleRate, 16);

  for(int i=0; i < (int) onsets.size(); i++)
  {
    n      = onsets[i].timeInSamples;
    x[n]   =  1.f;
    x[n+1] = -1.f;
  }
  writeToMonoWaveFile("D:/TmpAudio/Onsets.wav", x, numFrames, sampleRate, 16);

  if( x!= NULL )
    delete x;

	return 0;
}

