#include "SampleMapGenerator.h"

// a quick and dirty implementation to create a directory which may not yet exist (but its parent
// is assumed to exist) -  move to somewhere in RSLib, when when finished and cleaned:

#ifdef _MSC_VER
#include <direct.h>  // for msvc - otherwise maybe dirent.h
#else
#include <dirent.h>
//#include <sys/types.h>
//#include <sys/stat.h>
#endif

void rsCreateDirectoryIfNonExistent(const RSLib::rsString &path)
{
  char *cString = path.getAsZeroTerminatedString();

#ifdef _MSC_VER
  if( _mkdir( cString ) != 0 )
#else
  //if( mkdir(cString, 0777) != 0 )
  if( mkdir(cString) != 0 )
#endif
  {
    if( errno == EEXIST )
    {
      // directory exists already - not an error
    }
    else if( errno == ENOENT )
    {
      RSLib::rsError("Parent directory does not exist.");
        // _mkdir can only create one dircetory at a time, so if we want to use it to create a
        // sub-sub-sub... directory, we need to call this function recursively
    }
    else
    {
      RSLib::rsError( "Problem creating directory\n" );
        // errno may indicate some other condition like missing write access and stuff - complete
        // the implementation of that later
    }
  }

  delete[] cString;
  int dummy = 0;
}

// move to RSLib, write UnitTests:
void rsWriteToWaveFile(const RSLib::rsString &path, rsAudioBuffer &buffer, int sampleRate,
                       int numBits)
{
  // create an interleaved buffer of floats for writing:
  float *tmp = new float[buffer.getSize()];
  rsConvertBuffer(buffer.dataFlat, tmp, buffer.getSize());
  rsInterleave(tmp, buffer.numFrames, buffer.numChannels);

  // write to file and cleanup:
  rsOutputWaveFile file(path, sampleRate, numBits, buffer.numChannels);
  file.write(tmp, buffer.getSize());

  delete[] tmp;
}

// Construction/destruction:

SampleMapGenerator::SampleMapGenerator()
{
  outputDirectory = rsGetCurrentApplicationDirectory();
  setName("UnnamedSampleMap");

  for(int i = 0; i < 128; i++)
    gains[i] = ambGains[i] = 1.0;

  sampleRate = 44100.0;
  numBits    = 16;
  loKey      = 0;
  hiKey      = 127;
  ambience   = false;

  ft.setBlockSize(fftSize);
}

SampleMapGenerator::~SampleMapGenerator()
{

}

// Setup:

void SampleMapGenerator::setOutputDirectory(const RSLib::rsString& newDirectory)
{
  outputDirectory = newDirectory;
}

void SampleMapGenerator::setName(const RSLib::rsString& newName)
{
  name = newName;
  sampleSubDirectory = name + "\\";
  //sampleSubDirectory = "Samples" + name + "\\";
}

void SampleMapGenerator::setKeyRangeToRender(int newLowKey, int newHighKey)
{
  loKey = newLowKey;
  hiKey = newHighKey;
}

// Inquiry:

RSLib::rsString SampleMapGenerator::getSampleFileName(int key)
{
  key = rsClipToRange(key, loKey, hiKey);
  return name + "K"
    + RSLib::rsString(key)    // use toStringWithLeadingZeros(key, 3)
    + ".wav";
}

RSLib::rsString SampleMapGenerator::getSampleRelativePath(int key)
{
  return sampleSubDirectory + getSampleFileName(key);
}

RSLib::rsString SampleMapGenerator::getAmbienceRelativePath(int key)
{
  return sampleSubDirectory + "Amb" + getSampleFileName(key);
}

// Sample-map generation:

void SampleMapGenerator::generateSampleMap(bool printProgress)
{
  RSLib::rsString sampleDir = outputDirectory + sampleSubDirectory;
  rsCreateDirectoryIfNonExistent(sampleDir);

  generateAllSamples(printProgress); // temporarily commented

  // normalize gain factors such that the maximum gain factor is 0.5 (arbitrary - maybe choose
  // a better value later - match straightliner's presets gains)
  rsNormalize(gains, 128, 0.5);

  generateSoundFontFile();
}

RSLib::rsString SampleMapGenerator::getRegionString(int root, int lo, int hi)
{
  RSLib::rsString s;
  s += "<region>";
  s += " sample="          + getSampleRelativePath(root);
  s += " pitch_keycenter=" + RSLib::rsString(root);
  s += " lokey="           + RSLib::rsString(lo);
  s += " hikey="           + RSLib::rsString(hi);
  s += " volume="          + rsDoubleToString2(rsAmp2dB(gains[root])) + "\n";
  return s;
}

void SampleMapGenerator::generateSoundFontFile()
{
  RSLib::rsString sfz;

  // some global performance settings for all regions to make it more playable:
  sfz += "<group>";
  sfz += " amp_veltrack=90";
  sfz += " ampeg_release=0.5";
  sfz += " ampeg_attack=0.01";
  sfz += " ampeg_vel2attack=-0.01";
  sfz += " fil_type=lpf_1p";
  sfz += " cutoff=1200";
  sfz += " fil_keytrack=100";
  sfz += " fil_veltrack=9000";
  sfz += " bend_up=1200";
  sfz += " bend_down=-1200";
  //sfz += " pitchlfo_freq=7";  // this tiny unnoticable amount of vibrato is a workaround for
  //sfz += " pitchlfo_depth=1"; // sfz+'s failure to playback a sample as-is without aliasing
  sfz += "\n";

  // each rendered sample is used for a one-key wide region except for the lowest and highest
  // rendered sample - these have regions that extend to the left and to right of the keyboard
  // respectively:
  sfz += getRegionString(loKey, 0, loKey);
  for(int key = loKey+1; key < hiKey; key++)
    sfz += getRegionString(key, key, key);
  sfz += getRegionString(hiKey, hiKey, 127);

  // add the ambience group, if desired (move into function):
  if( ambience == true )
  {
    sfz += "\n";
    sfz += "<group>";
    sfz += " amp_veltrack=90";
    sfz += " ampeg_attack=0.05";
    sfz += " ampeg_decay=3.0";   // maybe make this dependent on the signal
    sfz += " ampeg_sustain=0.0";
    sfz += " ampeg_release=3.0";
    sfz += " fil_type=lpf_1p";
    sfz += " cutoff=600";
    sfz += " fil_keytrack=100";
    sfz += " fil_veltrack=3000 ";
    //sfz += " pitchlfo_freq=7";
    //sfz += " pitchlfo_depth=1";
    sfz += "\n";

    for(int key = loKey; key <= hiKey; key++)
    {
      if( key % 12 == 9 )
      {
        // \todo: make sure to extend the ambience region to the whole keyboard

        sfz += "<region>";
        sfz += " sample="          + getAmbienceRelativePath(key);
        sfz += " pitch_keycenter=" + RSLib::rsString(key);
        sfz += " lokey="           + RSLib::rsString(key);
        sfz += " hikey="           + RSLib::rsString(key+11);
        sfz += " volume="          + rsDoubleToString2(rsAmp2dB(ambGains[key])-12);
        sfz += " pan=-100";
        sfz += " loop_mode=loop_continuous";
        sfz += "\n";

        sfz += "<region>";
        sfz += " sample="          + getAmbienceRelativePath(key);
        sfz += " pitch_keycenter=" + RSLib::rsString(key);
        sfz += " lokey="           + RSLib::rsString(key);
        sfz += " hikey="           + RSLib::rsString(key+11);
        sfz += " volume="          + rsDoubleToString2(rsAmp2dB(ambGains[key])-12);
        sfz += " pan=100";
        sfz += " loop_mode=loop_continuous";
        sfz += " offset="          + RSLib::rsString(fftSize/2);
        sfz += "\n";

        // get rid of the code-duplication
      }
    }
  }

  // write the string to an .sfz file:
  RSLib::rsString sfzPath = outputDirectory + name + ".sfz";
  RSLib::rsFile sfzFile(sfzPath);
  sfzFile.writeStringToFile(sfz);
}

void SampleMapGenerator::generateAllSamples(bool printProgress)
{
  int numKeys = hiKey - loKey + 1;
  for(int k = loKey; k <= hiKey; k++)
  {
    // todo: instead of k++, use k+=keyIncrement which can be set by the user

    if( printProgress == true )
      printf("%s %d %s %d %s", "Generating sample ", k-loKey+1, " of ", numKeys, "\n");
    generateSampleForKey(k);
  }
  if( printProgress == true )
    printf("%s", "All samples done\n");
}

void SampleMapGenerator::generateSampleForKey(int key)
{
  writeUnnormalizedSampleForKeyToBuffer(key);

  // normalize buffer and store normalization gain:
  double max = rsMaxAbs(buffer.dataFlat, buffer.getSize());
  gains[key] = max;
  rsScale(buffer.dataFlat, buffer.getSize(), 1.0/max);

  // write buffer to wave-file:
  RSLib::rsString path = outputDirectory + getSampleRelativePath(key);
  rsWriteToWaveFile(path, buffer, (int) sampleRate, numBits);



  // move into function:
  if( ambience == true )
  {

    if( key % 12 == 9 ) // generate ambience-sample for every A
    {
      double *mag = new double[fftSize/2]; // FFT magnitudes
      double *phs = new double[fftSize/2]; // FFT phases

      // copy signal into temporary FFT buffer:
      rsAudioBuffer tmpBuf;
      tmpBuf.setSize(1, fftSize);
      int n;
      if( buffer.numChannels == 2 )
      {
        for(n = 0; n < RAPT::rsMin(fftSize, buffer.numFrames); n++)
          tmpBuf.data[0][n] = SQRT2_INV * (buffer.data[0][n] + buffer.data[1][n]);
      }
      else
      {
        for(n = 0; n < RAPT::rsMin(fftSize, buffer.numFrames); n++)
          tmpBuf.data[0][n] = buffer.data[0][n];
      }
      for(n = RAPT::rsMin(fftSize, buffer.numFrames); n < fftSize; n++)
          tmpBuf.data[0][n] = 0.0;
      if( fftSize < buffer.numFrames )
      {
        // fade out instead of truncating:
        int numFadeSamples = 1024;
        RAPT::rsFadeOut(tmpBuf.data[0], tmpBuf.numFrames-numFadeSamples-1, tmpBuf.numFrames-1);
          // check, if this is right
      }


      // transform, randomize phases and transform back:
      ft.getRealSignalMagnitudesAndPhases(tmpBuf.data[0], mag, phs);


      /*
      // todo: maybe apply a smoothing filter to the magnitude spectrum
      // moving average of length L:
      int L   = 10;
      double scl = 1.0 / (double) L;
      for(int k = L+1; k < fftSize/2-L; k++)
      {
        double accu = mag[k];
        for(int i = 0; i < L; i++)
          accu += mag[k-i] + mag[k+i];
        mag[k] = scl * accu;
      }
      */

      RAPT::rsRandomUniform(0.0, 1.0, 0);
      for(int n = 1; n < fftSize/2; n++)  // phs[0] has special meaning - start at 1
        phs[n] = RAPT::rsRandomUniform(0.0, 2*PI);
      ft.getRealSignalFromMagnitudesAndPhases(mag, phs, tmpBuf.data[0]);

      // normalize ambience buffer and store normalization gain:
      max = rsMaxAbs(tmpBuf.dataFlat, tmpBuf.getSize());
      ambGains[key] = max;
      rsScale(tmpBuf.dataFlat, tmpBuf.getSize(), 1.0/max);

      // write ambience buffer to wave-file:
      RSLib::rsString path = outputDirectory + getAmbienceRelativePath(key);
      rsWriteToWaveFile(path, tmpBuf, (int) sampleRate, numBits);

      delete[] mag;
      delete[] phs;
    }

    // \todo: this ambience signal will be used as is for the left channel and with
    // offset of length/2 for the right channel (it should be looped and the decay/release should
    // be set approriately (it should probably depend on the release-time of the sample)
  }
}

//=================================================================================================

SampleMapGeneratorModal::SampleMapGeneratorModal()
{
  truncationLevel = -80.0;
  for(int i = 0; i < 128; i++)
    keyParameters[i] = NULL;
}

SampleMapGeneratorModal::~SampleMapGeneratorModal()
{
  for(int i = 0; i < 128; i++)
    delete keyParameters[i];
}


void SampleMapGeneratorModal::setModalParametersForKey(int key,
  const rsModalBankParametersD& newParameters)
{
  delete keyParameters[key];
  keyParameters[key] = new rsModalBankParametersD(newParameters);
}


void SampleMapGeneratorModal::writeUnnormalizedSampleForKeyToBuffer(int key)
{
  rsModalBankParametersD p;

  if( keyParameters[key] != NULL )
    p = *(keyParameters[key]);      // use parameters from array as-is
  else
  {
    // find neighbours where parameters are defined:
    int keyL = key;
    int keyR = key;
    while( keyParameters[keyL] == NULL && keyL >= 0 )
      keyL--;
    while( keyParameters[keyR] == NULL && keyR <= 127 )
      keyR++;

    // inter- or extrapolate:
    if( keyL == -1 && keyR == 128 )
    {
      RAPT::rsError("No parameters for inter- or extrapolation available.");
    }
    if( keyL == -1 )
      p = *(keyParameters[keyR]);  // extrapolate leftward
    else if( keyR == 128 )
      p = *(keyParameters[keyL]);  // extrapolate rightward
    else
    {
      // interpolate:
      double delta = double (key-keyL) / double (keyR-keyL);
      p = interpolateParameters(*keyParameters[keyL], *keyParameters[keyR], delta);
    }
  }

  // set  fundamental frequency and truncate all arrays in p to below sampleRate/2:
  p.frequency = rsPitchToFreq((double)key);
  p = removeModesAbove(p, sampleRate/2);



  double fadeCycles = 5.0; // # pitch-cycles for fade out

  // setup the modal synthesizer:
  modalFilterBank.setSampleRate(sampleRate);
  modalFilterBank.setModalParameters(p.f, p.g, p.a, p.d, p.p);
  modalFilterBank.setReferenceFrequency(p.frequency);
  modalFilterBank.setReferenceAttack(p.attack);
  modalFilterBank.setReferenceDecay(p.decay);

  // compute number of samples to generate and fade-out time:
  int    numSamples     = (int) (modalFilterBank.getLength(truncationLevel)*sampleRate);
  double fadeOutTime    = fadeCycles / p.frequency;
  int    numFadeSamples = (int) (sampleRate * fadeOutTime);
  numSamples += numFadeSamples;

  // synthesize (and post-process) the sound and write the result into the buffer:
  modalFilterBank.reset();
  buffer.setSize(1, numSamples);
  buffer.data[0][0] = modalFilterBank.getSample(1.0);
  for(int n = 1; n < numSamples; n++)
    buffer.data[0][n] = modalFilterBank.getSample(0.0);
  rsFadeOut(buffer.data[0], numSamples-numFadeSamples-1, numSamples-1);




  // under construction - for stereoization (so far, this is all shit):

  /*
  rsRandomUniform(0.0, 1.0, 1);
  for(int m = 0; m < p.f.dim; m++)
  {
    double amount = 1.0;
    double rnd    = rsRandomUniform(-1.0, 1.0);
    //p.p[m] = p.f[m] + PI/2 + amount * rnd * (double) m / (p.f.dim-1);
    p.p[m] = p.f[m] + PI/2 + amount * rnd * PI/2;
  }
  */

  //rsRandomUniform(0.0, 1.0, 1);

  /*
  double fundamentalShiftSt = 0.5;  // shift of the fundametal in the side-signal in semitones
  double fundamentalDetune  = rsPitchOffsetToFreqFactor(fundamentalShiftSt);
  double fundamentalShiftHz = (fundamentalDetune-1) * p.frequency;
  double f0 = p.frequency;
  double df = fundamentalShiftHz;
  rsRandomUniform(0.0, 1.0, 1);

  //for(int m = 0; m < p.f.dim-1; m+=2)
  //{
  //  double rnd = rsRandomUniform(0.0, 1.0);
  //  p.f[m]   = p.f[m]   * (f0*p.f[m]   + rnd * df) / (f0*p.f[m]  );
  //  p.f[m+1] = p.f[m+1] * (f0*p.f[m+1] - rnd * df) / (f0*p.f[m+1]);
  //}
  for(int m = 0; m < p.f.dim; m++)
  {
    double amount = 0.1;
    double rnd    = rsRandomUniform(-1.0, 1.0);
    p.f[m] = p.f[m] + amount * rnd * (double) m / (p.f.dim-1);
  }

  // chek, if normalized average mode-distance is unity:
  static const int N = 101;
  double dbg[N];
  rsCopyBuffer(p.f.v, dbg, N);
  double tmp = 0;
  for(int i = 0; i < N-1; i++)
  {
    double d = dbg[i+1] - dbg[i];
    tmp += d;
  }
  */

  /*
  // render side-signal with modified modal parameters:
  modalFilterBank.setModalParameters(p.f, p.g, p.a, p.d, p.p);
  modalFilterBank.resetModalFilters();
  buffer.data[1][0] = modalFilterBank.getSample(1.0);
  for(int n = 1; n < numSamples; n++)
    buffer.data[1][n] = modalFilterBank.getSample(0.0);
  rsFadeOut(buffer.data[1], numSamples-numFadeSamples-1, numSamples-1);

  // convert from M/S to L/R:
  for(int n = 0; n < numSamples; n++)
  {
    double M = buffer.data[0][n];
    double S = buffer.data[1][n];
    double L = ONE_OVER_SQRT2 * (M+S);
    double R = ONE_OVER_SQRT2 * (M-S);
    buffer.data[0][n] = L;
    buffer.data[1][n] = R;
  }
  */

}

// move to RSLib:
double rsWeightedGeometricMean(double y1, double y2, double w1)
{
  double w2 = 1-w1;
  return pow(y1, w1) * pow(y2, w2);

}

rsModalBankParametersD SampleMapGeneratorModal::interpolateParameters(
  const rsModalBankParametersD &p1, const rsModalBankParametersD &p2, double proportion)
{
  rsModalBankParametersD r = p1;  // result, assigment just to reserve memory
  double w1 = 1-proportion;      // weight for parameter set 1

  r.frequency = rsWeightedGeometricMean(p1.frequency, p2.frequency, w1);
  r.gain      = rsWeightedGeometricMean(p1.gain,      p2.gain,      w1);
  r.attack    = rsWeightedGeometricMean(p1.attack,    p2.attack,    w1);
  r.decay     = rsWeightedGeometricMean(p1.decay,     p2.decay,     w1);
  for(size_t i = 0; i < r.g.size(); i++)
  {
    r.f[i] = rsWeightedGeometricMean(p1.f[i], p2.f[i], w1);
    r.g[i] = rsWeightedGeometricMean(p1.g[i], p2.g[i], w1);
    r.a[i] = rsWeightedGeometricMean(p1.a[i], p2.a[i], w1);
    r.d[i] = rsWeightedGeometricMean(p1.d[i], p2.d[i], w1);

    r.p[i] = w1*p1.p[i] + (1-w1)*p2.p[i];
     // weighted arithmetic mean - todo: use wrapped-interpolation (see sine-modeling code)
  }

  return r;
}

rsVectorDbl truncateVector(const rsVectorDbl v, int numElementsToRetain)
{
  return rsVectorDbl(numElementsToRetain, v.v);
}

rsModalBankParametersD SampleMapGeneratorModal::removeModesAbove(
  const rsModalBankParametersD &p, double cutoff)
{
  rsModalBankParametersD r = p;

  int cutoffIndex = 0;
  for(size_t i = 0; i < p.f.size(); i++)
  {
    if( p.frequency * p.f[i] < cutoff )
      cutoffIndex++;
  }

  /*
  // truncation via constructor of rsVectorDbl:
  r.f = rsVectorDbl(cutoffIndex, r.f.v);
  r.g = rsVectorDbl(cutoffIndex, r.g.v);
  r.a = rsVectorDbl(cutoffIndex, r.a.v);
  r.d = rsVectorDbl(cutoffIndex, r.d.v);
  r.p = rsVectorDbl(cutoffIndex, r.p.v);
  */

  // new - should also truncate (i hope)
  r.f.resize(cutoffIndex);
  r.g.resize(cutoffIndex);
  r.a.resize(cutoffIndex);
  r.d.resize(cutoffIndex);
  r.p.resize(cutoffIndex);
  //RAPT::rsAssertFalse; // check, if code above truncates vectors as desired - yes, works

  return r;
}
