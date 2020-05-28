#ifndef RAPT_FILEWRITING_H  // this is not a rapt class!!
#define RAPT_FILEWRITING_H

//#include "../RaptLibraryCode/RaptInstantiations.h"
#include "rosic/rosic.h"



bool writeImageToFilePPM(const char* path, unsigned char* buf, int w, int h);



bool writeImageToFilePPM(const RAPT::rsImage<rsPixelRGB>& image, const char* path);

/** Writes the passed monochrome image into a .ppm file */
bool writeImageToFilePPM(const RAPT::rsImage<float>& image, const char* path);
// make version that takes an rsImageF4 (4 floats per pixel representing RGBA)

/** Uses the 3 passed images (assumed to be of the same dimensions) as red, green and blue 
channel. */
bool writeImageToFilePPM(const RAPT::rsImage<float>& red, const RAPT::rsImage<float>& green,
  const RAPT::rsImage<float>& blue, const char* path);

/** Writes a naively (by pixel duplication) scaled up version of the given image into a file. This 
is useful for taking a close look at what the rendering algorithms do. */
bool writeScaledImageToFilePPM(RAPT::rsImage<float>& image, const char* path, int scale);




void writeToMonoWaveFile(std::string path, float *signal, int numFrames, int sampleRate,
  int numBits);

/** Writes a complex signal to a stereo wavefile using the real part as left channel and the 
imaginary part as right channel. */
void writeToWaveFile(std::string path, const std::vector<std::complex<double>>& signal, 
  int sampleRate = 44100);
// maybe rename to rsWavWrite



#endif