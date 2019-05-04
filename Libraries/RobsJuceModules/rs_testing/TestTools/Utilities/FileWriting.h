#ifndef RAPT_FILEWRITING_H
#define RAPT_FILEWRITING_H

//#include "../RaptLibraryCode/RaptInstantiations.h"
#include "rosic/rosic.h"

/** Writes the passed monochrome image into a .ppm file */
void writeImageToFilePPM(const RAPT::rsImage<float>& image, const char* path);
// make version that takes an rsImageF4 (4 floats per pixel representing RGBA)

void writeScaledImageToFilePPM(RAPT::rsImage<float>& image, const char* path, int scale);

void writeToMonoWaveFile(std::string path, float *signal, int numFrames, int sampleRate,
  int numBits);



#endif