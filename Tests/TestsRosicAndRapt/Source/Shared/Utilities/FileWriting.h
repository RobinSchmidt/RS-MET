#ifndef RAPT_FILEWRITING_H
#define RAPT_FILEWRITING_H

#include "../RaptLibraryCode/RaptInstantiations.h"
#include "rosic/rosic.h"

/** Writes the passed monochrome image into a .ppm file */
void writeImageToFilePPM(const rsImageF& image, const char* path);


void writeToMonoWaveFile(std::string path, float *signal, int numFrames, int sampleRate,
  int numBits);



#endif