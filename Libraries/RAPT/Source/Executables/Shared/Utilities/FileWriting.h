#ifndef RAPT_FILEWRITING_H
#define RAPT_FILEWRITING_H

#include "../../../Modules/RAPT.h"
using namespace RAPT;

/** Writes the passed monochrome image into a .ppm file */
void writeImageToFilePPM(const Image<float>& image, const char* path);



#endif