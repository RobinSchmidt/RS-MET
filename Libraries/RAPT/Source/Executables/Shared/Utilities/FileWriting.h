#ifndef RAPT_FILEWRITING_H
#define RAPT_FILEWRITING_H

#include "../RaptLibraryCode/RaptTypedefInstantiations.h"

//#include "../../../Modules/RAPT.h"
//using namespace RAPT;

/** Writes the passed monochrome image into a .ppm file */
void writeImageToFilePPM(const ImageF& image, const char* path);



#endif