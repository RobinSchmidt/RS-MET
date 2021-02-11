#ifndef rosic_FileInputOutput_h
#define rosic_FileInputOutput_h

namespace rosic
{

/** This file defines functions for handling audio file I/O. */

/** Writes a mono signal into a wavefile (single precision version). */
void writeToMonoWaveFile(const char* path, float* signal, int numFrames, int sampleRate,
  int numBits = 16);

/** Writes a mono signal into a wavefile (double precision version). */
void writeToMonoWaveFile(const char* path, double* signal, int numFrames, int sampleRate,
  int numBits = 16);

// convenience function that uses std::string instead of char*
void writeToMonoWaveFile(const std::string& path, double* signal, int numFrames, int sampleRate,
  int numBits = 16);


/** Writes a stereo signal into a wavefile. */
void writeToStereoWaveFile(const char* path, double* left, double* right, int numFrames,
  int sampleRate, int numBits = 16);

/** Reads audio-data from a wavefile and returns a pointer to the data. When you access the data
as an array, the first index is channel-number and the second index is the sample-number. Note
that the returned pointer may be NULL, if loading fails. If it is not NULL, it is the callers
responsibility to eventually free the allcoated memory (inner and outer pointers, via delete[]).
The reference parameters will be assigned to their correct values inside this function, so you
know what you've got when the function returns. */
double** readFromWaveFile(const char* path, int& numChannels, int& numFrames, int& sampleRate);

/** Writes the passed string into the file given by 'path' */
void writeStringToFile(const char* path, const char* stringToWrite);

/** Writes the values stored in x and y into a whitespace delimited ascii data file. */
void writeDataToFile(const char* path, int numValues, double* x, double* y1, double* y2 = NULL,
  double* y3 = NULL, double* y4 = NULL, double* y5 = NULL);

//void writeDataToFile(char* path, double *data, int pixelWidth, int pixelHeight);

void writeDataToFile(const char* path, double* data, int pixelWidth, int pixelHeight,
  double* xAxis, double* yAxis);



} // end namespace rosic

#endif // rosic_FileInputOutput_h
