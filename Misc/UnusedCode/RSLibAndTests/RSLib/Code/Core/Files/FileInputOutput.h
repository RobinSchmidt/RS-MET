#ifndef RS_FILEINPUTOUTPUT_H
#define RS_FILEINPUTOUTPUT_H

namespace RSLib
{

  /**

  This file contains convenience functions for simplified handling of file I/O without having to 
  instantiate an rsFile object.

  \todo: support 32 bit float wavefiles, 24 bit fixed, maybe even 64 bit float

  // see:
  http://soundfile.sapp.org/doc/WaveFormat/
  http://www-mmsp.ece.mcgill.ca/documents/AudioFormats/WAVE/WAVE.html
  http://wavefilegem.com/how_wave_files_work.html

  */

  /** Writes a mono signal into a wavefile (single precision version). */
  RSLib_API void writeToMonoWaveFile(const char* path, float *signal, int numFrames, 
    int sampleRate, int numBits);

  /** Writes a mono signal into a wavefile (double precision version). */
  RSLib_API void writeToMonoWaveFile(const char* path, double *signal, int numFrames, 
    int sampleRate, int numBits);

  /** Writes a stereo signal into a wavefile. */
  RSLib_API void writeToStereoWaveFile(const char* path, double *left, double *right, 
    int numFrames, int sampleRate, int numBits);
  
  // todo: templatize the function that have double and single precision versions
  // use an rsString as 1st parameter, prepend rs

  /** Reads audio-data from a wavefile and returns a pointer to the data. When you access the data 
  as an array, the first index is channel-number and the second index is the sample-number. Note 
  that the returned pointer may be NULL, if loading fails. If it is not NULL, it is the callers 
  responsibility to eventually free the allcoated memory (inner and outer pointers, via delete[]). 
  The reference parameters will be assigned to their correct values inside this function, so you 
  know what you've got when the function returns. The caller is responsible for freeing the memory
  associated with the pointer. */
  RSLib_API double** readFromWaveFile(const char* path, int& numChannels, int& numFrames, 
    int& sampleRate);

  /** Reads just the 1st channel of a wavefile. */
  RSLib_API double* readMonoWaveFile(const char* path, int& numFrames, int& sampleRate);


  /** Writes the passed string into the file given by 'path' */
  RSLib_API void writeStringToFile(const char* path, char* stringToWrite);

  /** Writes the values stored in x and y into a whitespace delimited ascii data file. */
  RSLib_API void writeDataToFile(const char* path, int numValues, double *x, double *y1, 
    double *y2 = NULL, double *y3 = NULL, double *y4 = NULL, double *y5 = NULL);

  //void writeDataToFile(char* path, double *data, int pixelWidth, int pixelHeight);

  RSLib_API void writeDataToFile(const char* path, double *data, int pixelWidth, int pixelHeight, 
    double *xAxis, double *yAxis);

}

#endif
