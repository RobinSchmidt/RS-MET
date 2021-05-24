//#include "rosic_FileInputOutput.h"
//using namespace rosic;

#include "../_third_party/soundtouch/WavFile.cpp"

void rosic::writeToMonoWaveFile(const char* path, float *signal, int numFrames, int sampleRate,
                                int numBits)
{
  WavOutFile file(path, sampleRate, numBits, 1);
  file.write(signal, numFrames);
}

void rosic::writeToMonoWaveFile(const char* path, double *signal, int numFrames, int sampleRate,
                                int numBits)
{
  WavOutFile file(path, sampleRate, numBits, 1);
  float *tmp = new float[numFrames];
  for(int n=0; n<numFrames; n++)
    tmp[n] = (float) signal[n];
  file.write(tmp, numFrames);
  delete[] tmp;
}

void rosic::writeToMonoWaveFile(const std::string& path, double* signal, int numFrames,
  int sampleRate, int numBits)
{
  writeToMonoWaveFile(path.c_str(), signal, numFrames, sampleRate, numBits);
}

void rosic::writeToStereoWaveFile(
  const char* path, double *left, double *right, int numFrames, int sampleRate, int numBits)
{
  WavOutFile file(path, sampleRate, numBits, 2);
  float *tmp = new float[2*numFrames];
  for(int n=0; n<numFrames; n++)
  {
    tmp[2*n]   = (float) left[n];
    tmp[2*n+1] = (float) right[n];
  }
  file.write(tmp, 2*numFrames);
  delete[] tmp;
}
void rosic::writeToStereoWaveFile(
  const char* path, float* left, float* right, int numFrames, int sampleRate, int numBits)
{
  WavOutFile file(path, sampleRate, numBits, 2);
  float *tmp = new float[2*numFrames];
  for(int n=0; n<numFrames; n++)
  {
    tmp[2*n]   = left[n];
    tmp[2*n+1] = right[n];
  }
  file.write(tmp, 2*numFrames);
  delete[] tmp;
}
// code is almost the same as for the double-precision version - maybe templatize and instantiate
// for float and double to get rid of duplication

float** rosic::readFloatFromWaveFile(
  const char* path, int& numChannels, int& numFrames, int& sampleRate)
{
  try
  {
    WavInFile file(path);
    numChannels = file.getNumChannels();
    numFrames   = file.getNumSamples();
    sampleRate  = file.getSampleRate();
    if( numChannels <= 0 || numFrames <= 0 )
      return nullptr;

    float *fBuffer = new float[numChannels*numFrames];
    file.read(fBuffer, numChannels*numFrames);
    if( numChannels > 1 )
      RAPT::rsArrayTools::deInterleave(fBuffer, numFrames, numChannels);

    float **pointers = new float*[numChannels];
    for(int c = 0; c < numChannels; c++)
      pointers[c] = &fBuffer[c*numFrames];

    return pointers;
  }
  catch(...)
  {
    numChannels = numFrames = sampleRate = 0; return nullptr;
  }

  // ToDo: 
  // -check, if the 1st or 2nd call to "new" return a nullptr (indicating allocation failure)
  //  -if the 1st: return nullptr
  //  -if the 2nd: clean up 1st and return nullptr
}

double** rosic::readFromWaveFile(
  const char* path, int& numChannels, int& numFrames, int& sampleRate)
{
  // todo: use the float version, avoid code duplication
  try
  {
    WavInFile file(path);
    numChannels = file.getNumChannels();
    numFrames   = file.getNumSamples();
    sampleRate  = file.getSampleRate();
    if( numChannels <= 0 || numFrames <= 0 )
      return NULL;

    float *fBuffer = new float[numChannels*numFrames];
    file.read(fBuffer, numChannels*numFrames);
    if( numChannels > 1 )
      RAPT::rsArrayTools::deInterleave(fBuffer, numFrames, numChannels);

    double *dBuffer = new double[numChannels*numFrames];
    for(int n = 0; n < numChannels*numFrames; n++)
      dBuffer[n] = (double) fBuffer[n];
    double **pointers = new double*[numChannels];
    for(int c = 0; c < numChannels; c++)
      pointers[c] = &dBuffer[c*numFrames];

    delete[] fBuffer;
    return pointers;
  }
  catch(...)
  {
    numChannels = 0;
    numFrames   = 0;
    sampleRate  = 0;
    return NULL;
  }
}

void rosic::writeStringToFile(const char* path, const char* stringToWrite)
{
  int length = (int) strlen(stringToWrite);
  FILE *f = fopen(path, "w");
  if( f != NULL )
  {
    fwrite(stringToWrite, sizeof(char), length, f);
    fclose(f);
  }
}

bool rosic::rsWriteStringToFile(const char* path, const char* str)
{
  // ToDo: provide optional argument for the length which defaults to 0 in which case we use strlen
  // and if it's not 0, we use the passed value and avoid strlen, but maybe we should do
  // rsAssert(strlen(str) == length).

  FILE* f = fopen(path, "w");
  if(f != NULL){
    fwrite(str, 1, strlen(str), f);
    fclose(f); 
    return true; }
  else {
    RAPT::rsError("Unable to open file");
    return false; }

  // https://www.tutorialspoint.com/cprogramming/c_file_io.htm
}

char* rosic::rsReadStringFromFile(const char *filename)
{
  // Code adapted from here:
  // https://stackoverflow.com/questions/3463426/in-c-how-should-i-read-a-text-file-and-print-all-strings
  // Note that it is pointed out there that the technique of inquiring the length via fseek, ftell 
  // may not work for binary files. Here, we are dealing with text files, so it should hopefully 
  // be ok. The unit tests pass, so far.

  char *buffer = NULL;
  size_t string_size, read_size;
  FILE *handler = fopen(filename, "r");
  if(handler)
  {
    fseek(handler, 0, SEEK_END);  // seek the last byte of the file
    string_size = ftell(handler); // offset from the first to the last byte, filesize
    rewind(handler);              // go back to the start of the file
    buffer = (char*) malloc(sizeof(char) * (string_size + 1) );    // allocate a string
    read_size = fread(buffer, sizeof(char), string_size, handler); // read it all in one go
    buffer[read_size] = '\0';     // put a \0 in the last position

    // The original code had this check:
    //if(string_size != read_size)
    //{
    //  // Something went wrong, free memory and set the buffer to NULL
    //  free(buffer);
    //  buffer = NULL;
    //}
    // But we actually do run into this branch even if all goes well, apparently due to CR/LF 
    // line-ending stuff. In a unit test, the file had 7 lines and 116 characters in total with 
    // CR/LF line endings, but the function only reads 109 characters. The extra 7 characters 
    // apparently are the additional line ending symbols. ToDo: implement unit tests that write 
    // and then read files with all 3 possible variants of line endings (CR, LF, CR+LF). Maybe we 
    // need to write the files in binary mode for this, i.e. use fopen(path, "wb").

    fclose(handler);    // close the file
  }
  return buffer;
}

void rosic::writeDataToFile(const char* path, int numValues, double *x, double *y1, double *y2,
                            double *y3, double *y4, double *y5)
{
  int charsPerNumber = 25;
  int numFunctions   = 1;
  if( y2 != 0 ) numFunctions++;
  if( y3 != 0 ) numFunctions++;
  if( y4 != 0 ) numFunctions++;
  if( y5 != 0 ) numFunctions++;

  int charsPerLine   = (2*numFunctions)*(charsPerNumber+2);
  int  numChars = numValues*charsPerLine;
  char  *str    = new char[numChars];
  int  pos      = 0;
  for(int i=0; i<numValues; i++)
  {
    pos += sprintf(str+pos, "%.15f %s", x[i], " ");
    pos += sprintf(str+pos, "%.15f %s", y1[i], " \n");
  }

  if( y2 != NULL )
  {
    pos += sprintf(str+pos, "%s", "\n");
    pos += sprintf(str+pos, "%s", "\n");
    for(int i=0; i<numValues; i++)
    {
      pos += sprintf(str+pos, "%.15f %s", x[i], " ");
      pos += sprintf(str+pos, "%.15f %s", y2[i], " \n");
    }
  }
  if( y3 != NULL )
  {
    pos += sprintf(str+pos, "%s", "\n");
    pos += sprintf(str+pos, "%s", "\n");
    for(int i=0; i<numValues; i++)
    {
      pos += sprintf(str+pos, "%.15f %s", x[i], " ");
      pos += sprintf(str+pos, "%.15f %s", y3[i], " \n");
    }
  }
  if( y4 != NULL )
  {
    pos += sprintf(str+pos, "%s", "\n");
    pos += sprintf(str+pos, "%s", "\n");
    for(int i=0; i<numValues; i++)
    {
      pos += sprintf(str+pos, "%.15f %s", x[i], " ");
      pos += sprintf(str+pos, "%.15f %s", y4[i], " \n");
    }
  }
  if( y5 != NULL )
  {
    pos += sprintf(str+pos, "%s", "\n");
    pos += sprintf(str+pos, "%s", "\n");
    for(int i=0; i<numValues; i++)
    {
      pos += sprintf(str+pos, "%.15f %s", x[i], " ");
      pos += sprintf(str+pos, "%.15f %s", y5[i], " \n");
    }
  }

  writeStringToFile(path, str);
  delete[] str;

  // ToDo: get rid of the duplication
}

/*
void rosic::writeDataToFile(char* path, double *data, int pixelWidth, int pixelHeight)
{
  int size = pixelWidth*pixelHeight;
  unsigned char r, g, b;

  FILE *f = fopen(path, "w");
  if( f != NULL )
  {
    //fwrite(stringToWrite, sizeof(char), length, f);
    for(int ix=0; ix<pixelWidth; ix++)
    {
      for(int iy=0; iy<pixelHeight; iy++)
      {
        //double z  = data[pixelWidth*ix+iy];
        //double z  = data[ix+pixelHeight*iy];
        //z = clip(z, 0.0, 1.0);
        //r = b = g = (unsigned char) 255*(1.0-z);

        if(      ix == 0            && iy == 0             ) { r = 255; g = 0;   b = 0;     }  // bottom-left: red
        else if( ix == pixelWidth-1 && iy == pixelHeight-1 ) { r = 255; g = 255; b = 0;     }  // top-right: yellow
        else if( ix == 0            && iy == pixelHeight-1 ) { r = 0;   g = 255; b = 0;     }  // top-left: green
        else if( ix == pixelWidth-1 && iy == 0             ) { r = 0;   g = 255; b = 255;   }  // bottom-right: cyan
        else                                                 { r = g = b = 255*((ix+iy)%2); }  // rest: checkerboard
        fwrite(&r, sizeof(unsigned char), 1, f);
        fwrite(&g, sizeof(unsigned char), 1, f);
        fwrite(&b, sizeof(unsigned char), 1, f);
      }
    }
    fclose(f);
  }
}
*/


void rosic::writeDataToFile(const char* path, double *data, int pixelWidth, int pixelHeight,
                            double *xAxis, double *yAxis)
{
  //float zDbg[50][50]; // for debug

  FILE *f = fopen(path, "w");
  int pos = 0;
  if( f != NULL )
  {
    // number of x-rows:
    float tmp = (float) pixelWidth;
    pos += (int) fwrite(&tmp, sizeof(float), 1, f);

    // the x-axis:
    for(int ix=0; ix<pixelWidth; ix++)
    {
      tmp = (float) xAxis[ix];

      rassert(tmp<10);

      pos += (int) fwrite(&tmp, sizeof(float), 1, f);
    }

    // the matix-data:
    for(int iy=0; iy<pixelHeight; iy++)
    {
      tmp = (float) yAxis[iy];
      rassert(tmp<10);
      pos += (int) fwrite(&tmp, sizeof(float), 1, f);
      for(int ix=0; ix<pixelWidth; ix++)
      {
        float z  = (float) data[pixelWidth*iy+ix];
        //zDbg[ix][iy] = z;
        rassert(z<10);
        pos += (int) fwrite(&z, sizeof(float), 1, f);
      }
    }
    fclose(f);
  }
}
