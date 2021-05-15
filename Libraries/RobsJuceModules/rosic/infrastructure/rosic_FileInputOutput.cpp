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

void rosic::writeToStereoWaveFile(const char* path, double *left, double *right, int numFrames,
                                  int sampleRate, int numBits)
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
