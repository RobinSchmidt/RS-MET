using namespace std;

namespace RSLib
{

  void writeToMonoWaveFile(const char* path, float *signal, int numFrames, int sampleRate, 
    int numBits)
  {
    rsOutputWaveFile file(path, sampleRate, numBits, 1);
    file.write(signal, numFrames);
  }

  void writeToMonoWaveFile(const char* path, double *signal, int numFrames, int sampleRate, 
    int numBits)
  {
    rsOutputWaveFile file(path, sampleRate, numBits, 1);
    float *tmp = new float[numFrames];
    for(int n = 0; n < numFrames; n++)
      tmp[n] = (float) signal[n];
    file.write(tmp, numFrames);
    delete[] tmp;
  }

  void writeToStereoWaveFile(const char* path, double *left, double *right, int numFrames, 
    int sampleRate, int numBits)
  {
    rsOutputWaveFile file(path, sampleRate, numBits, 2);
    float *tmp = new float[2*numFrames];
    for(int n = 0; n < numFrames; n++)
    {
      tmp[2*n]   = (float) left[n];
      tmp[2*n+1] = (float) right[n];
    }
    file.write(tmp, 2*numFrames);
    delete[] tmp;
  }

  double** readFromWaveFile(const char* path, int& numChannels, int& numFrames, 
    int& sampleRate)
  {
    try
    {
      rsInputWaveFile file(path);
      numChannels = file.getNumChannels();
      numFrames   = file.getNumSampleFrames();
      sampleRate  = file.getSampleRate();
      if( numChannels <= 0 || numFrames <= 0 )
        return NULL;

      float *fBuffer = new float[numChannels*numFrames];
      file.readAndConvertToFloat(fBuffer, numChannels*numFrames);
      if( numChannels > 1 )
        rsDeInterleave(fBuffer, numFrames, numChannels);

      double **pointers = new double*[numChannels];
      for(int c = 0; c < numChannels; c++)
      {
        pointers[c] = new double[numFrames];
        for(int i = 0; i < numFrames; i++)
          pointers[c][i] = (double) fBuffer[c*numFrames+i];
      }

      delete fBuffer;
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

  double* readMonoWaveFile(const char* path, int& numFrames, int& sampleRate)
  {
    int numChannels;
    double **px = readFromWaveFile(path, numChannels, numFrames, sampleRate);
    double *x = new double[numFrames];
    rsCopyBuffer(px[0], x, numFrames);
    for(int n = 0; n < numChannels; n++)
      delete px[n];
    delete[] px;
    return x;
  }


  void writeStringToFile(const char* path, char* stringToWrite)
  {
    int length = (int)strlen(stringToWrite);
    FILE *f = fopen(path, "w");
    if( f != NULL )
    {
      fwrite(stringToWrite, sizeof(char), length, f);
      fclose(f);
    }
  }

  void writeDataToFile(const char* path, int numValues, double *x, double *y1, double *y2,
    double *y3, double *y4, double *y5)
  {
    int charsPerNumber = 25;
    int numFunctions   = 1;
    if( y2 != 0 ) numFunctions++;
    if( y3 != 0 ) numFunctions++;
    if( y4 != 0 ) numFunctions++;
    if( y5 != 0 ) numFunctions++;

    int  charsPerLine = (2*numFunctions)*(charsPerNumber+2);
    int  numChars     = numValues*charsPerLine;
    char *str         = new char[numChars];
    int  pos          = 0;
    for(int i = 0; i < numValues; i++)
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
    delete str;
  }

  void writeDataToFile(const char* path, double *data, int pixelWidth, int pixelHeight,
    double *xAxis, double *yAxis)
  {
    float zDbg[50][50];

    FILE *f = fopen(path, "w");
    int pos = 0;
    if( f != NULL )
    {
      // number of x-rows:
      float tmp = (float) pixelWidth;
      pos += (int)fwrite(&tmp, sizeof(float), 1, f);

      // the x-axis:
      for(int ix = 0; ix < pixelWidth; ix++)
      {
        tmp = (float) xAxis[ix];
        rsAssert(tmp < 10);
        pos += (int)fwrite(&tmp, sizeof(float), 1, f);
      }

      // the matrix-data:
      for(int iy = 0; iy < pixelHeight; iy++)
      {
        tmp = (float) yAxis[iy];
        rsAssert(tmp < 10);
        pos += (int)fwrite(&tmp, sizeof(float), 1, f);
        for(int ix = 0; ix < pixelWidth; ix++)
        {
          float z  = (float) data[pixelWidth*iy+ix];
          zDbg[ix][iy] = z;
          rsAssert(z<10);
          pos += (int)fwrite(&z, sizeof(float), 1, f);
        }
      }
      fclose(f);
    }
  }

}
