//#include "rosic_WindowDesigner.h"
//using namespace rosic;

void WindowDesigner::getWindow(double *window, int length, int type)
{
  RAPT::rsArrayTools::fillWithValue(window, length, 1.0);
  applyWindow(window, length, type);
}

void WindowDesigner::applyWindow(double *buffer, int length, int type)
{
  switch( type )
  {
  case BLACKMAN:
    {
      for(int n=0; n<length; n++)
        buffer[n] *= blackmanWindow(n, length);
    }
    break;
  case COSINE_SQUARED:
    {
      for(int n=0; n<length; n++)
        buffer[n] *= cosinePowerWindow(n, length, 2.0);
    }
    break;
  case HAMMING:
    {
      for(int n=0; n<length; n++)
        buffer[n] *= hammingWindow(n, length);
    }
    break;
  case HANN:
    {
      for(int n=0; n<length; n++)
        buffer[n] *= hannWindow(n, length);
    }
    break;
  }
}

void WindowDesigner::getCosinePowerWindow(double *window, int length, double power)
{
  for(int n=0; n<length; n++)
    window[n] = cosinePowerWindow(n, length, power);
}

void WindowDesigner::getKaiserWindow(double *window, int length, double beta)
{
  for(int n=0; n<length; n++)
    window[n] = kaiserWindow(n, length, beta);
}
