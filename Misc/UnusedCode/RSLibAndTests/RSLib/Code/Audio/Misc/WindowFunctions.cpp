using namespace RSLib;

double RSLib::rsCosineSquaredWindow(double x, double length)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  double y = cos(PI*x/length);
  return y*y;
}

double RSLib::rsRaisedCosineWindow(double x, double length, double p)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  double c = cos(2*PI*x/length); // compute cosine shape
  c = 0.5*(c+1);                 // normalize to range 0...1
  return p + (1-p)*c;            // "crossfade" between unity and cosine
}

double RSLib::rsExactBlackmanWindow(double x, double length, double p)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  double c1 = cos(2*PI*x/length);   // compute cosine shape
  double c2 = 2*c1*c1-1;            // == cos(4*PI*x/length)

  const double a0 = 7938.0/18608.0; // 0.42659  ~ 0.42
  const double a1 = 9240.0/18608.0; // 0.49656  ~ 0.5
  const double a2 = 1430.0/18608.0; // 0.076849 ~ 0.08;

  return a0 + a1*c1 + a2*c2;
}



void RSLib::rsMakeBlackmanWindow(double *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.42 - 0.5*cos( 2.0*PI*n / (double) (length-1))
    + 0.08*cos(4.0*PI*n / (double) (length-1)) ;
}

void RSLib::rsMakeCosinePowerWindow(double *window, int length, double power)
{
  for(int n = 0; n < length; n++)
    window[n] = pow( sin(PI * (double) n / (double) length), power );
}

void RSLib::rsMakeHammingWindow(double *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.54 -  0.46 * cos(2.0*PI*n / (double) (length-1)) ;
}

void RSLib::rsMakeHanningWindow(double *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.5 * ( 1.0 - cos(2.0*PI*n / (double) (length-1)) );
}

void RSLib::rsMakeRectangularWindow(double *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 1.0;
}

double RSLib::rsWindowedSinc(double x, double length, double stretch)
{
  return rsNormalizedSinc(x/stretch) * rsCosineSquaredWindow(x, length);
}
