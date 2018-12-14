template<class T>
T rsWindowFunction::rsCosineSquaredWindow(T x, T length)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  T y = cos(PI*x/length);
  return y*y;
}

template<class T>
T rsWindowFunction::rsRaisedCosineWindow(T x, T length, T p)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  T c = cos(2*PI*x/length); // compute cosine shape
  c = 0.5*(c+1);                 // normalize to range 0...1
  return p + (1-p)*c;            // "crossfade" between unity and cosine
}

template<class T>
T rsWindowFunction::rsExactBlackmanWindow(T x, T length, T p)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  T c1 = cos(2*PI*x/length);   // compute cosine shape
  T c2 = 2*c1*c1-1;            // == cos(4*PI*x/length)

  const T a0 = 7938.0/18608.0; // 0.42659  ~ 0.42
  const T a1 = 9240.0/18608.0; // 0.49656  ~ 0.5
  const T a2 = 1430.0/18608.0; // 0.076849 ~ 0.08;

  return a0 + a1*c1 + a2*c2;
}

template<class T>
void rsWindowFunction::rsBlackmanWindow(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.42 - 0.5*cos( 2.0*PI*n / (T) (length-1))
    + 0.08*cos(4.0*PI*n / (T) (length-1)) ;
}

template<class T>
void rsWindowFunction::rsCosinePowerWindow(T *window, int length, T power)
{
  for(int n = 0; n < length; n++)
    window[n] = pow( sin(PI * (T) n / (T) length), power );
}

template<class T>
void rsWindowFunction::rsHammingWindow(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.54 -  0.46 * cos(2.0*PI*n / (T) (length-1)) ;
}

template<class T>
void rsWindowFunction::rsHanningWindow(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.5 * ( 1.0 - cos(2.0*PI*n / (T) (length-1)) );
}
// i think, this is the HanningWindowZZ? ...make also a HanningWindowNN

template<class T>
void rsWindowFunction::rsHanningWindowZN(T *w, int N)
{
  T s = 2*PI/N; // for a window that ends with zero: w[N-1]=0, this would be s=2*PI/(N-1)
  for(int n = 0; n < N; n++)
    w[n] = 0.5*(1-cos(s*n));
}


template<class T>
void rsWindowFunction::rsRectangularWindow(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 1.0;
}

template<class T>
T rsWindowFunction::rsWindowedSinc(T x, T length, T stretch)
{
  return rsNormalizedSinc(x/stretch) * rsCosineSquaredWindow(x, length);
}
