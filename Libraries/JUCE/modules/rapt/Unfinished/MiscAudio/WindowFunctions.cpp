template<class T>
void rsWindowFunction::createWindow(T* w, int N, int type, bool normalizeMean, T param)
{
  switch( type )
  {
  case RECTANGULAR_WINDOW:  rectangular(w, N); break;
  case HANNING_WINDOW:      hanning(    w, N); break;
  case HANNING_WINDOW_ZN:   hanningZN(  w, N); break;
  case HAMMING_WINDOW:      hamming(    w, N); break;
  case BLACKMAN_WINDOW:     blackman(   w, N); break;
  case BLACKMAN_HARRIS:     blackmanHarris(w, N); break;
  case BLACKMAN_NUTALL:     blackmanNutall(w, N); break;

    // more types to come...
  default: rectangular(w, N);
  }

  if(normalizeMean) {
    T m = rsArray::mean(w, N);
    rsArray::scale(w, N, T(1)/m);
  }
}

template<class T>
T rsWindowFunction::cosineSquared(T x, T length)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  T y = cos(PI*x/length);
  return y*y;
}

template<class T>
T rsWindowFunction::raisedCosine(T x, T length, T p)
{
  if( rsAbs(x) > 0.5*length )
    return 0.0;
  T c = cos(2*PI*x/length); // compute cosine shape
  c = 0.5*(c+1);                 // normalize to range 0...1
  return p + (1-p)*c;            // "crossfade" between unity and cosine
}

template<class T>
T rsWindowFunction::exactBlackman(T x, T length, T p)
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
void rsWindowFunction::blackman(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.42 - 0.5*cos( 2.0*PI*n / (T) (length-1))
    + 0.08*cos(4.0*PI*n / (T) (length-1)) ;
}

template<class T>
void rsWindowFunction::blackmanHarris(T* w, int N)
{
  T a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;
  for(int n = 0; n < N; n++) {
    T s = T(2*PI*n)/T(N-1);
    w[n] = a0 - a1*cos(s) + a2*cos(2*s) - a3*cos(3*s);
  }
  // https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window
}

template<class T>
void rsWindowFunction::blackmanNutall(T* w, int N)
{
  T a0 = 0.3635819, a1 = 0.4891775, a2 = 0.1365995, a3 = 0.0106411;
  for(int n = 0; n < N; n++) {
    T s = T(2*PI*n)/T(N-1);
    w[n] = a0 - a1*cos(s) + a2*cos(2*s) - a3*cos(3*s);
  }
  //https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Nuttall_window
}

// maybe make windows that minimize the sidlobe level with a given number of cosine terms
// Blackman-Nutall looks close to equiripple sidelobes so this would be close to an optimal
// 4-term window...could be used for windowed sinc filters and spectrum analysis




template<class T>
void rsWindowFunction::cosinePower(T *window, int length, T power)
{
  for(int n = 0; n < length; n++)
    window[n] = pow( sin(PI * (T) n / (T) length), power );
}

template<class T>
void rsWindowFunction::hamming(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.54 -  0.46 * cos(2.0*PI*n / (T) (length-1)) ;
}

template<class T>
void rsWindowFunction::hanning(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 0.5 * ( 1.0 - cos(2.0*PI*n / (T) (length-1)) );
}
// i think, this is the HanningWindowZZ? ...make also a HanningWindowNN

template<class T>
void rsWindowFunction::hanningZN(T *w, int N)
{
  T s = 2*PI/N; // for a window that ends with zero: w[N-1]=0, this would be s=2*PI/(N-1)
  for(int n = 0; n < N; n++)
    w[n] = 0.5*(1-cos(s*n));
}


template<class T>
void rsWindowFunction::rectangular(T *window, int length)
{
  for(int n = 0; n < length; n++)
    window[n] = 1.0;
}

template<class T>
T rsWindowFunction::windowedSinc(T x, T length, T stretch)
{
  return rsNormalizedSinc(x/stretch) * cosineSquared(x, length);
}
