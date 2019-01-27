template<class T>
void rsWindowFunction::createWindow(T* w, int N, int type, bool normalizeMean, T p)
{
  switch( type )
  {
  case RECTANGULAR_WINDOW:  rectangular(w, N); break;
  case TRIANGULAR_WINDOW:   triangular( w, N); break;
    // insert other polynomial windows here...

  case HANNING_WINDOW:      hanning(    w, N); break;
  case HANNING_WINDOW_ZN:   hanningZN(  w, N); break;
  case HAMMING_WINDOW:      hamming(    w, N); break;

  case BLACKMAN_WINDOW:     blackman(      w, N); break;
  case BLACKMAN_HARRIS:     blackmanHarris(w, N); break;
  case BLACKMAN_NUTALL:     blackmanNutall(w, N); break;
  case NUTALL:              nutall(        w, N); break;

  case TRUNCATED_GAUSSIAN:  truncatedGaussian(w, N, p); break; // p is the sigma

    // more types to come...
  default: rectangular(w, N);
  }

  if(normalizeMean) {
    T m = rsArray::mean(w, N);
    rsArray::scale(w, N, T(1)/m);
  }
}


template<class T>
T rsWindowFunction::getMainLobeWidth(int type, T param)
{
  switch( type )
  {
  case RECTANGULAR_WINDOW:  return T(2);

  //case TRIANGULAR_WINDOW:   return T(4);  // not sure yet -> plot and verify
  case HANNING_WINDOW:      return T(4);
  case HANNING_WINDOW_ZN:   return T(4);
  case HAMMING_WINDOW:      return T(4);

  case BLACKMAN_WINDOW:     return T(6);

  case BLACKMAN_HARRIS:     return T(8);
  case BLACKMAN_NUTALL:     return T(8);

  default: 
  {
    rsAssert(false); // not yet implemented for the given window type
    return T(0);
  }
  }
}

template<class T>
T rsWindowFunction::getSideLobeLevel(int type, T param)
{
  switch( type )
  {
  case RECTANGULAR_WINDOW:  return T(-13.2);

  case TRIANGULAR_WINDOW:   return T(-26.5);

  case HANNING_WINDOW:      return T(-31.5);
  case HANNING_WINDOW_ZN:   return T(-31.5);
  case HAMMING_WINDOW:      return T(-42.7);

  case BLACKMAN_WINDOW:     return T(-58);

  case BLACKMAN_HARRIS:     return T(-92);
  case BLACKMAN_NUTALL:     return T(-96.8);

  default: 
  {
    rsAssert(false); // not yet implemented for the given window type
    return T(0);
  }
  }
  // todo: verify these values in plots
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

template<class T>
void rsWindowFunction::nutall(T* w, int N)
{
  T a0 = 0.355768, a1 = 0.487396, a2 = 0.144232, a3 = 0.012604;
  for(int n = 0; n < N; n++) {
    T s = T(2*PI*n)/T(N-1);  // factor out: c = T(2*PI)/T(N-1), do here: s = c*n
    w[n] = a0 - a1*cos(s) + a2*cos(2*s) - a3*cos(3*s);
  }
  // https://en.wikipedia.org/wiki/Window_function#Nuttall_window,_continuous_first_derivative
}

template<class T>
void rsWindowFunction::flatTop(T* w, int N)
{
  T a0 = 0.21557895, a1 = 0.41663158, a2 = 0.277263158, a3 = 0.083578947, a4 = 0.006947368;
  for(int n = 0; n < N; n++) {
    T s = T(2*PI*n)/T(N-1);
    w[n] = a0 - a1*cos(s) + a2*cos(2*s) - a3*cos(3*s) + a4*cos(4*s);
  }
  //https://en.wikipedia.org/wiki/Window_function#Flat_top_window
}

template<class T>
void rsWindowFunction::truncatedGaussian(T* w, int N, T sigma)
{
  T m = T(0.5)*(N-1);  // mu, midpoint, center
  for(int n = 0; n < N; n++) {
    T t = (n-m)/(sigma*m);
    w[n] = exp(T(-0.5) * t*t);
  }
  // https://en.wikipedia.org/wiki/Window_function#Gaussian_window 
}

// todo: confined gaussian: 
// https://en.wikipedia.org/wiki/Window_function#Confined_Gaussian_window



// maybe make windows that minimize the sidelobe level with a given number of cosine terms
// Blackman-Nutall looks close to equiripple sidelobes so this would be close to an optimal
// 4-term window...could be used for windowed sinc filters and spectrum analysis
// but maybe for filters, equiripple is less desirable than a falloff for the tails




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
void rsWindowFunction::triangular(T *w, int N)
{
  T m = T(0.5)*T(N-1);
  T L = N;        // could also be N-1 or N+1, see wikipedia article - maybe let the user choose
  T s = T(2)/L;
  for(int n = 0; n < N; n++)
    w[n] = 1 - rsAbs((n-m)*s);
  // https://en.wikipedia.org/wiki/Window_function#Triangular_window
}

template<class T>
T rsWindowFunction::windowedSinc(T x, T length, T stretch)
{
  return rsNormalizedSinc(x/stretch) * cosineSquared(x, length);
}

// ideas:

// implement minimax optimized windows that have the minimum (maximum) sidelobe level
// for a 2-term window, use Hamming as starting point, for 3-term start with blackman and 
// 4-term with blackman-nutall -> this can probably done in SciPy using the optimizer
// povide them here under the names MiniMax2, MiniMax3, etc.
// or MinSideLobe2, etc. - or maybe we should start numbering at 0 - a rectangular window is
// 0th order, hanning/hamming windows are 1st order, blackman is 2nd order, blackman-harris 3rd, 
// etc. ...maybe we can have a similar hirarchy for polynomial windows: rectangular is 0th order,
// triangular is 1st order, parabolic (aka Welch) 2nd, etc.
// hmm - here, it looks like that the flat-top window is also equiripple, so maybe the minimax 
// optimized windows all will end up being flat-top? (bcs i think minimax and equiripple conditions 
// imply each other - verify that):
// https://www.mathworks.com/help/signal/ug/generalized-cosine-windows.html

// or maybe just use the Dolph-Chebychev window:
// https://en.wikipedia.org/wiki/Window_function#Dolph%E2%80%93Chebyshev_window
// https://ccrma.stanford.edu/~jos/sasp/Dolph_Chebyshev_Window.html

// here is a paper:
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.77.5098&rep=rep1&type=pdf

// this has c-code - but should probably be used only for a reference prototype since the 
// computation has complexity O(N^2)
// http://practicalcryptography.com/miscellaneous/machine-learning/implementing-dolph-chebyshev-window/

// but the Dolph-Chebychev window is not a sum-of-cosines type, so maybe we should obtain coeffs
// for optimized sum-of-cosine types anyway as they can be used for efficient windowed-sinc 
// interpolation via trig-recursions

// can we also somehow optimize the mainlobe? maybe, if we put a constraint on the maximum 
// sidelobe-level, i.e. find the window that has the smallest mainlobe width given a maximum
// sidelobe-level -> a constrained optimization problem? maybe use the min-sidelobe windows as
// starting points. if it works out, maybe provide the resulting windows for several selected 
// sidelobe levels, like -20, -30, -40, ... - but: which order (i.e. number of cosine terms) should 
// we use for that? that should probably depend on the desired sidelobe level maybe make and 
// interact in the notebook where the user can enter the desired parameters