template<class T>
void rsWindowFunction::createWindow(T* w, int N, WindowType type, bool normalizeMean, T p)
{
  typedef WindowType WT;
  switch( type )
  {
  case WT::rectangular:    rectangular(w, N); break;
  case WT::triangularNN:   triangular( w, N); break;
    // insert other polynomial windows here...
  case WT::hanningZZ:      hanning(       w, N); break; // ZZ windows are useless in practice!
  case WT::hanningZN:      hanningZN(     w, N); break;
  case WT::hamming:        hamming(       w, N); break;
  case WT::blackman:       blackman(      w, N); break;
  case WT::blackmanHarris: blackmanHarris(w, N); break;
  case WT::blackmanNutall: blackmanNutall(w, N); break;
  case WT::nutall:         nutall(        w, N); break;

  // parameterized windows:
  case WT::truncatedGaussian: truncatedGaussian(w, N, p); break; // p is the sigma
  case WT::dolphChebychev:    dolphChebychev(   w, N, p); break; // p is sidelobe attenuation

  default: 
  {
    rsError("not yet implemented");    // more types to come...
    rectangular(w, N);
  }
  }

  if(normalizeMean) {
    T m = rsArrayTools::mean(w, N);
    rsArrayTools::scale(w, N, T(1)/m);
  }
}


template<class T>
T rsWindowFunction::getMainLobeWidth(WindowType type, T p)
{
  typedef WindowType WT;
  switch( type )
  {
  case WT::rectangular:  return T(2);

  //case triangularNN:   return T(4);  // not sure yet -> plot and verify
  case WT::hanningZZ:      return T(4);
  case WT::hanningZN:      return T(4);
  case WT::hamming:        return T(4);

  case WT::blackman:       return T(6);

  case WT::blackmanHarris: return T(8);
  case WT::blackmanNutall: return T(8);


  case WT::dolphChebychev: return rsAbs(p) / T(10);
  // verify this formula! this is a very coarse ad-hoc approximation based on the observation that 
  // with p=60 we get a roughly Blackman-like and with p=40 a roughly Hamming-like window



  default: 
  {
    rsAssert(false); // not yet implemented for the given window type
    return T(0);
  }
  }
}

template<class T>
T rsWindowFunction::getSideLobeLevel(WindowType type, T param)
{
  typedef WindowType WT;
  switch( type )
  {
  // piecewise polynomial:
  case WT::rectangular:    return T(-13.3);
  case WT::triangularNN:   return T(-26.6);

  // 2 term cosine sum :
  case WT::hanningZZ:      return T(-31.5);
  case WT::hanningZN:      return T(-31.5);
  case WT::hamming:        return T(-42.4);

  // 3 term cosine sum:
  case WT::blackman:       return T(-58.1);

  // 4 term cosine sum:
  case WT::blackmanHarris: return T(-92.1);
  case WT::blackmanNutall: return T(-93.8);
  case WT::nutall:         return T(-93.3);

  default: 
  {
    rsAssert(false); // not yet implemented for the given window type
    return T(0);
  }
  }

  // The values have been obtained by creating a length 64 window with normalized DC gain (mean of
  // the window values = 1) and doing a length 16384 FFT and reading off the value of the 2nd 
  // largest maximum and rounding it to 1 digit after the point. 
  // ToDo: 
  // -automate that process by writing a function getSideLobeLevel(T* w, int N, int fftSize)
  // -maybe give two figures after the point (maybe use a larger FFT size to get more precise 
  //  readouts)
  // -figure out, if the values change, if we use other window-lengths (they shouldn't)
  // -implement it for the missing window types
}

template<class T>
void rsWindowFunction::cosineSum(T* w, int N, T* c, int K, bool normalizeMean)
{
  for(int n = 0; n < N; n++) {
    T z  = T(2*PI*n)/T(N);
    w[n] = c[0];
    for(int k = 1; k < K; k++)
      w[n] += c[k] * cos(k*z);  // (2), Eq. 123
  }
  if(normalizeMean) {
    T m = rsArrayTools::mean(w, N);
    rsArrayTools::scale(w, N, T(1)/m);
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
void rsWindowFunction::salFlatTopFast3(T* w, int N, bool norm)
{
  T c[3] = { 0.26526, -0.5, 0.23474 };
  cosineSum(w, N, c, 3, norm);
}

template<class T>
void rsWindowFunction::salFlatTopFast4(T* w, int N, bool norm)
{
  T c[4] = { 0.21706, -0.42103, 0.28294, -0.07897 };
  cosineSum(w, N, c, 4, norm);
}

template<class T>
void rsWindowFunction::salFlatTopFast5(T* w, int N, bool norm)
{
  T c[5] = { 0.1881, -0.36923, 0.28702, -0.13077, 0.02488 };
  cosineSum(w, N, c, 5, norm);
}

template<class T>
void rsWindowFunction::salFlatTopMin3(T* w, int N, bool norm)
{
  T c[3] = { 0.28235, -0.52105, 0.19659 };
  cosineSum(w, N, c, 3, norm);
}

template<class T>
void rsWindowFunction::salFlatTopMin4(T* w, int N, bool norm)
{
  T c[4] = { 0.241906, -0.460841, 0.255381, -0.041872 };
  cosineSum(w, N, c, 4, norm);
}

template<class T>
void rsWindowFunction::salFlatTopMin5(T* w, int N, bool norm)
{
  T c[5] = { 0.209671, -0.407331, 0.281225, -0.092669, 0.0091036 };
  cosineSum(w, N, c, 5, norm);
}

template<class T>
void rsWindowFunction::hrsFlatTop70(T* w, int N, bool norm)
{
  T c[4] = { 1, -1.90796, 1.07349, -0.18199 };
  cosineSum(w, N, c, 4, norm);
}

template<class T>
void rsWindowFunction::hrsFlatTop95(T* w, int N, bool norm)
{
  T c[5] = { 1, -1.9383379, 1.3045202, -0.4028270, 0.0350665 };
  cosineSum(w, N, c, 5, norm);
}

template<class T>
void rsWindowFunction::hrsFlatTop90D(T* w, int N, bool norm)
{
  T c[5] = { 1, -1.942604, 1.340318, -0.440811, 0.043097 };
  cosineSum(w, N, c, 5, norm);
}

template<class T>
void rsWindowFunction::hrsFlatTop116D(T* w, int N, bool norm)
{
  T c[6] = { 1, -1.9575375, 1.4780705, -0.6367431, 0.1228389, -0.0066288 };
  cosineSum(w, N, c, 6, norm);
}

template<class T>
void rsWindowFunction::hrsFlatTop144D(T* w, int N, bool norm)
{
  T c[7] = { 1, -1.96760033, 1.57983607, -0.81123644, 0.22583558, -0.02773848, 0.00090360 };
  cosineSum(w, N, c, 7, norm);

}

template<class T>
void rsWindowFunction::hrsFlatTop169D(T* w, int N, bool norm)
{
  T c[8] = { 1, -1.97441842, 1.65409888, -0.95788186, 0.33673420, -0.06364621, 0.00521942, 
    -0.00010599 };
  cosineSum(w, N, c, 8, norm);
}

template<class T>
void rsWindowFunction::hrsFlatTop196D(T* w, int N, bool norm)
{
  T c[9] = { 1, -1.979280420, 1.710288951, -1.081629853, 0.448734314, -0.112376628, 0.015122992,
    -0.000871252, 0.000011896 };
  cosineSum(w, N, c, 9, norm);

}

template<class T>
void rsWindowFunction::hrsFlatTop223D(T* w, int N, bool norm)
{
  T c[10] = { 1, -1.98298997309, 1.75556083063, -1.19037717712, 0.56155440797, -0.17296769663, 
    0.03233247087, -0.00324954578, 0.00013801040, -0.00000132725 };
  cosineSum(w, N, c, 10, norm);
}

template<class T>
void rsWindowFunction::hrsFlatTop248D(T* w, int N, bool norm)
{
  T c[11] = { 1, -1.985844164102, 1.791176438506, -1.282075284005, 0.667777530266, -0.240160796576,
    0.056656381764, -0.008134974479, 0.000624544650, -0.000019808998, 0.000000132974 };
  cosineSum(w, N, c, 11, norm);
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
static void rsWindowFunction::dolphChebychev(T* w, int M, T atten)
{
  int order = M-1;
  T beta = cosh(acosh(pow(T(10), (rsAbs(atten)/T(20)))) / order);
    // the rsAbs here allows for the attenuation to be given with or without the minus sign

  // Compute the complex spectrum of the window:
  std::vector<std::complex<T>> p(M);    // heap allocation
  for(int k = 0; k < M; k++) {
    T x = beta * cos(k*PI/M);
    p[k] = rsPolynomial<T>::chebychevDirect(x, order); }

  // Compute window by FFT (shouldn't it be an IFFT? maybe that just gives rise to a shift?):
  using Trafo = rsFourierTransformerBluestein<T>;
  int shift;
  if(rsIsOdd(M)) {
    Trafo::fft(&p[0], M, false);       // heap allocation
    shift = (M+1) / 2; }
  else {
    std::complex<T> j(T(0), T(1));
    for(int k = 0; k < M; k++)
      p[k] *= exp(j*(k*PI/M));         // additional modulation required for even lengths
    Trafo::fft(&p[0], M, false);       // heap allocation
    shift = (M/2) + 1; }

  // Apply a circular shift to the window and store it in the output array:
  for(int k = 0; k < M; k++)
    w[k] = p[(k+shift)%M].real();

  rsArrayTools::normalizeMean(w, M);
  // maybe make normalization optional - but if we don't normalize at all, what will we get? maybe 
  // we should normalize the peak if "normalize mean" is not desired?
  // seems like with 60 dB attenuation, we get a mean of 1000...figure out!
}
// This implementation follows the one from scipy
// https://github.com/scipy/scipy/blob/v0.19.0/scipy/signal/windows.py#L1293-L1416
// i'm not quite sure, why they use the forward FFT and not an inverse FFT - the difference is 
// probably just a phase-shift and does not matter in this case because we shift the result later 
// anyway?

template<class T>
T rsWindowFunction::windowedSinc(T x, T length, T stretch)
{
  return rsNormalizedSinc(x/stretch) * cosineSquared(x, length);
}

/*

The formulas for the ZZ, NZ, ZN, NN variants of the Hann window are:

  0.5 * (1 - cos(2*PI*n     /  N));     // ZN
  0.5 * (1 - cos(2*PI*n     / (N-1)));  // ZZ
  0.5 * (1 - cos(2*PI*(n+1) /  N));     // NZ
  0.5 * (1 - cos(2*PI*(n+1) / (N+1)));  // NN

and similarly for other window types (ZZ is probably useless, but for the sake of completeness - or 
maybe it's useful for plots for comparison purposes but not for practical usage). Maybe use the 
general formula:

  0.5 * (1 - cos(2*PI*(n+k1) / (N+k2)));

where k1 = 0 or 1, k2 = -1 or 0 or +1 - use this for all cosine-sum windows


*/

/*

ToDo:
-implement "exact" and "optimal" variants of Hamming (a0=25/46, a1=1-a0) and Blackman (see 
 Wikipedia)

// Ideas:



// implement minimax optimized windows that have the minimum (maximum) sidelobe level
// for a 2-term window, use Hamming as starting point, for 3-term start with blackman and 
// 4-term with blackman-nutall -> this can probably done in SciPy using the optimizer
// povide them here under the names MiniMax2, MiniMax3, etc.
// or MinSideLobe2, etc. - or maybe we should start numbering at 0 - a rectangular window is
// 0th order, hanning/hamming windows are 1st order, blackman is 2nd order, blackman-harris 3rd, 
// etc. ...maybe we can have a similar hirarchy for polynomial windows: rectangular is 0th order,
// triangular is 1st order, parabolic (aka Welch) 2nd, etc. - for the even order windows, also 
// provide their inverses ...or maybe have two indexes that determine the smoothness at the center
// 0 and the endpoints +-1. like polyWindow23 is 2nd order smooth at the center and 3rd order 
// smooth at +-1
// drag over the cosine-sum windows from the prototypes that satisfy smoothness constraints
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
// we use for that? that should probably depend on the desired sidelobe level. maybe make an 
// interact in the notebook where the user can enter the desired parameters

// add closed form formulas for the window spectra, where such formulas are available - for example
// for a sum-of-cosines window type, we get a corresponing sum-of-sincs in the frequency domain


// see also the paper:
// "Spectrum and spectral density estimation by the Discrete Fourier transform (DFT), including a 
// comprehensive list of window functions and some new flat-top windows"
// http://edoc.mpg.de/395068
// it has *many* windows - especially the flat-top ones may be interesting

*/

