#ifndef RAPT_WINDOWFUNCTIONS_H
#define RAPT_WINDOWFUNCTIONS_H

/** A class to create various window functions that are useful for spectral analysis and FIR filter
design. Some functions fill an array with values of the window function, other functions allow to 
evaluate the continuous time window function to be evaluated at arbitrary inputs. The functions are 
all static and the class only serves for putting them all under the same umbrella. 

References:
 (1) https://en.wikipedia.org/wiki/Window_function
 (2) http://edoc.mpg.de/395068 
     "Spectrum and spectral density estimation by the Discrete Fourier transform (DFT), including 
      a comprehensive list of window functions and some new flat-top windows"  */

class rsWindowFunction
{

public:

  /** The available window types. Some have a ZN, NN, ZZ qualifiers. ZN means, that the first 
  sample of the window is (Z)ero and the last sample ins (N)onzero. Mostly, zero values at the ends
  are undesirable because they artificially make the window shorter than it has to be. In some 
  cases, however, a zero sample at the start or end may be needed to make the window satisfy other
  conditions such as adding up to a constant when being overlapped with shifted versions of 
  itself. For example, a hanningZN window of length N, when overlapped with itself by N/2 (i.e. the
  hop-size is N/2) sums to unity at all times. With a NN version, one could perhaps use a hop size 
  of N/2+1 and with a ZZ version N/2-1 - but that may be an inconvenient hop size */
  enum class WindowType
  {
    // polynomial:
    rectangular,
    triangularNN, // add triangularZN,ZZ
    // add parabolic/welch

    // cosine sum, 2 terms:
    hanningZZ, // ...actually "ZZ" versions are pretty useless, i think
    hanningZN, // start at (Z)ero and end ends (N)onzero - sums to constant with overlap 1
    // hanningNN
    hamming,

    // cosine sum, 3 terms:
    blackman,  // ZZ

    // cosine sum, 4 terms:
    blackmanHarris, // NN
    blackmanNutall, // NN
    nutall,         // ZZ

    dolphChebychev,

    /*
    //wkpdFlatTop

    // Salvatore flat-top windows, see (2):
    salFlatTopFast3,   // 3 terms, fast sidelobe decay, ZN
    salFlatTopFast4,   // 4 terms, fast sidelobe decay, ZN
    salFlatTopFast5,   // 5 terms, fast sidelobe decay, ZN
    salFlatTopMin3,    // 3 terms, minimum sidelobe level, NN, asymmetrical
    salFlatTopMin4,    // 4 terms, minimum sidelobe level, NN, asymmetrical
    salFlatTopMin5,    // 5 terms, minimum sidelobe level, NN, asymmetrical

    // Heinzel/Rüdiger/Schilling flat-top windows, see (2):
    hrsFlatTop70, 
    hrsFlatTop95, 
    hrsFlatTop90D, 
    hrsFlatTop116D, 
    hrsFlatTop144D,
    hrsFlatTop169D, 
    hrsFlatTop196D, 
    hrsFlatTop223D, 
    hrsFlatTop248D,
    */




    truncatedGaussian
  };



  /** Writes window function values into the array w of length N. The type should be one of the 
  values in enum windowTypes. If normalizeMean is true, the window values will be scaled such that 
  they have a mean value of unity. This will give the window function unit gain at DC which is 
  often desirable. Some windows have an adjustable parameter - for these, the value of this 
  parameter is passed in param. The meaning of the parameter may vary from one window to 
  another. */
  template<class T>
  static void createWindow(T* w, int N, WindowType type, bool normalizeMean, T param = 0);

  /** Returns the width of the main lobe of given window in frequency bins. Those values are not
  exact but just rule-of-thumb values obtained from reading off the spectrum of the respective 
  window. For windows that don't have any parameters, you should nevertheless pass a dummy 
  value, so the compiler can figure out the template type. ...later it may be used to compute 
  mainlobe widths for parametrized windows but this is not yet implemented */
  template<class T>
  static T getMainLobeWidth(WindowType type, T param, int length);

  /** Returns the level (in decibels) of the highest sidelobe of the given window (just as the 
  mainlobe width, these values are also just rule-of-thumb values). The mainlobe is supposed to be
  normalized to 0dB and the returned value is a negative number. */
  template<class T>
  static T getSideLobeLevel(WindowType type, T param);

  // maybe have also a getSideLobeRollOff function

  /** Creates a window as a sum of cosine functions:
  w[n] = sum_k c[k] cos(k * 2*pi*n/N),     n = 0,...,N-1; k = 0,...,K-1
  where: w: window, N: length, c: coeffs, K: numTerms, n: sample index, k: term index. */
  template<class T>
  static void cosineSum(T* window, int length, T* coeffs, int numTerms, bool normalizeMean);


  /** Returns the value of a zero-centered cosine-squared shaped window of given length, having
  nonzero values in the range -length/2...+length/2 */
  template<class T>
  static T cosineSquared(T x, T length);
    // turn this into a rsCosinePowerWindow function

  /** Returns the value of a zero-centered cosine shaped window on a platform of given length,
  having nonzero values in the range -length/2...+length/2. The p parameter defines the height of
  the platform. The window can be seen as a compromise between the Hann window and the
  rectangular window. For p = 0.0 the Hann window (which is equivalent to the cosine-squared
  window) will result. For p = 1.0, the rectangular window will result. For a value of p = 0.08,
  the Hamming window will result. A value of p = 0.07672 = 0.53836-0.46164 will give equiripple in
  the sidelobes according to:
  http://en.wikipedia.org/wiki/Window_function#Hamming_window
  \todo verify this experimentally */
  template<class T>
  static T raisedCosine(T x, T length, T p = 0.0);

  /** Returns a value of the "exact" Blackman window. Compared to the unqualified "Blackman"
  window, this exact version uses coefficients for the terms that place zeros at the third and 
  fourth sidelobes, thereby reducing the sidelobe levels further at the expense of the sidelobe 
  rolloff. The third parameter is just a dummy to make the function suitable for usage with 
  function pointers that generally may point to parametrized windows.
  References:
  http://en.wikipedia.org/wiki/Window_function#Blackman_windows  */
  template<class T>
  static T exactBlackman(T x, T length, T dummy = 0.0);

  /** Blackman window.
  mainlobe width: ~6
  sidelobe rejection: ~58 dB
  sidelobe rolloff: yes (figure out numerical value)  */
  template<class T>
  static void blackman(T *window, int length);

  /** Blackman-Harris window.
  mainlobe width: ~8, 
  sidelobe rejection: ~92 dB */
  template<class T>
  static void blackmanHarris(T *window, int length);

  /** Blackman-Nutall window. 
  mainlobe width: ~8 (slighty narrower than Blackman-Harris)
  sidelobe rejection: ~96.8 dB */
  template<class T>
  static void blackmanNutall(T *window, int length);

  /** Nutall window. */
  template<class T>
  static void nutall(T *window, int length);

  /** Flat top window. This has a flat top in the frequency domain, making it suitable for 
  estimation of spectral amplitudes (when the frequency is off the bin-center, its amplitude is 
  still well represented by the bin-center). */
  template<class T>
  static void flatTop(T *window, int length);

  /** Salvatore 3-term flat-top window with fast sidelobe decay. see (2). */
  template<class T>
  static void salFlatTopFast3(T *window, int length, bool normalizeMean = true);

  /** Salvatore 4-term flat-top window with fast sidelobe decay. see (2). */
  template<class T>
  static void salFlatTopFast4(T *window, int length, bool normalizeMean = true);

  /** Salvatore 4-term flat-top window with fast sidelobe decay. see (2). */
  template<class T>
  static void salFlatTopFast5(T *window, int length, bool normalizeMean = true);

  /** Salvatore 3-term flat-top window with minimum sidelobe level. see (2). */
  template<class T>
  static void salFlatTopMin3(T *window, int length, bool normalizeMean = true);

  /** Salvatore 4-term flat-top window with minimum sidelobe level. see (2). */
  template<class T>
  static void salFlatTopMin4(T *window, int length, bool normalizeMean = true);

  /** Salvatore 5-term flat-top window with minimum sidelobe level. see (2). */
  template<class T>
  static void salFlatTopMin5(T *window, int length, bool normalizeMean = true);

  /** Heinzel/Rüdiger/Schilling flat-top window, 3 cosine terms, ~70dB sidelobes. see (2) */
  template<class T>
  static void hrsFlatTop70(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop95(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop90D(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop116D(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop144D(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop169D(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop196D(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop223D(T *window, int length, bool normalizeMean = true);

  template<class T>
  static void hrsFlatTop248D(T *window, int length, bool normalizeMean = true);



  /** Truncated Gaussian window. The magnitude response of this has also an (approximate) gaussian 
  shape which translates to a parabola for the respective dB values. This makes it suitable for 
  frequency estimation by parabolic interpolation. */
  template<class T>
  static void truncatedGaussian(T *window, int length, T sigma);


  /** Fills the window-array with a cosine power window. */
  template<class T>
  static void cosinePower(T *window, int length, T power = 2.0);




  /** Fills the window-array with a Hamming-window. 
  mainlobe width: 4 
  sidelobe rejection: ~42.7 dB  */
  template<class T>
  static void hamming(T *window, int length);

  /** Fills the window-array with a Hanning-window. 
  mainlobe width: 4 
  sidelobe rejection: ~31.5 dB  */
  template<class T>
  static void hanning(T *window, int length);
  // todo: rename to specify the values at endpoints as done in rsHanningWindowZN

  /** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero
  value in w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again.
  That means, the window has a period length of N. Such a window is suitable for applications
  where it is important that suitably overlapped windows sum up to a constant, like when identity
  resynthesis is required. */
  template<class T>
  static void hanningZN(T *window, int length);
  // maybe have a version NZ, ZZ, NN

  
  /** Fills the array with a rectangular window (i.e. just all ones).
  mainlobe width: 2 
  sidelobe rejection: ~13.2 dB  
  sidelodbe rolloff: 6 dB/oct  (verify)  */
  template<class T>
  static void rectangular(T *window, int length);

  /** Triangular window.
  mainlobe width: 4 (verify)
  sidelobe rejection: 26.5 dB
  sidelobe rolloff: 12 dB/oct   */
  template<class T>
  static void triangular(T *window, int length);

  /** Dolph-Chebychev window with a given attenuation of the sidelobes in dB (it can be given with
  or without the minus sign). This window has equiripple sidebands and achieves the lowest sidelobe
  level for a given mainlobe width - or the narrowest mainlobe width for a given sidelobe level. */
  template<class T>
  static void dolphChebychev(T *window, int length, T attenuation);
  // allocates temporary heap-memory

  // maybe provide a function dolphChebychevByWidth in which the parameter gives the mainlobe width

  // todo: should compute the mainlobe-width from the sidelobe attenuation (length may not be 
  // needed when the width is normalized somehow, we'll see)
  template<class T>
  static T dolphChebychevMainLobeWidth(int length, T attenuation);


  /** Returns the value of a cosine-squared windowed (normalized) sinc function. It has nonzero
  values in the range -length/2...+length/2 and zero crossings at integer multiples of the
  "stretch" parameter. */
  template<class T>
  static T windowedSinc(T x, T length, T stretch);
  //...hmm...this does not really belong here - maybe move to FIR filter or interpolation class

};

// todo:
// -implement Hann-Poisson window - it has no sidelobes(!!!) - that may be a useful feature for the 
//  sinusoidal analysis:
//  https://en.wikipedia.org/wiki/Window_function#Hann%E2%80%93Poisson_window
// -Generalize this to make a generic sum-of-cosines window where the Hanning- and Hamming windows
//  are a special case. This class includes also Blackman, Nutall, etc. windows. We may be able to
//  come up with other window shapes that supress sidebands even more than the existing standard
//  windows.
// -Maybe have a parameter that switches between both ends zero (like octave),
//  both ends nonzero (like matlab), and start zero and end nonzero (periodic, suitable for 
//  phase-vocoder analysis/resnthesis)...or start nonzero and end zero
//  -for spectral analysis (without resynthesis), it seems best to have both ends nonzero in order 
//   to not artificially shorten the window
// -For reassingment, compute also (optionally) a derivative window wd and a time-ramped window wr
//  these can be optional arguments that default to nullptr

// requirements on windows:
// for spectrum analysis:
// -a narrower mainlobe allows for better frequency discrimination
// -a flat top in the frequency domain will give smaller amplitude error for frequencies that
//  fall between the analysis bins (see SRS flat top window)
// -a decreasing sidelobe level may be less interesting as feature since it doesn't make much of
//  a difference, if leakage error comes from nearby or more distant bins - it will be just "noise"
//  or "error" all the same
//
// for windowed-sinc filters and -interpolation:
// -low sidelobes are most important - their height ultimately determines the ripple and leakage
// -mainlobe width is less important since the transition width can be narrowed arbitrarily by
//  using a longer filter/interpolator
// -a sidelobe rolloff means less aliasing into the low frequency range in case of windowed sinc
//  interpolators - aliasing is confined to higher frequency ranges where it is less objectionable

#endif
