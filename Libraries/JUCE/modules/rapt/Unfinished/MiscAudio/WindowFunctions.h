#ifndef RAPT_WINDOWFUNCTIONS_H
#define RAPT_WINDOWFUNCTIONS_H

// todo: wrap into a class rsWindowFunctions...get rid of the "make" prefixes - just use
// rsBlahWindow

//template<class T>
class rsWindowFunction
{

public:

  /** Returns the value of a zero-centered cosine-squared shaped window of given length, having
  nonzero values in the range -length/2...+length/2 */
  template<class T>
  static T rsCosineSquaredWindow(T x, T length);
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
  static T rsRaisedCosineWindow(T x, T length, T p = 0.0);

  /** Retruns a value of the "exact" Blackman window. Compared to the unqualified "Blackman"
  window, this exact version uses coefficients for the terms that place place zeros at the third
  and fourth sidelobes, thereby reducing the sidelobe levels further. The third parameter is just a
  dummy to make the function suitable for usage with function pointers that generally may point to
  parametrized windows.
  References:
  http://en.wikipedia.org/wiki/Window_function#Blackman_windows  */
  template<class T>
  static T rsExactBlackmanWindow(T x, T length, T dummy = 0.0);

  /** Fills the window-array with a Hamming-window. */
  template<class T>
  static void rsBlackmanWindow(T *window, int length);

  /** Fills the window-array with a cosine power window. */
  template<class T>
  static void rsCosinePowerWindow(T *window, int length, T power = 2.0);

  /** Fills the window-array with a Hamming-window. */
  template<class T>
  static void rsHammingWindow(T *window, int length);

  /** Fills the window-array with a Hanning-window. */
  template<class T>
  static void rsHanningWindow(T *window, int length);
  // todo: rename to specify the values at endpoints as done in rsHanningWindowZN

  /** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero
  value in w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again.
  That means, the window has a period length of N. Such a window is suitable for applications
  where it is important that suitably overlapped windows sum up to a constant, like when identity
  resynthesis is required. */
  template<class T>
  static void rsHanningWindowZN(T *window, int length);
  // maybe have a version NZ, ZZ, NN

  /** Fills the window-array with a rectangular window. ...ahem trivial */
  template<class T>
  static void rsRectangularWindow(T *window, int length);

  /** Returns the value of a cosine-squared windowed (normalized) sinc function. It has nonzero
  values in the range -length/2...+length/2 and zero crossings at integer multiples of the
  "stretch" parameter. */
  template<class T>
  static T rsWindowedSinc(T x, T length, T stretch);
  //...hmm...this does not really belong here - maybe move to and FIR filter or interpolation class

};

// todo:
// -Generalize this to make a generic sum-of-cosines window where the Hanning- and Hamming windows
//  are a special case. This class includes also Blackman, Nutall, etc. windows. We may be able to
//  come up with other window shapes that supress sidebands even more than the existing stadard
//  windows.
// -Maybe have a parameter that switches between both ends zero (like octave),
//  both ends nonzero (like matlab), and start zero and end nonzero (periodic, suitable for 
//  phase-vocoder analysis/resnthesis)...or start nonzero and end zero
//  for spectral analysis (without resynthesis), it seems best to have both ends nonzero in order 
//  to not artificially shorten the window
// -For reassingment, compute also (optionally) a derivative window wd and a time-ramped window wr
//  these can be optional arguments that default to nullptr

// requirements on windows:
// for spectrum analysis:
// -a narrower mainlobe allows for better frequency discrimination
// -a flat top in the frequency domain will give smaller amplitude error for frequencies that
//  fall between the analysis bins (see SRS flat top window)
// -a decreasing sidelobe level may be less interesting as feature since it doesn't make much of
//  a difference, if leakage error comes form nearby or more distant bins
//
// for windowed-sinc filters and -interpolation:
// -low sidelobes are most important - their height ultimately determines the ripple and leakage
// -mainlobe width is less important since the transition width can be narrowed arbitrarily by
//  using a longer filter/interpolator
// -a decreasing sidelobe level means less aliasing into the low frequency range in case of 
//  windowed sinc interpolators - aliasing is confined to higher frequency ranges where it is
//  less objectionable

#endif
