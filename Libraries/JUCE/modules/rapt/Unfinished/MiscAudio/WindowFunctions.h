#ifndef RS_WINDOWFUNCTIONS_H
#define RS_WINDOWFUNCTIONS_H

namespace RSLib
{

  /** Returns the value of a zero-centered cosine-squared shaped window of given length, having 
  nonzero values in the range -length/2...+length/2 */
  RSLib_API double rsCosineSquaredWindow(double x, double length);
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
  RSLib_API double rsRaisedCosineWindow(double x, double length, double p = 0.0);


  /** Retruns a value of the "exact" Blackman window. Compared to the unqualified "Blackman" 
  window, this exact version uses coefficients for the terms that place place zeros at the third 
  and fourth sidelobes, thereby reducing the sidelobe levels further. The third parameter is just a 
  dummy to make the function suitable for usage with function pointers that generally may point to
  parametrized windows.
  References:
  http://en.wikipedia.org/wiki/Window_function#Blackman_windows  */
  RSLib_API double rsExactBlackmanWindow(double x, double length, double dummy = 0.0);


  /** Fills the window-array with a Hamming-window. */
  RSLib_API void rsMakeBlackmanWindow(double *window, int length);

  /** Fills the window-array with a cosine power window. */
  RSLib_API void rsMakeCosinePowerWindow(double *window, int length, double power = 2.0);

  /** Fills the window-array with a Hamming-window. */
  RSLib_API void rsMakeHammingWindow(double *window, int length);

  /** Fills the window-array with a Hanning-window. */
  RSLib_API void rsMakeHanningWindow(double *window, int length);

  /** Fills the window-array with a rectangular window. ...ahem trivial */
  RSLib_API void rsMakeRectangularWindow(double *window, int length);

  /** Returns the value of a cosine-squared windowed (normalized) sinc function. It has nonzero 
  values in the range -length/2...+length/2 and zero crossings at integer multiples of the 
  "stretch" parameter. */
  RSLib_API double rsWindowedSinc(double x, double length, double stretch);



  // todo:
  // Generalize this to make a generic sum-of-cosines window where the Hanning- and Hamming windows
  // are a special case. This class includes also Blackman, Nutall, etc. windows. We may be able to
  // come up with other window shapes that supress sidebands even more than the existing stadard
  // windows.
  // Maybe have a parameter that switches between both ends zero (like octave),
  // both ends nonzero (like matlab), and start zero and end nonzero (peridoic, suitable for 
  // phase-vocoder analysis/resnthesis)
  // For reassingment, compute also (optionally) a derivative window wd and a time-ramped window wr
  // these can be optinal arguments that default to nullptr

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



}

#endif
