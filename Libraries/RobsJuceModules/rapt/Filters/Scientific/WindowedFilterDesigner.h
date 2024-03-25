#ifndef RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED
#define RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED


/** This class has functions to design various filters using the windowing method. The impulse 
response will be written into the array "h" which must be "numTaps" long. You also need to 
specify which window function shall be used via the "type" parameter as one of the types 
defined in the  rsWindowFunction::WindowType enum class. */

class rsWindowedFilterDesigner
{

public:

  /** Creates a Hilbert transformer filter. An ideal Hilbert transformer shifts the phases of all 
  frequencies by 90°. ...TBC... */
  template<class T>
  static void hilbert(T* h, int numTaps, rsWindowFunction::WindowType type);
  // ToDo: explain properties of even and odd length Hilbert filters

  /** A function to design smoothed Hilbert filters. They represent a Hilbert filter with an 
  additional MA filter applied. Such smoothed Hilbert filters will always have odd lengths, so 
  numTaps must be odd. If the nominal length is even, the smoothed length will be 1 sample longer 
  because we use a two-sample MA with kernel [0.5 0.5] for smoothing which increases the length by 
  one. If the nominal length is odd, the smoothed length will be 2 samples longer because we use 
  3-sample MA with kernel [0.25 0.5 0.25]. For odd nominal lengths, such smoothing will remove the 
  Nyquist ripple artifacts present in the approximated Hilbert trafo of a saw wave obtained by the 
  filter. For even nominal length, the smoothing will bring the orignally half-integer delay that 
  the filter introduces back to a full integer delay, making it easier to align the Hilbert trafo 
  with the original signal.  ....TBC...*/
  template<class T>
  static void hilbertSmoothed(T* h, int numTaps, rsWindowFunction::WindowType type, 
    bool evenNominalLength);




};



#endif