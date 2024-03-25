#ifndef RAPT_HILBERTFILTER_H_INCLUDED
#define RAPT_HILBERTFILTER_H_INCLUDED



//=================================================================================================
// This file is still very much under ocnstruction and needs to be cleaned up




//=================================================================================================

/** Under construction. */

template<class TSig, class TPar>
class rsHilbertFilter
{

public:



  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  //rsHilbertFilter() {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the maximal (nominal) length of the impulse response. The actual length can be upt to
  two samples longer due to smoothing, if smoothing is active. @see setSmoothing. When you intend
  to operate the filter in a realtime context, it's advisable to call this function with the
  maximum length you want to support some time soon after construction or inside some sort of 
  prepareToPlay function (IIRC, that is what is was called in VST2) to avoid memory re-allocations
  in calls to setNominalLength. */
  void setMaxNominalLength(int newMaxLength) { convolver.setMaxLength(newMaxLength+2); }
  // +2 because smoothing can increase actual length by at most 2 sample with respect to nominal
  // length.

  /** Sets the nominal length of the Hilber filter impulse response, i.e. the number of taps of the
  raw Hilbert filter. The actual filter may be up to 2 taps longer when smoothing is activated. */
  void setNominalLength(int newLength) { nominalLength = newLength; setDirty();  }

  /** Sets the window function to be used to window the impulse response. The default window is a 
  Blackman window, so if you never call this function, a Blackman window will be used. */
  void setWindow(rsWindowFunction::WindowType newWindow) { window = newWindow; setDirty(); }
  // I'm not yet sure, if Blackman is a good default. We'll see...

  /** Activates smoothing. This will modify the Hilbert filter's impulse response in such a way 
  that the end result is that of a Hilbert filter with a moving average filter after it */
  void setSmoothing(bool shouldSmooth) { smooth = shouldSmooth; setDirty(); }
  // ToDo: make unit test that check, if the smoothing result is indeed the same as when literally
  // applying the MA filter. Maybe also test, if literally using an MA filter rather than smooting
  // has numerical precision advantages. Maybe provide more advanced smoothing modes for even 
  // lengths based on higher order interpolation


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  //TPar getDelay() const {  }
  // This will be a half-integer for even nominal length without smoothing


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes one output sample for a given input sample at a time. */
  inline TSig getSample(TSig in) 
  { 
    if(isDirty())
      updateCoeffs();
    return convolver.getSample(in); 
  }

  /** Resets the filter state. */
  void reset() { convolver.reset(); }



protected:


  /** Called from the setters to indicate tha coefficient recalculation is necesarry. */
  void setDirty() { dirty = true; }

  /** Called from getSample to trigger coefficient recalculation on the audio thread, if needed. */
  bool isDirty() const { return dirty; }

  /** Updates the filter coefficients according to the settings and sets the dirty flag to 
  false. */
  void updateCoeffs();


  // Data members:
  rsConvolverNaive<TSig, TPar> convolver;
  rsWindowFunction::WindowType window = rsWindowFunction::WindowType::blackman;
  int nominalLength  = 101;
  bool smooth = false;
  std::atomic<bool> dirty  = true;

};

template<class TSig, class TPar>
void rsHilbertFilter<TSig, TPar>::updateCoeffs()
{
  using WFD = rsWindowedFilterDesigner;
  int M = nominalLength;
  if(smooth)
  {
    if(rsIsEven(M))
    {
      M += 1;
      convolver.setLength(M);
      WFD::hilbertSmoothed(convolver.getCoeffPointer(), convolver.getLength(), window, true);
    }
    else
    {
      M += 2;
      convolver.setLength(M);
      WFD::hilbertSmoothed(convolver.getCoeffPointer(), convolver.getLength(), window, false);
    }
  }
  else
  {
    convolver.setLength(M);
    WFD::hilbert(convolver.getCoeffPointer(), convolver.getLength(), window);
  }
  dirty = false;
}

//=================================================================================================

/** Under construction

Combines an FIR Hilbert filter with a suitable compensation delay for the direct signal to produce
a complex signal with synchronized, delayed real and imaginary parts.

. ...TBC...  */

template<class TSig, class TPar>
class rsComplexifier
{

public:

  void setMaxLength(int newMaxLength)  // maybe rename to setMaxNominalLength
  {
    hilbert.setMaxNominalLength(newMaxLength);

    delay.setCapacity(newMaxLength);     // Maybe rename function to setMaxDelay
    // actually half of that should be enough. maybe plus 1
  }

  void setLength(int newLength) 
  { 
    hilbert.setNominalLength(newLength); 
    delay.setLength(newLength/2);   // call hilbert.getDelay
  }
  // Needs unit tests for even and odd lengths.

  void setSmoothing(bool shouldSmooth) { hilbert.setSmoothing(shouldSmooth); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes one output sample frame at a time. The first parameter is an in/out parameter for 
  the real part which will be delayed on output. The second parameter is an out parameter for the
  imaginary part. */
  inline void processSampleFrame(TSig* inOutRe, TSig* outIm)
  {
    *outIm   = hilbert.getSample(*inOutRe);  // Imaginary part
    *inOutRe = delay.getSample(*inOutRe);    // Real part (delayed)
  }

  void reset()
  {
    hilbert.reset();
    delay.reset();
  }

protected:

  rsHilbertFilter<TSig, TPar> hilbert;
  rsDelayBuffer<TSig>         delay;

};



#endif