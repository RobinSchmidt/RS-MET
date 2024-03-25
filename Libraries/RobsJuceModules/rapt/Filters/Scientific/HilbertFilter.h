#ifndef RAPT_HILBERTFILTER_H_INCLUDED
#define RAPT_HILBERTFILTER_H_INCLUDED



//=================================================================================================
// This file is still very much under ocnstruction and needs to be cleaned up


/** Under construction.... */

template<class TSig, class TPar>
class rsConvolverNaive
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  //rsConvolverNaive();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */


  void setMaxLength(int newMaxLength) { h.reserve(newMaxLength); buf.reserve(newMaxLength); }


  /** Sets up the length of the impulse response. Note that after calling setLength, the content of
  impulse response is undefined. The intention for this function is to be called in a sequence 
  like:

    c.setLength(newLength);
    computeCoeffs(c.getCoeffPointer(), c.getLength());

  where "c" is some convolver object and computeCoeffs is some filter design routine. c.getLength()
  should return "newLength". The idea is that we want client code to be able to set up the impulse 
  response without having to have an additional buffer for it. It can invoke the design routing 
  directly on our buffer here. */
  void setLength(int newLength);
  // This may allocate if the newLength is greater than our current capacity
  // Maybe have a boolean parameter "init"
  // Document the intention better


  /** Sets up the impulse response to be used. This can be used by client code, if it already has 
  the impulse response in some buffer. */
  void setImpulseResponse(TPar* newImpulseResponse, int newLength);
  // This may allocate if the newLength is greater than our current capacity


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  int getLength() const { return length; }

  TPar* getCoeffPointer() { return &h[0]; }
  // Idea: If client code wants to invoke a filter design routine directly on our member array here
  // instead of first computing the impulse response in its own buffer and then calling 
  // setImpulseResponse, it can do so by doing:
  //
  //   computeCoeffs(c.getCoeffPointer(), c.getLength())
  //
  // when "c" is the convolver object



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes one output sample for a given input sample at a time. */
  inline TSig getSample(TSig in);

  /** Resets the filter state. */
  void reset();

protected:


  std::vector<TPar> h;    // Impulse response
  std::vector<TSig> buf;  // Circular buffer for input samples
  int length = 0;         // Number of filter taps
  int tapIn  = 0;         // Keeps track of position in circular buffer

  //template<class U, class V> friend class rsHilbertFilter;

};

template<class TSig, class TPar>
void rsConvolverNaive<TSig, TPar>::setLength(int newLength)
{
  length = newLength;
  h.resize(length);
  buf.resize(length);
}

template<class TSig, class TPar>
void rsConvolverNaive<TSig, TPar>::setImpulseResponse(TPar* hNew, int newLength)
{
  setLength(newLength);
  rsCopyToVector(hNew, length, h);
}

template<class TSig, class TPar>
TSig rsConvolverNaive<TSig, TPar>::getSample(TSig x)
{
  // Accept new input:
  buf[tapIn] = x;

  // Compute output:
  TSig sum = TSig(0);
  for(int i = 0; i <= tapIn; i++)
    sum += h[i] * buf[tapIn-i];
  for(int i = tapIn+1; i < length; i++)
    sum += h[i] * buf[tapIn-i+length];
  
  // Update write index:
  tapIn++;
  if(tapIn >= length)
    tapIn -= length;

  // Return result:
  return sum;
}
// Needs test
// Optimize this: use a buffer length that is a power of 2 wrap around via masking...maybe

template<class TSig, class TPar>
void rsConvolverNaive<TSig, TPar>::reset()
{
  rsFill(buf, TSig(0));
  tapIn = 0;
}


// see rosic::ConvolverBruteForce, FiniteImpulseResponseDesigner::getHilbertTransformerResponse


//=================================================================================================

/** A free function that designs an FIR Hilbert filter using the windowing method. The impulse 
response will be written into the array "h" which must be "numTaps" long. You also need to specify 
which window function shall be used via the "type" parameter as one of the types defined in the 
rsWindowFunction::WindowType enum class. */
template<class T>
void makeHilbertFilter(T* h, int numTaps, RAPT::rsWindowFunction::WindowType type)
{
  // Create the window:
  RAPT::rsWindowFunction::createWindow(h, numTaps, type, false);

  // Multiply in the Hilbert-filter weights:
  int m = numTaps/2;                          // Middle tap
  if(rsIsOdd(numTaps))
  {
    for(int k = m % 2; k < numTaps; k += 2)   // k starts at 0 if m is even and at 1 if m is odd
      h[k] = T(0);
    for(int k = 1; k <= m; k += 2)
    {
      T hk = T(2) / T(k*PI);
      h[m+k] *= +hk;
      h[m-k] *= -hk;
    }
  }
  else
  {
    for(int k = 0; k < m; k++)
    {
      T t  = T(k) + T(0.5); 
      T hk = T(1) / (t*PI);
      h[m+k]   *= +hk;
      h[m-k-1] *= -hk;
    }
  }

  // See:
  // https://www.kvraudio.com/forum/viewtopic.php?t=608320
  // https://en.wikipedia.org/wiki/Hilbert_transform#Discrete_Hilbert_transform
  // https://www.dsprelated.com/freebooks/sasp/Hilbert_Transform_Design_Example.html
  // https://www.intechopen.com/chapters/39362
  //
  // ToDo:
  // -Compare the results of this routine with those of some reference implementations from octave 
  //  or numpy/scipy
  // -Move this as static member into a class rsFiniteImpulseResponseDesigner
}


// A function to design smoothed Hilbert filters. They represent a Hilbert filter with an 
// additional MA filter applied. Such smoothed Hilbert filters will always have odd lengths. If 
// the nominal length is even, the smoothed length will be 1 sample longer because we use a 
// two-sample MA with kernel [0.5 0.5] for smoothing which increases the length by one. If the 
// nominal length is odd, the smoothed length will be 2 samples longer because we use 3-sample MA 
// with kernel [0.25 0.5 0.25]. For odd nominal lengths, such smoothing will remove the Nyquist 
// ripple artifacts present in the approximated Hilbert trafo of a saw wave obtained by the filter.
// For even nominal length, the smoothing will bring the orignally half-integer delay that the 
// filter introduces back to a full integer delay, making it easier to align the Hilbert trafo 
// with the original signal.  ....TBC...
template<class T>
void makeSmoothOddHilbertFilter(T* h, int numTaps, 
  rsWindowFunction::WindowType type, bool evenNominalLength)
{
  rsAssert(rsIsOdd(numTaps));

  int M = numTaps;
  if(evenNominalLength)
  {
    makeHilbertFilter(h, M-1, type);
    h[M-1] = T(0);
    rsArrayTools::movingAverage2ptForward(h, M, h);
  }
  else
  {
    makeHilbertFilter(&h[1], M-2, type);
    h[0]   = 0;
    h[M-1] = 0;
    rsArrayTools::weightedAverage3pt(h, M, h, T(0.25), T(0.5), T(0.25));
  }
}
// ToDo: 
// -In the case of even nominal length, the main purposes of the MA smoothing is to bring the delay
//  from a half-integer to a full integer. Maybe that can be done in better ways by interpolating
//  the impulse response with polynomial interpolators. The 2-sample MA is basically a linear 
//  interpolator that reads out the raw impulse response at half-integer poistions. We can do 
//  better than linear! Doing so may give more desirable magnitude responses, i.e. less lowpassing.
//  Although, that lowpassing might not be such a bad thing in the context of instantaneous 
//  envelope detection.



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
  /** \name Design */

  //void computeCoeffs(TPar* h, int numTaps, RAPT::rsWindowFunction::WindowType type);

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
  int M = nominalLength;
  if(smooth)
  {
    if(rsIsEven(M))
    {
      M += 1;
      convolver.setLength(M);
      makeSmoothOddHilbertFilter(convolver.getCoeffPointer(), convolver.getLength(), 
        window, true);
    }
    else
    {
      M += 2;
      convolver.setLength(M);
      makeSmoothOddHilbertFilter(convolver.getCoeffPointer(), convolver.getLength(), 
        window, false);
    }
  }
  else
  {
    convolver.setLength(M);
    makeHilbertFilter(convolver.getCoeffPointer(), convolver.getLength(), window);
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



// ToDo:
// -Implement a convenience class that combines a Hilbert filter with a suitable compensation
//  delay to produce the delayed complex analytic signal. Maybe call is rsComplexifier or 
//  rsQuadratureFilter (but that clashes with the IIR-variant)





#endif