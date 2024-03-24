#ifndef RAPT_HILBERTFILTER_H_INCLUDED
#define RAPT_HILBERTFILTER_H_INCLUDED



//=================================================================================================


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


  void setMaxLength(int newMaxLength) { convolver.setMaxLength(newMaxLength); }
  // should be called before going into realtime operation to allocate enough memory fo the 
  // buffers

  void setLength(int newLength) { convolver.setLength(newLength); setDirty(); }

  void setWindow(rsWindowFunction::WindowType newWindow) { window = newWindow; setDirty(); }




  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  bool isDirty() const { return dirty; }


  //TPar getDelay() const {  }
  // This can be a half-integer


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

  void updateCoeffs();

  void setDirty() { dirty = true; }


  rsConvolverNaive<TSig, TPar> convolver;


  rsWindowFunction::WindowType window = rsWindowFunction::WindowType::blackman;
  // I'm not yet sure, if Blackman is a good default. We'll see...

  bool dirty = true;  // Maybe use std::atomic<bool>

};

template<class TSig, class TPar>
void rsHilbertFilter<TSig, TPar>::updateCoeffs()
{
  makeHilbertFilter(convolver.getCoeffPointer(), convolver.getLength(), window);
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

  void setMaxLength(int newMaxLength)
  {
    hilbert.setMaxLength(newMaxLength);
    delay.setCapacity(newMaxLength);     // Maybe rename function to setMaxDelay
  }

  void setLength(int newLength) 
  { 
    hilbert.setLength(newLength); 
    delay.setLength(newLength/2); 
  }
  // Needs unit tests for even and odd lengths.

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

protected:

  rsHilbertFilter<TSig, TPar> hilbert;
  rsDelayBuffer<TSig>         delay;

};



// ToDo:
// -Implement a convenience class that combines a Hilbert filter with a suitable compensation
//  delay to produce the delayed complex analytic signal. Maybe call is rsComplexifier or 
//  rsQuadratureFilter (but that clashes with the IIR-variant)





#endif