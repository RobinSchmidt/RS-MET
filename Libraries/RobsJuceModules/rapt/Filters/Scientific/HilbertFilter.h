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

  /** Sets up the impulse response to be used. */
  void setImpulseResponse(TPar* newImpulseResponse, int newLength);
  // This may allocate if the newLength is greater than our current capacity


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

};

template<class TSig, class TPar>
void rsConvolverNaive<TSig, TPar>::setImpulseResponse(TPar* hNew, int newLength)
{
  length = newLength;
  h.resize(length);
  buf.resize(length);
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
  rsFill(buf, T(0));
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


/** Under construction. Not yet usable */

template<class TSig, class TPar>
class rsHilbertFilter : public rsConvolverNaive<TSig, TPar>
{

public:



  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  //rsHilbertFilter() {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */



  //-----------------------------------------------------------------------------------------------
  /** \name Design */

  void computeCoeffs(TPar* h, int numTaps, RAPT::rsWindowFunction::WindowType type);





protected:

  // ToDo: embedd an object of class rsConvolverNaive rather than subclassing:
  //rsConvolverNaive<TSig, TPar> convolver;
  // ..this will need a bit of delegation but it's cleaner API-wise.

};


#endif