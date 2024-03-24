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