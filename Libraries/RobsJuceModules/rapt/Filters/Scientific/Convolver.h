#ifndef RAPT_CONVOLVER_H_INCLUDED
#define RAPT_CONVOLVER_H_INCLUDED


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
  response without having to have an additional buffer for it. It can invoke the design routine
  directly on our buffer here. */
  void setLength(int newLength);
  // This may allocate if the newLength is greater than our current capacity
  // Maybe have a boolean parameter "init"
  // Document the intention better


  /** Sets up the impulse response to be used. This can be used by client code, if it already has 
  the desired impulse response in some buffer. */
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

// see rosic::ConvolverBruteForce


#endif