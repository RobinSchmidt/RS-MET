#ifndef RAPT_HILBERTFILTER_H_INCLUDED
#define RAPT_HILBERTFILTER_H_INCLUDED



//=================================================================================================


/** Under construction. Not yet usable */

template<class TSig, class TPar>
class rsConvolverNaive
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  rsConvolverNaive();


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline TSig getSample(TSig x);


protected:

  std::vector<TPar> h;    // Impulse response
  std::vector<TSig> buf;  // Circular buffer for input samples
  rsUint32 length = 0;    // Number of filter taps (nominal, i.e. zero-valued taps also count)
  rsUint32 tapIn  = 0;    // Keeps track of position in circular buffer


};

template<class TSig, class TPar>
TSig rsConvolverNaive<TSig, TPar>::getSample(TSig x)
{
  // Accept new input:
  buf[tapIn] = x;

  // Compute output:
  TSig sum = T(0);
  for(rsUint32 i = 0; i <= tapIn; i++)
    sum += h[i] * buf[tapIn-i];
  for(rsUint32 i = tapIn+1; i < length; i++)
    sum += h[i] * buf[tapIn-i+length];
  
  // Update write index:
  tapIn++;
  if(tapIn >= length)
    tapIn -= length;

  // Return result:
  return sum;
}
// Needs test
// Optimize this: use s buffer length that is a power of 2 wrap around via masking...maybe


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
  rsHilbertFilter();


protected:



};


#endif