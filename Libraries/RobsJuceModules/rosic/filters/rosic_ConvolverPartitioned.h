#ifndef rosic_ConvolverPartitioned_h
#define rosic_ConvolverPartitioned_h

namespace rosic
{

/** This class implements a convolution via partioning the impulse response into several blocks 
and convolving the input signal with each of these blocks separately. In order to achieve a zero 
input/output delay, convolution with the first block will be realized by direct convolution and 
later blocks will be computed by FFT/IFFT. The algorithms used here partitions the impulse response
in a way which minimizes the the total computation per sample but does not ensure uniform CPU load. 
Best efficiency in relation to the impulse response length is achieved when the length is a power 
of two, the most wasteful case occurs at a length of a power of two plus one.

\todo:
-use a more proper size for the last fftConvolver and delay it's output appropriately
-maybe try to achieve more uniform load by changing the algo...  */

class ConvolverPartitioned
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  //ConvolverPartitioned();

  /** Destructor. */
  //~ConvolverPartitioned();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the impulse response to be used. */
  void setImpulseResponse(double* newImpulseResponse, int newLength);

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a single convolved output-sample. */
  INLINE double getSample(double in);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Sets the impulse response to an impulse (1 0 0 0....). */
  void clearImpulseResponse();

  /** Sets the buffer for the previous input samples to zero. */
  void clearInputBuffers();

  //===============================================================================================

protected:

  static const int directConvolutionLength = 64;
  // \todo: check whether 64 is the optimal value, maybe use some value that depends on the 
  // overhanging length of the impulse-response

  ConvolverBruteForce directConvolver;     // a single direct convolver for the first block...
  std::vector<ConvolverFFT> fftConvolvers; // and a bunch of FFT convolvers of increasing length
  int M = 0;                               // length of the impulse response
};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE double ConvolverPartitioned::getSample(double in)
{
  double tmp = directConvolver.getSample(in);
  for(size_t c = 0; c < fftConvolvers.size(); c++)
    tmp += fftConvolvers[c].getSample(in);
  return tmp;
}

}

#endif 
