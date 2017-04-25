#ifndef rosic_ConvolverFFT_h
#define rosic_ConvolverFFT_h

// rosic-indcludes:
#include "../transforms/rosic_FourierTransformerRadix2.h"
#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  /**

  This class implements a convolution via multiplying the spectra of input signal blocks
  with the spectrum of an impulse response. It introduces an input/output delay of
  nextPowerOfTwo(L) where L is the length of the impulse response in samples.

  */

  class ConvolverFFT
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ConvolverFFT();

    /** Destructor. */
    ~ConvolverFFT();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the impulse response to be used. */
    void setImpulseResponse(double *newImpulseResponse, int newLength);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single convolved output-sample. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the impulse response to an impulse (1 0 0 0....). */
    void clearImpulseResponse();

    /** Sets the input-and output buffers to all zeros. */
    void clearBuffers();

    /** Sets the input buffer to all zeros. */
    void clearInputBuffer();

    /** Sets the output buffer 1 to all zeros. */
    void clearOutputBuffer1();

    /** Sets the output buffer 2 to all zeros. */
    void clearOutputBuffer2();

    //=============================================================================================

  protected:

    /** (Re)-allocates memory for the internal buffers if necessarry. */
    void allocateBuffers(int newImpulseResponseLength);

    doubleA  *x;             // buffer for past input samples
    doubleA  *y1, *y2;       // 2 buffers for the output samples (overlap/add)
    doubleA  *h;             // array containing the impulse response (actually redundant with H)
    Complex  *X;             // DFT of the currently processed input block
    Complex  *H;             // DFT of the impulse response (length M)
    int       L;             // actual length of the impulse response
    int       M;             // length of the zero padded impulse response: M = nextPowerOfTwo(2*L)
    int       writeCounter;  // sample-counter for the circular input buffer
    int       readCounter1;  // sample-counter for the first output buffer
    int       readCounter2;  // sample-counter for the second output buffer
    bool      useOutBuffer2; // flag to indicate that xOut2 should be used for output
    MutexLock mutex;         // mutex-lock for accessing the buffers in a thread safe way

    FourierTransformerRadix2 forwardTransformer, inverseTransformer;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double ConvolverFFT::getSample(double in)
  {
    if( x == NULL || h == NULL )
    {
      DEBUG_BREAK;
      return 0.0;
    }

    double out;

    mutex.lock();

    // accept an input sample in the circular input buffer:
    x[writeCounter] = in;

    // if necessarry, compute a new output buffer and update the circular output buffer:
    if( readCounter1 >= M ||  readCounter2 >= M )  //
    {
      double *y;
      if( useOutBuffer2 == true )
      {
        y            = y2;
        readCounter2 = 0;
      }
      else
      {
        y            = y1;
        readCounter1 = 0;
      }

      int k;                        // position in y buffer
      int p = writeCounter-(M/2);   // position in circular input buffer
      if( p < 0 )
        p += M;

      // fill work buffer with input block:
      for(k=0; k<M/2; k++)
      {
        if( p >= M )
          p -= M;
        y[k] = x[p];
        p++;
      }
      for(k=M/2; k<M; k++)
        y[k] = 0.0; // zero padding

      // do the transform / spectral multiplication / inverse transform:
      forwardTransformer.transformRealSignal(y, X);
      X[0].re *= H[0].re; // we must treat the first complex spectral value separately because it's
      X[0].im *= H[0].im; // a pair of real numbers (and does not represent a complex number)
      for(k=1; k<M/2; k++)
        X[k] *= H[k];
      inverseTransformer.transformSymmetricSpectrum(X, y);

      useOutBuffer2 = !useOutBuffer2; // use the other buffer next time
    }

    // compute one output sample by by overlap/adding the 2 output buffers:
    out = y1[readCounter1] + y2[readCounter2];

    readCounter1++;
    readCounter2++;
    writeCounter++;
    if( writeCounter >= M )
      writeCounter = 0;

    mutex.unlock();

    return out;
  }

} // end namespace rosic

#endif // rosic_ConvolverFFT_h
