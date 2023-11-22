#ifndef rosic_OverlapAddProcessor_h
#define rosic_OverlapAddProcessor_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

//=================================================================================================
// class OverlapAddProcessor:

/** This class serves a baseclass for all effects that need overlap-add processing. It
encapsulates all the necesarry buffering, and provides the usual per-sample input/output
function. Subclasses just need to override processBlock to do their actual processing.
Furthermore, this baseclass may also apply a window function to the input and/or the output
buffers.   */

class OverlapAddProcessor
{

public:

  //-----------------------------------------------------------------------------------------------
  // \Lifetime

  /** Constructor. You must pass the maximum blocksize here and may optionally pass a maximum
  overlap- and zero-padding factor. For spectral processing with a cosine^2 window on the input,
  an overlap and zero-padding of 2 is usually a good choice. */
  OverlapAddProcessor(int maxBlockSize, int maxOverlapFactor = 4, int maxPaddingFactor = 4);

  /** Destructor. */
  virtual ~OverlapAddProcessor();

  //-----------------------------------------------------------------------------------------------
  // \Setup

  /** Sets the size of the blocks that are taken from the input signal - should be a power of
  two. */
  void setInputBlockSize(int newSize);

  /** Sets the factor by which the hopsize is smaller than the blocksize - should be a power of
  two. */
  void setOverlapFactor(int newFactor);

  /** Sets the factor by which the output buffers are larger than the input buffers, for example
  due to zero-padding to avoid circular convolution in spectral manipulations - should be a power
  of two. */
  void setPaddingFactor(int newFactor);

  /** Selects whether the window function should be applied to the inputput buffers. */
  void setUseInputWindow(bool shouldBeUsed)
  {
    useInputWindow = shouldBeUsed; calculateCompensationGain();
  }

/** Selects whether the window function should be applied to the output buffers. */
  void setUseOutputWindow(bool shouldBeUsed)
  {
    useOutputWindow = shouldBeUsed; calculateCompensationGain();
  }

  /** We use a cos^n window here, where n is some integer that defines the power of the cosine
  that should be used - this power to be used can be set up with this function. Generally, it is
  desirable that the overlapped window functions add up to a constant. For example, when you
  apply the window only the input and use a cos^2 window and an overlap-factor of 2, the
  overlapping windows will add up to unity. If you apply it also to the output signal, you will
  need an overlap-factor of 4 or alternatively, use a cos^1 window, etc. */
  void setWindowPower(int newPower);

  //-----------------------------------------------------------------------------------------------
  // \Inquiry

  //-----------------------------------------------------------------------------------------------
  // \Processing

  /** Accepts and returns one sample at a time. */
  INLINE double getSample(double in);

  /** Resets the object into initial state with respect to the received signal. It doesn't reset
  the parameters. */
  void reset() { clearBuffers(); }
  // Maybe rename to resetState. But that would be bad for consistency of the library.

  //-----------------------------------------------------------------------------------------------
  // \Misc

  /** Clears the content of all buffers. */
  void clearBuffers();

  //===============================================================================================

protected:

  /** This function is the one, you should override in your subclass to do the actual processing.
  The baseclass implementation here does nothing. It is a template method in the sense of the
  "template method" design pattern. */
  virtual void processBlock(double* /*block*/, int /*blockSize*/) {}
  // ToDo: maybe make purely virtual

  /** Prepares the 'tmp' member variable such that is contains the (possibly windowed and padded)
  block to be processed. */
  void prepareBlockForProcessing();

  /** Post-processes a block stored in 'tmp' that is assumed to have just been undergone the
  actual processing that was done by the subclass'es overriden implementation of processBlock().
  This post-processing will consist of copying the data into the appropriate y-buffer and
  possibly applying the output window. */
  void postProcessBlock();
  // ToDo: rename to finishBlockProcessing or addProcessedBlockToOutput

  /** Initializes the internal state (clears buffer, resets read/write pointers, generates
  window, etc.) - all the stuff that is necessary when the blockSize, overlap or zero-padding
  setting changes. */
  void initInternalState();

  /** Initializes the write pointer in the circular input buffer. */
  void initWritePointer();

  /** Initializes the read pointers in the output buffers. */
  void initReadPointers();

  /** Fills the w-buffer with the window function. */
  void generateWindowFunction();

  /** Calculates the gain compensation factor to account for overlapping windows. */
  void calculateCompensationGain();

  // to be moved to protected later:
  /*
  static const int maxBlockSize      = 32;
  static const int maxOverlapFactor  = 4;
  static const int maxPaddingFactor  = 4;
  double x[maxBlockSize];                                                     // circular input buffer
  double tmp[maxBlockSize*maxPaddingFactor];                                  // temporary buffer for processing
  double w[maxBlockSize];                                                     // window function
  int    readPositions[maxPaddingFactor*maxOverlapFactor];                    // read positions in the y-buffers
  double y[maxPaddingFactor*maxOverlapFactor][maxPaddingFactor*maxBlockSize]; // output buffers
  */

  int    maxBlockSize, maxOverlapFactor, maxPaddingFactor;
  double* x;                    // circular input buffer
  double* tmp;                  // temporary buffer for processing
  double* w;                    // window function
  int* readPositions;           // read positions in the y-buffers
  double** y;                   // output buffers
  double gain;
  int    blockSize, overlapFactor, paddingFactor, hopSize;
  bool   useInputWindow, useOutputWindow;
  int    writePosition;                                     // write position in the x-buffer
  int    nextOutBuffer;                                     // next output buffer to write into
  int    windowPower;


};

//-----------------------------------------------------------------------------------------------
// inlined functions:

INLINE double OverlapAddProcessor::getSample(double in)
{
  // check whether a new block should be established and processed, if yes, do it:
  if(writePosition % hopSize == 0)
  {
    prepareBlockForProcessing();
    processBlock(tmp, blockSize*paddingFactor);  // Calls to overriden function of subclass.
    postProcessBlock();
  }

  // accept the input sample in the circular input buffer:
  x[writePosition] = in;

  // overlap/add by accumulating the output buffers:
  double accu = 0.0;
  for(int b=0; b<overlapFactor*paddingFactor; b++)
    accu += y[b][readPositions[b]];

  // increment the read/write pointers:
  writePosition++;
  if(writePosition >= maxBlockSize)
    writePosition = 0;
  for(int b=0; b<overlapFactor*paddingFactor; b++)
  {
    readPositions[b]++; // reset to 0 is done in postProcessBlock()
    if(readPositions[b] >= paddingFactor*blockSize)
      readPositions[b] = 0;
  }

  return gain*accu;
}



// ToDo: 
// -Perhaps we need a mutex-lock to ensure consistency of the internal state...but no, this should
//  be done on a higher level





} // end namespace rosic

#endif // rosic_OverlapAddProcessor_h
