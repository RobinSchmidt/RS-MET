#ifndef rosic_SpectralProcessor_h
#define rosic_SpectralProcessor_h

namespace rosic
{

/** WARNING: This code is not yet realtime-ready because processBlock allocates. ToDo: keep
the complex buffer as member variable!

This class serves a baseclass for effects that are based on spectral manipulations. */

class SpectralProcessor : public OverlapAddProcessor
{

public:

  //-----------------------------------------------------------------------------------------------
  // \Lifetime

  /** Constructor. You must pass the maximum blocksize here and may optionally pass a maximum
  overlap- and zero-padding factor. For spectral processing with a cosine^2 window on the input,
  an overlap and zero-padding of 2 is usually a good choice. */
  SpectralProcessor(int maxBlockSize, int maxOverlapFactor = 4, int maxPaddingFactor = 4);

  /** Destructor. */
  ~SpectralProcessor();


  //-----------------------------------------------------------------------------------------------
  // \Inquiry

  /** Returns the maximum size of the spectrum */
  int getMaxSpectrumSize() const { return maxBlockSize * maxPaddingFactor / 2; }
   // we divide by two because we use only the positive frequencies

  // write function getSpectrumSize. it should return blockSize * paddingFactor / 2;




  //maxBlockSize * maxPaddingFactor / 2;  

  //===============================================================================================



protected:

  //-----------------------------------------------------------------------------------------------
  // \Processing

  /** Overriden callback from OverlapAddProcessor. We take the block, transform it to the frequency
  domain via FFT and pass it to the virtual processSpectrum() method that the subclass is
  supposed to override. After the subclass has finished processing the spectrum, we transform it 
  back to the time domain via IFFT. Then our baseclass OverlapAddProcessor takes over to take care 
  of the overlap/add business. */
  virtual void processBlock(double* block, int blockSize) override;

  /** This function is to be overriden by your subclass to process the spectrum in some way. The
  passed 'spectrumSize' will be N/2 where N is the FFT-size - thus, the passed spectrum contains
  only the non-redundant (positive) frequencies with bin-indices 0...N/2-1. The purely real
  coefficients for DC and the Nyquist-frequency will be contained in the real and imaginary parts
  of spectrum[0] respectively. */
  virtual void processSpectrum(Complex* /*spectrum*/, int /*spectrumSize*/) {}
  // ToDo: This method allocates! Get rid of that to make the code realtime ready!


  FourierTransformerRadix2 transformer; // Embedded object to take care of the FFT and IFFT

  int     maxSpectrumSize;
  Complex *spectrum = nullptr;

};

} // end namespace rosic

#endif // rosic_SpectralProcessor_h
