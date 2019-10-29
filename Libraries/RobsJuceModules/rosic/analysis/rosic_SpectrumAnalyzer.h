#ifndef rosic_SpectrumAnalyzer_h
#define rosic_SpectrumAnalyzer_h

//// rosic-indcludes:
//#include "rosic_EnvelopeFollower.h"
//#include "rosic_SignalMeasures.h"
//#include "../transforms/rosic_FourierTransformerRadix2.h"

namespace rosic
{

  /**

  This is a meaurement device which measures certain parameters of an incoming audio signal. It is
  mainly intended for visualization purposes.

  \todo use block-processing, introduce a flag to suspend processing

  */

  class SpectrumAnalyzer
  {

  public:

    enum estimationMethods
    {
      FFT = 1,
      FILTER_BANK,
      OSCILLATOR_BANK,
      LINEAR_PREDICTION
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SpectrumAnalyzer();   

    /** Destructor. */
    ~SpectrumAnalyzer(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Switches the analysis into Mid/Side mode. */
    void setMidSideMode(bool shouldUseMidSideMode) { midSideMode = shouldUseMidSideMode; }

    /** Sets up the block-size (which equals the FFT-size). */
    void setBlockSize(int newBlockSize);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the block-size (which equals the FFT-size). */
    int getBlockSize() const { return blockSize; }

    /** Returns the number of signal channels. */
    int getNumChannels() const { return numChannels; }

    /** Returns the number of non-redundant bins (i.e. FFT-size / 2) */
    int getNumNonRedundantBins() const { return fftSize/2; }

    /** Returns a pointer to an array with the bin-frequencies - the first index indicates the 
    channel, the second index indicates the bin. */
    double* getBinFrequencies() { return &(frequencyBuffer[0]); }

    /** Returns a pointer to the current spectra - the first index indicates the channel, the 
    second index indicates the bin. */
    double** getCurrentSpectra() { return magnitudePointer; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Measures a stereo-ouput frame. */
    INLINE void measureSampleFrameStereo(double* inL, double* inR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Clears all the internal buffers. */
    void clearBuffers();

    /** Causes and update of all the display-related buffers. */
    void updateDisplayBuffers();

    //=============================================================================================

  protected:

    /** Calculates the window function and writes it into the array. */
    void makeWindow();  

    static const int maxNumChannels = 2;
    static const int maxBlockSize   = 65536;
    //static const int maxBlockSize   = 2048;

    double inBuffer[maxNumChannels][maxBlockSize];
    double magnitudeBuffer[maxNumChannels][maxBlockSize];
    double tmpBuffer[maxNumChannels][maxBlockSize];
    double frequencyBuffer[maxBlockSize];
    double windowBuffer[maxBlockSize];


    double sampleRate;
    int    numChannels;
    int    blockSize;
    int    fftSize;
    int    sampleCounter;
    bool   midSideMode;
    RAPT::rsWindowFunction::WindowType windowType;


    double*  magnitudePointer[maxNumChannels];
    double** magnitudePointer2;

    double*  frequencyPointer;


    FourierTransformerRadix2 fourierTransformer;

  };

  //---------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void SpectrumAnalyzer::measureSampleFrameStereo(double* inL,  double* inR)
  {
    // wraparound circular buffer if necesarry:
    while( sampleCounter >= maxBlockSize )
      sampleCounter -= maxBlockSize;

    // establish the 2 input channel signals:
    double tmp1, tmp2;
    if( midSideMode == true )
    {
      tmp1 = SQRT2_INV * (*inL + *inR);
      tmp2 = SQRT2_INV * (*inL - *inR);
    }
    else
    {
      tmp1 = *inL;
      tmp2 = *inR;
    }

    // buffer the incoming samples:
    inBuffer[0][sampleCounter] = tmp1;
    inBuffer[1][sampleCounter] = tmp2;

    // increment sampleCounter:
    sampleCounter++;
  }

} // end namespace rosic

#endif // rosic_SpectrumAnalyzer_h
