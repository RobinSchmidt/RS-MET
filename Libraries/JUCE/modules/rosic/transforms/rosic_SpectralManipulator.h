#ifndef rosic_SpectralManipulator_h
#define rosic_SpectralManipulator_h

//// rosic-indcludes:
//#include "../math/rosic_ComplexFunctions.h"

namespace rosic
{

  /**

  This class can apply various manipulations of an FFT spectrum. Most of the manipulations are 
  available as static member functions such that you don't need to instantiate an objectbin order
  to use this class. On the other hand, the non-static manipulations may require some background
  information about the data being processed such as the sample-rate etc. which is stored in member
  variables and thus requires an object to be instatiated.

  */

  class SpectralManipulator  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SpectralManipulator();  

    /** Destructor. */
    ~SpectralManipulator(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample rate at which the spectrum is to be interpreted. */
    void setSampleRate(double newSampleRate); 

    //---------------------------------------------------------------------------------------------
    // static magnitude-spectrum manipulation functions:

    /** Increases the ratio between high and low magnitudes by applying a bin-wise power function
    with an exponent given by the 'power' argument. A value of 1 is neutral, values > 1 emphasize 
    differences in magnitude, values < 1 even out differences in magnitude, a value of 0 whitens 
    the magnitude spectrum and negative values invert the spectrum. */
    static void applyContrast(double *magnitudes, int length, double power);

    /** Sets a slope for the spectrum in dB/oct - negative values will cause the tables to get a 
    lowpassish character, positive value lead to a highpassish character. 
    \todo: introduce a parameter 'cutoffBin' to introduce a knee at which the slope starts to 
    become effective */
    static void applySlope(double *magnitudes, int length, double slope /*,int cutoffBin = -1*/);

    /** Applies a sharp (brickwall) lowpass-filter to the spectrum keeping only magnitudes up to 
    (and including) 'highestBinToKeep' and zeroing out all others. */
    static void applyBrickwallLowpass(double *magnitudes, int length, int highestBinToKeep);

    /** Applies a sharp (brickwall) highpass-filter to the spectrum keeping only magnitudes below
    and including 'lowestBinToKeep' and zeroing out all others. */
    static void applyBrickwallHighpass(double *magnitudes, int length, int lowestBinToKeep);

    /** Adjusts the balance between even and odd harmonics. */
    static void applyEvenOddBalance(double *magnitudes, int length, double balance);

    /** !!! NOT YET IMPLEMENTED !!!
    Multiplies each magnitude bin with a random value between 0...1. The 'amount' parameter
    performs a linear fade between 1 (amount=0) and the random number (amount=0). */
    static void applyMagnitudeRandomization(double *magnitudes, int length, double amount, 
      int seed);

    //.....multiply, envelope (forward/backward attack-release filter), vocoder, warp, 
    // scale, shift, addSpectra, applyMagnitudeRandomization

    //---------------------------------------------------------------------------------------------
    // static phase-spectrum manipulation functions:

    /** Applies a phase offset between even and odd harmonics in degrees. Half of this value will 
    be applied with positive sign to the odd harmonics and with negative sign to the even harmonics
    (leaving DC alone). */
    static void applyEvenOddPhaseShift(double *phases, int length, double shiftInDegrees);

    /** Applies a phase-shift of the the harmonics between left and right channel in degrees. The
    left channel's harmonics will be shifted half this value into the negative direction and the
    right channel's harmonics by the same amount into the positive direction. */
    static void applyStereoPhaseShift(double *phasesL, double *phasesR, int length, 
      double shiftInDegrees);

    /** Applies a stereo phase-shift between even and odd harmonics. Half of the value is applied 
    with positive sign to the left channel's odd and right channel's even harmonics and with 
    negative sign to the left channel's even and right channel's odd harmonics. */
    static void applyEvenOddStereoPhaseShift(double *phasesL, double *phasesR, int length, 
      double shiftInDegrees);

    // applyPhaseShift, applyPhaseRandomization, 

    //=============================================================================================

  protected:

    double sampleRate;

  };

} // end namespace rosic

#endif // rosic_SpectralManipulator_h
