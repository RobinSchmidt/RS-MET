#ifndef rosic_WaveTable_h
#define rosic_WaveTable_h

// rosic-indcludes:
//#include "../datastructures/rosic_Array.h"
#include "../rendering/rosic_WaveformRenderer.h"
#include "../others/rosic_SlewRateLimiter.h"
#include "../transforms/rosic_FourierTransformerRadix2.h"

namespace rosic
{

  //===============================================================================================

  /**

  This is a class for generating and storing a single-cycle-waveform in a lookup-table and 
  retrieving values form it at arbitrary positions by means of interpolation. It doesn't do any
  anti-aliasing, so it's recommended only for use in low-frequency oscillators (LFOs).

  \todo: factor out a WaveTableManipulator class

  */

  class WaveTable
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WaveTable();          

    /** Destructor. */
    ~WaveTable();         

    //---------------------------------------------------------------------------------------------
    // parameter-settings:

    /** Switches automatic updating of the waveform buffers when a parameter changes on/off. The 
    update of the buffer is costly so when you change more than one parameter at a time, it is 
    advisable to turn the automatic updating off, change all parameters as desired and then update
    the buffer only once via updateBuffers(). */
    void setAutoUpdateOnParameterChange(bool shouldUpdateAutomatically)
    { autoUpdate = shouldUpdateAutomatically; }

    // the following waveform manipulations are applied to the waveform in the order of their 
    // declaration

    // time domain manipulations:

    /** Sets a time warping coeffcient 'a' which maps the phase pointer according to the formula:
    ..... - this is a mapping of the domain 0...1 into the range 0...1
    which is the identity mapping for a=0, introduces concavity for a<0 and convexity for a>0.
    The coefficient should be in the range ]-1...+1[ (not including the extremes). */
    void setFullWaveWarp(double newWarpCoefficient)
    { fullWaveWarp = newWarpCoefficient; if(autoUpdate) updateBuffers(); }

    /** Sets the coefficient 'b' in the formula: .... for the time warping of the phase. */
    void setHalfWaveWarp(double newWarpCoefficient)
    { halfWaveWarp = newWarpCoefficient; if(autoUpdate) updateBuffers(); }

    /** Time-reverses the waveform. */
    void setTimeReverse(bool shouldBeReversed)
    { timeReverse = shouldBeReversed; if(autoUpdate) updateBuffers(); }

    // magnitude spectrum manipulations:

    /** Sets the spectral contrast - this contrast value acts as an exponent on the magnitude of
    each harmonic, thus, a value of 1 is neutral, values > 1 emphasize differences in magnitude,
    values < 1 even out differences in magnitude, a value of 0 whitens the magnitude spectrum
    and negative values invert the spectrum. Contrast is applied before all the other spectral
    modifiers (such as slope, even/odd ratio, etc.). */
    void setSpectralContrast(double newContrast)
    { spectralContrast = newContrast; if(autoUpdate) updateBuffers(); }

    /** Sets a slope for the spectrum in dB/oct - negative values will cause a lowpassish 
    character, positive value lead to a highpassish character. */
    void setSpectralSlope(double newSlope)
    { spectralSlope = newSlope; if(autoUpdate) updateBuffers(); }

    /** Adjusts the ratio between even and odd harmonics. */
    void setEvenOddRatio(double newRatio)
    { evenOddRatio = newRatio; if(autoUpdate) updateBuffers(); }

    /** Applies a sharp (brickwall) lowpass-filter to the harmonics of the waveform. */
    void setHighestHarmonicToKeep(int newHighestHarmonicToKeep)
    { highestHarmonicToKeep = newHighestHarmonicToKeep; if(autoUpdate) updateBuffers(); }

    /** Applies a sharp (brickwall) highpass-filter to the harmonics of the waveform. */
    void setLowestHarmonicToKeep(int newLowestHarmonicToKeep)
    { lowestHarmonicToKeep = newLowestHarmonicToKeep; if(autoUpdate) updateBuffers(); }

    // phase spectrum manipulations:

    /** Sets a scaling factor for the harmonis phases. */
    void setPhaseScale(double newPhaseScale)
    { phaseScale = newPhaseScale; if(autoUpdate) updateBuffers(); }

    /** Sets a constant shifting offset for the phases of the harmonics. */
    void setPhaseShift(double newPhaseShift)
    { phaseShift = newPhaseShift; if(autoUpdate) updateBuffers(); }

    /** Sets a phase offset between even and odd harmonics */
    void setEvenOddPhaseShift(double newPhaseShift)
    { evenOddPhaseShift = newPhaseShift; if(autoUpdate) updateBuffers(); }

    /** Sets the phase-shift of the the harmonics between left and right channel in degrees. The
    left channel's harmoincs will be shifted half this value into the negative direction and the
    right channel's harmonics by the same amount into the positive direction. */
    void setStereoPhaseShift(double newPhaseShift)
    { stereoPhaseShift = newPhaseShift; if(autoUpdate) updateBuffers(); }

    /** Sets an additional phase-shift between even and odd harmonics, half of which is applied
    with negative sign to the left channel and with positive sign to the right channel. */
    void setEvenOddStereoPhaseShift(double newPhaseShift)
    { evenOddStereoPhaseShift = newPhaseShift; if(autoUpdate) updateBuffers(); }

    // smoothing/filtering manipulations:

    /** Sets the harmonic number to which the comb filter is tuned. */
    void setCombHarmonic(double newCombHarmonic)
    { combHarmonic = newCombHarmonic; if(autoUpdate) updateBuffers(); }

    /** Sets the amount of comb-filtering in percent -100...100. */
    void setCombAmount(double newCombAmount)
    { combAmount = newCombAmount; if(autoUpdate) updateBuffers(); }

    /** Sets the attack time-constant of the slewrate limiter. The time unit is milliseconds under 
    the assumption that the waveform is played back at a frequency of 1 Hz. */
    void setSmoothAttack(double newAttack)
    { slewRateLimiter.setAttackTime(newAttack); if(autoUpdate) updateBuffers(); }

    /** Sets the release time-constant of the slewrate limiter. The time unit is milliseconds under 
    the assumption that the waveform is played back at a frequency of 1 Hz. */
    void setSmoothRelease(double newRelease)
    { slewRateLimiter.setReleaseTime(newRelease); if(autoUpdate) updateBuffers(); }


    // \todo: (iterated) integration/differencing (with DC-remove and normalization)

    // range manipulations:

    /** Inverts the polarity of the waveform. */
    void setPolarityInversion(bool shouldBeInverted)
    { polarityInvert = shouldBeInverted; if(autoUpdate) updateBuffers(); }

    /** Removes the mean (aka DC offset) from the waveform. */
    void setRemoveMean(bool shouldBeRemoved)
    { meanRemove = shouldBeRemoved; if(autoUpdate) updateBuffers(); }

    /** Normalizes the waveform such that the mximum absolute value is unity. */
    void setNormalize(bool shouldNormalize)
    { normalize = shouldNormalize; if(autoUpdate) updateBuffers(); }

    /** Fits the waveform into the range specified by setRangeMin(), setRangeMax() - if this is
    turned off, the range specified there will be ignored. */
    void setFitToRange(bool shouldFit)
    { fitToRange = shouldFit; if(autoUpdate) updateBuffers(); }

    /** Sets the minimum value of the waveform (i.e. the maximum downward excursion). Only relevant 
    when setFitToRange was (or will be) called with 'true' as argument. */
    void setRangeMin(double newMin)
    { rangeMin = newMin; if(autoUpdate) updateBuffers(); }

    /** Sets the maximum value of the waveform (i.e. the maximum upward excursion). Only relevant 
    when setFitToRange was (or will be) called with 'true' as argument. */
    void setRangeMax(double newMax)
    { rangeMax = newMax; if(autoUpdate) updateBuffers(); }


    // start-phase (should be the last in chain, although it would actually belong into time-domain
    // manipulations:

    /** Sets the start-phase (in degrees) of the (left channel) waveform - the right channel 
    waveform may have n additional offset. */
    void setStartPhase(double newStartPhase)
    { startPhase = newStartPhase; if(autoUpdate) updateBuffers(); }

    /** Sets the phase offset for the right channel with respect to the left channel. */
    void setRightChannelPhaseOffset(double newOffset)
    { rightPhaseOffset = newOffset; if(autoUpdate) updateBuffers(); }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the length of the wavetable (excluding the extra-samples for the interpolator). */
    int getTableLength() const { return tableLength; }

    /** Writes the resulting waveform (with all manipulations applied) into the passed buffers for 
    left and right channel (waveforms are potentially stereo). */
    void getWaveform(double *targetBufferL, double *targetBufferR, int targetBufferLength);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Returns the value at position 'integerPart+fractionalPart' with linear interpolation - this 
    function may be preferred over getValueAt(double phaseIndex) when you want to calculate the 
    integer and fractional part of the phase-index yourself. */
    INLINE double getValueAt(int integerPart, double fractionalPart);

    /** Returns the value at position 'phaseIndex' of table 'tableIndex' with linear 
    interpolation - this function computes the integer and fractional part of the phaseIndex
    internally. */
    INLINE double getValueAt(double phaseIndex);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Updates waveBufferL/R by copying the prototypeBufferL/R into them and applying the desired
    manipulations. */
    void updateBuffers();

    //=============================================================================================

    WaveformRenderer waveformRenderer;

    // to be moved to protected area:
    static const int tableLength = 2048;  // length of the internal buffer

    double prototypeBufferL[tableLength];
    double prototypeBufferR[tableLength];
    double waveBufferL[tableLength+1];     // waveform buffer for left channel
    double waveBufferR[tableLength+1];     // waveform buffer for right channel

  protected:

    // functions to fill table with the built-in waveforms (these functions are
    // called from setWaveform(int newWaveform):
    void fillWithZeros();

    /** Returns a value of the prototype wave at some (possibly non-integer) phase index measured
    in samples (range: 0...prototypeWaveNumSamples). It uses linear interpolation */
    double getPrototypeValueAt(int channel, double phaseIndex);

    /** Warps a phase index measured in samples according to the warping map (input (and output)
    range: 0...prototypeWaveNumSamples). */
    double warpPhaseIndex(double unwarpedIndex);

    // manipulation parameters:
    double startPhase, rightPhaseOffset;
    double fullWaveWarp, halfWaveWarp;
    double spectralContrast, spectralSlope, evenOddRatio;
    double phaseScale, phaseShift, evenOddPhaseShift, stereoPhaseShift, evenOddStereoPhaseShift;
    double combHarmonic, combAmount;
    double rangeMin, rangeMax;
    int    lowestHarmonicToKeep, highestHarmonicToKeep;
    bool   timeReverse;
    bool   polarityInvert, meanRemove, normalize, fitToRange;
    bool   autoUpdate;

    // embedded objects to perform synthesis and manipulations:
    SlewRateLimiter          slewRateLimiter;
    FourierTransformerRadix2 fourierTransformer;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions: 

  INLINE double WaveTable::getValueAt(int i, double f)
  {
    return (1.0-f)*waveBufferL[i] + f*waveBufferL[i+1];
  }

  INLINE double WaveTable::getValueAt(double phaseIndex)
  {
    int    intIndex = floorInt(phaseIndex);
    double frac     = phaseIndex  - (double) intIndex;
    return getValueAt(intIndex, frac);
  }

} // end namespace rosic

#endif // rosic_WaveTable_h
