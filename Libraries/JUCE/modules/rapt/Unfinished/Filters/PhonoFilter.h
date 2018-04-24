#ifndef RS_PHONOFILTER_H
#define RS_PHONOFILTER_H

namespace RSLib
{

  /**

  This is an implementation of a filter that models the RIAA equalization curve for phonograph 
  records. This curve is a (kind of) highpass filter (actually more a shelving filter that 
  boosts high and attenuates low frequencies), called "pre-emphasis" filter. This pre-emphasis 
  filter applied to the audio material before engraving it into the vinyl for various physical 
  reasons. On playback, and the inverse curve (the "de-emphasis" filter) is applied. You can 
  select between pre-/de-emphasis by using setMode().

  The analog transfer function is modeled by matching the magnitude response of two digital 1st 
  order filters to the corresponding 1st order stages in the analog prototype at 3 frequencies, 
  namely: 0 (DC), 1 kHz and sampleRate/3. Then, the 2 digital 1st order sections are consolidated 
  into a single biquad.

  \todo: factor out a general biquad baseclass

  References: 
  (1) http://en.wikipedia.org/wiki/RIAA_equalization
  (2) http://www.hagtech.com/pdf/riaa.pdf 

  */

  class RSLib_API rsPhonoFilter
  {

  public:

    /** This is an enumeration of the available filter modes. */
    enum modes
    {
      PRE_EMPHASIS = 0,
      DE_EMPHASIS
    };


    /** \name Construction/Destruction */

    /** Constructor. */
    rsPhonoFilter();   


    /** \name Setup */

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Chooses the filter mode. See the enumeration for available modes. */
    void setMode(int newMode);


    /** \name Audio Processing */

    /** Calculates a single filtered output-sample. */
    RS_INLINE double getSample(double in);

    /** Processes a block of samples *in and *out are pointers to the input and output block (which
    may or may not be distinct - in the latter case, the content is overwritten). */
    void processBlock(float *in, float *out, int blockSize);


    /** \name Misc */

    /** Resets the internal filter buffers to zero. */
    void reset();

    /** Returns the magnitude response of this filter at the given frequency. */
    double getMagnitudeAt(double frequency);

    /** Returns the magnitude of the analog prototype at some given physical frequency in Hz. 
    Useful for comparing the magnitude response of the digital model to the analog filter. */ 
    static double prototypeMagnitudeAt(double frequency);

  protected:

    /** \name Internal Functions */

    /** Calculates filter coefficients from filter parameters. */
    void calcCoeffs();  

    // magnitude responses of the the two stages of the analog prototype (normalized and 
    // unnormalized with respect to having unity gain at 1 kHz):
    static double unnormalizedMagnitude1(double frequency);
    static double unnormalizedMagnitude2(double frequency);
    static double prototypeMagnitude1(double frequency);
    static double prototypeMagnitude2(double frequency);


    /** \name Data */

    // buffering:
    double x1, x2, y1, y2;

    // filter coefficients:
    double b0, b1, b2; // feedforward coeffs
    double a1, a2;     // feedback coeff

    // filter parameters:
    int    mode;
    double sampleRate; 

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  RS_INLINE double rsPhonoFilter::getSample(double in)
  {
    double y = b0*in + b1*x1 + b2*x2 - a1*y1 - a2*y2;
    y2 = y1;
    y1 = y;
    x2 = x1;
    x1 = in;
    return y;
    // \todo factor out a biquad class, maybe use DF2 instead of DF1 (measure performance)
  }

}

#endif
