#ifndef rosic_StandardWaveformRenderer_h
#define rosic_StandardWaveformRenderer_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This is a class for generating a standard single-cycle-waveform by means of various algorithms.

  */

  class StandardWaveformRenderer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    StandardWaveformRenderer();               

    //---------------------------------------------------------------------------------------------
    // parameter-settings:

    /** Chooses one of the standard waveform @see: rosic::standardWaveforms. */
    void setWaveform(int newWaveform) { waveform = newWaveform; }

    //---------------------------------------------------------------------------------------------
    // waveform rendering:

    /** Renders the waveform into the passed buffer. */
    void renderWaveform(double *targetBuffer, int length);

    /** Renders a sine waveform into the passed buffer. */
    static void renderSineWaveform(    double *buffer, int length);

    /** Renders a sawtooth waveform into the passed buffer. */
    static void renderSawWaveform(     double *buffer, int length);

    /** Renders a square waveform into the passed buffer. */
    static void renderSquareWaveform(  double *buffer, int length);

    /** Renders a triangle waveform into the passed buffer. */
    static void renderTriangleWaveform(double *buffer, int length);

    //=============================================================================================

  protected:

    int waveform;

  };

} // end namespace rosic

#endif 
