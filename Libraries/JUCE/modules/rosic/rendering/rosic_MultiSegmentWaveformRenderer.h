#ifndef rosic_MultiSegmentWaveformRenderer_h
#define rosic_MultiSegmentWaveformRenderer_h

// rosic-indcludes:
#include "../modulators/rosic_BreakpointModulator.h"

namespace rosic
{

  /**

  This is a class for generating a waveform by means of a BreakPointModulator.

  */

  class MultiSegmentWaveformRenderer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    MultiSegmentWaveformRenderer();               

    //---------------------------------------------------------------------------------------------
    // parameter-settings:

    /** Selects whether only the loop portion of the envelope should be used for the waveform. */
    //setUseOnlyLoop(bool shouldUseOnlyLoop);

    //---------------------------------------------------------------------------------------------
    // waveform rendering:

    /** Renders the waveform into the passed buffer. */
    void renderWaveform(double *targetBuffer, int length);

    //=============================================================================================

    BreakpointModulator breakpointModulator;

  protected:

    //bool useOnlyLoop;

  };

} // end namespace rosic

#endif 
