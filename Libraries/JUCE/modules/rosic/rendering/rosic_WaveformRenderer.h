#ifndef rosic_WaveformRenderer_h
#define rosic_WaveformRenderer_h

// rosic-indcludes:
#include "rosic_AlgorithmicWaveformRenderer.h"
#include "rosic_MultiSegmentWaveformRenderer.h"
#include "rosic_StandardWaveformRenderer.h"
#include "rosic_WaveformBuffer.h"

namespace rosic
{

  /**

  This is a class for generating a single-cycle-waveform by different methods.

  \todo: maybe define an abstract WaveformRenderer baseclass for all the different methods with 
  purely virtual renderWaveform and create the concrete renderers on demand (instead of keeping 
  them around all the time)

  */

  class WaveformRenderer
  {

  public:

    /** Enumeration of the modes that can be used to create the waveforms. */
    enum modes
    {
      STANDARD_WAVEFORM,
      AUDIO_FILE,
      ALGORITHM,
      MULTI_SEGMENT,
      TIME_DOMAIN_FORMULA,
      ADDITIVE_FORMULA,

      NUM_MODES
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WaveformRenderer();          

    /** Destructor. */
    ~WaveformRenderer();         

    //---------------------------------------------------------------------------------------------
    // parameter-settings:

    /** Chooses one of the modes for waveform generation. @see: modes. */
    void setMode(int newMode) { mode = newMode; }

    //---------------------------------------------------------------------------------------------
    // waveform rendering:

    /** Renders the waveform into the passed buffer. */
    void renderWaveForm(double *targetBuffer, int length);

    //=============================================================================================

    StandardWaveformRenderer     standardRenderer;
    WaveformBuffer               waveBuffer;
    AlgorithmicWaveformRenderer  algorithmicRenderer;
    MultiSegmentWaveformRenderer multiSegmentRenderer;

  protected:

    int  mode;
    bool dcRemove, normalize, fitToUnitRange;

  };

} // end namespace rosic

#endif 
