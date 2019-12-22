#ifndef rosic_Module_h
#define rosic_Module_h

namespace rosic
{
  //===============================================================================================
  // class Module:

  /**

  This class defines the interface (I/O functions, etc.) for audio-processing modules that
  take a stereo-pair of input signals and produce a stereo-pair of output signals. Each concrete
  effect module should be subclassed from this Module baseclass and implement these I/O
  functions in a suitable manner. The class facilitates the use of the modules in a
  (semi) modular framework such as Quadrifex.

  */

  class Module
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Module() {}

    /** Destructor. */
    virtual ~Module() {}

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Override this to set up the sample-rate. */
    virtual void setSampleRate(double newSampleRate) = 0;

    /** Override this to set up the tempo, if your effect needs this. */
    virtual void setTempoInBPM(double /*newTempo*/) {}

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of modulatable parameters. Override this in your subclass if you want
    to allow some parameter to be modulated over a modulation matrix (or similar routing
    facility). */
    virtual int getNumModulatableParameters() { return modulatableParameters.getNumElements(); }

    ///** Returns the name of one of the modulatable parameters as a c-string. */
    //virtual char* getModulatableParameterName(int index);

    /** Returns the name of one of the modulatable parameters. */
    virtual std::string getModulatableParameterName(int index);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Override this to render one output sample-frame from and input sample-frame at a time. The
    function is supposed to be wrapped into calls to acquireLock/releaseLock to make it
    thread-safe.*/
    virtual void processSampleFrame(double *inOutL, double *inOutR) = 0;

    //---------------------------------------------------------------------------------------------
    // others:

    /** Override this to reset the internal state (buffers, oscillator-phases and such) of the
    effect. */
    virtual void reset() = 0;

    /** Override this to re-trigger oscillators. */
    virtual void trigger() {}

    //=============================================================================================

  protected:

    rsDynamicArray<ModulatableParameter> modulatableParameters;

  };

} // end namespace rosic

#endif // rosic_Module_h
