#ifndef rosic_AlgorithmicWaveformRenderer_h
#define rosic_AlgorithmicWaveformRenderer_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /** Returns a function-pointer to the function that generates a particular waveform, based on the
  index of the waveform. For example, if you pass rosic::SAW as argument, you will get a pointer to 
  the function that generates a sawtooth waveform. @see rosic::standardWaveforms */
  UnaryFunctionPointer getWaveformFunction(int waveformIndex);


  //===============================================================================================
  // class ModulationWaveformRenderer

  /**

  This is a baseclass for the different modulation waveform renderers.

  // \todo: make subclasses: (true) FrequenyModulationRenderer, UpperSidebandModulationRenderer, 
  // LowerSidebandModulationRenderer (-> write a hilbertTransform function for buffers)

  */

  class ModulationWaveformRenderer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulationWaveformRenderer();            

    //---------------------------------------------------------------------------------------------
    // parameter-settings:

    /** Selects the waveform for the carrier. @see rosic::standardWaveforms */
    void setCarrierWaveform(int newWaveform) { cWave = newWaveform; }

    /** Selects the waveform for the modulator. @see rosic::standardWaveforms */
    void setModulatorWaveform(int newWaveform) { mWave = newWaveform; }

    /** Sets the relative frequency for the carrier - this is the 'c' in the c/m ratio. */
    void setCarrierRelativeFrequency(int newRelativeFrequency) { cFreq = newRelativeFrequency; }

    /** Sets the relative frequency for the modulator - this is the 'm' in the c/m ratio. */
    void setModulatorRelativeFrequency(int newRelativeFrequency) { mFreq = newRelativeFrequency; }

    /** Sets the start phase for the carrier in degrees. */
    void setCarrierPhase(int newPhase) { cPhase = newPhase; }

    /** Sets the start phase for the modulator in degrees. */
    void setModulatorPhase(int newPhase) { mPhase = newPhase; }

    /** Sets the modulation index. */  
    void setModulationIndex(double newIndex) { modIndex = newIndex; }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Assigns the passed function pointers to the selected waveform functions, wc and wm are 
    assigned to the carrier's and modulator's radian frequency and cp and am are assigned to the
    carrier's and modulator's start phase. */
    void assignRenderingParameters(UnaryFunctionPointer &cFun, UnaryFunctionPointer &mFun, 
      double &wc, double &wm, double &pc, double &pm, int tableLength);

    //=============================================================================================

  protected:

    int    cWave, mWave;
    int    cFreq, mFreq;    // c and m in (c/m)-ratio ->ratio of carrier and modulator frequencies
    double cPhase, mPhase;
    double modIndex;

  };


  //===============================================================================================
  // class AmplitudeModulationWaveformRenderer:

  class AmplitudeModulationWaveformRenderer : public ModulationWaveformRenderer
  {

  public:

    /** Renders the waveform into the passed buffer. */
    void renderWaveform(double *targetBuffer, int length)
    {
      double (*cFun) (double); double (*mFun) (double); double wc, wm, pc, pm; 
      assignRenderingParameters(cFun, mFun, wc, wm, pc, pm, length);
      for(int i=0; i<length; i++)
        targetBuffer[i] = cFun(wc*i+pc) * (1.0+modIndex*mFun(wm*i+pm));
    }

  };


  //===============================================================================================
  // class PhaseModulationWaveformRenderer:

  class PhaseModulationWaveformRenderer : public ModulationWaveformRenderer
  {

  public:

    /** Renders the waveform into the passed buffer. */
    void renderWaveform(double *targetBuffer, int length)
    {
      double (*cFun) (double); double (*mFun) (double); double wc, wm, pc, pm; 
      assignRenderingParameters(cFun, mFun, wc, wm, pc, pm, length);
      for(int i=0; i<length; i++)
        targetBuffer[i] = cFun(wc*i+pc+modIndex*mFun(wm*i+pm));
    }

  };


  //===============================================================================================
  // class RingModulationWaveformRenderer:

  class RingModulationWaveformRenderer : public ModulationWaveformRenderer
  {

  public:

    /** Renders the waveform into the passed buffer. */
    void renderWaveform(double *targetBuffer, int length)
    {
      double (*cFun) (double); double (*mFun) (double); double wc, wm, pc, pm; 
      assignRenderingParameters(cFun, mFun, wc, wm, pc, pm, length);
      for(int i=0; i<length; i++)
        targetBuffer[i] = cFun(wc*i+pc) * mFun(wm*i+pm);
    }

  };


  //===============================================================================================
  // class AlgorithmicWaveformRenderer

  /**

  This is a class for generating a single-cycle-waveform by means of various algorithms.

  */

  class AlgorithmicWaveformRenderer
  {

  public:

    /** Enumeration of the algorithms that can be used to create the waveforms. */
    enum algorithms
    {
      AMPLITUDE_MODULATION,
      PHASE_MODULATION,
      RING_MODULATION,
      WINDOW_FUNCTION,
      OCTAVES_OF_PRIMES,  // maybe better done in the frequency domian

      //.... more to come

      NUM_ALGORITHMS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    AlgorithmicWaveformRenderer();          

    /** Destructor. */
    ~AlgorithmicWaveformRenderer();         

    //---------------------------------------------------------------------------------------------
    // parameter-settings:

    /** Chooses one of the algorithms for wavefrom generation. @see: algorithms. */
    void setAlgorithm(int newAlgorithm) { algorithm = newAlgorithm; }

    //---------------------------------------------------------------------------------------------
    // waveform rendering:

    /** Renders the waveform into the passed buffer. */
    void renderWaveform(double *targetBuffer, int length);

    //=============================================================================================

    AmplitudeModulationWaveformRenderer amplitudeModulationRenderer;
    PhaseModulationWaveformRenderer     phaseModulationRenderer;
    RingModulationWaveformRenderer      ringModulationRenderer;

  protected:

    int algorithm;

  };

} // end namespace rosic

#endif 
