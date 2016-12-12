#ifndef jura_PhaseScope_h
#define jura_PhaseScope_h
  

/** Implements the buffering for a phasescope analyzer (maybe move this class to RAPT). */

class JUCE_API PhaseScopeBuffer
{

public:

  /** Constructor. */
  PhaseScopeBuffer();

  /** Destructor. */
  virtual ~PhaseScopeBuffer();

  /** Sets the sample rate. */
  void setSampleRate(double newSampleRate);

  /** Sets the time it takes for "color" to decay away. */
  void setDecayTime(double newDecayTime);

  /** Sets the size of the matrix into which we buffer the incoming samples. This should correspond 
  to the pixel size of the display. */
  void setSize(int newWidth, int newHeight);

  /** Accepts one input sample frame for buffering. */
  void bufferSampleFrame(double left, double right);

  /** Resets the internal buffer to all zeros. */
  void reset();

protected:

  void allocateBuffer();

  void freeBuffer();

  void updateDecayFactor();

  double sampleRate;
  double decayTime;
  float decayFactor;
  float *bufferFlat;
  float **buffer;
  int width, height;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeBuffer)
};

//=================================================================================================

/** Implements a phasescope analyzer */

class JUCE_API PhaseScope : public jura::AudioModule
{

public:

  PhaseScope(CriticalSection *lockToUse);

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

protected:

  PhaseScopeBuffer phaseScopeBuffer;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScope)
};

//=================================================================================================

/** Implements the GUI editor for the phase scope  */

class JUCE_API PhaseScopeEditor : public AudioModuleEditor
{


public:

  PhaseScopeEditor(jura::PhaseScope *newPhaseScopeToEdit);

  virtual void resized() override;


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeEditor)
};



#endif 