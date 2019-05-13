#pragma once




// maybe rename to AnlogOscArray, MorphOscArray (when we later allow waveshape-morphing)
class JUCE_API BlepOscArrayModule : public jura::AudioModuleWithMidiIn
{

public:

  BlepOscArrayModule(CriticalSection *lockToUse);

  //AudioModuleEditor* createEditor(int type) override;
  // override this later to create a custom editor - for the moment, we use the generic editor

  virtual void setSampleRate(double newSampleRate) override
  {
    oscArrayCore.setSampleRate(newSampleRate);
  }

  virtual void noteOn(int noteNumber, int velocity) override
  {
    oscArrayCore.setFrequency(RAPT::rsPitchToFreq(double(noteNumber)));
    oscArrayCore.reset();
  }

  virtual void processBlock(double **buf, int numChannels, int numSamples) override
  {
    jassert(numChannels == 2);
    for(int n = 0; n < numSamples; n++)
    {
      //buf[0][n] = buf[1][n] = oscArrayCore.getSample();
      oscArrayCore.getSampleFrameStereo(&buf[0][n], &buf[1][n]);
    }
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    //oscArrayCore.getSampleFrameStereo(left, right);
    *left = *right = oscArrayCore.getSample();
  }

  virtual void reset() override
  {
    oscArrayCore.reset();
  }

protected:

  void createParameters();

  rosic::rsOscArrayPolyBlep1 oscArrayCore; 
  // maybe use a pointer to allow wrapping the AudioModule aorund existing cores (like in a synth)

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(BlepOscArrayModule)
};

//=================================================================================================

class JUCE_API BlepOscArrayEditor : public jura::AudioModuleEditor
{

public:

  BlepOscArrayEditor(BlepOscArrayModule* oscArrayToEdit);

  virtual void resized() override;

protected:

  BlepOscArrayModule* oscArrayModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(BlepOscArrayEditor)
};