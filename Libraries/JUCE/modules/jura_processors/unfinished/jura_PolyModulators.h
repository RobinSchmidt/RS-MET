
/** A polyphonic modulator section. The user can plug in an arbitrary number of modulator 
modules. */

class JUCE_API PolyModulatorsAudioModule : public jura::AudioModule
{

public:

protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolyModulatorsAudioModule)
};

//=================================================================================================

class JUCE_API PolyModulatorsEditor : public jura::AudioModuleEditor
{

public:

protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolyModulatorsEditor)
};