#ifndef rosof_PitchShifterModuleEditor_h
#define rosof_PitchShifterModuleEditor_h

#include "rosof_PitchShifterAudioModule.h"
#include "../rosof_AudioModuleEditor.h"

namespace rosof
{

  class PitchShifterModuleEditor : public AudioModuleEditor, public RComboBoxObserver
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    PitchShifterModuleEditor(CriticalSection *newPlugInLock, PitchShifterAudioModule* newPitchShifterAudioModule);

    //---------------------------------------------------------------------------------------------
    // callbacks:

    virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
    virtual void resized();
    virtual void updateWidgetsAccordingToState();

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Makes currently required widgets visible and currently not required widgets invisible. */
    virtual void updateWidgetVisibility();

    PitchShifterAudioModule *pitchShifterModuleToEdit;

    RSlider *coarseSlider, *fineSlider, *grainLengthInMillisecondsSlider, *grainLengthInCyclesSlider,
      *grainLengthInBeatsSlider, *feedbackSlider, *dryWetSlider;
    RComboBox *grainLengthUnitComboBox;
    RButton *invertButton, *reverseButton, *antiAliasButton; // *formantPreserveButton, *monoButton;

  };

}

#endif
