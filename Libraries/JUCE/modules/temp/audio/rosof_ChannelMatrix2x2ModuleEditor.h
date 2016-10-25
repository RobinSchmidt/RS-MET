#ifndef rosof_ChannelMatrix2x2ModuleEditor_h
#define rosof_ChannelMatrix2x2ModuleEditor_h

#include "rosof_ChannelMatrix2x2AudioModule.h"
#include "../rosof_AudioModuleEditor.h"

namespace rosof
{

  class ChannelMatrix2x2ModuleEditor : public AudioModuleEditor, public RTextEntryFieldObserver
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    ChannelMatrix2x2ModuleEditor(CriticalSection *newPlugInLock, ChannelMatrix2x2AudioModule* newChannelMatrix2x2AudioModule);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // callbacks:

    virtual void textChanged(RTextEntryField *rTextEntryFieldThatHasChanged);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    virtual void resized();

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    virtual void updateWidgetsAccordingToState();

    ChannelMatrix2x2AudioModule *channelMatrix2x2AudioModule;

    // the widgets:
    RTextField      *labelLeftToLeft, *labelRightToLeft, *labelLeftToRight, *labelRightToRight, *leftEquationLabel, *rightEquationLabel;  
    RTextEntryField *editLabelLeftToLeft, *editLabelRightToLeft, *editLabelLeftToRight, *editLabelRightToRight;

  };

}

#endif
