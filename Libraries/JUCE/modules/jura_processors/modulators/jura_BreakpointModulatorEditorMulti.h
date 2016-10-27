#ifndef rosof_BreakpointModulatorEditorMulti_h
#define rosof_BreakpointModulatorEditorMulti_h

#include "../../rosic_displays/rosof_ModulatorCurveEditorMulti.h"
#include "rosof_BreakpointModulatorEditor.h"

namespace rosof
{

  /**

  This class is...

  \todo: make the RadioGroup-stuff for the Edit-buttons work again

  */

  class BreakpointModulatorEditorMulti : public BreakpointModulatorEditor
  {

  public:


    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    BreakpointModulatorEditorMulti(CriticalSection *newPlugInLock, BreakpointModulatorAudioModule* newBreakpointModulatorAudioModule); 

    /** Destructor. */
    //virtual ~BreakpointModulatorEditorMulti(); 

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Adds a rosic::Breakpointmodulator object to the array of edited modulators. */
    //virtual void addModulatorToEdit(rosic::BreakpointModulator* newModulatorToEdit);
    virtual void addModulatorToEdit(BreakpointModulatorAudioModule* newModulatorToEdit, 
      bool createWidgets = true);

    /** Selects one of the rosic::Breakpointmodulator objects for editing. */
    virtual void selectModulatorToEdit(int index);

    /** Sets the juce::Label in which the descriptions for the widgets will appear. */
    //virtual void setDescriptionField(RLabel* newDescriptionField);

    /** De-selects the currently selected breakpoint (if any) and updates the GUI accordingly. */
    virtual void deSelectBreakpoint();

    /** Sets the label for the modulator that appears top-left to the preset-filenam field. */
    virtual void setModulatorLabel(int index, const juce::String& newLabel);

    /** Sets up the ColourScheme for one of the widget-sets. */
    virtual void setChildColourScheme(int index, const EditorColourScheme& newEditorColourScheme, 
      const WidgetColourScheme& newWidgetColourScheme);
    //virtual void setWidgetSetColourScheme(int index, const WidgetColourScheme& newColourScheme);

    //---------------------------------------------------------------------------------------------
    // callbacks:

    virtual void rButtonClicked(RButton *buttonThatWasClicked);
    virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
    virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);

    virtual void copyColourSettingsFrom(const ColourSchemeComponent *componentToCopyFrom);
    virtual void paint(Graphics &g);
    virtual void resized();

    //---------------------------------------------------------------------------------------------
    // others:

    /** Updates the sliders, buttons, etc. according to the state of the rosic::Modulator object 
    which is being edited. It makes sense to de-select any possibly selected after restoring a 
    state, so this might be done here optionally too. */
    virtual void updateWidgetsAccordingToState(bool deSelectBreakpoint);

    /** Calls updateWidgetsAccordingToState(bool) with true as argument (we need to implement this 
    purely virtual function here). */
    virtual void updateWidgetsAccordingToState();

    //---------------------------------------------------------------------------------------------
    // public data members:

    // widgets:
    juce::Array<BreakpointModulatorGlobalEditor*, CriticalSection> globalEditors;
    //juce::Array<StateLoadSaveWidgetSet*,       CriticalSection> stateWidgetSets;

    // the plot editor:
    ModulatorCurveEditorMulti* breakpointEditorMulti;
    //CoordinateSystemZoomer*    breakpointZoomer;

    //===============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Creates the widgets for a new modulator - called from addModulatorToEdit. */
    virtual void createWidgetsForNewModulator(BreakpointModulatorAudioModule* newModulator);

    /** Automatically adjusts the x-axis plot-range according to the current content. */    
    virtual void autoAdjustPlotRangeX();

    /** Automatically adjusts the y-axis plot-range according to the current content. */
    virtual void autoAdjustPlotRangeY();

    juce::Array<BreakpointModulatorAudioModule*, CriticalSection> modulatorModules;
    int editedModulatorIndex; // -1, if none
    //Rectangle leftSectionRectangle;
    juce::OwnedArray<Rectangle<int> > leftSectionRectangles;

    //juce::Array<BreakpointModulatorWidgetSet*> widgetSets;


  };

}

#endif  