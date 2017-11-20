#ifndef RPolyphonicInstrumentEditor_h
#define RPolyphonicInstrumentEditor_h

//#include "RAudioProcessor.h"
#include "RPlugInEngineEditor.h"
#include "../rosic_state_management/PolyphonicInstrumentXmlTools.h"

#include "../../rojue/components/widgets/RButton.h"
#include "../../rojue/components/widgets/RSlider.h"
#include "../../rojue/misc/TuningFileManager.h"
#include "../../rojue/components/rosic_editors/PresetRemembererEditor.h"
//#include "../../../rojue/components/rosic_editors/MultiModeFilterEditor.h"
//#include "../../../rojue/components/rosic_editors/RAudioProcessorOscSectionEditor.h"
//#include "../../../rojue/components/rosic_editors/BreakpointModulatorEditorSmall.h"

class RPolyphonicInstrumentEditor : public RPlugInEngineEditor, public LabelListener, 
  public RSliderListener, public ChangeBroadcaster, public TuningFileManager
{

public:

  RPolyphonicInstrumentEditor(rosic::PolyphonicInstrument* newInstrumentToEdit = NULL);
  /**< Constructor. */

  ~RPolyphonicInstrumentEditor();
  /**< Destructor. */

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void buttonClicked(Button *buttonThatWasClicked);
  /**< Implements the purely virtual buttonClicked()-method of the ButtonListener base-class. */

  virtual void changeListenerCallback(void *objectThatHasChanged);
  /**< Implements the purely virtual changeListenerCallback()-method of the 
  ChangeListener base-class. */

  virtual void  labelTextChanged(Label *labelThatHasChanged);
  /**< Implements the purely virtual labelTextChanged()-method of the 
  LabelListener base-class. */

  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  /**< Implements the purely virtual sliderValueChanged()-method of the RSliderListener 
  base-class. */

  //-----------------------------------------------------------------------------------------------
  // state-management:

  virtual XmlElement* getStateAsXml(
    const String& stateName = String(T("RAudioProcessorState"))) const;
  /**< Overrides the method inherited from RobsEditorBase. */

  virtual bool setStateFromXml(const XmlElement& xmlState);
  /**< Overrides the method inherited from RobsEditorBase. */

  //-----------------------------------------------------------------------------------------------
  // others:

  virtual void setInstrumentToEdit(rosic::PolyphonicInstrument* newInstrumentToEdit);
  /**< Attaches this editor to the actual plugin which is to be edited. */

  virtual void resized();
  /**< Overrides the resized-method of the RobsEditorBase base-class. */

protected:

  virtual void updateWidgetsAccordingToState();
  /**< Overrides the inherited method. */

  //virtual void updateTuningWidgets();


  rosic::PolyphonicInstrument* instrumentEngine;

  RSlider *levelSlider, *levelByVelSlider, *numVoicesSlider, *compSlider;

  RLabel     *tuningLabel;
  RTextField *tuningFileNameLabel;  // for the currently loaded tuning-file
  RButton    *tuningMinusButton;
  RButton    *tuningPlusButton;
  RButton    *tuningLoadButton;
  RSlider    *masterTuneSlider;
  RSlider    *wheelRangeSlider;
  RButton    *glideButton;
  RSlider    *glideTimeSlider;

};

#endif
