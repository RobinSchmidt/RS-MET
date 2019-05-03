#ifndef jura_ColourSchemeSetupDialog_h
#define jura_ColourSchemeSetupDialog_h

//#include "../../rojue/components/dialogs/rojue_RDialogBox.h"
//#include "../../rojue/components/editors/rojue_EditorWithStateFile.h"
//using namespace rojue;

/** This class implements a dialog for setting up the colour-scheme of some Editor object. */

class JUCE_API ColourSchemeSetupDialog : virtual public RDialogBox, 
  virtual public EditorWithStateFile, public RComboBoxObserver, public RSliderListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. The calling component should pass itself here to allow us to display ourselves 
  with the same colorscheme. The second parameter specifies the number of hue-offset that the 
  editor defines (which in turn determines the number of hue-offset sliders that are shown. */
  ColourSchemeSetupDialog(ColourSchemeComponent *owner, int numHueOffsets);

  /** Destructor. Deletes the temporarily stored XmlElement. */
  virtual ~ColourSchemeSetupDialog();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton      *buttonThatWasClicked);
  virtual void rComboBoxChanged(RComboBox  *comboBoxThatHasChanged);
  virtual void rSliderValueChanged(RSlider *rSliderThatHasChanged);
  //virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(juce::ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();

  /** Override to update our colors according to the owner. */
  virtual void setVisible(bool shouldBeVisible);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean);
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  /** Updates the widgets according to the state. */
  virtual void updateWidgetsAccordingToState();

protected:

  // widgets:
  RNamedComboBox *editorAppearanceComboBox, *widgetAppearanceComboBox, 
    *plotAppearanceComboBox, *defaultsComboBox;
  RSlider *hueSlider, *saturationSlider, *gammaSlider; //, *contrastSlider, *gammaSlider;

  juce::Array<RSlider*> hueOffsetSliders;

   // maybe: hueSpreadSlider

  ColourSchemeComponent *ownerComponent;
  XmlElement            *xmlColorsOnOpening;
    // remembers the colorscheme settings as they were when opening the dialog to facilitate 
    // reverting on 'Cancel'

  juce_UseDebuggingNewOperator;
};

#endif
