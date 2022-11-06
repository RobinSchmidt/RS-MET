//#include "rosof_ColourSchemeSetupDialog.h"
//using namespace rosof;

ColourSchemeSetupDialog::ColourSchemeSetupDialog(ColourSchemeComponent *owner, int numHueOffsets)
{
  setHeadlineText(String("Color Setup"));

  okButton->setDescription(String("Apply settings and return"));
  cancelButton->setDescription(String("Return without applying new settings"));

  addWidget( editorAppearanceComboBox =
    new RNamedComboBox(juce::String("Editors"), juce::String("Editors:")) );
  editorAppearanceComboBox->setDescription(
    juce::String("Selects the general appearance of the editors"));
  editorAppearanceComboBox->addItem(ColourScheme::DARK_ON_BRIGHT, "Dark on bright");
  editorAppearanceComboBox->addItem(ColourScheme::BRIGHT_ON_DARK, "Bright on dark");
  editorAppearanceComboBox->selectItemByIndex(0, false, false);
  editorAppearanceComboBox->registerComboBoxObserver(this);

  addWidget( widgetAppearanceComboBox =
    new RNamedComboBox(juce::String("Widgets"), juce::String("Widgets:")) );
  widgetAppearanceComboBox->setDescription(
    juce::String("Selects the general appearance of the widgets"));
  widgetAppearanceComboBox->addItem(ColourScheme::DARK_ON_BRIGHT, "Dark on bright");
  widgetAppearanceComboBox->addItem(ColourScheme::BRIGHT_ON_DARK, "Bright on dark");
  widgetAppearanceComboBox->selectItemByIndex(0, false, false);
  widgetAppearanceComboBox->registerComboBoxObserver(this);

  addWidget( plotAppearanceComboBox =
    new RNamedComboBox(juce::String("Plots"), juce::String("Plots:")) );
  plotAppearanceComboBox->setDescription(
    juce::String("Selects the general appearance of the plots"));
  plotAppearanceComboBox->addItem(ColourScheme::DARK_ON_BRIGHT, "Dark on bright");
  plotAppearanceComboBox->addItem(ColourScheme::BRIGHT_ON_DARK, "Bright on dark");
  plotAppearanceComboBox->selectItemByIndex(0, false, false);
  plotAppearanceComboBox->registerComboBoxObserver(this);

  addWidget( saturationSlider = new RSlider(String("Saturation")) );
  saturationSlider->setSliderName(String("Saturation"));
  saturationSlider->setDescription(String("Adjusts the overall saturation of the colors"));
  saturationSlider->setRange(0.0, 1.0, 0.01, 0.0);
  saturationSlider->addListener(this);

  addWidget( gammaSlider = new RSlider(String("Gamma")) );
  gammaSlider->setSliderName(String("Gamma"));
  gammaSlider->setDescription(String("Adjusts the gamma for all colors"));
  gammaSlider->setRange(0.25, 4.0, 0.01, 1.0);
  gammaSlider->setScaling(Parameter::EXPONENTIAL);
  gammaSlider->addListener(this);

  addWidget( hueSlider = new RSlider(String("Hue")) );
  if( numHueOffsets > 0 )
  {
    hueSlider->setSliderName(String("Central Hue"));
    hueSlider->setDescription(String("Adjusts the central hue"));
  }
  else
  {
    hueSlider->setSliderName(String("Hue"));
    hueSlider->setDescription(String("Adjusts the hue."));
  }
  hueSlider->setRange(0.0, 1.0, 0.01, 0.0);
  hueSlider->addListener(this);

  int height = 200;
  RSlider *hueOffsetSlider;
  for(int i=0; i<numHueOffsets; i++)
  {
    String name = String("HueOffset") + String(i+1);
    addWidget( hueOffsetSlider = new RSlider(name) );
    hueOffsetSlider->setSliderName(String("Hue Offset ") + String(i+1));
    hueOffsetSlider->setDescription(
      String("Adjusts the hue offset for a certain part of the GUI"));
    hueOffsetSlider->setRange(0.0, 1.0, 0.01, 0.0);
    hueOffsetSlider->addListener(this);
    hueOffsetSliders.add(hueOffsetSlider);
    height += 14;
  }

  StateFileManager::setActiveDirectory(
    getSupportDirectory() + File::getSeparatorString() + String("ColorSchemes") );

  ownerComponent     = owner;
  xmlColorsOnOpening = NULL;

  setSize(300, height);

  // make more space for hue-offset sliders
}

ColourSchemeSetupDialog::~ColourSchemeSetupDialog()
{
  hueOffsetSliders.clear();
  delete xmlColorsOnOpening;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void ColourSchemeSetupDialog::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == cancelButton )
    setColourSchemeFromXml(*xmlColorsOnOpening);
  RDialogBox::rButtonClicked(buttonThatWasClicked);
}

void ColourSchemeSetupDialog::rComboBoxChanged(RComboBox* comboBoxThatHasChanged)
{
  if( comboBoxThatHasChanged == editorAppearanceComboBox )
    setEditorAppearance(editorAppearanceComboBox->getSelectedItemIdentifier());
  else if( comboBoxThatHasChanged == widgetAppearanceComboBox )
    setWidgetAppearance(widgetAppearanceComboBox->getSelectedItemIdentifier());
  else if( comboBoxThatHasChanged == plotAppearanceComboBox )
    setPlotAppearance(plotAppearanceComboBox->getSelectedItemIdentifier());

  sendChangeNotification();
}

void ColourSchemeSetupDialog::rSliderValueChanged(RSlider *rSliderThatHasChanged)
{
  if( rSliderThatHasChanged == hueSlider )
    setCentralHue( (float) hueSlider->getValue() );
  else if( rSliderThatHasChanged == saturationSlider )
    setSaturationMultiplier( (float) saturationSlider->getValue() );
  else if( rSliderThatHasChanged == gammaSlider )
    setBrightnessGamma( (float) gammaSlider->getValue() );
  else
  {
    for(int i=0; i<hueOffsetSliders.size(); i++)
    {
      if( rSliderThatHasChanged == hueOffsetSliders[i] )
        setHueOffset(i, (float) hueOffsetSliders[i]->getValue());
    }
  }

  sendChangeNotification();
}

void ColourSchemeSetupDialog::resized()
{
  RDialogBox::resized();
  EditorWithStateFile::resized();

  int x = 0;
  int y = stateWidgetSet->getBottom();;
  int w = getWidth();
  //int h = getHeight();

  editorAppearanceComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  widgetAppearanceComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  plotAppearanceComboBox->setBounds(  x+4, y+4, w-8, 16);
  y += 20;
  gammaSlider->setBounds(             x+4, y+4, w-8, 16);
  y += 20;
  saturationSlider->setBounds(        x+4, y+4, w-8, 16);
  y += 20;
  hueSlider->setBounds(               x+4, y+4, w-8, 16);

  y = hueSlider->getBottom()-2;
  for(int i=0; i<hueOffsetSliders.size(); i++)
  {
    hueOffsetSliders[i]->setBounds(x+4, y, w-8, 16);
    y += 14;
  }

  // make all the name-label widths in the comboboxes equal (using the maximum):
  w = jmax(editorAppearanceComboBox->getNameLabelWidth(),
    widgetAppearanceComboBox->getNameLabelWidth());
  w = jmax(w, plotAppearanceComboBox->getNameLabelWidth() );
  editorAppearanceComboBox->setNameLabelWidth(w);
  widgetAppearanceComboBox->setNameLabelWidth(w);
  plotAppearanceComboBox->setNameLabelWidth(w);
}

void ColourSchemeSetupDialog::setVisible(bool shouldBeVisible)
{
  if( shouldBeVisible == true )
  {
    copyColourSettingsFrom(ownerComponent);
    if( xmlColorsOnOpening != NULL )
      delete xmlColorsOnOpening;
    xmlColorsOnOpening = getColourSchemeAsXml();
    EditorWithStateFile::setStateName(String());
    updateWidgetsAccordingToState();
  }
  RDialogBox::setVisible(shouldBeVisible);
}

void ColourSchemeSetupDialog::changeListenerCallback(juce::ChangeBroadcaster *objectThatHasChanged)
{
  EditorWithStateFile::changeListenerCallback(objectThatHasChanged);
  if( objectThatHasChanged == stateWidgetSet )
  {
    updateWidgetsAccordingToState();
    sendChangeNotification();
  }
}

void ColourSchemeSetupDialog::updateWidgetsAccordingToState()
{
  editorAppearanceComboBox->selectItemByIndex(getEditorColourScheme().getAppearance(), false, false);
  widgetAppearanceComboBox->selectItemByIndex(getWidgetColourScheme().getAppearance(), false, false);
  plotAppearanceComboBox->selectItemByIndex(  getPlotColourScheme().getAppearance(),   false, false);

  saturationSlider->setValue(getEditorColourScheme().getSaturationMultiplier(), false);
  gammaSlider->setValue(     getEditorColourScheme().getBrightnessGamma(),      false);
  hueSlider->setValue(       getEditorColourScheme().getCentralHue(),           false);

  for(int i=0; i<hueOffsetSliders.size(); i++)
    hueOffsetSliders[i]->setValue(getEditorColourScheme().getHueOffset(i), false);
}

void ColourSchemeSetupDialog::setStateFromXml(const XmlElement& xmlState,
  const juce::String& stateName, bool markAsClean)
{
  ColourSchemeComponent::setColourSchemeFromXml(xmlState);
}

XmlElement* ColourSchemeSetupDialog::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  return ColourSchemeComponent::getColourSchemeAsXml();
}



/*

ToDo:
-Actually, it does't make much sense to have comboboxes to switch between the 2 options 
 dark-on-bright and broght-on-dark. Just use a button for each. Maybe just call it "Dark" which
 stands for "Bright on Dark", i.e. it means the background color as that is the dominant one.
-Maybe we could have brightness and contrast sliders as well. Maybe a power rule and some sort of
 sigmoid rule can be used? But which should be first? Or maybe use power -> sigmoid -> power
 ...or wait - no - power is the gamma. brightness would be an additive constant. maybe do 
 add -> contrast -> add, i.e. a sigmoid with pre- and post shift

*/