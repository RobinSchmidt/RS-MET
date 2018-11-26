
ContainerModuleEditor::ContainerModuleEditor(LibertyAudioModule *newLiberty, 
  romos::Module* newModuleToEdit)
  : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{

}

//-------------------------------------------------------------------------------------------------

ParameterModuleEditor::ParameterModuleEditor(LibertyAudioModule *newLiberty, 
  romos::Module* newModuleToEdit)
  : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  minValueLabel = new RTextField("Min");
  minValueLabel->setDescription("Minimum value of the parameter");
  addWidget(minValueLabel, true, true);

  minValueField = new LibertyTextEntryField("MinValue");
  minValueField->setDescription(minValueLabel->getDescription());
  minValueField->registerTextEntryFieldObserver(this);
  addWidget(minValueField, true, true);


  maxValueLabel = new RTextField("Max");
  maxValueLabel->setDescription("Maximum value of the parameter");
  addWidget(maxValueLabel, true, true);

  maxValueField = new LibertyTextEntryField("MaxValue");
  maxValueField->setDescription(maxValueLabel->getDescription());
  maxValueField->registerTextEntryFieldObserver(this);
  addWidget(maxValueField, true, true);


  valueSlider = new LibertySlider("Value");
  valueSlider->setSliderName("Value");
  valueSlider->setRange(0.0, 1.0, 0.0, 0.5, true);
  valueSlider->setStringConversionFunction(&valueToString5);
  valueSlider->setDescription("Current value of the parameter");
  valueSlider->addListener(this);
  addWidget(valueSlider, true, true);


  helpTextLabel = new RTextField("Help:");
  helpTextLabel->setDescription("Help text for the parameter");
  addWidget(helpTextLabel, true, true);

  helpTextField = new LibertyTextEntryField("HelpText"); 
  helpTextField->setDescription(helpTextLabel->getDescription());
  helpTextField->registerTextEntryFieldObserver(this);
  addWidget(helpTextField, true, true);


  parameterSetupLabel = new RTextField("Parameter Setup");
  parameterSetupLabel->setDescription("General setup for the parameter");
  parameterSetupLabel->setJustification(juce::Justification::centred);
  addWidget(parameterSetupLabel, true, true);

  valueField = new LibertyLabeledTextEntryField("Value"); 
  valueField->setDescription("Current value of the parameter");
  valueField->setLabelText("Value:");
  valueField->getTextEntryField()->registerTextEntryFieldObserver(this);
  addWidget(valueField, true, true);

  defaultField = new LibertyLabeledTextEntryField("DefaultValue"); 
  defaultField->setDescription("Default value of the parameter");
  defaultField->setLabelText("Default:");
  defaultField->getTextEntryField()->registerTextEntryFieldObserver(this);
  addWidget(defaultField, true, true);

  unitField = new LibertyLabeledTextEntryField("Unit"); 
  unitField->setDescription("Physical unit of the parameter");
  unitField->getTextEntryField()->registerTextEntryFieldObserver(this);
  unitField->setLabelText("Unit:");
  addWidget(unitField, true, true);

  scalingComboBox = new LibertyNamedComboBox("Scaling");
  scalingComboBox->setComboBoxName("Scaling:");
  scalingComboBox->setDescription("Scaling behavior of the parameter");
  scalingComboBox->addItem(romos::ParameterModule::LINEAR_MAPPING,      "Linear",      true, false);
  scalingComboBox->addItem(romos::ParameterModule::EXPONENTIAL_MAPPING, "Exponential", true, false);
  addWidget(scalingComboBox, true, true);
  scalingComboBox->registerComboBoxObserver(this);


  setToDefaultButton = new RClickButton("Use");
  setToDefaultButton->setDescription("Set parameter to default value");
  setToDefaultButton->addRButtonListener(this);
  addWidget(setToDefaultButton, true, true);

  //RTextField            *currentLabel, *defaultLabel, *unitLabel;
  //LibertyTextEntryField *currentField, *defaultField, *unitField;
  //LibertyNamedComboBox  *scalingComboBox;
  //RClickButton          *enterValueButton, *setToDefaultButton;


  //setToMinButton = new RClickButton(juce::String(("Use")));
  //setToMinButton->setDescription(juce::String(("Set parameter to min value")));
  //setToMinButton->addRButtonListener(this);
  //addWidget(setToMinButton, true, true);

  //defaultValueField = new LibertyTextEntryField(juce::String(("DefaultValue")));
  //defaultValueField->setDescription(juce::String(("Default value of the parameter")));
  //defaultValueField->registerTextEntryFieldObserver(this);
  //addWidget(defaultValueField, true, true);

  //defaultValueLabel = new RTextField(juce::String(("Default:")));
  //addWidget(defaultValueLabel, true, true);


  //setToMaxButton = new RClickButton(juce::String(("Use")));
  //setToMaxButton->setDescription(juce::String(("Set parameter to max value")));
  //setToMaxButton->addRButtonListener(this);
  //addWidget(setToMaxButton, true, true);

  //mappingComboBox = new LibertyNamedComboBox(juce::String(("Scaling")));
  //mappingComboBox->setDescription(juce::String(("Scaling behavior of the parameter")));
  //mappingComboBox->addItem(romos::ParameterModule::LINEAR_MAPPING,      "Linear",      true, false);
  //mappingComboBox->addItem(romos::ParameterModule::EXPONENTIAL_MAPPING, "Exponential", true, false);
  //addWidget(mappingComboBox, true, true);
  //mappingComboBox->registerComboBoxObserver(this);

  updateWidgetsFromModuleState(); 
  // ah - we should perhaps have GUI parameters in the parameter module...
}

void ParameterModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  ScopedLock scopedLock(*plugInLock);

  //if( buttonThatWasClicked == setToMinButton )
  //  valueSlider->setValue(valueSlider->getMinimum(), true);
  if( buttonThatWasClicked == setToDefaultButton )
    valueSlider->setToDefaultValue(true);
  //else if( buttonThatWasClicked == setToMaxButton )
  //  valueSlider->setValue(valueSlider->getMaximum(), true);
  else
    ModulePropertiesEditor::rButtonClicked(buttonThatWasClicked);
}

void ParameterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  ModulePropertiesEditor::resized();

  int x, y, w, h;
  x  = 0;
  y  = getHeadlineBottom()+4;
  w  = getWidth();
  h  = 80;

  topSectionRect.setBounds(x, y, w, h);
  y += h-1;
  h  = getHeight() - y - 24;
  setupRect.setBounds(x, y, w/2, h);
  x += w/2 - 1;
  w  = getWidth() - x - 1;
  controlSetupRect.setBounds(x, y, w, h);

  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(topSectionRect);
  guiLayoutRectangles.add(setupRect);
  guiLayoutRectangles.add(controlSetupRect);

  x = topSectionRect.getX();
  y = topSectionRect.getY();

  int ww = 60;  // widget width
  int wh = 20;

  minValueField->setBounds(x+4, y+4, ww, wh);
  x = topSectionRect.getRight() - ww;
  maxValueField->setBounds(x-4, y+4, ww, wh);
  x = minValueField->getRight();
  minValueLabel->setBounds(x, y+4, minValueLabel->getWidthToFitText(), wh);
  ww = maxValueLabel->getWidthToFitText();
  x  = maxValueField->getX() - ww;
  maxValueLabel->setBounds(x+2*RWidget::outlineThickness, y+4, ww, wh);


  x  = topSectionRect.getX();
  y  = minValueLabel->getBottom() - RWidget::outlineThickness;
  wh = 24;
  w  = topSectionRect.getWidth();
  valueSlider->setBounds(x+4, y, w-8, wh);

  y = valueSlider->getBottom();
  helpTextLabel->setBounds(x+4, y+8, helpTextLabel->getWidthToFitText(), 16);
  x = helpTextLabel->getRight();
  helpTextField->setBounds(x-2, y+8, topSectionRect.getWidth()-x-4+2, 16);




  x = setupRect.getX();
  y = setupRect.getY();
  w = setupRect.getWidth();
  h = setupRect.getHeight();

  parameterSetupLabel->setBounds(x, y+4, w, 16);

  y      = parameterSetupLabel->getBottom();
  ww     = w - 40;
  wh     = 16;
  int dy = wh - RWidget::outlineThickness;
  int lw = 60;  // label width 
  valueField->setBounds(x+4, y+4, ww, wh);
  y += dy;
  defaultField->setBounds(x+4, y+4, ww, wh);
  y += dy;
  scalingComboBox->setBounds(x+4, y+4, ww, wh);
  //y += dy;
  //unitField->setBounds(x+4, y+4, ww, wh);

  valueField->setLabelWidth(lw);
  defaultField->setLabelWidth(lw);
  unitField->setLabelWidth(lw);
  scalingComboBox->setNameLabelWidth(lw);

  x  = defaultField->getRight() - RWidget::outlineThickness;
  ww = setupRect.getWidth() - x - 4; 
  y  = defaultField->getY();
  //setToDefaultButton->setBounds(x, y, ww, wh);

  //x = topSectionRect.getX();
  //y = helpTextField->getBottom();

  //topSectionRect, setupRect, controlSetupRect;

  //x  = 32;
  //y  = getHeadlineBottom()+4;
  //w  = getWidth() - x - 8;
  //h  = 24;

  //valueSlider->setBounds(x, y, w, h);
  //x = valueSlider->getX()      - RWidget::outlineThickness;
  //y = valueSlider->getBottom() - RWidget::outlineThickness;

  //int fw = 90; // field-width for min/max/default fields
  //x += RWidget::outlineThickness;
  //h  = 20;
  //minValueLabel->setBounds(x-28, y, 28, h); // dont use magic numbers (28) - make a function RTextField::getOptimalWidth()
  //minValueField->setBounds(x,    y, fw, h);

  //x = (valueSlider->getX() + valueSlider->getRight()) / 2 - fw/2;
  //defaultValueField->setBounds(x,    y, fw, h);
  //defaultValueLabel->setBounds(x-52, y, 52, h);

  //x = valueSlider->getRight() - fw;
  //maxValueLabel->setBounds(x-34, y, 34, h);
  //maxValueField->setBounds(x,    y, fw, h);

  //int bw = 40;  // button width
  //y = minValueField->getBottom() - RWidget::outlineThickness;
  //x = minValueField->getX() + minValueField->getWidth()/2 - bw/2;
  //setToMinButton->setBounds(x, y, bw, h);
  //x = defaultValueField->getX() + defaultValueField->getWidth()/2 - bw/2;
  //setToDefaultButton->setBounds(x, y, bw, h);
  //x = maxValueField->getX() + maxValueField->getWidth()/2 - bw/2;
  //setToMaxButton->setBounds(x, y, bw, h);

  //x = 4;
  //w = getWidth()/4 - 8;
  //y = setToMaxButton->getBottom() + 8;
  //mappingComboBox->setBounds(x, y, w, h);
}

void ParameterModuleEditor::updateWidgetsFromModuleState()
{
  ScopedLock scopedLock(*plugInLock);
  romos::ParameterModule* pm = (romos::ParameterModule*) moduleToEdit;
  if(pm != nullptr)
  {
    double minValue     = pm->getMinValue();
    double maxValue     = pm->getMaxValue();
    double defaultValue = pm->getDefaultValue();
    double value        = pm->getValue();
    double quantization = 0.0;  // preliminary
    valueSlider->setRange(minValue, maxValue, quantization, defaultValue, false);
    valueSlider->setValue(value, false, false);

    int scaling = pm->getMappingFunction();
    if(scaling == romos::ParameterModule::EXPONENTIAL_MAPPING)
      valueSlider->setScaling(Parameter::EXPONENTIAL);
    else
      valueSlider->setScaling(Parameter::LINEAR);

    ModulePropertiesEditor::updateWidgetsFromModuleState();
  }


}

//-------------------------------------------------------------------------------------------------

TopLevelModuleEditor::TopLevelModuleEditor(LibertyAudioModule *newLiberty, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{

}


//-------------------------------------------------------------------------------------------------

VoiceKillerModuleEditor::VoiceKillerModuleEditor(LibertyAudioModule *newLiberty, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{
  /*
  thresholdField = new RLabeledTextEntryField(juce::String(("Threshold:")), juce::String(("0.0001")));
  thresholdField->setDescription(juce::String(("Amplitude threshold below which voice gets killed")));
  //thresholdField->setDescriptionField(descriptionField);
  addWidget(thresholdField, true, true);


  timeOutField = new RLabeledTextEntryField(juce::String(("TimeOut:")), juce::String(("0.01")));
  timeOutField->setDescription(juce::String(("Time until voice gets killed after amplitude falls below threshold")));
  //timeOutField->setDescriptionField(descriptionField);
  addWidget(timeOutField, true, true);
  */

  ScopedLock scopedLock(*plugInLock);

  thresholdSlider = new LibertySlider("Threshold");
  thresholdSlider->setSliderName("Threshold");
  thresholdSlider->setRange(-180.0, -40.0, 1.0, -100.0, true);
  thresholdSlider->setStringConversionFunction(&decibelsToStringWithUnit);
  thresholdSlider->setDescription("Amplitude threshold below which voice gets killed");
  addWidget(thresholdSlider, true, true);
  thresholdSlider->addListener(this);


  timeOutSlider = new LibertySlider("TimeOut");
  timeOutSlider->setSliderName("TimeOut");
  timeOutSlider->setRange(0.01, 1.0, 0.01, 0.01, true);
  timeOutSlider->setStringConversionFunction(&valueToString2);
  timeOutSlider->setScaling(Parameter::EXPONENTIAL);
  timeOutSlider->setDescription("Time until voice gets killed after amplitude falls below threshold");
  addWidget(timeOutSlider, true, true);
  timeOutSlider->addListener(this);

  updateWidgetsFromModuleState();
}

void VoiceKillerModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  ModulePropertiesEditor::resized();

  int x, y, w, h, dy;

  x  = 4;
  y  = getHeadlineBottom()+4;
  w  = getWidth()/4 - 8;
  h  = 16;
  dy = h - RWidget::outlineThickness;

  thresholdSlider->setBounds(x, y, w, h);
  y += dy;
  timeOutSlider->setBounds(x, y, w, h);
}

//-------------------------------------------------------------------------------------------------

WhiteNoiseModuleEditor::WhiteNoiseModuleEditor(LibertyAudioModule *newLiberty, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  seedSlider = new LibertySlider("Seed");
  seedSlider->setSliderName("Seed");
  seedSlider->setRange(0.0, 1000.0, 1.0, 0.0, true);
  seedSlider->setStringConversionFunction(&valueToString);
  seedSlider->setDescription("Seed for the pseudo-random number generator");
  //seedSlider->setDescriptionField(descriptionField);
  addWidget(seedSlider, true, true);
  seedSlider->addListener(this);

  updateWidgetsFromModuleState();
}

void WhiteNoiseModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  ModulePropertiesEditor::resized();

  int x, y, w, h;
  x  = 4;
  y  = getHeadlineBottom()+4;
  w  = getWidth()/4 - 8;
  h  = 16;

  seedSlider->setBounds(x, y, w, h);
}

//-------------------------------------------------------------------------------------------------

BiquadDesignerModuleEditor::BiquadDesignerModuleEditor(LibertyAudioModule *newLiberty, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  typedef romos::BiquadDesigner BQD;

  modeComboBox = new LibertyNamedComboBox("Mode");
  modeComboBox->addItem(BQD::BYPASS,                        "Bypass",                      true, false);
  modeComboBox->addItem(BQD::LOWPASS_6_BILINEAR,            "Lowpass, 6 dB/oct, BLT",      true, false);
  modeComboBox->addItem(BQD::HIGHPASS_6_BILINEAR,           "Highpass, 6 dB/oct, BLT",     true, false);
  modeComboBox->addItem(BQD::LOW_SHELF_1_BILINEAR,          "Low Shelf, 1st order, BLT",   true, false);
  modeComboBox->addItem(BQD::HIGH_SHELF_1_BILINEAR,         "High Shelf, 1st order, BLT",  true, false);
  modeComboBox->addItem(BQD::ALLPASS_1_BILINEAR,            "Allpass, 1st order, BLT",     true, false);
  modeComboBox->addItem(BQD::LOWPASS_12_BILINEAR,           "Lowpass, 12 dB/oct, BLT",     true, false);
  modeComboBox->addItem(BQD::HIGHPASS_12_BILINEAR,          "Highpass, 12 dB/oct, BLT",    true, false);
  modeComboBox->addItem(BQD::BANDPASS_CONST_SKIRT_BILINEAR, "Bandpass, const. skirt, BLT", true, false);
  modeComboBox->addItem(BQD::BANDPASS_CONST_PEAK_BILINEAR,  "Bandpass, const. peak, BLT",  true, false);
  modeComboBox->addItem(BQD::BANDREJECT_BILINEAR,           "Bandreject, BLT",             true, false);
  modeComboBox->addItem(BQD::PEAK_BILINEAR,                 "Peak, BLT",                   true, false);
  modeComboBox->addItem(BQD::LOW_SHELF_2_BILINEAR,          "Low Shelf, 2nd order, BLT",   true, false);
  modeComboBox->addItem(BQD::HIGH_SHELF_2_BILINEAR,         "High Shelf, 2nd order, BLT",  true, false);
  modeComboBox->addItem(BQD::ALLPASS_2_BILINEAR,            "Allpass, 2nd order, BLT",     true, false);
  modeComboBox->setDescription("Mode of the filter to be designed");
  modeComboBox->setComboBoxName("Mode:");
  modeComboBox->setNameLabelWidth(44);
  modeComboBox->registerComboBoxObserver(this);
  addWidget(modeComboBox, true, true);

  updateWidgetsFromModuleState();
}

void BiquadDesignerModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  ModulePropertiesEditor::resized();

  int x, y, w, h;
  x  = 4;
  y  = getHeadlineBottom()+4;
  w  = getWidth()/2 - 8;
  h  = 16;

  modeComboBox->setBounds(x, y, w, h);
}

//-------------------------------------------------------------------------------------------------

LibertyLadderFilterModuleEditor::LibertyLadderFilterModuleEditor(LibertyAudioModule *newLiberty, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  int labelWidth = 72;

  typedef romos::LadderFilter LDR;

  filterModeComboBox = new LibertyNamedComboBox("Mode");
  filterModeComboBox->addItem(LDR::FLAT,    "Flat",                  true, false);
  filterModeComboBox->addItem(LDR::LP_6,    "Lowpass, 6 dB/oct",     true, false);
  filterModeComboBox->addItem(LDR::LP_12,   "Lowpass, 12 dB/oct",    true, false);
  filterModeComboBox->addItem(LDR::LP_18,   "Lowpass, 18 dB/oct",    true, false);
  filterModeComboBox->addItem(LDR::LP_24,   "Lowpass, 24 dB/oct",    true, false);
  filterModeComboBox->addItem(LDR::HP_6,    "Highpass, 6 dB/oct",    true, false);
  filterModeComboBox->addItem(LDR::HP_12,   "Highpass, 12 dB/oct",   true, false);
  filterModeComboBox->addItem(LDR::HP_18,   "Highpass, 18 dB/oct",   true, false);
  filterModeComboBox->addItem(LDR::HP_24,   "Highpass, 24 dB/oct",   true, false);
  filterModeComboBox->addItem(LDR::BP_6_6,  "Bandpass, 6/6 dB/oct",  true, false);
  filterModeComboBox->addItem(LDR::BP_6_12, "Bandpass, 6/12 dB/oct", true, false);
  filterModeComboBox->addItem(LDR::BP_12_6, "Bandpass, 12/6 dB/oct", true, false);
  filterModeComboBox->addItem(LDR::BP_6_18, "Bandpass, 6/18 dB/oct", true, false);
  filterModeComboBox->addItem(LDR::BP_18_6, "Bandpass, 18/6 dB/oct", true, false);
  filterModeComboBox->setDescription("Mode of the filter");
  filterModeComboBox->setComboBoxName("Mode:");
  filterModeComboBox->setNameLabelWidth(labelWidth);
  filterModeComboBox->registerComboBoxObserver(this);
  addWidget(filterModeComboBox, true, true);

  saturationModeComboBox = new LibertyNamedComboBox("SaturationMode");
  saturationModeComboBox->addItem(LDR::NO_SATURATION, "No Saturation", true, false);
  saturationModeComboBox->addItem(LDR::LAST_STAGE,    "Last Stage",    true, false);
  saturationModeComboBox->addItem(LDR::FEEDBACK,      "Feedback",      true, false);
  saturationModeComboBox->addItem(LDR::EACH_STAGE,    "Each Stage",    true, false);
  saturationModeComboBox->setDescription("Point(s) in the filter where saturation is applied");
  saturationModeComboBox->setComboBoxName("Saturation:");
  saturationModeComboBox->setNameLabelWidth(labelWidth);
  saturationModeComboBox->registerComboBoxObserver(this);
  addWidget(saturationModeComboBox, true, true);

  updateWidgetsFromModuleState();
}

void LibertyLadderFilterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  ModulePropertiesEditor::resized();

  int x, y, w, h;
  x  = 4;
  y  = getHeadlineBottom()+4;
  w  = getWidth()/2 - 8;
  h  = 16;
  int dy = h + 4;

  filterModeComboBox->setBounds(x, y, w, h);
  y += dy;
  saturationModeComboBox->setBounds(x, y, w, h);
}

//-------------------------------------------------------------------------------------------------

LibertyFormulaModuleEditor::LibertyFormulaModuleEditor(LibertyAudioModule *newLiberty,
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newLiberty, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  formula1In1OutModule = dynamic_cast<FormulaModule_1_1*> (newModuleToEdit);

  formulaField = new LibertyTextEntryField("Formula");
  formulaField->registerTextEntryFieldObserver(this);
  addWidget(formulaField, true, true);

  updateWidgetsFromModuleState();
}

void LibertyFormulaModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  ModulePropertiesEditor::resized();

  int x, y, w, h;
  x  = 4;
  y  = getHeadlineBottom()+4;
  w  = getWidth()/2 - 8;
  h  = 16;
  int dy = h + 4;

  formulaField->setBounds(x, y, w, h);
  y += dy;
}

void LibertyFormulaModuleEditor::textChanged(RTextEntryField *entryField)
{
  if(entryField == formulaField) {
    std::string newFormula = entryField->getText().toStdString();
    if(formula1In1OutModule->isFormulaValid(newFormula))    {
      formula1In1OutModule->setFormula(newFormula);
      formulaField->markTextAsInvalid(false);
    }
    else
      formulaField->markTextAsInvalid(true);
  }
  else
    ModulePropertiesEditor::textChanged(entryField);
}

void LibertyFormulaModuleEditor::updateWidgetsFromModuleState()
{
  std::string formula = formula1In1OutModule->getFormula();
  formulaField->setText(formula);
}

// todo: 
// -save/recall formula (maybe use rosic::KeyValueMap for that - the romos::Module baseclass
//  should define member functions getState/setState that return/take a KeyValueMap
// -make field bigger and multiline

// Maybe:
// -allow use to set names of in/out variables
// -allow user to set name of the formula field (like "y=sin(x)", for example)

// Allow multiple inputs and outputs:
// -have two additional text entry fields "Inputs", "Outputs" that list the I/O variables
//  for example: Inputs: Freq, Q, Outputs: a1, k for ladder coeff computation