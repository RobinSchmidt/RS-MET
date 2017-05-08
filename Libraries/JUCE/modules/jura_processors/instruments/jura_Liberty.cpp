//=================================================================================================
// class LibertyInterfaceState:

LibertyInterfaceState::LibertyInterfaceState()
{
  activePanel = STRUCTURE_PANEL;
}

//=================================================================================================
// class LibertyAudioModule:

LibertyAudioModule::LibertyAudioModule(CriticalSection *newPlugInLock, 
  romos::Liberty *modularSynthToWrap)
//: PolyphonicInstrumentAudioModule(newPlugInLock, modularSynthToWrap)
: AudioModule(newPlugInLock)
{
  jassert(modularSynthToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedLiberty       = modularSynthToWrap;
  //underlyingRosicInstrument = modularSynthToWrap;
  setModuleName(juce::String(("Liberty")));
  setActiveDirectory(getApplicationDirectory() + juce::String(("/LibertyPresets")) );
  macroDirectory = getApplicationDirectory() + juce::String(("/LibertyMacros")) ;
}

//-------------------------------------------------------------------------------------------------
// persistence:

void LibertyAudioModule::writeModuleTypeSpecificStateDataToXml(romos::Module *module, 
  XmlElement* xmlState)
{
  romos::ParameterMixIn *m = dynamic_cast<romos::ParameterMixIn*> (module);
  if( m != NULL )
  {
    for(int i = 0; i < m->getNumParameters(); i++)
      xmlState->setAttribute(rosicToJuce(m->getParameterName(i)), 
        rosicToJuce(m->getParameterValue(i)));
  }
}

void LibertyAudioModule::restoreModuleTypeSpecificStateDataFromXml(romos::Module *module, 
  const XmlElement& xmlState)
{
  romos::ParameterMixIn *m = dynamic_cast<romos::ParameterMixIn*> (module);
  if( m != NULL )
  {
    if( module->getTypeIdentifier() == romos::ModuleTypeRegistry::PARAMETER )
    {
      // we need special treatment of the parameter module - minValue/maxValue shlould be set 
      // simultaneuously, because otherwise, min/max might not be recalled correctly, for example 
      // when newMin > oldMax, the module will refuse to use the new minimum, etc. moreover, 
      // min/max must be set up before attempting to set the default and current value:
      double min = xmlState.getDoubleAttribute(juce::String(("MinValue")), 0.0);
      double max = xmlState.getDoubleAttribute(juce::String(("MaxValue")), 1.0);

      juce::String mappingString = xmlState.getStringAttribute(juce::String(("MaxValue")), 
        juce::String(("Linear")));
      int mappingIndex = romos::ParameterModule::LINEAR_MAPPING;
      if( mappingString == juce::String(("Exponential")) )
        mappingIndex = romos::ParameterModule::EXPONENTIAL_MAPPING;

      ((romos::ParameterModule*) m)->setMinMaxAndMapping(min, max, mappingIndex); 
    }
    for(int i = 0; i < m->getNumParameters(); i++)
    {
      juce::String name         = rosicToJuce(m->getParameterName(i));
      juce::String defaultValue = rosicToJuce(m->getParameterDefaultValue(i));
      juce::String storedValue  = xmlState.getStringAttribute(name, defaultValue);
      m->setParameter(i, juceToRosic(storedValue));
    }
  }
}

XmlElement* LibertyAudioModule::getModuleStateAsXml(romos::Module *module, 
  bool withExternalConnections)
{
  //ScopedLock scopedLock(*plugInLock); // this function is static

  XmlElement *xmlState = new XmlElement(rosicToJuce(module->getTypeName()));

  xmlState->setAttribute(juce::String(("Name")), rosicToJuce(module->getName())    );
  xmlState->setAttribute(juce::String(("X")),    juce::String(module->getPositionX())  );
  xmlState->setAttribute(juce::String(("Y")),    juce::String(module->getPositionY())  );
  xmlState->setAttribute(juce::String(("Poly")), juce::String(module->isPolyphonic())  );

  unsigned int i;
  romos::ModuleContainer *container = dynamic_cast<romos::ModuleContainer*> (module);
  if( container != NULL )
  {
    for(i=0; i<container->getNumChildModules(); i++)
      xmlState->addChildElement( getModuleStateAsXml(container->getChildModule(i), true) );
  }
  else
    writeModuleTypeSpecificStateDataToXml(module, xmlState);

  // store audio connections:

  // input modules have invisible "ghost" connections to the source that goes into the 
  // respective container - we don't want to store these:
  if( !module->isInputModule() )  
  {
    if( withExternalConnections == true )
    {
      romos::AudioConnection ac; // connection under investigation
      juce::String connectionString;
      std::vector<romos::AudioConnection> incomingConnections = 
        module->getIncomingAudioConnections();
      for(i = 0; i < incomingConnections.size(); i++)
      {
        ac  = incomingConnections[i];
        connectionString +=   juce::String(("A_")) 
          + juce::String(ac.getSourceModule()->getIndexWithinParentModule()) + juce::String("_") 
          + juce::String(ac.getSourceOutputIndex())                          + juce::String("_") 
          + juce::String(ac.getTargetInputIndex())                           + juce::String(", ");
      }

      // \todo add event connections to the string..

      if( connectionString != juce::String::empty )
      {
        connectionString = connectionString.trimCharactersAtEnd(juce::String(", "));
        xmlState->setAttribute(juce::String(("IncomingConnections")), connectionString);
      }
    }
  }

  return xmlState;
}

void LibertyAudioModule::createAndSetupEmbeddedModulesFromXml(const XmlElement& xmlState, 
  romos::Module *module)
{
  if( module == NULL )
  {
    jassertfalse;
    return;
  }

  module->setModuleName( juceToRosic(xmlState.getStringAttribute(juce::String(("Name")), 
    juce::String(("NoName"))))         );
  if( !module->isTopLevelModule() )
    module->setPolyphonic( xmlState.getBoolAttribute(juce::String(("Poly")), false) );

  // this call is needed only to recall the position of the top-level module's I/O modules 
  // - the other modules are already inserted at their right positions inside the recursive call:
  //module->setPositionXY( xmlState.getIntAttribute(juce::String(("X")), 0), 
  //                       xmlState.getIntAttribute(juce::String(("Y")), 0), false);
  //...mmhh. does not work

  // module->setPositionXY not needed because we assign x, y in the loop and the toplevel module 
  // should always reamain at (0,0) 

  if(  module->getTypeIdentifier() == romos::ModuleTypeRegistry::CONTAINER 
    || module->getTypeIdentifier() == romos::ModuleTypeRegistry::TOP_LEVEL_MODULE )
  {
    romos::ModuleContainer *container = dynamic_cast<romos::ModuleContainer*> (module);
    for(int i=0; i<xmlState.getNumChildElements(); i++)
    {
      XmlElement* childState = xmlState.getChildElement(i);
      rosic::String moduleTypeName = juceToRosic(childState->getTagName());
      int typeIdentifier = romos::ModuleTypeRegistry::getSoleInstance()
        ->getModuleIdentifierFromTypeString(moduleTypeName);

      if( typeIdentifier != romos::ModuleTypeRegistry::UNKNOWN_MODULE_TYPE )
      {
        if(  module->getTypeIdentifier() == romos::ModuleTypeRegistry::TOP_LEVEL_MODULE 
          && romos::ModuleTypeRegistry::isIdentifierInputOrOutput(typeIdentifier) )
        {
          // do nothing when this is the top-level module and the to-be-added child is an I/O module
        }
        else
        {
          rosic::String moduleName = juceToRosic(childState
            ->getStringAttribute(juce::String("Name")));
          int x = childState->getIntAttribute(juce::String("X"), 0);
          int y = childState->getIntAttribute(juce::String("Y"), 0);
          bool poly = childState->getBoolAttribute(juce::String("Poly"), false);
          romos::Module *newModule = container->addChildModule(typeIdentifier, moduleName, 
            x, y, poly, false);
          createAndSetupEmbeddedModulesFromXml(*childState, newModule);
        }
      }
    }
    container->sortChildModuleArray();
  }
  else
  {
    restoreModuleTypeSpecificStateDataFromXml(module, xmlState);
  }  
}

void LibertyAudioModule::createConnectionsFromXml(const XmlElement& xmlState, romos::Module *module)
{
  if( module == NULL )
  {
    jassertfalse;
    return;
  }

  juce::String connectionString = xmlState.getStringAttribute(juce::String("IncomingConnections"), 
    juce::String::empty);

  if( connectionString != juce::String::empty )
  {
    connectionString = connectionString.removeCharacters(juce::String((" ")));
    juce::String remainingString = connectionString + juce::String((","));
    juce::String currentString;
    juce::String tmpString1, tmpString2;
    int smi, spi, tpi; // source module- and pin-index, target pin index
    while( remainingString != juce::String::empty )
    {
      currentString   = remainingString.upToFirstOccurrenceOf((","), false, false);
      remainingString = remainingString.fromFirstOccurrenceOf((","), false, false);
      tmpString1      = currentString;
      int connectionKind;

      if( tmpString1.startsWithChar('A') )
        connectionKind = romos::AUDIO;
      else if( tmpString1.startsWithChar('E') )
        connectionKind = romos::EVENT;
      else
      {
        jassertfalse;  // unknown connection type
        connectionKind = romos::AUDIO;
      }
   
      tmpString1 = tmpString1.fromFirstOccurrenceOf(("_"), false, false); 
      tmpString2 = tmpString1.upToFirstOccurrenceOf(("_"), false, false);
      smi        = tmpString2.getIntValue();
      tmpString1 = tmpString1.fromFirstOccurrenceOf(("_"), false, false); 
      tmpString2 = tmpString1.upToFirstOccurrenceOf(("_"), false, false);
      spi        = tmpString2.getIntValue();
      tmpString1 = tmpString1.fromFirstOccurrenceOf(("_"), false, false); 
      tmpString2 = tmpString1.upToFirstOccurrenceOf(("_"), false, false);
      tpi        = tmpString2.getIntValue();

      romos::ModuleContainer *parentModule = module->getParentModule();
      if( parentModule != NULL ) // maybe NULL when this is called from ModularBlockDiagramPanel::openContainerSaveDialog
      {
        romos::Module *sourceModule = parentModule->getChildModule(smi);
        if( connectionKind == romos::AUDIO )
          parentModule->addAudioConnection(sourceModule, spi, module, tpi);
        else if( connectionKind == romos::EVENT )
        {
          // \todo add the event connection here
        }
      }
    }
  }

  // recursion:
  romos::ModuleContainer *container = dynamic_cast<romos::ModuleContainer*> (module);
  if( container != NULL )
  {
    for(unsigned int i=0; i<container->getNumChildModules(); i++)
    {
      XmlElement    *childState  = xmlState.getChildElement(i);
      romos::Module *childModule = container->getChildModule(i);
      if( childState != NULL )
        createConnectionsFromXml(*childState, childModule);
      else
        romos::triggerRuntimeError("Number of embedded modules does not match number of child elements \
                                   in LibertyAudioModule::createConnectionsFromXml");
        //jassertfalse; // number of embedded modules does not match number of child elements
    }
  }
}

void LibertyAudioModule::setModuleStateFromXml(const XmlElement& xmlState, romos::Module *module)
{
  createAndSetupEmbeddedModulesFromXml(xmlState, module);
  createConnectionsFromXml(xmlState, module);    
}

XmlElement* LibertyAudioModule::getStateAsXml(const juce::String &stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  //XmlElement* xmlState = PolyphonicInstrumentAudioModule::getStateAsXml(stateName, false);
  XmlElement* xmlState = AudioModule::getStateAsXml(stateName, false);
  romos::Module *topLevelMasterVoice = wrappedLiberty->getTopLevelModule();
  xmlState->addChildElement( getModuleStateAsXml(topLevelMasterVoice, true) );
  return xmlState;
}
    
void LibertyAudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  wrappedLiberty->reset();

  romos::TopLevelModule *topLevelModule = wrappedLiberty->getTopLevelModule();
  topLevelModule->deleteAllChildModules();       
  topLevelModule->disconnectAudioOutputModules();

  XmlElement* topLevelModuleState = xmlState.getChildByName(juce::String(("TopLevelModule")));
  if( topLevelModuleState != NULL )
    setModuleStateFromXml(*topLevelModuleState, topLevelModule);
  restoreTopLevelInOutStates(*topLevelModuleState);
  //PolyphonicInstrumentAudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

void LibertyAudioModule::restoreTopLevelInOutStates(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*lock);
  romos::ModuleContainer *topLevelMasterVoice = wrappedLiberty->getTopLevelModule();

  XmlElement *childElement;
  int numInsRestored  = 0;
  int numOutsRestored = 0;

  int x, y;
  for(int i=0; i<xmlState.getNumChildElements(); i++)
  {
    childElement = xmlState.getChildElement(i);
    if( childElement->hasTagName(juce::String(("AudioInput"))) )
    {
      x = childElement->getIntAttribute(juce::String(("X")), 2);
      y = childElement->getIntAttribute(juce::String(("Y")), 2+2*numInsRestored);
      topLevelMasterVoice->getAudioInputModule(numInsRestored)->setPositionXY(
        childElement->getIntAttribute(juce::String(("X")), 2), 
        childElement->getIntAttribute(juce::String(("Y")), 2+2*numInsRestored), false);
      numInsRestored++;
    }
    else if( childElement->hasTagName(juce::String(("AudioOutput"))) )
    {
      topLevelMasterVoice->getAudioOutputModule(numOutsRestored)->setPositionXY(
        childElement->getIntAttribute(juce::String(("X")), 2), 
        childElement->getIntAttribute(juce::String(("Y")), 2+2*numOutsRestored), false);
      numOutsRestored++;
    }
    if( numInsRestored >= 2 && numOutsRestored >= 2 )
      return;
  }
}

//-------------------------------------------------------------------------------------------------
// event handling:
    
void LibertyAudioModule::noteOn(int noteNumber, int velocity)
{
  wrappedLiberty->noteOn(noteNumber, velocity);
}
   
void LibertyAudioModule::noteOff(int noteNumber)
{
  wrappedLiberty->noteOff(noteNumber);
}


//=================================================================================================

//=================================================================================================
// class ModulePropertiesEditor:

ModulePropertiesEditor::ModulePropertiesEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit)
{
  plugInLock   = newPlugInLock;
  moduleToEdit = newModuleToEdit;
  ScopedLock scopedLock(*plugInLock);

  setHeadlineStyle(Editor::SUB_HEADLINE);
  setHeadlineText(rosicToJuce(moduleToEdit->getName()));

  moduleTypeLabel = new RTextField(juce::String("Type:"));
  moduleTypeLabel->setDescription(juce::String("Type of the module"));
  addWidget(moduleTypeLabel, true, true);

  moduleTypeField = new RTextField(rosicToJuce(moduleToEdit->getTypeName()));
  moduleTypeField->setJustification(Justification::centredLeft);
  moduleTypeField->setDescription(moduleTypeLabel->getDescription());
  addWidget(moduleTypeField, true, true);


  polyButton = new RButton(juce::String("Poly"));
  polyButton->setDescription(juce::String("Switch between polyphonic/monophonic mode"));
  polyButton->setToggleState(moduleToEdit->isPolyphonic(), false);
  addWidget(polyButton, true, true);
}

//ModulePropertiesEditor::~ModulePropertiesEditor()
//{
//
//}

void ModulePropertiesEditor::rSliderValueChanged(RSlider* rSlider)
{
  ScopedLock scopedLock(*plugInLock);
  widgetChanged(rSlider);
}

void ModulePropertiesEditor::rComboBoxChanged(RComboBox* comboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  widgetChanged(comboBoxThatHasChanged);
}

void ModulePropertiesEditor::textChanged(RTextEntryField *rTextEntryFieldThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  widgetChanged(rTextEntryFieldThatHasChanged);
}

void ModulePropertiesEditor::rButtonClicked(RButton *buttonThatWasClicked)
{  
  ScopedLock scopedLock(*plugInLock);
  widgetChanged(buttonThatWasClicked);
}

void ModulePropertiesEditor::widgetChanged(RWidget *widgetThatHasChanged)
{  
  ScopedLock scopedLock(*plugInLock);
  romos::ParameterMixIn *m = dynamic_cast<romos::ParameterMixIn*> (moduleToEdit);
  if( m !=  NULL )
  {
    LibertyModuleWidget *lmw = dynamic_cast<LibertyModuleWidget*> (widgetThatHasChanged);
    if( lmw != NULL )
    {
      rosic::String name  = juceToRosic(lmw->getWidgetParameterName());
      rosic::String value = juceToRosic(widgetThatHasChanged->getStateAsString());
      m->setParameter(name, value);
    }
  }
  updateWidgetsFromModuleState(); // to update the other widgets, implements widget-interdependence
}

void ModulePropertiesEditor::resized()
{  
  ScopedLock scopedLock(*plugInLock);
  int x, y, w, h;

  x = 0;
  //y = getHeadlineBottom()+4;
  y = getHeight()-20;
  h = 16;

  moduleTypeLabel->setBounds(x, y, 36, 16); 
  x = moduleTypeLabel->getRight();
  w = getWidth()-x;
  moduleTypeField->setBounds(x, y, w, 16);

  w = 32;
  x = getWidth() - w - 4;
  polyButton->setBounds(x, y, w, 16);

}

void ModulePropertiesEditor::updateWidgetsFromModuleState()
{
  ScopedLock scopedLock(*plugInLock);
  romos::ParameterMixIn *m = dynamic_cast<romos::ParameterMixIn*> (moduleToEdit);
  if( m != NULL )
  {
    for(int i = 0; i < widgets.size(); i++)
    {
      RWidget *widget = widgets[i]; // for debug
      LibertyModuleWidget *lmw = dynamic_cast<LibertyModuleWidget*> (widgets[i]);

      if( lmw != NULL )
      {
        rosic::String name  = juceToRosic(lmw->getWidgetParameterName());
        juce::String  value = rosicToJuce(m->getParameterValue(name));
        widgets[i]->setStateFromString(value, false);
      }
    }
  }
}

//=================================================================================================
//=================================================================================================
// subclasses for particular module types:

ContainerModuleEditor::ContainerModuleEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit)
  : ModulePropertiesEditor(newPlugInLock, newModuleToEdit)
{

}

//-------------------------------------------------------------------------------------------------

ParameterModuleEditor::ParameterModuleEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit)
  : ModulePropertiesEditor(newPlugInLock, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  minValueLabel = new RTextField(juce::String(("Min")));
  minValueLabel->setDescription(juce::String(("Minimum value of the parameter")));
  addWidget(minValueLabel, true, true);

  minValueField = new LibertyTextEntryField(juce::String(("MinValue")));
  minValueField->setDescription(minValueLabel->getDescription());
  minValueField->registerTextEntryFieldObserver(this);
  addWidget(minValueField, true, true);


  maxValueLabel = new RTextField(juce::String(("Max")));
  maxValueLabel->setDescription(juce::String(("Maximum value of the parameter")));
  addWidget(maxValueLabel, true, true);

  maxValueField = new LibertyTextEntryField(juce::String(("MaxValue")));
  maxValueField->setDescription(maxValueLabel->getDescription());
  maxValueField->registerTextEntryFieldObserver(this);
  addWidget(maxValueField, true, true);


  valueSlider = new LibertySlider(juce::String(("Value")));
  valueSlider->setSliderName(juce::String(("Value")));
  valueSlider->setRange(0.0, 1.0, 0.0, 0.5, true);
  valueSlider->setStringConversionFunction(&valueToString5);
  valueSlider->setDescription(juce::String(("Current value of the parameter")));
  valueSlider->addListener(this);
  addWidget(valueSlider, true, true);


  helpTextLabel = new RTextField(juce::String(("Help:")));
  helpTextLabel->setDescription(juce::String(("Help text for the parameter")));
  addWidget(helpTextLabel, true, true);

  helpTextField = new LibertyTextEntryField(juce::String(("HelpText"))); 
  helpTextField->setDescription(helpTextLabel->getDescription());
  helpTextField->registerTextEntryFieldObserver(this);
  addWidget(helpTextField, true, true);


  parameterSetupLabel = new RTextField(juce::String(("Parameter Setup")));
  parameterSetupLabel->setDescription(juce::String(("General setup for the parameter")));
  parameterSetupLabel->setJustification(juce::Justification::centred);
  addWidget(parameterSetupLabel, true, true);

  valueField = new LibertyLabeledTextEntryField(juce::String(("Value"))); 
  valueField->setDescription(juce::String(("Current value of the parameter")));
  valueField->setLabelText(juce::String(("Value:")));
  valueField->getTextEntryField()->registerTextEntryFieldObserver(this);
  addWidget(valueField, true, true);

  defaultField = new LibertyLabeledTextEntryField(juce::String(("DefaultValue"))); 
  defaultField->setDescription(juce::String(("Default value of the parameter")));
  defaultField->setLabelText(juce::String(("Default:")));
  defaultField->getTextEntryField()->registerTextEntryFieldObserver(this);
  addWidget(defaultField, true, true);

  unitField = new LibertyLabeledTextEntryField(juce::String(("Unit"))); 
  unitField->setDescription(juce::String(("Physical unit of the parameter")));
  unitField->getTextEntryField()->registerTextEntryFieldObserver(this);
  unitField->setLabelText(juce::String(("Unit:")));
  addWidget(unitField, true, true);

  scalingComboBox = new LibertyNamedComboBox(juce::String(("Scaling")));
  scalingComboBox->setComboBoxName(juce::String(("Scaling:")));
  scalingComboBox->setDescription(juce::String(("Scaling behavior of the parameter")));
  scalingComboBox->addItem(romos::ParameterModule::LINEAR_MAPPING,      "Linear",      true, false);
  scalingComboBox->addItem(romos::ParameterModule::EXPONENTIAL_MAPPING, "Exponential", true, false);
  addWidget(scalingComboBox, true, true);
  scalingComboBox->registerComboBoxObserver(this);


  setToDefaultButton = new RClickButton(juce::String(("Use")));
  setToDefaultButton->setDescription(juce::String(("Set parameter to default value")));
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
  ModulePropertiesEditor::updateWidgetsFromModuleState();


  romos::ParameterModule* pm = (romos::ParameterModule*) moduleToEdit;

  double minValue     = pm->getMinValue();
  double maxValue     = pm->getMaxValue();
  double defaultValue = pm->getDefaultValue();
  double value        = pm->getValue();
  double quantization = 0.0;  // preliminary
  valueSlider->setRange(minValue, maxValue, quantization, defaultValue, false);
  valueSlider->setValue(value, false, false);

  int scaling = pm->getMappingFunction();
  if( scaling == romos::ParameterModule::EXPONENTIAL_MAPPING )
    valueSlider->setScaling(Parameter::EXPONENTIAL);
  else
    valueSlider->setScaling(Parameter::LINEAR);
}

//-------------------------------------------------------------------------------------------------

TopLevelModuleEditor::TopLevelModuleEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newPlugInLock, newModuleToEdit)
{

}


//-------------------------------------------------------------------------------------------------

VoiceKillerModuleEditor::VoiceKillerModuleEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newPlugInLock, newModuleToEdit)
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

  thresholdSlider = new LibertySlider(juce::String(("Threshold")));
  thresholdSlider->setSliderName(juce::String(("Threshold")));
  thresholdSlider->setRange(-180.0, -40.0, 1.0, -100.0, true);
  thresholdSlider->setStringConversionFunction(&decibelsToStringWithUnit);
  thresholdSlider->setDescription(juce::String(("Amplitude threshold below which voice gets killed")));
  addWidget(thresholdSlider, true, true);
  thresholdSlider->addListener(this);


  timeOutSlider = new LibertySlider(juce::String(("TimeOut")));
  timeOutSlider->setSliderName(juce::String(("TimeOut")));
  timeOutSlider->setRange(0.01, 1.0, 0.01, 0.01, true);
  timeOutSlider->setStringConversionFunction(&valueToString2);
  timeOutSlider->setScaling(Parameter::EXPONENTIAL);
  timeOutSlider->setDescription(juce::String(("Time until voice gets killed after amplitude falls below threshold")));
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

WhiteNoiseModuleEditor::WhiteNoiseModuleEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newPlugInLock, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  seedSlider = new LibertySlider(juce::String(("Seed")));
  seedSlider->setSliderName(juce::String(("Seed")));
  seedSlider->setRange(0.0, 1000.0, 1.0, 0.0, true);
  seedSlider->setStringConversionFunction(&valueToString);
  seedSlider->setDescription(juce::String(("Seed for the pseudo-random number generator")));
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

BiquadDesignerModuleEditor::BiquadDesignerModuleEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newPlugInLock, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  modeComboBox = new LibertyNamedComboBox(juce::String(("Mode")));
  modeComboBox->addItem(romos::BiquadDesigner::BYPASS,                        "Bypass",                     true, false);
  modeComboBox->addItem(romos::BiquadDesigner::LOWPASS_6_BILINEAR,            "Lowpass, 6 dB/oct, BLT",     true, false);
  modeComboBox->addItem(romos::BiquadDesigner::HIGHPASS_6_BILINEAR,           "Highpass, 6 dB/oct, BLT",    true, false);
  modeComboBox->addItem(romos::BiquadDesigner::LOW_SHELF_1_BILINEAR,          "Low Shelf, 1st order, BLT",  true, false);
  modeComboBox->addItem(romos::BiquadDesigner::HIGH_SHELF_1_BILINEAR,         "High Shelf, 1st order, BLT", true, false);
  modeComboBox->addItem(romos::BiquadDesigner::ALLPASS_1_BILINEAR,            "Allpass, 1st order, BLT",    true, false);
  modeComboBox->addItem(romos::BiquadDesigner::LOWPASS_12_BILINEAR,           "Lowpass, 12 dB/oct, BLT",     true, false);
  modeComboBox->addItem(romos::BiquadDesigner::HIGHPASS_12_BILINEAR,          "Highpass, 12 dB/oct, BLT",    true, false);
  modeComboBox->addItem(romos::BiquadDesigner::BANDPASS_CONST_SKIRT_BILINEAR, "Bandpass, const. skirt, BLT", true, false);
  modeComboBox->addItem(romos::BiquadDesigner::BANDPASS_CONST_PEAK_BILINEAR,  "Bandpass, const. peak, BLT",  true, false);
  modeComboBox->addItem(romos::BiquadDesigner::BANDREJECT_BILINEAR,           "Bandreject, BLT",             true, false);
  modeComboBox->addItem(romos::BiquadDesigner::PEAK_BILINEAR,                 "Peak, BLT",                   true, false);
  modeComboBox->addItem(romos::BiquadDesigner::LOW_SHELF_2_BILINEAR,          "Low Shelf, 2nd order, BLT",  true, false);
  modeComboBox->addItem(romos::BiquadDesigner::HIGH_SHELF_2_BILINEAR,         "High Shelf, 2nd order, BLT", true, false);
  modeComboBox->addItem(romos::BiquadDesigner::ALLPASS_2_BILINEAR,            "Allpass, 2nd order, BLT",    true, false);
  modeComboBox->setDescription(juce::String(("Mode of the filter to be designed")));
  modeComboBox->setComboBoxName(juce::String(("Mode:")));
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

LibertyLadderFilterModuleEditor::LibertyLadderFilterModuleEditor(CriticalSection *newPlugInLock, 
  romos::Module* newModuleToEdit) : ModulePropertiesEditor(newPlugInLock, newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);

  int labelWidth = 72;

  filterModeComboBox = new LibertyNamedComboBox(juce::String(("Mode")));
  filterModeComboBox->addItem(romos::LadderFilter::FLAT,    "Flat",                  true, false);
  filterModeComboBox->addItem(romos::LadderFilter::LP_6,    "Lowpass, 6 dB/oct",     true, false);
  filterModeComboBox->addItem(romos::LadderFilter::LP_12,   "Lowpass, 12 dB/oct",    true, false);
  filterModeComboBox->addItem(romos::LadderFilter::LP_18,   "Lowpass, 18 dB/oct",    true, false);
  filterModeComboBox->addItem(romos::LadderFilter::LP_24,   "Lowpass, 24 dB/oct",    true, false);
  filterModeComboBox->addItem(romos::LadderFilter::HP_6,    "Highpass, 6 dB/oct",    true, false);
  filterModeComboBox->addItem(romos::LadderFilter::HP_12,   "Highpass, 12 dB/oct",   true, false);
  filterModeComboBox->addItem(romos::LadderFilter::HP_18,   "Highpass, 18 dB/oct",   true, false);
  filterModeComboBox->addItem(romos::LadderFilter::HP_24,   "Highpass, 24 dB/oct",   true, false);
  filterModeComboBox->addItem(romos::LadderFilter::BP_6_6,  "Bandpass, 6/6 dB/oct",  true, false);
  filterModeComboBox->addItem(romos::LadderFilter::BP_6_12, "Bandpass, 6/12 dB/oct", true, false);
  filterModeComboBox->addItem(romos::LadderFilter::BP_12_6, "Bandpass, 12/6 dB/oct", true, false);
  filterModeComboBox->addItem(romos::LadderFilter::BP_6_18, "Bandpass, 6/18 dB/oct", true, false);
  filterModeComboBox->addItem(romos::LadderFilter::BP_18_6, "Bandpass, 18/6 dB/oct", true, false);
  filterModeComboBox->setDescription(juce::String(("Mode of the filter")));
  filterModeComboBox->setComboBoxName(juce::String(("Mode:")));
  filterModeComboBox->setNameLabelWidth(labelWidth);
  filterModeComboBox->registerComboBoxObserver(this);
  addWidget(filterModeComboBox, true, true);

  saturationModeComboBox = new LibertyNamedComboBox(juce::String(("SaturationMode")));
  saturationModeComboBox->addItem(romos::LadderFilter::NO_SATURATION, "No Saturation", true, false);
  saturationModeComboBox->addItem(romos::LadderFilter::LAST_STAGE,    "Last Stage",    true, false);
  saturationModeComboBox->addItem(romos::LadderFilter::FEEDBACK,      "Feedback",      true, false);
  saturationModeComboBox->addItem(romos::LadderFilter::EACH_STAGE,    "Each Stage",    true, false);
  saturationModeComboBox->setDescription(juce::String(("Point(s) in the filter where saturation is applied")));
  saturationModeComboBox->setComboBoxName(juce::String(("Saturation:")));
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
