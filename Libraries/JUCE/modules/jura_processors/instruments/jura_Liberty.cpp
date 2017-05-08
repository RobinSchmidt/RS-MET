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

