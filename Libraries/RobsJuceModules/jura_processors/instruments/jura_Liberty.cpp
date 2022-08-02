using namespace romos;

//=================================================================================================
// class LibertyInterfaceState:

LibertyInterfaceState::LibertyInterfaceState()
{
  activePanel = STRUCTURE_PANEL;
}

//=================================================================================================
// class LibertyAudioModule:
/*
LibertyAudioModule::LibertyAudioModule(CriticalSection *newPlugInLock, 
  romos::Liberty *modularSynthToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, modularSynthToWrap)
//: AudioModule(newPlugInLock)
{
  jassert(modularSynthToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedLiberty = modularSynthToWrap;
  underlyingRosicInstrument = modularSynthToWrap;
  init();
}
*/
LibertyAudioModule::LibertyAudioModule(CriticalSection *newPlugInLock) 
  : PolyphonicInstrumentAudioModule(newPlugInLock)
{
  wrappedLiberty = new romos::Liberty;
  //underlyingRosicInstrument = wrappedLiberty;
  wrappedLibertyIsOwned = true;
  init();
}

void LibertyAudioModule::init()
{
  setModuleTypeName("Liberty");
  macroDirectory = getPresetDirectory() + "/Macros";
}

LibertyAudioModule::~LibertyAudioModule()
{
  if(wrappedLibertyIsOwned)
    delete wrappedLiberty;
}

AudioModuleEditor* LibertyAudioModule::createEditor(int type)
{
  return new jura::LibertyEditor(lock, this); // get rid of passing the lock
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
    if ( module->isParameterModule() )  // to be used later
    //if( module->getTypeIdentifierOld() == romos::ModuleTypeRegistry::PARAMETER )
    {
      // we need special treatment of the parameter module - minValue/maxValue shlould be set 
      // simultaneuously, because otherwise, min/max might not be recalled correctly, for example 
      // when newMin > oldMax, the module will refuse to use the new minimum, etc. moreover, 
      // min/max must be set up before attempting to set the default and current value:
      double min = xmlState.getDoubleAttribute("MinValue", 0.0);
      double max = xmlState.getDoubleAttribute("MaxValue", 1.0);

      juce::String mappingString = xmlState.getStringAttribute("MaxValue", "Linear");
      int mappingIndex = romos::ParameterModule::LINEAR_MAPPING;
      if( mappingString == "Exponential" )
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

  std::map<std::string, std::string> stateMap = module->getState();
  if(!stateMap.empty())
    addAttributesFromMap(*xmlState, stateMap);
  else
    jassertfalse;

  unsigned int i;
  romos::ContainerModule *container = dynamic_cast<romos::ContainerModule*> (module);
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

      if( connectionString != juce::String() )
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
  if( module == nullptr ) { jassertfalse; return; }

  std::map<std::string, std::string> moduleState = getAttributesAsMap(xmlState);
  module->setState(moduleState);

  if(  module->isContainerModule() || module->isTopLevelModule() ) {
    romos::ContainerModule *container = dynamic_cast<romos::ContainerModule*> (module);
    for(int i=0; i<xmlState.getNumChildElements(); i++) {
      XmlElement* childState = xmlState.getChildElement(i);
      rosic::rsString moduleTypeName = juceToRosic(childState->getTagName()); // get rid of that intermediate format

      int typeIdentifier = romos::moduleFactory.getModuleId(moduleTypeName.asStdString());
      if( typeIdentifier != -1 ) {

        // verify that this is useless and then delete this old code:
        //if(  module->isTopLevelModule() 
        //  && romos::ModuleTypeRegistry::isIdentifierInputOrOutput(typeIdentifier) )
        //if( module->isTopLevelModule() && (module->isInputModule() || module->isOutputModule()) )
           // noo - this is wrong - we are not interested in whethere the "module" is I/O but rather
           // the child to be added should be I/O

        if( module->isTopLevelModule() && (moduleTypeName == "AudioInput" || moduleTypeName == "AudioOutput") ) {
          // do nothing when this is the top-level module and the to-be-added child is an I/O module
        }
        else {
          romos::Module *newModule = container->addChildModule(
            moduleTypeName.asStdString(), "", 0, 0, false, false);
          createAndSetupEmbeddedModulesFromXml(*childState, newModule);
        }
      }
    }
    container->sortChildModuleArray();
  }
  else
    restoreModuleTypeSpecificStateDataFromXml(module, xmlState);
}

void LibertyAudioModule::createConnectionsFromXml(const XmlElement& xmlState, romos::Module *module)
{
  if( module == NULL )
  {
    jassertfalse;
    return;
  }

  juce::String connectionString = xmlState.getStringAttribute("IncomingConnections", juce::String());

  if( connectionString != juce::String() )
  {
    connectionString = connectionString.removeCharacters(" ");
    juce::String remainingString = connectionString + juce::String(",");
    juce::String currentString;
    juce::String tmpString1, tmpString2;
    int smi, spi, tpi; // source module- and pin-index, target pin index
    while( remainingString != juce::String() )
    {
      currentString   = remainingString.upToFirstOccurrenceOf(",", false, false);
      remainingString = remainingString.fromFirstOccurrenceOf(",", false, false);
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
   
      tmpString1 = tmpString1.fromFirstOccurrenceOf("_", false, false); 
      tmpString2 = tmpString1.upToFirstOccurrenceOf("_", false, false);
      smi        = tmpString2.getIntValue();
      tmpString1 = tmpString1.fromFirstOccurrenceOf("_", false, false); 
      tmpString2 = tmpString1.upToFirstOccurrenceOf("_", false, false);
      spi        = tmpString2.getIntValue();
      tmpString1 = tmpString1.fromFirstOccurrenceOf("_", false, false); 
      tmpString2 = tmpString1.upToFirstOccurrenceOf("_", false, false);
      tpi        = tmpString2.getIntValue();

      romos::ContainerModule *parentModule = module->getParentModule();
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
  romos::ContainerModule *container = dynamic_cast<romos::ContainerModule*> (module);
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

  sendActionMessage("StateRecall"); //  GUI should delete all pointers to the modules

  wrappedLiberty->reset();
  romos::TopLevelModule *topLevelModule = wrappedLiberty->getTopLevelModule();
  topLevelModule->deleteAllChildModules();       
  topLevelModule->disconnectAudioOutputModules();

  XmlElement* topLevelModuleState = xmlState.getChildByName("TopLevelModule");
  if( topLevelModuleState != nullptr )
    setModuleStateFromXml(*topLevelModuleState, topLevelModule);
  restoreTopLevelInOutStates(*topLevelModuleState);
  //PolyphonicInstrumentAudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

void LibertyAudioModule::restoreTopLevelInOutStates(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*lock);
  romos::ContainerModule *topLevelMasterVoice = wrappedLiberty->getTopLevelModule();

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

//ModulePropertiesEditor::ModulePropertiesEditor(CriticalSection *newPlugInLock, 
//  romos::Module* newModuleToEdit)

ModulePropertiesEditor::ModulePropertiesEditor(LibertyAudioModule *newLiberty, romos::Module* newModuleToEdit)
{

  //plugInLock   = newPlugInLock;  // old

  // new:
  libertyModule = newLiberty;
  plugInLock    = libertyModule->getCriticalSection();;

  ScopedLock scopedLock(*plugInLock);

  moduleToEdit = newModuleToEdit; 
  libertyModule->addActionListener(this);

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

ModulePropertiesEditor::~ModulePropertiesEditor()
{
  libertyModule->removeActionListener(this);
}

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
      rosic::rsString name  = juceToRosic(lmw->getWidgetParameterName());
      rosic::rsString value = juceToRosic(widgetThatHasChanged->getStateAsString());
      m->setParameter(name, value);
    }
  }
  updateWidgetsFromModuleState(); // to update the other widgets, implements widget-interdependence
}

void ModulePropertiesEditor::actionListenerCallback(const String& message)
{
  if(message == "StateRecall")
    moduleToEdit = nullptr; // invalidate pointer
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
        rosic::rsString name  = juceToRosic(lmw->getWidgetParameterName());
        juce::String  value = rosicToJuce(m->getParameterValue(name));
        widgets[i]->setStateFromString(value, false);
      }
    }
  }
}

//=================================================================================================
//=================================================================================================

LibertyInterfaceComponent::LibertyInterfaceComponent(
  LibertyInterfaceMediator *interfaceMediatorToUse)
{
  mediator = interfaceMediatorToUse;
}

LibertyInterfaceMediator* LibertyInterfaceComponent::getInterfaceMediator() const 
{ 
  return dynamic_cast<LibertyInterfaceMediator*> (mediator); 
}

//=================================================================================================

LibertyInterfaceMediator::LibertyInterfaceMediator(CriticalSection *newPlugInLock, 
  LibertyAudioModule* newLibertyModuleToEdit)
{
  plugInLock = newPlugInLock;     
  ScopedLock scopedLock(*plugInLock);
  modularSynthModuleToEdit = newLibertyModuleToEdit;
  topLevelModule           = modularSynthModuleToEdit->wrappedLiberty->getTopLevelModule();
  containerShownInDiagram  = topLevelModule;
  moduleToShowEditorFor    = containerShownInDiagram;
}

LibertyInterfaceMediator::~LibertyInterfaceMediator()
{
  //int dummy = 0;
}

void LibertyInterfaceMediator::setContainerToShowInDiagram(romos::ContainerModule* containerToShow)
{
  ScopedLock scopedLock(*plugInLock);
  containerShownInDiagram = containerToShow;
  //moduleToShowEditorFor = containerShownInDiagram; // superfluous?
  setModuleToShowEditorFor(containerShownInDiagram);
  Mediator::sendNotificationToColleagues(NULL, 
    LibertyInterfaceComponent::CONTAINER_SHOWN_IN_DIAGRAM);  // NULL is preliminary
}

void LibertyInterfaceMediator::setModuleToShowEditorFor(romos::Module *module)
{
  ScopedLock scopedLock(*plugInLock);
  moduleToShowEditorFor = module;
  Mediator::sendNotificationToColleagues(NULL, 
    LibertyInterfaceComponent::MODULE_TO_SHOW_EDITOR_FOR);  // NULL is preliminary
}

//=================================================================================================
// class ModularStructureTreeView

ModularStructureTreeView::ModularStructureTreeView(LibertyInterfaceMediator *interfaceMediatorToUse)
  : LibertyInterfaceComponent(interfaceMediatorToUse)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  //mediator->registerColleague(this);
  setMediator(interfaceMediatorToUse); // will register "this" as colleague

  setRootNode( new RTreeViewNode(juce::String("TopLevelModule")) );
  rootNode->setDeleteChildNodesOnDestruction(true);
  rootNode->setUserData(getInterfaceMediator()->topLevelModule);
  rebuildTree();
}

ModularStructureTreeView::~ModularStructureTreeView()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  //getInterfaceMediator()->deRegisterInterfaceComponent(this); // should be done in baseclass destructor

  delete rootNode;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void ModularStructureTreeView::handleMediatorNotification(
  MediatedColleague *originatingColleague, int messageCode, rsMessageData* messageData)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));


  if( messageCode == NUM_CHILDREN || messageCode == MODULE_NAME 
    || messageCode == CONTAINER_SHOWN_IN_DIAGRAM ) 
    // actually, the last one is overkill - we only need to rebuild if the top-level module 
    // changes - i.e. a patch is loaded
  {
    rebuildTree();
    // \todo often, rebuilding the whole tree is overkill and it may suffice to just add or remove 
    // a single node - we need to somehow distinguis these cases....
  }

  if( messageCode == MODULE_TO_SHOW_EDITOR_FOR )
    updateNodeHighlighting();

}

void ModularStructureTreeView::nodeClicked(RTreeViewNode *nodeThatWasClicked, 
  const MouseEvent &mouseEvent, int clickPosition)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( clickPosition == RTreeView::ON_TEXT )
  {
    romos::Module *clickedNodeModule =  
      static_cast<romos::Module*> (nodeThatWasClicked->getUserData());
    if( dynamic_cast<romos::ContainerModule*> (clickedNodeModule) != NULL )
      getInterfaceMediator()->setContainerToShowInDiagram( 
        static_cast<romos::ContainerModule*> (nodeThatWasClicked->getUserData()) );
    else
    {  
      getInterfaceMediator()->setContainerToShowInDiagram(clickedNodeModule->getParentModule());
      //if( romos::ModuleTypeRegistry::hasModuleTypeEditor(clickedNodeModule->getTypeIdentifierOld()) )
      if( clickedNodeModule->hasEditor() )
        getInterfaceMediator()->setModuleToShowEditorFor(clickedNodeModule);
    }
  }
  else if( clickPosition == RTreeView::ON_PLUSMINUS )
    RTreeView::nodeClicked(nodeThatWasClicked, mouseEvent, clickPosition);  // opens/closes the node
}

//-------------------------------------------------------------------------------------------------
// others:

void ModularStructureTreeView::rebuildTree()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  rootNode->deleteChildNodesRecursively();
  rootNode->setNodeText(rosicToJuce(getInterfaceMediator()->topLevelModule->getName()));
  for(unsigned int i=0; i<getInterfaceMediator()->topLevelModule->getNumChildModules(); i++)
  {
    // \todo add I/O modules ... or maybe not...
    createAndHangInSubTree(rootNode, getInterfaceMediator()->topLevelModule->getChildModule(i));
  }
  updateScrollBarBoundsAndVisibility();
  updateNodeHighlighting();
  repaint();
}

void ModularStructureTreeView::createAndHangInSubTree(RTreeViewNode *parentNodeToUse, 
  romos::Module *moduleToCreateSubTreeFor)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  // check if node is container - if not do nothing and return - we don't want leaf nodes in the tree
  //romos::ModuleContainer *containerModule = dynamic_cast<romos::ModuleContainer*> (moduleToCreateSubTreeFor);
  //if( containerModule == NULL )
  //  return;

  //if( romos::ModuleTypeRegistry::hasModuleTypeEditor(moduleToCreateSubTreeFor->getTypeIdentifierOld()) )
  if( moduleToCreateSubTreeFor->hasEditor() ) // maybe optionally show all modules in the tree
  //if(true)  // preliminary
  //if( moduleToCreateSubTreeFor->hasTreeNode() )
  {
    juce::String name = juce::String( moduleToCreateSubTreeFor->getName().c_str() );
    RTreeViewNode *newNode = new RTreeViewNode(name);
    newNode->setUserData(moduleToCreateSubTreeFor);
    parentNodeToUse->addChildNode(newNode);


    // recursion to take care of the new node's child nodes:
    romos::ContainerModule *containerModule = 
      dynamic_cast<romos::ContainerModule*> (moduleToCreateSubTreeFor);
    if( containerModule != NULL )
    {
      for(unsigned int i=0; i<containerModule->getNumChildModules(); i++) // conatinerModule was already checked for NULL
        createAndHangInSubTree(newNode, containerModule->getChildModule(i));
    }
  }
}

void ModularStructureTreeView::updateNodeHighlighting()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  rootNode->setAllNodesUnticked();
  //RTreeViewNode *selectedNode = rootNode->findNodeByData(getInterfaceMediator()->getContainerShownInDiagram());
  RTreeViewNode *selectedNode 
    = rootNode->findNodeByData(getInterfaceMediator()->getModuleToShowEditorFor());
  jassert( selectedNode != NULL );  // one module/node should always be focused...
  if( selectedNode != NULL )
    selectedNode->setTicked(true);
  repaint();
}

//=================================================================================================
// class ModulePropertiesEditorHolder:

ModulePropertiesEditorHolder::ModulePropertiesEditorHolder(
  LibertyInterfaceMediator *interfaceMediatorToUse)
  : LibertyInterfaceComponent(interfaceMediatorToUse)
{  
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock)); 
  // getInterfaceMediator() returns already valid mediator because we pass interfaceMediatorToUse to our basclass constructor

  setMediator(interfaceMediatorToUse); // will register "this" as colleague


  currentEditor = NULL; // init to NULL required because the next call first deletes the old pointer
  createPropertiesEditorForSelectedModule();
}

ModulePropertiesEditorHolder::~ModulePropertiesEditorHolder()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void ModulePropertiesEditorHolder::handleMediatorNotification(
  MediatedColleague *originatingColleague, int messageCode, rsMessageData* messageData)
{
  // \todo maybe include a switch on the messageCode later - we may not want to re-create the 
  // editor on all kinds of messages

  createPropertiesEditorForSelectedModule();
}

void ModulePropertiesEditorHolder::paint(Graphics &g)
{
  // overriden with empty function to avoid painting a gradient that will be invisible anyway
}

void ModulePropertiesEditorHolder::resized()
{
  currentEditor->setBounds(0, 0, getWidth(), getHeight());
}

//-------------------------------------------------------------------------------------------------
// others:

void ModulePropertiesEditorHolder::createPropertiesEditorForSelectedModule()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  romos::Module *moduleToShowEditorFor = getInterfaceMediator()->getModuleToShowEditorFor();
  if( moduleToShowEditorFor == NULL )
    moduleToShowEditorFor = getInterfaceMediator()->getContainerShownInDiagram(); // preliminary

  removeChildColourSchemeComponent(currentEditor, true);

  /*
  // old:
  // this switch statement sucks - use std::map or something - maybe when this object is created, 
  // create a map that uses the module-id as key and the creator function (pointer) as value
  switch( moduleToShowEditorFor->getTypeIdentifierOld() )
  {
  case romos::ModuleTypeRegistry::PARAMETER:
    currentEditor = new ParameterModuleEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);    break;
  case romos::ModuleTypeRegistry::CONTAINER:   
    currentEditor = new ContainerModuleEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);    break;
  case romos::ModuleTypeRegistry::TOP_LEVEL_MODULE:   
    currentEditor = new TopLevelModuleEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);    break;
  case romos::ModuleTypeRegistry::WHITE_NOISE:   
    currentEditor = new WhiteNoiseModuleEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);    break;
  case romos::ModuleTypeRegistry::BIQUAD_DESIGNER:
    currentEditor = new BiquadDesignerModuleEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);    break;
  case romos::ModuleTypeRegistry::LADDER_FILTER:
    currentEditor = new LibertyLadderFilterModuleEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);    break;
  case romos::ModuleTypeRegistry::VOICE_KILLER:   
    currentEditor = new VoiceKillerModuleEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);    break;
  default:
  {
    //jassertfalse;  // for every module-type, there should be case
    currentEditor = new ModulePropertiesEditor(getInterfaceMediator()->modularSynthModuleToEdit, 
      moduleToShowEditorFor);
  }
  }
  */

  // new:
  // abbreviations for convenience:
  LibertyAudioModule* lbrtyMd = getInterfaceMediator()->modularSynthModuleToEdit;
  romos::Module* mdl = moduleToShowEditorFor;
  std::string type = mdl->getTypeName();
  if(     type == "Parameter")      currentEditor = new ParameterModuleEditor(lbrtyMd, mdl);
  else if(type == "Container")      currentEditor = new ContainerModuleEditor(lbrtyMd, mdl);
  else if(type == "TopLevelModule") currentEditor = new TopLevelModuleEditor(lbrtyMd, mdl);
  else if(type == "WhiteNoise")     currentEditor = new WhiteNoiseModuleEditor(lbrtyMd, mdl);
  else if(type == "BiquadDesigner") currentEditor = new BiquadDesignerModuleEditor(lbrtyMd, mdl);
  else if(type == "LadderFilter")   currentEditor = new LibertyLadderFilterModuleEditor(lbrtyMd, mdl);
  else if(type == "VoiceKiller")    currentEditor = new VoiceKillerModuleEditor(lbrtyMd, mdl);
  //else if(type == "Formula_1_1")    currentEditor = new LibertyFormulaModuleEditor(lbrtyMd, mdl);
  //else if(type == "Formula_N_1")    currentEditor = new LibertyFormula_N_1ModuleEditor(lbrtyMd, mdl);
  else if(type == "Formula")        currentEditor = new LibertyFormula_N_MModuleEditor(lbrtyMd, mdl);
  else                              currentEditor = new ModulePropertiesEditor(lbrtyMd, mdl); // generic
  // todo: optimize away all these string-comparisons
  // maybe make a map from type-id to creator-function

  currentEditor->setDescriptionField(descriptionField, true );
  addChildColourSchemeComponent(currentEditor, true, true);
  resized();  // will set the bounds of the child
}

//=================================================================================================
//=================================================================================================
// class ModularBlockDiagramPanel:

ModularBlockDiagramPanel::ModularBlockDiagramPanel(LibertyInterfaceMediator *interfaceMediatorToUse)
  : LibertyInterfaceComponent(interfaceMediatorToUse)
{  
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock)); 
  // getInterfaceMediator() returns already valid mediator because we pass interfaceMediatorToUse 
  // to our basclass constructor

  setMediator(interfaceMediatorToUse); // will register "this" as colleague

  availableModulesTreeView = new RTreeView();
  availableModulesTreeView->setOpenOrCloseNodesOnClick(true);
  availableModulesTreeView->setOpaque(true);
  availableModulesTreeView->setAlwaysOnTop(true);
  availableModulesTreeView->setDrawRootNode(false);
  fillAvailableModulesTreeView();
  availableModulesTreeView->registerTreeViewObserver(this);

  actOnSelectionMenu = new RPopUpMenu(this);  // can we pass NULL too?
  actOnSelectionMenu->setDismissOnFocusLoss(false); // because it immediately seems to loose focus after opening for some reason (?)
  actOnSelectionMenu->registerPopUpMenuObserver(this);
  addWidget(actOnSelectionMenu, false, false);

  nameEntryField = new RTextEntryField();
  nameEntryField->registerTextEntryFieldObserver(this);
  addWidget(nameEntryField, true, false);


  mouseDownX        = 0; 
  mouseDownY        = 0;
  selectionOffsetX  = 0;
  selectionOffsetY  = 0;
  availableWidth    = 0;
  availableHeight   = 0;
  //gridStyle       = DOTTED_GRID;
  //gridStyle       = GRID_LINES;
  gridStyle         = NO_GRID;

  // define metric:
  m  = 2;                    // margin between text and outlines ...use later 2 here - maybe
  t  = 2;                    // thickness of outlines
  s  = 2;                    // stickout for the pins

  bigFontHeight    = bigFont->getFontHeight();    
  normalFontHeight = normalFont->getFontHeight();
  smallFontHeight  = smallFont->getFontHeight();

  pinDistance      = smallFontHeight + m;
  arrowLength      = 12;  
  arrowHeadLength  = 8;

  // todo: use colors from colorscheme
  wireColour      = Colours::black;
  highlightColour = Colours::red;
  // maybe override ColorSchemeComponent::colorSchmeChanged and set the colors there

  pinHighlighter = new RectangleComponent(highlightColour, highlightColour, 0);
  pinHighlighter->setInterceptsMouseClicks(false, false);
  addChildComponent(pinHighlighter);

  lassoComponent = new RectangleComponent(highlightColour.withMultipliedAlpha(0.0625f), 
    highlightColour.withMultipliedAlpha(0.25f), 2);
  addChildComponent(lassoComponent);
}

ModularBlockDiagramPanel::~ModularBlockDiagramPanel()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  //getInterfaceMediator()->deRegisterInterfaceComponent(this);  // done in baseclass

  delete actOnSelectionMenu;
  delete availableModulesTreeView;
  delete treeRootNode;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void ModularBlockDiagramPanel::setAvailabeSizeForCanvas(int w, int h)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  availableWidth  = w;
  availableHeight = h;
  updateCanvasBounds(0, 0, w, h);
}

void ModularBlockDiagramPanel::updateCanvasBounds(int x, int  y, int w, int h)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  x = jmin(x, getMinXInPixels());
  y = jmin(y, getMinYInPixels()); 
  w = jmax(w, getMaxXInPixels())-x;
  h = jmax(h, getMaxYInPixels())-y;
  setBounds(x, y, w, h);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

int ModularBlockDiagramPanel::getMaxXInPixels() const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  int result = 0;
  int tmp;
  romos::Module *child;
  for(unsigned int i=0; i<getInterfaceMediator()->getContainerShownInDiagram()->getNumChildModules(); i++)
  {
    child = getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i);
    tmp = inPixels(child->getPositionX()) + getRequiredWidthForModuleInPixels(child, true);
    if( tmp > result )
      result = tmp;
  }
  return result;
}

int ModularBlockDiagramPanel::getMaxYInPixels() const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  int result = 0;
  int tmp;
  romos::Module *child;

  for(unsigned int i=0; i<getInterfaceMediator()->getContainerShownInDiagram()->getNumChildModules(); i++)
  {
    child = getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i);
    tmp = inPixels(child->getPositionY()) - getOffsetY(child) + getRequiredHeightForModuleInPixels(child);
    // maybe use getModuleRectangle here
    if( tmp > result )
      result = tmp;
  }
  return result;
}

int ModularBlockDiagramPanel::getMinXInPixels() const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  int result = 0;
  romos::Module *child;

  for(unsigned int i=0; i<getInterfaceMediator()->getContainerShownInDiagram()->getNumChildModules(); i++)
  {
    child = getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i);
    if( child->getPositionX() < result )
      result = child->getPositionX();
  }
  result = inPixels(result);
  if( result < 0 )
    result -= s; // the stickout should also be seen
  return result; 
}

int ModularBlockDiagramPanel::getMinYInPixels() const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  int result = 0;
  int y;
  romos::Module *child;
  for(unsigned int i=0; i<getInterfaceMediator()->getContainerShownInDiagram()->getNumChildModules(); i++)
  {
    child = getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i);
    y     = inPixels(child->getPositionY()) - getOffsetY(child);
    if( y  < result )
      result = y;
  }
  return result;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void ModularBlockDiagramPanel::mouseExit(const MouseEvent &e)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  ColourSchemeComponent::mouseExit(e);
  pinHighlighter->setVisible(false);
}

void ModularBlockDiagramPanel::mouseMove(const MouseEvent &e)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  romos::ContainerModule *container       = getInterfaceMediator()->getContainerShownInDiagram();
  romos::Module          *module          = getModuleAtPixels(e.x, e.y, true);
  juce::Rectangle<int>    moduleRectangle = getRectangleForModuleInPixels(module, true);

  int kindOfPin, directionOfPin, indexOfPin;
  bool isOverPin = getPinPropertiesAtPixels(e.x, e.y, module, moduleRectangle, kindOfPin, directionOfPin, indexOfPin);

  if( isOverPin )
  {
    pinHighlighter->setBounds( getPinBounds(kindOfPin, directionOfPin, indexOfPin, module, moduleRectangle) );
    pinHighlighter->setVisible(true);
  }
  else
    pinHighlighter->setVisible(false);
}

void ModularBlockDiagramPanel::mouseDown(const MouseEvent &e)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  availableModulesTreeView->setVisible(false);
  actOnSelectionMenu->setVisible(false);
  nameEntryField->setVisible(false);

  romos::ContainerModule *container = getInterfaceMediator()->getContainerShownInDiagram();

  mouseDownX = e.x; 
  mouseDownY = e.y;

  lassoComponent->setVisible(false);
  lassoComponent->setBounds(e.x, e.y, 0, 0);

  romos::Module          *module    = getModuleAtPixels(e.x, e.y, true);
  romos::AudioConnection connection = getConnectionAtPixels(e.x, e.y);




  if( module == NULL && connection.isNull() )
  {
    if( e.getNumberOfClicks() == 2 )
    {
      romos::ContainerModule *shownContainer  = getInterfaceMediator()->getContainerShownInDiagram();
      romos::ContainerModule *parentContainer = shownContainer->getParentModule();
      if( parentContainer != NULL )
      {
        getInterfaceMediator()->setContainerToShowInDiagram(parentContainer);
        getInterfaceMediator()->setModuleToShowEditorFor(shownContainer);
      }
    }
    else if( e.mods.isRightButtonDown() )
      openModuleInsertionMenu(getScreenX()+e.x, getScreenY()+e.y);
    else if( e.mods.isLeftButtonDown() )
    {
      if( !e.mods.isShiftDown() )
      {
        selectedModules.clear();
        selectedAudioConnections.clear();
      }
      lassoComponent->setVisible(true);
      repaint();
    }
  }
  else
  {
    if( e.mods.isRightButtonDown() )
    {
      openActOnSelectionMenu(getScreenX()+e.x, getScreenY()+e.y);
      //repaint();
      return;
    }

    // check if the click was on a pin:
    int kindOfPin, directionOfPin, indexOfPin;
    bool pinClicked = getPinPropertiesAtPixels(e.x, e.y, module, getRectangleForModuleInPixels(module, true), 
      kindOfPin, directionOfPin, indexOfPin);
    if( pinClicked )
    {
      // start drawing a new connection:
      if( kindOfPin == romos::AUDIO )
      {
        tmpAudioConnection.resetToNull(); 
        if( directionOfPin == romos::OUTGOING )
        {
          tmpAudioConnection.setSourceModule(module);
          tmpAudioConnection.setSourceOutputIndex(indexOfPin);
        }
        else 
        {
          tmpAudioConnection.setTargetModule(module);
          tmpAudioConnection.setTargetInputIndex(indexOfPin);
        }
      }
    }

    bool selected;
    bool connectionClicked = !connection.isNull();
    if( connectionClicked )
    {
      selected = rosic::containsElement(selectedAudioConnections, connection);
      //selected = selectedAudioConnections.hasElement(connection);
      if( e.mods.isShiftDown() )
      {
        if( selected )
          rosic::removeElementByValue(selectedAudioConnections, connection);
        //selectedAudioConnections.removeElementByValue(connection);
        else
          rosic::appendElement(selectedAudioConnections, connection);
        //selectedAudioConnections.appendElement(connection);
      }
      else
      {
        if( selected )
        {
          // keep selected (important for start dragging)
        }
        else
        {
          selectedModules.clear();
          selectedAudioConnections.clear();
          rosic::appendElement(selectedAudioConnections, connection);
          //selectedAudioConnections.appendElement(connection);
        }
      }
    }
    else
    {
      // click was not on the pin or connections, so it must have been on a module's body - change module selection:
      selected = isModuleSelected(module);
      if( e.mods.isShiftDown() )
      {
        if( selected )
          removeModuleWithConnectionsFromArray(module, selectedModules, selectedAudioConnections);
        else
          addModuleWithConnectionsToArray(module, selectedModules, selectedAudioConnections);
      }
      else if( e.getNumberOfClicks() == 2 )
      {
        if( selectedModules.size() == 1 )
        {
          //openModuleNameEntryField(); ...nah - maybe make this available via the PopUp - double clicks naviagte into the conatiner
          if( selectedModules[0]->isContainerModule() )
            getInterfaceMediator()->setContainerToShowInDiagram( ((ContainerModule*)selectedModules[0]) );
        }
      }
      else
      {
        if( selected )
        {
          // keep selected (important for start dragging)
        }
        else
        {
          //if( romos::ModuleTypeRegistry::hasModuleTypeEditor(module->getTypeIdentifierOld()) )
          if( module->hasEditor() )
            getInterfaceMediator()->setModuleToShowEditorFor(module); // will also select it via the callback that we'll receive
          else
            selectSingleModule(module);
        }
      }
    }
    repaint();
    return;
  }
}

/*
void ModularBlockDiagramPanel::mouseDoubleClick(const MouseEvent &e)
{

// todo: navigate up and down in the container hierarchy

int dummy = 0;

}
*/

void ModularBlockDiagramPanel::mouseDrag(const MouseEvent &e)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  mouseMove(e);
  if( !tmpAudioConnection.isNull() ) // we are drawing a new audio connection
  {

  }
  else if( lassoComponent->isVisible() )  // we are opening a lasso-selector
  {
    int x = mouseDownX;
    int y = mouseDownY;
    int w = e.x-mouseDownX;
    int h = e.y-mouseDownY;
    if( w < 0 )
    {
      w  = -w;
      x -= w;
    }
    if( h < 0 )
    {
      h  = -h;
      y -= h;
    }
    lassoComponent->setBounds(x, y, w, h);
    modulesInLasso.clear();
    audioConnectionsInLasso.clear();
    std::vector<romos::Module*> tmpArray = getModulesInRectangle(lassoComponent->getBounds());
    for(unsigned int i = 0; i < tmpArray.size(); i++)
      addModuleWithConnectionsToArray(tmpArray[i], modulesInLasso, audioConnectionsInLasso);
  }
  else // we are dragging a bunch of modules around
  {
    selectionOffsetX = snapPixelPositionToGrid(e.getDistanceFromDragStartX());
    selectionOffsetY = snapPixelPositionToGrid(e.getDistanceFromDragStartY());
  }
  repaint();
}

void ModularBlockDiagramPanel::mouseUp(const MouseEvent &e)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  if( !tmpAudioConnection.isNull() ) // we were drawing a new audio connection... 
  {
    romos::Module *module = getModuleAtPixels(e.x, e.y, true);
    if( module == NULL )           
    {
      tmpAudioConnection.resetToNull();   //...and dismiss it now
      repaint();
      return;
    }

    // check if the mouseUp was on a pin:
    int kindOfPin, directionOfPin, indexOfPin;
    bool wasOnPin = getPinPropertiesAtPixels(e.x, e.y, module, getRectangleForModuleInPixels(module, true), 
      kindOfPin, directionOfPin, indexOfPin);
    if( wasOnPin )
    {
      if( kindOfPin == romos::AUDIO )
      {
        if( directionOfPin == romos::INCOMING && tmpAudioConnection.getTargetModule() == NULL )
        {
          tmpAudioConnection.setTargetModule(module);
          tmpAudioConnection.setTargetInputIndex(indexOfPin);
          getInterfaceMediator()->getContainerShownInDiagram()->addAudioConnection(&tmpAudioConnection);
          tmpAudioConnection.resetToNull();
          repaint();
        }
        else if( directionOfPin == romos::OUTGOING && tmpAudioConnection.getSourceModule() == NULL )
        {
          tmpAudioConnection.setSourceModule(module);
          tmpAudioConnection.setSourceOutputIndex(indexOfPin);
          getInterfaceMediator()->getContainerShownInDiagram()->addAudioConnection(&tmpAudioConnection);
          tmpAudioConnection.resetToNull();
          repaint();
        }
      }
    }

    tmpAudioConnection.resetToNull();
    updateAudioConnectionArray();
    repaint();
    return;
  }
  else if( lassoComponent->isVisible() )  // we were openenig a lasso selector
  {
    modulesInLasso = getModulesInRectangle(lassoComponent->getBounds());
    for(unsigned int i = 0; i < modulesInLasso.size(); i++)
      addModuleWithConnectionsToArray(modulesInLasso[i], selectedModules, selectedAudioConnections);
    lassoComponent->setVisible(false);
    modulesInLasso.clear();
    audioConnectionsInLasso.clear();
    repaint();
  }
  else  // we were dragging a bunch of modules around
  {
    for(unsigned int i = 0; i < selectedModules.size(); i++)
    {
      selectedModules[i]->setPositionXY(selectedModules[i]->getPositionX() + inPinDistances(selectionOffsetX), 
        selectedModules[i]->getPositionY() + inPinDistances(selectionOffsetY));
    }
    selectionOffsetX = 0;
    selectionOffsetY = 0;
    updateCanvasBounds(0, 0, availableWidth, availableHeight);
  }
  //repaint();
}

void ModularBlockDiagramPanel::paint(Graphics &g)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  //g.fillAll(Colours::lavender);   // \todo obtain the background colour from the colorscheme 
  g.fillAll(getPlotColourScheme().topLeft);  
    // we should actually call drawBilinearGradient with the colors for the 4 corners

  drawGrid(g);
  drawDiagram(g);
}

void ModularBlockDiagramPanel::paintOverChildren(Graphics &g)
{
  //ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  // do nothing - avoid outline drawing from the baseclass
}

/*
void ModularBlockDiagramPanel::resized()
{
ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
ColourSchemeComponent::resized();
}
*/

/*
void ModularBlockDiagramPanel::updateDiagram()
{
ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
repaint();
}
*/

void ModularBlockDiagramPanel::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( menuThatHasChanged == actOnSelectionMenu )
  {
    RTreeViewNode *selectedItem = actOnSelectionMenu->getSelectedItem();
    if( selectedItem != NULL )
    {
      switch( selectedItem->getNodeIdentifier() )
      {
      case EDIT_NAME:        openModuleNameEntryField();      break;
      case SAVE_CONTAINER:   openContainerSaveDialog();       break;
      case EXPORT_TO_CODE:   openExportToCodeDialog();        break;
      case SET_POLYPHONIC:   setPolyphonyForSelection(true);  break;
      case SET_MONOPHONIC:   setPolyphonyForSelection(false); break;
      case DELETE_SELECTION: deleteSelection();               break;
      case CONTAINERIZE:     containerizeSelection();         break;
      case UNCONTAINERIZE:   unContainerizeSelection();       break;
        //case 3: minimizeNumberOfInputs();  break;
      default:
      {
        int dummy = 0;
        jassertfalse;
      }
      }
    }
    else
      jassertfalse;
    actOnSelectionMenu->dismiss();
  }
  else
    jassertfalse;
}

void ModularBlockDiagramPanel::textChanged(RTextEntryField *rTextEntryFieldThatHasChanged)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( selectedModules.size() == 1 )
  {
    char* newName = toZeroTerminatedString(nameEntryField->getText());
    selectedModules[0]->setModuleName(std::string(newName));
    nameEntryField->setVisible(false);
    grabKeyboardFocus(); // to move it away from the entry field
    delete newName;
    //getInterfaceMediator()->sendModuleChangeNotification(selectedModules[0], MODULE_NAME);
    notifyMediator(MODULE_NAME);
    //repaint();
  }
}

void ModularBlockDiagramPanel::treeNodeClicked(RTreeView *treeView, RTreeViewNode *nodeThatWasClicked, const MouseEvent &mouseEvent, 
  int clickPosition)
{
  WRITE_TO_LOGFILE("ModularBlockDiagramPanel::treeNodeClicked entered\n");
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( treeView == availableModulesTreeView && nodeThatWasClicked->isLeafNode() )
  {
    availableModulesTreeView->setVisible(false);
    availableModulesTreeView->removeFromDesktop();
    //if( nodeThatWasClicked->getNodeIdentifier() == LOAD_CONTAINER )
    if( nodeThatWasClicked->getNodeText() == "Load Container..." )
      openContainerLoadDialog();
    else
    {
      // old:
      //insertModule(nodeThatWasClicked->getNodeIdentifier(), inPinDistances(mouseDownX), inPinDistances(mouseDownY));

      //new:
      insertModule(nodeThatWasClicked->getNodeText(), inPinDistances(mouseDownX), inPinDistances(mouseDownY));
    }
  }
  WRITE_TO_LOGFILE("ModularBlockDiagramPanel::treeNodeClicked finished\n");
}

void ModularBlockDiagramPanel::treeNodeChanged(RTreeView *treeView, RTreeViewNode *nodeThatHasChanged)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));


}

void ModularBlockDiagramPanel::handleMediatorNotification(MediatedColleague *originatingColleague, 
  int messageCode, rsMessageData* messageData)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  selectedAudioConnections.clear();
  selectedModules.clear();
  updateAudioConnectionArray();
  // maybe clear more arrays?

  if( messageCode == MODULE_TO_SHOW_EDITOR_FOR )
  {
    romos::ContainerModule *containerShownInDiagram = getInterfaceMediator()->getContainerShownInDiagram();
    romos::Module          *moduleToShowEditorFor   = getInterfaceMediator()->getModuleToShowEditorFor();
    if( moduleToShowEditorFor != containerShownInDiagram )
      selectSingleModule(moduleToShowEditorFor);
  }


  updateCanvasBounds(0, 0, availableWidth, availableHeight);
  repaint(); // \todo: put in an if( modulePropertiesThatHasChanged == focusedModule ), maybe restrict on the causes also
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void ModularBlockDiagramPanel::openModuleInsertionMenu(int x, int y)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  availableModulesTreeView->setTopLeftPosition(x, y);
  availableModulesTreeView->addToDesktop(ComponentPeer::windowIsTemporary, 0);
  availableModulesTreeView->setVisible(true);
}

void ModularBlockDiagramPanel::openActOnSelectionMenu(int x, int y)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  actOnSelectionMenu->clear();
  if( selectedModules.size() == 1 )
  {
    actOnSelectionMenu->addItem(EDIT_NAME, ("Edit Name"), true, false);
    //if( selectedModules[0]->getTypeIdentifierOld() == romos::ModuleTypeRegistry::CONTAINER )
    if( selectedModules[0]->getTypeName() == "Container" )
    {
      actOnSelectionMenu->addItem(SAVE_CONTAINER, ("Save Container..."), true, false);
      actOnSelectionMenu->addItem(EXPORT_TO_CODE, ("Export to Code..."), true, false);
    }
  }


  actOnSelectionMenu->addItem(SET_POLYPHONIC,   ("Set Polyphonic"),         true, false);
  actOnSelectionMenu->addItem(SET_MONOPHONIC,   ("Set Monophonic"),         true, false);
  actOnSelectionMenu->addItem(DELETE_SELECTION, ("Delete Selection"),       true, false);
  actOnSelectionMenu->addItem(CONTAINERIZE,     ("Containerize Selection"), true, false);
  actOnSelectionMenu->addItem(UNCONTAINERIZE,   ("Unpack Container(s)"),    true, false);
  //actOnSelectionMenu->addItem(MINIMIZE_PINS,  ("Minimize Pins"),    true, false);

  actOnSelectionMenu->setBounds(100, 100, 150, 200); // \todo make the size adapt to the content
  actOnSelectionMenu->showAtMousePosition(false);
}

void ModularBlockDiagramPanel::openModuleNameEntryField()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( selectedModules.size() == 1 )
  {
    juce::String name = juce::String(selectedModules[0]->getName().c_str());
    nameEntryField->setText(name);
    int x, y, w, h;
    getRectangleForModuleInPixels(selectedModules[0], x, y, w, h, false);
    //nameEntryField->setBounds(x, y, w, getModuleTitleHeightInPixels(selectedModules[0]));
    nameEntryField->setBounds(x, y, jmax(w, 40), 2*t+2*m+bigFontHeight);

    //int type = selectedModules[0]->getTypeIdentifierOld();
    int typeId = selectedModules[0]->getTypeId();
    std::string typeName = selectedModules[0]->getTypeName();
    if(  !romos::moduleFactory.isModuleNameEditable(typeId) )
      return;  
    //if( getInterfaceMediator()->getContainerShownInDiagram()->getTypeIdentifierOld() == romos::ModuleTypeRegistry::TOP_LEVEL_MODULE )
    if( getInterfaceMediator()->getContainerShownInDiagram()->isTopLevelModule() )
    {
      //if( romos::ModuleTypeRegistry::isIdentifierInputOrOutput(type) )
      if( selectedModules[0]->isInputOrOutput() )
        return; // disallow editing of toplevel I/O module names
    }

    //if( type == romos::ModuleTypeRegistry::CONSTANT )
    if( typeName == "Constant" )
      nameEntryField->setPermittedCharacters(juce::String(("0123456789.-")));
    else
      nameEntryField->setPermittedCharacters(juce::String(("0123456789.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")));

    nameEntryField->setVisible(true);
    nameEntryField->grabKeyboardFocus();
  }
}

void ModularBlockDiagramPanel::openContainerLoadDialog()
{
  WRITE_TO_LOGFILE("ModularBlockDiagramPanel::openContainerLoadDialog entered\n");
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  juce::File   initialDirectory = getInterfaceMediator()->modularSynthModuleToEdit->macroDirectory;
  juce::String allowedPatterns  = juce::String(("*.xml"));
  FileChooser chooser(juce::String(("Load Container")), initialDirectory, allowedPatterns, true);
  if( chooser.browseForFileToOpen() )
  {
    juce::File fileToLoad = chooser.getResult();
    getInterfaceMediator()->modularSynthModuleToEdit->macroDirectory = fileToLoad.getParentDirectory();
    XmlElement *xmlState = getXmlFromFile(fileToLoad);
    WRITE_TO_LOGFILE("Module state was loaded\n");
    if( xmlState != NULL )
    {
      //romos::ModuleContainer *newModule = new romos::ModuleContainer(NULL);

      // deprecated:
      //romos::ModuleContainer *newModule = (ModuleContainer*) ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER);

      // later use:
      romos::ContainerModule *newModule = (ContainerModule*) moduleFactory.createModule("Container");

      WRITE_TO_LOGFILE("Module created\n");
      LibertyAudioModule::setModuleStateFromXml(*xmlState, newModule);
      WRITE_TO_LOGFILE("Module state was set\n");
      newModule->setPositionXY(inPinDistances(mouseDownX), inPinDistances(mouseDownY), false);
      WRITE_TO_LOGFILE("Module position was set\n");
      getInterfaceMediator()->getContainerShownInDiagram()->addChildModule(newModule, true);
      //getInterfaceMediator()->sendModuleChangeNotification(getInterfaceMediator()->getContainerShownInDiagram(), NUM_CHILDREN);
      notifyMediator(NUM_CHILDREN);
      delete xmlState;
    }
  }
  WRITE_TO_LOGFILE("ModularBlockDiagramPanel::openContainerLoadDialog finished\n");
}

void ModularBlockDiagramPanel::openContainerSaveDialog()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  //if( selectedModules.size() == 1 && selectedModules[0]->getTypeIdentifierOld() == ModuleTypeRegistry::CONTAINER )
  if( selectedModules.size() == 1 && selectedModules[0]->isContainerModule() ) 
  {
    //juce::File   initialDirectory = rojue::getApplicationDirectory(); // preliminary - we should have a member containerDirectory or sth
    juce::File   initialDirectory = getInterfaceMediator()->modularSynthModuleToEdit->macroDirectory;
    juce::String extension        = juce::String(("xml"));
    juce::String allowedPatterns  = juce::String(("*.xml"));
    FileChooser chooser(juce::String(("Save Container")), initialDirectory, allowedPatterns, true);
    if( chooser.browseForFileToSave(true) )
    {
      juce::File fileToSaveTo = chooser.getResult();
      getInterfaceMediator()->modularSynthModuleToEdit->macroDirectory = fileToSaveTo.getParentDirectory();
      if ( !fileToSaveTo.hasFileExtension(extension) )
        fileToSaveTo = fileToSaveTo.withFileExtension( extension ) ;
      XmlElement *xmlState = LibertyAudioModule::getModuleStateAsXml(selectedModules[0], false);
      if( xmlState != NULL )
      {
        saveXmlToFile(*xmlState, fileToSaveTo, true);
        delete xmlState;
      }
    }
  }
}

void ModularBlockDiagramPanel::openExportToCodeDialog()
{
  jassertfalse; // to be re-activated
  //ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  //if( selectedModules.size() == 1 && selectedModules[0]->getTypeIdentifier() == ModuleTypeRegistry::CONTAINER )
  //{
  //  juce::File   initialDirectory = getInterfaceMediator()->modularSynthModuleToEdit->macroDirectory;
  //  juce::String extension        = juce::String(("txt"));
  //  juce::String allowedPatterns  = juce::String(("*.txt"));
  //  FileChooser chooser(juce::String(("Export to Code")), initialDirectory, allowedPatterns, true);
  //  if( chooser.browseForFileToSave(true) )
  //  {
  //    juce::File fileToSaveTo = chooser.getResult();
  //    //if ( !fileToSaveTo.hasFileExtension(extension) )
  //    //  fileToSaveTo = fileToSaveTo.withFileExtension(extension) ;
  //    if( fileToSaveTo.existsAsFile() )
  //    {
  //      fileToSaveTo.deleteFile();
  //      fileToSaveTo.create();
  //    }
  //    juce::String codeString = rosicToJuce(ModuleBuildCodeGenerator::getCodeForModule(selectedModules[0]));
  //    fileToSaveTo.appendText(codeString);
  //  }
  //}
}



/*
void ModularBlockDiagramPanel::openPropertiesDialog()
{
ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
if( selectedModules.getNumElements() == 1 )
{
propertiesDialog->retrievePropertiesFromModule(selectedModules[0]);
propertiesDialog->setVisible(true);
}
}
*/

/*
// old:
void ModularBlockDiagramPanel::insertModule(int moduleIdentifer, int xInPinDistances, int yInPinDistances)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  getInterfaceMediator()->getContainerShownInDiagram()->addChildModule(moduleIdentifer, 
    rosic::rsString(),                                               // will fall back to a default name based on the module-type
    xInPinDistances, yInPinDistances, 
    getInterfaceMediator()->getContainerShownInDiagram()->isPolyphonic()  // use same polyphony as the outlying container 
    );

  //getInterfaceMediator()->sendModuleChangeNotification(getInterfaceMediator()->getContainerShownInDiagram(), NUM_CHILDREN);
  notifyMediator(NUM_CHILDREN);
}
*/

// new:
void ModularBlockDiagramPanel::insertModule(const juce::String& typeName, int x, int y)
{
  LibertyInterfaceMediator* med = getInterfaceMediator();
  ScopedLock scopedLock(*(med->plugInLock));
  med->getContainerShownInDiagram()->addChildModule(typeName.toStdString(), 
    "",                                                  // use a default name based on the module-type
    x, y, 
    med->getContainerShownInDiagram()->isPolyphonic());  // use same polyphony as the outlying container 
  notifyMediator(NUM_CHILDREN);
}

void ModularBlockDiagramPanel::selectSingleModule(romos::Module *moduleToSelect)
{  
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  selectedAudioConnections.clear();
  selectedModules.clear();
  addModuleWithConnectionsToArray(moduleToSelect, selectedModules, selectedAudioConnections);
}

void ModularBlockDiagramPanel::addModuleWithConnectionsToArray(romos::Module *moduleToAdd, 
  std::vector<romos::Module*> &moduleArray, 
  std::vector<romos::AudioConnection> &connectionArray)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  rosic::appendIfNotAlreadyThere(moduleArray, moduleToAdd);

  // wrap these 3-liners into a functions rosic::appendArrayElementsIfNotAlreadyThere(vector, vector):

  if( !moduleToAdd->isInputModule() )
  {
    std::vector<AudioConnection> incomingAudioConnections = moduleToAdd->getIncomingAudioConnections();
    for(unsigned int i = 0; i < incomingAudioConnections.size(); i++)
      rosic::appendIfNotAlreadyThere(connectionArray, incomingAudioConnections[i]);
  }

  std::vector<AudioConnection> outgoingAudioConnections = moduleToAdd->getOutgoingAudioConnections();
  for(unsigned int i = 0; i < outgoingAudioConnections.size(); i++)
    rosic::appendIfNotAlreadyThere(connectionArray, outgoingAudioConnections[i]);
}

void ModularBlockDiagramPanel::removeModuleWithConnectionsFromArray(romos::Module *moduleToRemove, 
  std::vector<romos::Module*> &moduleArray, 
  std::vector<romos::AudioConnection> &connectionArray)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  rosic::removeElementByValue(moduleArray, moduleToRemove);
  romos::AudioConnection ac;
  unsigned int i;

  /*
  for(i = 0; i < moduleToRemove->getNumIncomingAudioConnections(); i++)
  {
  ac = moduleToRemove->getIncomingAudioConnection(i);
  if( !rosic::containsElement(moduleArray, ac->getSourceModule()) )
  rosic::removeElementByValue(connectionArray, ac);
  }
  */


  std::vector<AudioConnection> incomingAudioConnections = moduleToRemove->getIncomingAudioConnections(); 
  for(i = 0; i < incomingAudioConnections.size(); i++)
  {
    ac = incomingAudioConnections[i];
    if( !rosic::containsElement(moduleArray, ac.getSourceModule()) )
      rosic::removeElementByValue(connectionArray, ac);
  }

  std::vector<AudioConnection> outgoingAudioConnections = moduleToRemove->getOutgoingAudioConnections(); 
  for(i = 0; i < outgoingAudioConnections.size(); i++)
  {
    ac = outgoingAudioConnections[i];
    if( !rosic::containsElement(moduleArray, ac.getTargetModule()) )
      rosic::removeElementByValue(connectionArray, ac);
  }
}

void ModularBlockDiagramPanel::setPolyphonyForSelection(bool shouldBePolyphonic)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock)); 
  getInterfaceMediator()->getContainerShownInDiagram()->setPolyphonyForModules(selectedModules, shouldBePolyphonic);
  notifyMediator(POLYPHONY);
}

void ModularBlockDiagramPanel::deleteSelection()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock)); 
  getInterfaceMediator()->getContainerShownInDiagram()->deleteAudioConnections(selectedAudioConnections);
  getInterfaceMediator()->getContainerShownInDiagram()->deleteModules(selectedModules);
  //getInterfaceMediator()->sendModuleChangeNotification(getInterfaceMediator()->getContainerShownInDiagram(), NUM_CHILDREN);

  getInterfaceMediator()->setModuleToShowEditorFor(getInterfaceMediator()->getContainerShownInDiagram());
  notifyMediator(NUM_CHILDREN);
}

void ModularBlockDiagramPanel::containerizeSelection()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  romos::ContainerModule *createdContainer = getInterfaceMediator()->getContainerShownInDiagram()->containerizeModules(selectedModules);
  updateAudioConnectionArray();
  selectedModules.clear();
  addModuleWithConnectionsToArray(createdContainer, selectedModules, selectedAudioConnections); // select the just created container
                                                                                                //getInterfaceMediator()->sendModuleChangeNotification(getInterfaceMediator()->getContainerShownInDiagram(), NUM_CHILDREN);

  getInterfaceMediator()->setModuleToShowEditorFor(getInterfaceMediator()->getContainerShownInDiagram());
  notifyMediator(NUM_CHILDREN);
}

void ModularBlockDiagramPanel::unContainerizeSelection()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  // keep track of what's being unpacked in order to later select it:
  std::vector<romos::Module*> unpackedModules;
  for(unsigned int i = 0; i < selectedModules.size(); i++)
  {
    ContainerModule *tmpContainer = dynamic_cast<ContainerModule*> (selectedModules[i]);
    if( tmpContainer != NULL )
      rosic::appendVector(unpackedModules, tmpContainer->getNonInOutChildModules());
    //unpackedModules.appendArray( tmpContainer->getNonInOutChildModules() );
  }

  getInterfaceMediator()->getContainerShownInDiagram()->unContainerizeModules(selectedModules);
  updateAudioConnectionArray();  
  selectedModules.clear();

  getInterfaceMediator()->setModuleToShowEditorFor(getInterfaceMediator()->getContainerShownInDiagram());
  notifyMediator(NUM_CHILDREN);  

  for(unsigned int i = 0; i < unpackedModules.size(); i++) // maybe wrap this loop into a function
    addModuleWithConnectionsToArray(unpackedModules[i], selectedModules, selectedAudioConnections);
}

void ModularBlockDiagramPanel::minimizeNumberOfInputs()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  for(unsigned int i = 0; i < selectedModules.size(); i++)
  {
    ContainerModule *tmpContainer = dynamic_cast<ContainerModule*> (selectedModules[i]);
    if( tmpContainer != NULL )
      tmpContainer->minimizeNumberOfAudioInputs();
  }
  //getInterfaceMediator()->sendModuleChangeNotification(getInterfaceMediator()->getContainerShownInDiagram(), NUM_CONNECTIONS);
  notifyMediator(NUM_CONNECTIONS);
}

void ModularBlockDiagramPanel::fillAvailableModulesTreeView()
{
  // ToDo: this code should be replaced with some code that loops through the module types in the
  // global moduleTypeRegistry object

  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  treeRootNode = new RTreeViewNode("Module", 0);
  treeRootNode->setDeleteChildNodesOnDestruction(true);
  // -> all the child-nodes that are subsequently created here via "new" (and added to the tree) will be deleted when the treeRootNode 
  // itself is deleted

  // \todo write a short description for each node to appear in the infoline

  RTreeViewNode *tmpNode1; //, *tmpNode2;  // we need one temporary node for each tree-level except the lowest ...or not?


  tmpNode1 = new RTreeViewNode("Load Container...", LOAD_CONTAINER);
  treeRootNode->addChildNode(tmpNode1);

  RTreeViewNode *insertModuleNode = new RTreeViewNode("Insert");
  treeRootNode->addChildNode(insertModuleNode);

  // new - activate soon:: 
  for(size_t i = 0; i < moduleFactory.getNumModuleTypes(); i++)
  {
    romos::ModuleTypeInfo* typeInfo = moduleFactory.getModuleTypeInfo(i);
    juce::String typeName = typeInfo->fullName;
    juce::String category = typeInfo->category;
    RTreeViewNode* categNode = insertModuleNode->findNodeByText(category);
    if(categNode == nullptr) {
      categNode = new RTreeViewNode(category);
      insertModuleNode->addChildNode(categNode);
    }
    categNode->addChildNode(new RTreeViewNode(typeName, typeInfo->id));
    int dummy = 0;
  }
  // this loop should replace (almost) all the code below


  /*
  //---------------------------------------------------------------------------
  // test modules:

  tmpNode1 = new RTreeViewNode("Test Modules");
  tmpNode1->addChildNode(new RTreeViewNode("Gain",            romos::ModuleTypeRegistry::TEST_GAIN));
  tmpNode1->addChildNode(new RTreeViewNode("SumDiff",         romos::ModuleTypeRegistry::TEST_SUM_DIFF));
  tmpNode1->addChildNode(new RTreeViewNode("WrappedSumDiff",  romos::ModuleTypeRegistry::TEST_WRAPPED_SUM_DIFF));
  tmpNode1->addChildNode(new RTreeViewNode("SummedDiffs",     romos::ModuleTypeRegistry::TEST_SUMMED_DIFFS));
  tmpNode1->addChildNode(new RTreeViewNode("MovingAverage",   romos::ModuleTypeRegistry::TEST_MOVING_AVERAGE));
  tmpNode1->addChildNode(new RTreeViewNode("LeakyIntegrator", romos::ModuleTypeRegistry::TEST_LEAKY_INTEGRATOR));
  tmpNode1->addChildNode(new RTreeViewNode("TestFilter1",     romos::ModuleTypeRegistry::TEST_FILTER1));
  tmpNode1->addChildNode(new RTreeViewNode("Biquad",          romos::ModuleTypeRegistry::TEST_BIQUAD));
  tmpNode1->addChildNode(new RTreeViewNode("AddedConstants",  romos::ModuleTypeRegistry::TEST_ADDED_CONSTANTS));
  tmpNode1->addChildNode(new RTreeViewNode("PinSorting",      romos::ModuleTypeRegistry::TEST_PIN_SORTING));
  tmpNode1->addChildNode(new RTreeViewNode("Blip",            romos::ModuleTypeRegistry::TEST_BLIP));
  tmpNode1->addChildNode(new RTreeViewNode("PolyBlipStereo",  romos::ModuleTypeRegistry::TEST_POLY_BLIP_STEREO));
  tmpNode1->addChildNode(new RTreeViewNode("NoiseFlute",      romos::ModuleTypeRegistry::TEST_NOISE_FLUTE));

  //tmpNode1->addChildNode(new RTreeViewNode(("Moog Filter"),      romos::ModuleTypeRegistry::EXAMPLE_MOOG_FILTER));
  //tmpNode1->addChildNode(new RTreeViewNode(("Containerize"),      romos::ModuleTypeRegistry::TEST_CONTAINERIZE));
  //tmpNode1->addChildNode(new RTreeViewNode(("UnContainerize"),    romos::ModuleTypeRegistry::TEST_UNCONTAINERIZE));
  //tmpNode1->addChildNode(new RTreeViewNode(("Minimize Inputs 1"), romos::ModuleTypeRegistry::TEST_MINIMIZE_INS1));

  insertModuleNode->addChildNode(tmpNode1);

  //---------------------------------------------------------------------------
  // special modules:

  tmpNode1 = new RTreeViewNode(("Infrastructural"));
  tmpNode1->addChildNode(new RTreeViewNode(("Parameter"),              romos::ModuleTypeRegistry::PARAMETER));
  tmpNode1->addChildNode(new RTreeViewNode(("Container"),              romos::ModuleTypeRegistry::CONTAINER));
  tmpNode1->addChildNode(new RTreeViewNode(("Audio Input"),            romos::ModuleTypeRegistry::AUDIO_INPUT));
  tmpNode1->addChildNode(new RTreeViewNode(("Audio Output"),           romos::ModuleTypeRegistry::AUDIO_OUTPUT));
  tmpNode1->addChildNode(new RTreeViewNode(("Voice Combiner"),         romos::ModuleTypeRegistry::VOICE_COMBINER));
  tmpNode1->addChildNode(new RTreeViewNode(("SampleRate"),             romos::ModuleTypeRegistry::SYSTEM_SAMPLE_RATE));
  //tmpNode1->addChildNode(new RTreeViewNode(("SamplePeriod"),           romos::ModuleTypeRegistry::SYSTEM_SAMPLE_PERIOD));

  // insert event I/O here
  insertModuleNode->addChildNode(tmpNode1);

  //---------------------------------------------------------------------------
  // event modules:

  tmpNode1 = new RTreeViewNode(("Events"));
  tmpNode1->addChildNode(new RTreeViewNode(("NoteGate"),              romos::ModuleTypeRegistry::NOTE_GATE));
  tmpNode1->addChildNode(new RTreeViewNode(("NoteOnTrigger"),         romos::ModuleTypeRegistry::NOTE_ON_TRIGGER));
  tmpNode1->addChildNode(new RTreeViewNode(("NoteOffTrigger"),        romos::ModuleTypeRegistry::NOTE_OFF_TRIGGER));
  tmpNode1->addChildNode(new RTreeViewNode(("VoiceKiller"),           romos::ModuleTypeRegistry::VOICE_KILLER));
  tmpNode1->addChildNode(new RTreeViewNode(("NoteFrequency"),         romos::ModuleTypeRegistry::NOTE_FREQUENCY));
  tmpNode1->addChildNode(new RTreeViewNode(("NoteVelocity"),          romos::ModuleTypeRegistry::NOTE_VELOCITY));
  insertModuleNode->addChildNode(tmpNode1);

  //---------------------------------------------------------------------------
  // signal generator modules:

  tmpNode1 = new RTreeViewNode(("Signal Generators"));
  tmpNode1->addChildNode(new RTreeViewNode(("Phasor"),                  romos::ModuleTypeRegistry::PHASOR));
  tmpNode1->addChildNode(new RTreeViewNode(("WhiteNoise"),              romos::ModuleTypeRegistry::WHITE_NOISE));
  tmpNode1->addChildNode(new RTreeViewNode(("BandlimitedImpulseTrain"), romos::ModuleTypeRegistry::BANDLIMITED_IMPULSE_TRAIN));
  tmpNode1->addChildNode(new RTreeViewNode(("BlitSaw"),                 romos::ModuleTypeRegistry::BLIT_SAW_OSCILLATOR)); 
  tmpNode1->addChildNode(new RTreeViewNode(("DualBlitSaw"),             romos::ModuleTypeRegistry::DUAL_BLIT_SAW_OSCILLATOR));

  insertModuleNode->addChildNode(tmpNode1);

  //---------------------------------------------------------------------------
  // modulator modules:

  tmpNode1 = new RTreeViewNode(("Modulators"));
  tmpNode1->addChildNode(new RTreeViewNode(("EnvelopeADSR"),   romos::ModuleTypeRegistry::ENVELOPE_ADSR));
  insertModuleNode->addChildNode(tmpNode1);


  //---------------------------------------------------------------------------
  // arithmetic modules:

  tmpNode1 = new RTreeViewNode(("Arithmetic"));

  // unary:
  tmpNode2 = new RTreeViewNode(("Unary"));
  tmpNode2->addChildNode(new RTreeViewNode(("Constant"),   romos::ModuleTypeRegistry::CONSTANT));
  tmpNode2->addChildNode(new RTreeViewNode(("UnaryMinus"), romos::ModuleTypeRegistry::UNARY_MINUS));
  tmpNode2->addChildNode(new RTreeViewNode(("Reciprocal"), romos::ModuleTypeRegistry::RECIPROCAL));
  //tmpNode2->addChildNode(new RTreeViewNode(("Formula"),  romos::ModuleTypeRegistry::UNARY_FORMULA));
  tmpNode1->addChildNode(tmpNode2);

  // binary:
  tmpNode2 = new RTreeViewNode(("Binary"));
  tmpNode2->addChildNode(new RTreeViewNode(("Adder"),      romos::ModuleTypeRegistry::ADDER));
  tmpNode2->addChildNode(new RTreeViewNode(("Subtractor"), romos::ModuleTypeRegistry::SUBTRACTOR));
  tmpNode2->addChildNode(new RTreeViewNode(("Multiplier"), romos::ModuleTypeRegistry::MULTIPLIER));
  tmpNode2->addChildNode(new RTreeViewNode(("Divider"),    romos::ModuleTypeRegistry::DIVIDER));
  tmpNode2->addChildNode(new RTreeViewNode(("AdderN"),     romos::ModuleTypeRegistry::ADDER_N));
  //tmpNode2->addChildNode(new RTreeViewNode(("Formula"),  romos::ModuleTypeRegistry::BINARY_FORMULA));
  tmpNode1->addChildNode(tmpNode2);

  // n-ary:
  //tmpNode2 = new RTreeViewNode(("Multiple Ins"));
  //tmpNode2->addChildNode(new RTreeViewNode(("Sum"),     romos::ModuleTypeRegistry::SUM));
  //tmpNode2->addChildNode(new RTreeViewNode(("Product"), romos::ModuleTypeRegistry::PRODUCT));
  //tmpNode2->addChildNode(new RTreeViewNode(("Formula"), romos::ModuleTypeRegistry::MULTI_IN_FORMULA));
  //tmpNode1->addChildNode(tmpNode2);

  // n ins, m outs:
  //tmpNode2 = new RTreeViewNode(("Multiple Ins and Outs"));
  //tmpNode2->addChildNode(new RTreeViewNode(("Matrix"),        romos::ModuleTypeRegistry::MATRIX));
  //tmpNode2->addChildNode(new RTreeViewNode(("Formula Array"), romos::ModuleTypeRegistry::FORMULA_ARRAY));
  //tmpNode1->addChildNode(tmpNode2);

  insertModuleNode->addChildNode(tmpNode1);

  //---------------------------------------------------------------------------
  // functional modules:

  tmpNode1 = new RTreeViewNode("Functions");
  tmpNode1->addChildNode(new RTreeViewNode("Clipper", romos::ModuleTypeRegistry::CLIPPER));
  tmpNode1->addChildNode(new RTreeViewNode("SinCos",  romos::ModuleTypeRegistry::SIN_COS));
  tmpNode1->addChildNode(new RTreeViewNode("TriSaw",  romos::ModuleTypeRegistry::TRISAW));
  insertModuleNode->addChildNode(tmpNode1);



  //---------------------------------------------------------------------------
  // delay modules:

  tmpNode1 = new RTreeViewNode(("Delays"));
  tmpNode1->addChildNode(new RTreeViewNode(("Unit Delay"),              romos::ModuleTypeRegistry::UNIT_DELAY));
  //tmpNode1->addChildNode(new RTreeViewNode(("Integer Delay"),           romos::ModuleTypeRegistry::INTEGER_DELAY));
  //tmpNode1->addChildNode(new RTreeViewNode(("Tapped Integer Delay"),    romos::ModuleTypeRegistry::TAPPED_INTEGER_DELAY));
  //tmpNode1->addChildNode(new RTreeViewNode(("Fractional Delay"),        romos::ModuleTypeRegistry::FRACTIONAL_DELAY));
  //tmpNode1->addChildNode(new RTreeViewNode(("Tapped Fractional Delay"), romos::ModuleTypeRegistry::TAPPED_FRACTIONAL_DELAY));
  insertModuleNode->addChildNode(tmpNode1);

  //---------------------------------------------------------------------------
  // filter modules:
  // \todo use another hierarch level here to distinguish several classes of filters

  tmpNode1 = new RTreeViewNode(("Filters"));
  tmpNode1->addChildNode(new RTreeViewNode(("FirstOrderLowpass"), romos::ModuleTypeRegistry::FIRST_ORDER_LOWPASS));
  tmpNode1->addChildNode(new RTreeViewNode(("FirstOrderFilter"),  romos::ModuleTypeRegistry::FIRST_ORDER_FILTER));
  tmpNode1->addChildNode(new RTreeViewNode(("Biquad"),            romos::ModuleTypeRegistry::BIQUAD));
  tmpNode1->addChildNode(new RTreeViewNode(("BiquadDesigner"),    romos::ModuleTypeRegistry::BIQUAD_DESIGNER));
  tmpNode1->addChildNode(new RTreeViewNode(("LadderFilter"),      romos::ModuleTypeRegistry::LADDER_FILTER));
  insertModuleNode->addChildNode(tmpNode1);
  */







  availableModulesTreeView->setRootNode(treeRootNode);
  availableModulesTreeView->setBounds(0, 0, 200, 400);
}

void ModularBlockDiagramPanel::updateAudioConnectionArray()
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  selectedAudioConnections.clear();
  allAudioConnections.clear();

  unsigned int i;
  romos::Module *m;
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumChildModules(); i++)
  {
    m = getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i);
    rosic::appendVector(allAudioConnections, m->getIncomingAudioConnections());

    //for(j = 0; j < m->getNumIncomingAudioConnections(); j++)
    //  rosic::appendElement(allAudioConnections, m->getIncomingAudioConnection(j));
    //allAudioConnections.appendElement(m->getIncomingAudioConnection(j));

  }

  /*
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumOutputPins(); i++)
  {
  m = getInterfaceMediator()->getContainerShownInDiagram()->getAudioOutputModule(i);
  for(j = 0; j < m->getNumIncomingAudioConnections(); j++)
  rosic::appendElement(allAudioConnections, m->getIncomingAudioConnection(j));
  //allAudioConnections.appendElement(m->getIncomingAudioConnection(j));
  }
  */
}

inline void setPixel(Graphics &g, int x, int y)
{
  // preliminary - try to find a more efficient way ...and move to GraphicsTools
  g.fillRect(x, y, 1, 1);
}

void ModularBlockDiagramPanel::drawGrid(Graphics &g)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  int numVerticalLines   = inPinDistances(getWidth())+1;
  int numHorizontalLines = inPinDistances(getHeight())+1;

  //g.setColour(Colours::black);  // preliminary
  g.setColour(getPlotColourScheme().coarseGrid);  

  if( gridStyle == GRID_LINES )
  {
    for(int x=1; x<numVerticalLines; x++)
      g.drawVerticalLine(inPixels(x), 0.f, (float) getHeight());
    for(int y=1; y<numHorizontalLines; y++)
      g.drawHorizontalLine(inPixels(y), 0.f, (float) getWidth());
  }
  else if( gridStyle == DOTTED_GRID )
  {
    for(int x=1; x<numVerticalLines; x++)
    {
      for(int y=1; y<numHorizontalLines; y++)
      {
        //g.setPixel(inPixels(x), inPixels(y)); // not available in juce 5.2.0 anymore
        //g.fillRect(inPixels(x), inPixels(y), 1, 1);
        setPixel(g, inPixels(x), inPixels(y));
      }
    }
  }
}

void ModularBlockDiagramPanel::drawDiagram(Graphics &g)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  romos::ContainerModule *mc = getInterfaceMediator()->getContainerShownInDiagram();
  unsigned int i;

  // draw the embedded modules:
  int test = mc->getNumInputPins();

  for(i = 0; i < mc->getNumInputPins(); i++)
    drawModule(g, mc->getAudioInputModule(i));
  for(i = 0; i < mc->getNumChildModules(); i++)
    drawModule(g, mc->getChildModule(i));
  for(i = 0; i < mc->getNumOutputPins(); i++)
    drawModule(g, mc->getAudioOutputModule(i));

  // draw the connections:
  for(i = 0; i < mc->getNumInputPins(); i++)
    drawIncomingConnectionsForModule(g, mc->getAudioInputModule(i));
  for(i = 0; i < mc->getNumChildModules(); i++)
    drawIncomingConnectionsForModule(g, mc->getChildModule(i));
  for(i = 0; i < mc->getNumOutputPins(); i++)
    drawIncomingConnectionsForModule(g, mc->getAudioOutputModule(i));

  // \todo draw event I/O modules ...maybe we can collapse these loops into one when Module provides functions like 
  // getNumContainedModules/getContainedModule - these would inlcude all the enclosed modules (child/I/O)

  // draw temporary connection:
  if( !tmpAudioConnection.isNull() )
  {
    Point<int> mouse = getMouseXYRelative();
    int x1 = snapPixelPositionToGrid(mouseDownX);   // write a function for inPixels(inPinDistances)), call it snapToGrid
    int y1 = snapPixelPositionToGrid(mouseDownY);
    int x2 = snapPixelPositionToGrid(mouse.getX());
    int y2 = snapPixelPositionToGrid(mouse.getY());
    g.drawLine((float) x1, (float) y1, (float) x2, (float) y2);
  }

  // draw selected connections in red (preliminary):
  g.setColour(highlightColour);
  for(i = 0; i < selectedAudioConnections.size(); i++)
    g.drawLine(getLineForConnection(selectedAudioConnections[i]));
  for(i = 0; i < audioConnectionsInLasso.size(); i++)
    g.drawLine(getLineForConnection(audioConnectionsInLasso[i]));
}



void ModularBlockDiagramPanel::drawModule(Graphics &g, romos::Module *moduleToDraw)
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  if( moduleToDraw == NULL )
    return;


  //if( moduleToDraw->getTypeIdentifier() == romos::ModuleTypeRegistry::AUDIO_OUTPUT )
  //  DEBUG_BREAK;

  int x, y, w, h;
  getRectangleForModuleInPixels(moduleToDraw, x, y, w, h, false);
  juce::Rectangle<int> moduleReactangle(x, y, w ,h);


  if( moduleToDraw->isPolyphonic() )
  {
    // wrap this into a function (drawPolyphonyShadow() or something):
    g.setColour(plotColourScheme.getCurveColourUniform(0).withAlpha(0.375f));
    int offset = 8;
    g.drawRect(x-offset, y+offset, w, h, t);

    setPixel(g, x-2, y+2);
    setPixel(g, x-2, y+3);
    setPixel(g, x-3, y+2);
    setPixel(g, x-3, y+3);

    setPixel(g, x-5, y+5);
    setPixel(g, x-5, y+6);
    setPixel(g, x-6, y+5);
    setPixel(g, x-6, y+6);

    setPixel(g, x+w-2, y+h+2);
    setPixel(g, x+w-2, y+h+3);
    setPixel(g, x+w-3, y+h+2);
    setPixel(g, x+w-3, y+h+3);

    setPixel(g, x+w-5, y+h+5);
    setPixel(g, x+w-5, y+h+6);
    setPixel(g, x+w-6, y+h+5);
    setPixel(g, x+w-6, y+h+6);

    setPixel(g, x-2+1, y+h+2-2);
    setPixel(g, x-2+1, y+h+3-2);
    setPixel(g, x-3+1, y+h+2-2);
    setPixel(g, x-3+1, y+h+3-2);

    setPixel(g, x-5+1, y+h+5-2);
    setPixel(g, x-5+1, y+h+6-2);
    setPixel(g, x-6+1, y+h+5-2);
    setPixel(g, x-6+1, y+h+6-2);

    /*
    // setPixel not available anymore in juce 5.2 - find replacement
    g.setPixel(x-2, y+2);
    g.setPixel(x-2, y+3);
    g.setPixel(x-3, y+2);
    g.setPixel(x-3, y+3);

    g.setPixel(x-5, y+5);
    g.setPixel(x-5, y+6);
    g.setPixel(x-6, y+5);
    g.setPixel(x-6, y+6);

    g.setPixel(x+w-2, y+h+2);
    g.setPixel(x+w-2, y+h+3);
    g.setPixel(x+w-3, y+h+2);
    g.setPixel(x+w-3, y+h+3);

    g.setPixel(x+w-5, y+h+5);
    g.setPixel(x+w-5, y+h+6);
    g.setPixel(x+w-6, y+h+5);
    g.setPixel(x+w-6, y+h+6);

    g.setPixel(x-2+1, y+h+2-2);
    g.setPixel(x-2+1, y+h+3-2);
    g.setPixel(x-3+1, y+h+2-2);
    g.setPixel(x-3+1, y+h+3-2);

    g.setPixel(x-5+1, y+h+5-2);
    g.setPixel(x-5+1, y+h+6-2);
    g.setPixel(x-6+1, y+h+5-2);
    g.setPixel(x-6+1, y+h+6-2);
    */
  }

  //g.setColour(Colours::lavender); // preliminary

  g.setColour(getPlotColourScheme().topLeft);  

  g.fillRect(x, y, w, h);

  g.setColour(plotColourScheme.getCurveColourUniform(0));
  g.drawRect(x, y, w, h, t);

  juce::String name = juce::String( moduleToDraw->getName().c_str() );

  if( moduleToDraw->hasHeader() )
  {
    g.drawRect(x, y, w, bigFontHeight+2*m+2*t, t);  // encloses the title
    drawBitmapFontText(g, x+t+m, y+t+m, name, bigFont, plotColourScheme.text, -1, Justification::topLeft);
  }
  else
  {
    int xArrow;
    int xText = x+t+s+m;
    if( dynamic_cast<romos::AudioInputModule*> (moduleToDraw) )
    {
      xArrow  = x+t+m;
      xText  += arrowLength + m;
      g.drawArrow(Line<float>((float) xArrow, (float) (y+h/2), (float) (xArrow+arrowLength), (float) (y+h/2)), 
        (float) (t+1), (float) (h-(2*t+2*m)), (float) arrowHeadLength); 
      // preliminary - write our own drawing function (no-anti-aliasing)
    }
    else if( dynamic_cast<romos::AudioOutputModule*> (moduleToDraw) )
    {
      xArrow = x+w-arrowLength-t-m;
      g.drawArrow(Line<float>((float) xArrow, (float) (y+h/2), (float) (xArrow+arrowLength), (float) (y+h/2)), 
        (float) (t+1), (float) (h-(2*t+2*m)), (float) arrowHeadLength); 
    }
    //drawBitmapFontText(g, xText, y+t+m-2, name, &boldFont10px, plotColourScheme.text, -1, Justification::topLeft);
    drawBitmapFontText(g, xText, y+t+m-1, name, bigFont, plotColourScheme.text, -1, Justification::topLeft);
    //drawBitmapFontText(g, xText, y+t+m-2+1, name, normalFont, plotColourScheme.text, -1, Justification::topLeft);
  }

  if( isModuleSelected(moduleToDraw) || rosic::containsElement(modulesInLasso, moduleToDraw) )
    //if( isModuleSelected(moduleToDraw) || modulesInLasso.hasElement(moduleToDraw) )
  {
    g.setColour(highlightColour);
    g.drawRect(x, y, w, h, 1);  // \todo red frame is preliminary, find a better selection indicator
  }

  if( !dynamic_cast<romos::AudioInputModule*> (moduleToDraw) )
    drawInputPins(g, moduleToDraw, moduleReactangle);

  if( !dynamic_cast<romos::AudioOutputModule*> (moduleToDraw) )
    drawOutputPins(g, moduleToDraw, moduleReactangle);

  //drawIncomingConnectionsForModule(g, moduleToDraw, moduleReactangle);
}

void ModularBlockDiagramPanel::drawInputPins(Graphics &g, romos::Module *module, 
  juce::Rectangle<int> moduleRectangle) 
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
    return;

  // the position where to draw the pin's text:
  int x = moduleRectangle.getX() + s + t + m;
  int y = moduleRectangle.getY() + getModuleTitleHeightInPixels(module) + m;

  // coordinates of the pin itself:
  float px = (float) (x - m - 2*s - t); 
  float py = (float) (y + 1);
  float pw = (float) (t+2*s);
  float ph = (float) (smallFontHeight - 2);

  unsigned int i;
  for(i=0; i<module->getNumInputPins(); i++)
  {
    juce::String pinName = juce::String( module->getPinName(AUDIO, INCOMING, i).getRawString() );

    if(!module->isInputPinConnected(i) && module->hasHeader())
      pinName += "=" + juce::String(module->getInputPinDefaultValue(i));

    drawBitmapFontText(g, x, y, pinName, smallFont, plotColourScheme.text, -1, 
      Justification::topLeft);
    py = (float) (y+1);
    g.fillRect(px, py, pw, ph);
    y += (smallFontHeight+m);
    //y += pinDistance;
  }
}

void ModularBlockDiagramPanel::drawOutputPins(Graphics &g, romos::Module *module, juce::Rectangle<int> moduleRectangle)
{  
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
    return;

  int x = moduleRectangle.getX();
  int y = moduleRectangle.getY();
  int w = moduleRectangle.getWidth();
  int h = moduleRectangle.getHeight();

  y += getModuleTitleHeightInPixels(module) + m;

  float px = (float) (x + w - t - s);
  float py = (float) (y + 1);
  float pw = (float) (t+2*s);
  float ph = (float) (smallFontHeight - 2);

  unsigned int i;
  for(i = 0; i < module->getNumOutputPins(); i++)
  {
    drawBitmapFontText(g, x+w-s-t-m, y, juce::String( module->getPinName(AUDIO, OUTGOING, i).getRawString() ), 
      smallFont, plotColourScheme.text, -1,  Justification::topRight);
    py = (float) (y + 1);
    g.fillRect(px, py, pw, ph);
    y += (smallFontHeight+m);
  }

  // \todo perhaps we can avoid some code duplication by factoring out commonalities in drawing all the different kinds of pins
}

void ModularBlockDiagramPanel::drawIncomingConnectionsForModule(Graphics &g, romos::Module *module) 
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  if( module == NULL )
    return;

  if( module->isInputModule() )  // input modules do not show their connections (they are invisible to the GUI)
    return;

  juce::Rectangle<int> moduleRectangle = getRectangleForModuleInPixels(module, false);

  //g.setColour(Colours::red); // only for test

  int xs, xt, ys, yt;  // center coordinates of source and target pin
  unsigned int i;
  romos::AudioConnection c;
  romos::Module *sourceModule;
  unsigned int sourcePinIndex, targetPinIndex;

  std::vector<romos::AudioConnection> incomingConnections = module->getIncomingAudioConnections();

  for(i = 0; i < incomingConnections.size(); i++)
  {
    c = incomingConnections[i];
    sourceModule   = c.getSourceModule();

    if( sourceModule == NULL ) // can happen for temporary connections
      continue;

    sourcePinIndex = c.getSourceOutputIndex();
    targetPinIndex = c.getTargetInputIndex();
    getPinCenterCoordinates(romos::AUDIO, romos::INCOMING,c.getTargetInputIndex(), module, moduleRectangle, xt, yt);
    juce::Rectangle<int> sourceRectangle = getRectangleForModuleInPixels(c.getSourceModule(), false);
    getPinCenterCoordinates(romos::AUDIO, romos::OUTGOING, c.getSourceOutputIndex(), c.getSourceModule(), sourceRectangle, xs, ys);
    Line<float> connectionLine((float) xs, (float) ys, (float) xt, (float) yt);

    if( ys == yt )
      g.drawHorizontalLine(ys, (float) xs, (float) xt);
    else
    {
      g.drawLine(connectionLine);
      // \todo draw the line in 3 segments - 2 horizontals and a vertical - but the decision where to put the vertical should actually
      // depend on external factors such as how many other connections there are which might be obscured...mmmhhh
      // probably better to let the user do this by letting connections have breakpoints ...they could then also provide a method
      // isPointOnConnection...
    }

    if( c.hasImplicitDelay() )
    {
      int xMid = roundToInt( 0.5*(xs+xt) );
      int yMid = roundToInt( 0.5*(ys+yt) );
      int w = bigFont->getTextPixelWidth(juce::String(("D")), bigFont->getDefaultKerning())+2*m+2*t;
      int h = bigFont->getFontHeight()+2*m+2*t;
      int xd = xMid - w/2;
      int yd = yMid - h/2;
      g.setColour(Colours::white);  // preliminary
      g.fillRect(xd, yd, w, h);
      g.setColour(Colours::black);  // preliminary
      g.drawRect(xd, yd, w, h, t);
      drawBitmapFontText(g, xd+t+m, yd+t+m, juce::String(("D")), bigFont, plotColourScheme.text, 
        bigFont->getDefaultKerning(), Justification::topLeft);
    }
  }
}

void ModularBlockDiagramPanel::getRectangleForModuleInPixels(romos::Module *module, int &x, int &y, int &w, int &h, bool includingPins) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
  {
    x = y = w = h = 0;
    return;
  }

  x = inPixels(module->getPositionX());
  y = inPixels(module->getPositionY()) - getOffsetY(module);
  if( isModuleSelected(module) )
  {
    x += selectionOffsetX;
    y += selectionOffsetY;
  }
  w = getRequiredWidthForModuleInPixels(module, includingPins);
  h = getRequiredHeightForModuleInPixels(module);
  if( includingPins )
    x -= s;  // left border goes one stickout "s" leftward
}

juce::Rectangle<int> ModularBlockDiagramPanel::getRectangleForModuleInPixels(romos::Module *module, bool includingPins) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  int x, y, w, h;
  getRectangleForModuleInPixels(module, x, y, w, h, includingPins);
  return juce::Rectangle<int>(x, y, w, h);
}

int ModularBlockDiagramPanel::getModuleTitleHeightInPixels(romos::Module *module) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
    return 0;
  if( module->hasHeader() )
    return 2*t + 2*m + bigFontHeight;
  else
    return t;
}

bool ModularBlockDiagramPanel::isModuleInsideRectangle(romos::Module *module, juce::Rectangle<int> rectangle, bool includingPins) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  juce::Rectangle<int> moduleRectangle = getRectangleForModuleInPixels(module, includingPins);
  return rectangle.intersects(moduleRectangle);
}

int ModularBlockDiagramPanel::getRequiredWidthForModuleInPixels(romos::Module *module, bool includingPins) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
    return 0;

  juce::String name = juce::String( module->getName().c_str() );
  int width = bigFont->getTextPixelWidth(name, bigFont->getDefaultKerning());

  if( dynamic_cast<romos::AudioInputModule*> (module) || dynamic_cast<romos::AudioOutputModule*> (module) )
    width += arrowLength + m;

  if( includingPins == false )
    width += 2*t + 2*m;
  else
    width += 2*t + 2*m + 2*s;  // pins stick out one stickout "s" to left and right, so width is 2*s more

  if( !module->hasHeader() )
    width += 2*s;

  return width;
  // \todo maybe refine the width-computation to take into account the names of the pins
}

int ModularBlockDiagramPanel::getRequiredHeightForModuleInPixels(romos::Module *module) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
    return 0;
  int numPins = jmax(module->getNumInputPins(), module->getNumOutputPins() );
  return getModuleTitleHeightInPixels(module) + numPins*(m+smallFontHeight) + m + t;
}

juce::Rectangle<int> ModularBlockDiagramPanel::getPinBounds(int kindOfPin, int direction, int pinIndex, romos::Module *module, 
  juce::Rectangle<int> moduleRectangle) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
    return juce::Rectangle<int>(0, 0, 0, 0);

  int x = moduleRectangle.getX();
  if( direction == romos::OUTGOING )
    x += moduleRectangle.getWidth() - (t+2*m);

  int y = moduleRectangle.getY();
  y += getModuleTitleHeightInPixels(module) + m + 1; // + (int) floor(0.5 * hs);
  y += pinIndex * pinDistance;
  if( kindOfPin == romos::EVENT )
  {
    if( direction == romos::INCOMING )
      y += module->getNumInputPins() * pinDistance;
    else
      y += module->getNumOutputPins() * pinDistance;
  }

  int w = t+2*m;
  int h = smallFontHeight-2;

  return juce::Rectangle<int>(x, y, w, h);
}

void ModularBlockDiagramPanel::getPinCenterCoordinates(int kindOfPin, int direction, int pinIndex, romos::Module *module, 
  juce::Rectangle<int> moduleRectangle, int &xPin, int &yPin) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  if( module == NULL )
  {
    xPin = yPin;
    return;
  }

  xPin = moduleRectangle.getX();
  if( direction == romos::OUTGOING )
    xPin += moduleRectangle.getWidth()-1;
  else
    xPin += 1;

  yPin = moduleRectangle.getY();
  yPin += getModuleTitleHeightInPixels(module) + m + (int) floor(0.5 * smallFontHeight);
  yPin += pinIndex * pinDistance;
  if( kindOfPin == romos::EVENT )
  {
    if( direction == romos::INCOMING )
      yPin += module->getNumInputPins() * pinDistance;
    else
      yPin += module->getNumOutputPins() * pinDistance;
  }
}

juce::Line<float> ModularBlockDiagramPanel::getLineForConnection(romos::AudioConnection connection) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  romos::Module *sourceModule = connection.getSourceModule();
  romos::Module *targetModule = connection.getTargetModule();

  int sourcePinIndex = connection.getSourceOutputIndex(); 
  int targetPinIndex = connection.getTargetInputIndex();

  juce::Rectangle<int> sourceRectangle = getRectangleForModuleInPixels(sourceModule, false);
  juce::Rectangle<int> targetRectangle = getRectangleForModuleInPixels(targetModule, false);

  if( sourceModule == NULL ) 
    return Line<float>(0.f, 0.f, 200.f, 200.f); // if this line appears on the screen, something went wrong

  int xs, xt, ys, yt; 
  getPinCenterCoordinates(romos::AUDIO, romos::OUTGOING, sourcePinIndex, sourceModule, sourceRectangle, xs, ys);
  getPinCenterCoordinates(romos::AUDIO, romos::INCOMING, targetPinIndex, targetModule, targetRectangle, xt, yt);

  return Line<float>((float) xs, (float) ys, (float) xt, (float) yt);
}

std::vector<romos::Module*> ModularBlockDiagramPanel::getModulesInRectangle(juce::Rectangle<int> rectangle) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  std::vector<romos::Module*> result;

  unsigned int i;
  romos::Module *module;
  /*
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumInputs(); i++)
  {
  module = getInterfaceMediator()->getContainerShownInDiagram()->getAudioInputModule(i);
  if( isModuleInsideRectangle(module, rectangle, false) )
  rosic::appendElement(result, module);
  //result.appendElement(module);
  }
  */
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumChildModules(); i++)
  {
    module = getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i);
    if( isModuleInsideRectangle(module, rectangle, false) )
      rosic::appendElement(result, module);
    //result.appendElement(module);
  }
  /*
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumOutputs(); i++)
  {
  module = getInterfaceMediator()->getContainerShownInDiagram()->getAudioOutputModule(i);
  if( isModuleInsideRectangle(module, rectangle, false) )
  rosic::appendElement(result, module);
  //result.appendElement(module);
  }
  */
  // todo we need a means to iterate over all embedded modules at once to avoid this kind of code-duplication
  // ->done ->commented code may be deleted ...hopefully - but check that it works before

  return result;
}

romos::Module* ModularBlockDiagramPanel::getModuleAtPixels(int x, int y, bool considerPins) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  // loop through the modules, return the first which is found to include the point:
  // todo: maybe use a local variable for getInterfaceMediator()->getContainerShownInDiagram() as abbreviation
  unsigned int i;
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumInputPins(); i++)
  {
    if( getRectangleForModuleInPixels(getInterfaceMediator()->getContainerShownInDiagram()->getAudioInputModule(i), considerPins).contains(x, y) )
      return getInterfaceMediator()->getContainerShownInDiagram()->getAudioInputModule(i);
  }
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumChildModules(); i++)
  {
    if( getRectangleForModuleInPixels(getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i), considerPins).contains(x, y) )
      return getInterfaceMediator()->getContainerShownInDiagram()->getChildModule(i);
  }
  for(i = 0; i < getInterfaceMediator()->getContainerShownInDiagram()->getNumOutputPins(); i++)
  {
    if( getRectangleForModuleInPixels(getInterfaceMediator()->getContainerShownInDiagram()->getAudioOutputModule(i), considerPins).contains(x, y) )
      return getInterfaceMediator()->getContainerShownInDiagram()->getAudioOutputModule(i);
  }

  return NULL;
}

romos::AudioConnection ModularBlockDiagramPanel::getConnectionAtPixels(int x, int y) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));
  romos::AudioConnection ac;
  juce::Line<float> connectionLine;
  float tolerance = 2.f;
  for(unsigned int i = 0; i < allAudioConnections.size(); i++)
  {
    ac = allAudioConnections.at(i);
    connectionLine = getLineForConnection(ac);

    Point<float> dummy;
    if( connectionLine.getDistanceFromPoint(Point<float>((float) x, (float) y), dummy) < tolerance )
      return ac;
    //if( connectionLine.getDistanceFromLine((float) x, (float) y) < tolerance )
    //  return ac;

  }
  return AudioConnection();
}

bool ModularBlockDiagramPanel::getPinPropertiesAtPixels(int x, int y, romos::Module* module, juce::Rectangle<int> moduleRectangle, 
  int &kindOfPin, int &directionOfPin, int &indexOfPin) const
{
  ScopedLock scopedLock(*(getInterfaceMediator()->plugInLock));

  kindOfPin      = -1;
  directionOfPin = -1;
  indexOfPin     = -1;

  if( module == NULL )
    return false;

  unsigned int i;
  int tolerance  = pinDistance/2;
  int yPin;
  int yFirstPin  = moduleRectangle.getY() + getModuleTitleHeightInPixels(module) + m + (int) floor(0.5 * smallFontHeight);

  if( x >= moduleRectangle.getX() && x <= moduleRectangle.getX() + tolerance ) // x is near left border
  {
    if( dynamic_cast<romos::AudioInputModule*> (module) != NULL )
      return false;

    yPin = yFirstPin;
    for(i = 0; i < module->getNumInputPins(); i++)
    {
      if( abs(y-yPin) <= tolerance )
      {
        kindOfPin      = romos::AUDIO;
        directionOfPin = romos::INCOMING;
        indexOfPin     = i;
        return true;
      }
      yPin += pinDistance;
    }
  }

  if( x <= moduleRectangle.getRight() && x >= moduleRectangle.getRight() - tolerance ) // x is near right border
  {
    if( dynamic_cast<romos::AudioOutputModule*> (module) != NULL )
      return false;

    yPin = yFirstPin;
    for(i=0; i<module->getNumOutputPins(); i++)
    {
      if( abs(y-yPin) <= tolerance )
      {
        kindOfPin      = romos::AUDIO;
        directionOfPin = romos::OUTGOING;
        indexOfPin     = i;
        return true;
      }
      yPin += pinDistance;
    }
  }

  return false;
}


//=========================================================================================================================================
// class LibertyEditor:

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

LibertyEditor::LibertyEditor(CriticalSection *newPlugInLock, LibertyAudioModule* newLibertyAudioModule) 
//: PolyphonicInstrumentEditor(newPlugInLock, newLibertyAudioModule)
  : AudioModuleEditor(newLibertyAudioModule)
{
  ScopedLock scopedLock(*lock); 
  setHeadlineStyle(MAIN_HEADLINE);

  jassert(newLibertyAudioModule != NULL ); // you must pass a valid module here

  modularSynthAudioModule = newLibertyAudioModule;

  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);


  interfaceMediator = new LibertyInterfaceMediator(newPlugInLock, newLibertyAudioModule);

  structureTreeView = new ModularStructureTreeView(interfaceMediator);
  structureTreeView->setOpenOrCloseNodesOnClick(true);
  structureTreeView->setDescription(("Patch structure represented as tree"));
  structureTreeView->setDescriptionField(descriptionField);
  addWidget(structureTreeView);

  moduleEditorHolder = new ModulePropertiesEditorHolder(interfaceMediator);
  moduleEditorHolder->setDescription(("Editor for currently selected module"));
  moduleEditorHolder->setDescriptionField(descriptionField, true);
  addChildColourSchemeComponent(moduleEditorHolder);

  blockDiagramPanel = new ModularBlockDiagramPanel(interfaceMediator);
  blockDiagramPanel->setDescription(("Block diagram representation of selected container"));
  blockDiagramPanel->setDescriptionField(descriptionField, true);
  diagramScrollContainer = new ComponentScrollContainer(blockDiagramPanel);
  addChildColourSchemeComponent(diagramScrollContainer);

  presetSectionPosition = BELOW_HEADLINE;
  //isTopLevelEditor      = true;
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);

  updateWidgetsAccordingToState();

  setSize(600, 400);
}

LibertyEditor::~LibertyEditor()
{
  // for the mediated components, we must take care to delete them  before the mediator is deleted, so we do it manually here:
  removeChildComponent(structureTreeView);
  removeChildComponent(moduleEditorHolder);
  removeChildComponent(diagramScrollContainer);
  delete structureTreeView;
  delete moduleEditorHolder;
  delete diagramScrollContainer;
  delete interfaceMediator;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

/*
void LibertyEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
ScopedLock scopedLock(*plugInLock);

updateWidgetsAccordingToState();
}
*/

void LibertyEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x, y, w, h;

  x = 0;
  y = stateWidgetSet->getBottom() - 8;
  w = getWidth()/4 ;
  h = 240;
  stateWidgetSet->setBounds(0, 4, getWidth()/4, 32);

  structureTreeView->setBounds(x, y, w, h);
  x += w - RWidget::outlineThickness; 
  w  = getWidth() - x ;
  moduleEditorHolder->setBounds(x, y, w, h);


  x = 0;
  y = moduleEditorHolder->getBottom();
  w = getWidth();
  h = infoField->getY() - y;
  diagramScrollContainer->setBounds(x, y, w, h);
  blockDiagramPanel->setAvailabeSizeForCanvas(w, h);
}

void LibertyEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*lock);
  structureTreeView->getInterfaceMediator()->setContainerToShowInDiagram(structureTreeView->getInterfaceMediator()->getTopLevelModule());
  // the mediator will take care to update all panels
}


/*

Bugs:
-The Poly switch on the GUI seems to have no effect


*/