rsModulationSetup::rsModulationSetup(AutomatableWidget* widgetToModulate, 
  MetaParameterManager* metaManagerTouse)
  : widget(widgetToModulate), metaManager(metaManagerTouse)
  , rsDeletionRequester(widgetToModulate)
{
  addWidget( modulationsLabel = new RTextField("Modulations") );
  modulationsLabel->setNoBackgroundAndOutline(true);
  modulationsLabel->setDescription(juce::String("Modulation setup"));

  addWidget( closeButton = new RClickButtonNotifyOnMouseUp(RButton::CLOSE) );
  closeButton->setDescription(juce::String("Closes the modulation setup window"));
  closeButton->setClickingTogglesState(false);
  closeButton->addRButtonListener(this);

  addWidget( addButton = new RButton("Add") );
  addButton->setDescription(juce::String("Adds a new modulation connection"));
  addButton->setClickingTogglesState(false);
  addButton->addRButtonListener(this);

  addWidget( removeButton = new RButton("Remove") );
  removeButton->setDescription(juce::String("Removes an existing modulation connection"));
  removeButton->setClickingTogglesState(false);
  removeButton->addRButtonListener(this);

  addWidget( clipMinField = new RLabeledTextEntryField("Min:") );
  clipMinField->setEntryFieldText(String(getClipMin()));
  clipMinField->setDescription(juce::String("Range minimum for modulated value"));
  clipMinField->setLabelWidth(36);
  clipMinField->getTextEntryField()->registerTextEntryFieldObserver(this);

  addWidget( clipMaxField = new RLabeledTextEntryField("Max:") );
  clipMaxField->setEntryFieldText(String(getClipMax()));
  clipMaxField->setDescription(juce::String("Range maximum for modulated value"));
  clipMaxField->setLabelWidth(36);
  clipMaxField->getTextEntryField()->registerTextEntryFieldObserver(this);

  updateConnectionWidgetsArray();

  //modManager = nullptr;
  //ModulatableParameter* mp = widgetToModulate->getModulatableParameter();
  //if(mp != nullptr)
  //  modManager = mp->getModulationManager();
}

rsModulationSetup::~rsModulationSetup()
{
  delete connectableSourcesPopUp;
  delete removableSourcesPopUp;
}

//void rsModulationSetup::paint(Graphics& g)
//{
//  ColourSchemeComponent::paint(g);
//}

void rsModulationSetup::resized()
{
  int d   = sliderDistance;
  int x   = d;
  int y   = d;
  int w   = getWidth();
  int h   = getHeight();
  int sh  = sliderHeight;
  int inc = sh+d;

  closeButton->setBounds(w-16, 0, 16, 16);
  modulationsLabel->setBounds(x, y, w-8-16, sh); y += inc; 

  for(int i = 0; i < size(connectionWidgets); i++) {
    connectionWidgets[i]->setBounds(x, y, w-8, sh); y += inc; }

  y = h - sh - d;
  addButton->setBounds(x, y, 40, 16);
  x = addButton->getRight() + d;
  removeButton->setBounds(x, y, 60, 16);

  y -= inc;
  x  = 0;
  clipMinField->setBounds(x+d, y, w/2-2*d, 16);
  x = w/2;
  clipMaxField->setBounds(x+d, y, w/2-2*d, 16);
}

void rsModulationSetup::rButtonClicked(RButton *button)
{
  if(button == closeButton)
    requestDeletion();
  else if(button == addButton)
    showConnectableSourcesPopUp();
  else if(button == removeButton)
    showRemovableSourcesPopUp();
  for(int i = 0; i < size(connectionWidgets); i++){
    if(button == connectionWidgets[i]->removeButton)
      removeConnection(i);
  }
}

void rsModulationSetup::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  int id;
  if(menuThatHasChanged == connectableSourcesPopUp)
  {
    id = connectableSourcesPopUp->getSelectedIdentifier();
    if(id > 0)
      addConnection(id-1);
  }
  else if(menuThatHasChanged == removableSourcesPopUp)
  {
    id = removableSourcesPopUp->getSelectedIdentifier();
    if(id > 0)
      removeConnection(id-1);
  }
}

void rsModulationSetup::textChanged(RTextEntryField *rTextEntryFieldThatHasChanged)
{
  //jassertfalse;
  // this doesn't work yet: String::getDoubleValue doesn't parse -inf correctly - it returns inf

  if(rTextEntryFieldThatHasChanged == clipMinField->getTextEntryField())
    setClipMin(toDouble(clipMinField->getTextEntryField()->getText()));
  else if(rTextEntryFieldThatHasChanged == clipMaxField->getTextEntryField())
    setClipMax(toDouble(clipMaxField->getTextEntryField()->getText()));
}

void rsModulationSetup::addConnection(int index)
{
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationSource*> sources = mp->getDisconnectedSources();
    mp->addModulationSource(sources[index]);
    addWidgetsForConnection(mp->getConnectionTo(sources[index]));
  }
}

void rsModulationSetup::removeConnection(int index)
{
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    removeWidgetsForConnection(index);
    std::vector<ModulationSource*> sources = mp->getConnectedSources();
    mp->removeModulationSource(sources[index]);
    updateSize();
  }
}

void rsModulationSetup::updateConnectionWidgetsArray()
{
  clearConnectionWidgets();
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationConnection*> connections = mp->getConnections();
    for(int i = 0; i < size(connections); i++)
    {
      if(!hasSlider(connections[i]->getDepthParameter()))
        addWidgetsForConnection(connections[i]);
    }
  }
  updateSize();
}

void rsModulationSetup::showConnectableSourcesPopUp()
{
  // create popup, if necessary:
  if(connectableSourcesPopUp == nullptr)
  {
    connectableSourcesPopUp = new RPopUpMenu(this); // maybe attach to the addButton instead of this?
    connectableSourcesPopUp->registerPopUpMenuObserver(this);
    connectableSourcesPopUp->setDismissOnFocusLoss(true);
  }

  // populate it:
  connectableSourcesPopUp->clear();
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationSource*> sources = mp->getDisconnectedSources();
    for(int i = 0; i < size(sources); i++)
    {
      juce::String name = sources[i]->getModulationSourceDisplayName();
      connectableSourcesPopUp->addItem(i+1, name);     // +1 bcs 0 is not allowed for the id
    }
  }

  // show it:
  int w = connectableSourcesPopUp->getRequiredWidth(true);
  int h = connectableSourcesPopUp->getRequiredHeight(true);
  connectableSourcesPopUp->show(true, RPopUpComponent::BELOW, w, h); // showModally = true
}

void rsModulationSetup::showRemovableSourcesPopUp()
{
  // lots of code duplication from showConnectableSourcesPopUp - can this be refactored?

  if(removableSourcesPopUp == nullptr)
  {
    removableSourcesPopUp = new RPopUpMenu(this);
    removableSourcesPopUp->registerPopUpMenuObserver(this);
    removableSourcesPopUp->setDismissOnFocusLoss(true);
  }

  removableSourcesPopUp->clear();
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationSource*> sources = mp->getConnectedSources();
    for(int i = 0; i < size(sources); i++)
    {
      juce::String name = sources[i]->getModulationSourceDisplayName();
      removableSourcesPopUp->addItem(i+1, name);
    }
  }

  int w = removableSourcesPopUp->getRequiredWidth(true);
  int h = removableSourcesPopUp->getRequiredHeight(true);
  removableSourcesPopUp->show(true, RPopUpComponent::BELOW, w, h);
}

bool rsModulationSetup::hasSlider(MetaControlledParameter* p)
{
  for(int i = 0; i < size(connectionWidgets); i++)
  {
    Parameter* ps = connectionWidgets[i]->depthSlider->getAssignedParameter();
    if(ps == p)
      return true;
  }
  return false;
}

void rsModulationSetup::addWidgetsForConnection(ModulationConnection* c)
{
  rsModulationConnectionWidget* w = new rsModulationConnectionWidget(c, this);
  connectionWidgets.push_back(w);
  w->depthSlider->assignParameter(c->getDepthParameter());
  w->depthSlider->setSliderName(c->getSource()->getModulationSourceDisplayName());
  w->removeButton->addRButtonListener(this);
  addWidget(w);
  updateSize();
}

void rsModulationSetup::removeWidgetsForConnection(int i)
{
  connectionWidgets[i]->removeButton->removeRButtonListener(this);
  deleteObject(connectionWidgets[i]);              // mark for later deletion
  removeWidget(connectionWidgets[i], true, false); // false, to not delete it immediately
  remove(connectionWidgets, i);
}

void rsModulationSetup::clearConnectionWidgets()
{
  for(int i = 0; i < size(connectionWidgets); i++)
  {  
    // maybe factor out these 3 lines into a function, they appear also in 
    // removeWidgetsForConnection:
    connectionWidgets[i]->removeButton->removeRButtonListener(this);
    deleteObject(connectionWidgets[i]);              // mark for later deletion
    removeWidget(connectionWidgets[i], true, false); // false, to not delete it immediately
  }
  connectionWidgets.clear();
}

void rsModulationSetup::updateSize()
{
  int width  = 250;  // maybe we should use the widget's width...but maybe not
  int height = 100;  // preliminary

  height  = (sliderHeight+sliderDistance) * size(connectionWidgets);
  height += 68;

  setSize(width, height); 
  resized(); // needed during development - might be redundant when finished
}

void rsModulationSetup::setClipMin(double newMin) 
{ 
  widget->getModulatableParameter()->setModulationRangeMin(newMin); 
}

void rsModulationSetup::setClipMax(double newMax) 
{ 
  widget->getModulatableParameter()->setModulationRangeMax(newMax); 
}

double rsModulationSetup::getClipMin()
{ 
  return widget->getModulatableParameter()->getModulationRangeMin(); 
}

double rsModulationSetup::getClipMax()
{ 
  return widget->getModulatableParameter()->getModulationRangeMax(); 
}

//=================================================================================================

AutomatableWidget::AutomatableWidget(RWidget *widgetToWrap)
{
  wrappedWidget = widgetToWrap;
}

AutomatableWidget::~AutomatableWidget()
{
  delete rightClickPopUp;
  delete modSetup;
}

bool AutomatableWidget::isPopUpOpen()
{
  return popUpIsOpen;

  // The code below doesn't work because apparently, when the user clicks on a widget while the
  // popup is open, the popup gets closed first and only after that, the mouseDown callback of the
  // widget is received, so isPopUpOpen would always return false in the mouseDown callback. The
  // desired behavior is that one right-click on the widget opens the popup and a second click
  // closes it. This behavior now requires that we maintain a popUpIsOpen flag here and that flag
  // should be set in the mouseDown method of the widget. Without that, a second right-click would
  // make the menu disappear for a fraction of a second and immediately reappear.

  //if(rightClickPopUp == nullptr)
  //  return false;
  //else
  //  return rightClickPopUp->isOpen();
}

void AutomatableWidget::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  if(menuThatHasChanged != rightClickPopUp)
    return;

  int selectedIdentifier = rightClickPopUp->getSelectedIdentifier();

  if(selectedIdentifier == MODULATION_SETUP)
  {
    showModulationSetup();
    return;
  }

  MetaControlledParameter* mcp = getMetaControlledParameter();
  if(mcp != nullptr)
  {
    switch(selectedIdentifier)
    {
    case META_ATTACH: mcp->attachToMetaParameter(
      (int)wrappedWidget->openModalNumberEntryField(mcp->getMetaParameterIndex()) ); break;
    case META_DETACH: mcp->detachFromMetaParameter();    break;
    }
  }

  AutomatableParameter* ap = getAutomatableParameter();
  if(ap != nullptr)
  {
    if(selectedIdentifier != MIDI_LEARN)
      ap->switchIntoMidiLearnMode(false); // turn off, if currently on and something else is selected
    switch(selectedIdentifier)
    {
    case MIDI_LEARN:  ap->switchIntoMidiLearnMode();             break;
    case MIDI_ASSIGN:
    {
      int result = (int)wrappedWidget->openModalNumberEntryField(ap->getAssignedMidiController());
      result = (int)clip(result, 0, 127);
      ap->assignMidiController(result);
    } break;
    case MIDI_MIN:    ap->setLowerAutomationLimit(ap->getValue());   break;
    case MIDI_MAX:    ap->setUpperAutomationLimit(ap->getValue());   break;
    case MIDI_REVERT: ap->revertToDefaults(false, false, false);     break;
    }
  }
}

void AutomatableWidget::rPopUpMenuDismissed(RPopUpMenu* menuThatwasDismissed)
{
  popUpIsOpen = false;
}

void AutomatableWidget::updatePopUpMenu()
{
  if(rightClickPopUp == nullptr)
  {
    // popup used the 1st time - we need to create it:
    rightClickPopUp = new RPopUpMenu(wrappedWidget);
    rightClickPopUp->registerPopUpMenuObserver(this);
    rightClickPopUp->setDismissOnFocusLoss(true);
    wrappedWidget->addChildWidget(rightClickPopUp, false, false);
  }
  rightClickPopUp->clear();
  addPopUpMenuItems();
}

void AutomatableWidget::addPopUpMenuItems()
{
  addPopUpMetaItems();
  addPopUpMidiItems();
  addPopUpModulationItems();
}

void AutomatableWidget::addPopUpMidiItems()
{
  AutomatableParameter* ap = getAutomatableParameter();
  if(ap != nullptr)
  {
    // prepare some strings for the popup menu:
    int cc = ap->getAssignedMidiController();
    String ccString;
    if(cc > -1)
      ccString = "(currently CC" + String(cc) + ")";
    else
      ccString = "(currently none)";

    int defaultCc = ap->getDefaultMidiController();
    String defaultString;
    if(defaultCc > -1)
      defaultString = "CC" + String(defaultCc);
    else
      defaultString = "none";
    String minString = wrappedWidget->stringConversionFunction(ap->getLowerAutomationLimit());
    String maxString = wrappedWidget->stringConversionFunction(ap->getUpperAutomationLimit());

    rightClickPopUp->addItem(MIDI_LEARN,  "MIDI learn " + ccString);
    rightClickPopUp->addItem(MIDI_ASSIGN, "MIDI assign");
    rightClickPopUp->addItem(MIDI_MIN,    "use value as lower limit (currently " + minString + String(")"));
    rightClickPopUp->addItem(MIDI_MAX,    "use value as upper limit (currently " + maxString + String(")"));
    rightClickPopUp->addItem(MIDI_REVERT, "revert MIDI mapping to defaults");
  }
}

void AutomatableWidget::addPopUpMetaItems()
{
  MetaControlledParameter* mcp = getMetaControlledParameter();
  if(mcp != nullptr)
  {
    int mi = mcp->getMetaParameterIndex();
    String miString;
    if(mi > -1)
      miString = "(currently " + String(mi) + ": " + mcp->getMetaParameterName() + ")";
    else
      miString = "(currently none)";

    rightClickPopUp->addItem(META_ATTACH, "Meta attach " + miString);
    rightClickPopUp->addItem(META_DETACH, "Meta detach");
  }
}

void AutomatableWidget::addPopUpModulationItems()
{
  ModulatableParameter* mp = getModulatableParameter();
  if(mp != nullptr)
    rightClickPopUp->addItem(MODULATION_SETUP, "Modulation setup");
}

void AutomatableWidget::openRightClickPopupMenu()
{
  updatePopUpMenu();
  int w = jmax(wrappedWidget->getWidth(), rightClickPopUp->getRequiredWidth(true));
  int h = jmin(200,                       rightClickPopUp->getRequiredHeight(true));
  //rightClickPopUp->show(false, RPopUpComponent::BELOW, w, h); // showModally = false
  rightClickPopUp->show(true, RPopUpComponent::BELOW, w, h); // showModally = true
  // If we don't show it modally (1st parameter = true), it will be immediately dismissed
  // after opening (so it appears as if it doesn't open at all). We could avoid it by calling
  // setDismissOnFocusLoss(false) in our constructor, but then it will stay open all the time
  // until we choose some option.

  popUpIsOpen = true;
}

void AutomatableWidget::closePopUp()
{
  if(rightClickPopUp != nullptr)
    rightClickPopUp->dismiss();
  popUpIsOpen = false;
}

void AutomatableWidget::showModulationSetup()
{
  //int ww = wrappedWidget->getWidth();        // widget width
  int wh = wrappedWidget->getHeight();       // widget height
  int x  = wrappedWidget->getScreenX();
  int y  = wrappedWidget->getScreenY() + wh; // preliminary

  if(modSetup == nullptr)
    modSetup = new rsModulationSetup(this, getMetaParameterManager());
  modSetup->setTopLeftPosition(x, y);
  modSetup->addToDesktop(ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);
  modSetup->setVisible(true);
  modSetup->toFront(true);
}

void AutomatableWidget::deleteObject(rsDeletionRequester* objectToDelete)
{
  if(objectToDelete == modSetup)
    modSetup->setVisible(false);
  else
    jassertfalse;
}

AutomatableParameter* AutomatableWidget::getAutomatableParameter()
{
  return dynamic_cast<AutomatableParameter*> (wrappedWidget->assignedParameter);
}

MetaControlledParameter* AutomatableWidget::getMetaControlledParameter()
{
  return dynamic_cast<MetaControlledParameter*> (wrappedWidget->assignedParameter);
}

ModulatableParameter* AutomatableWidget::getModulatableParameter()
{
  return dynamic_cast<ModulatableParameter*> (wrappedWidget->assignedParameter);
}

MetaParameterManager* AutomatableWidget::getMetaParameterManager()
{
  MetaControlledParameter* mcp = 
    dynamic_cast<MetaControlledParameter*> (wrappedWidget->assignedParameter);
  if(mcp)
    return mcp->getMetaParameterManager();
  else
    return nullptr;
}

//=================================================================================================

AutomatableSlider::AutomatableSlider()
  : AutomatableWidget(this)
{

}

void AutomatableSlider::mouseDown(const MouseEvent& e)
{
  if(e.mods.isRightButtonDown())
  {
    openRightClickPopupMenu();
    //if(!isPopUpOpen())
    //  openRightClickPopupMenu();
    //else
    //  closePopUp();
  }
  else
  {
    RSlider::mouseDown(e);

    // doesn't work:
    //if(!isPopUpOpen())
    //  RSlider::mouseDown(e);
    //else
    //  closePopUp();
  }
}

void AutomatableSlider::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  if( menuThatHasChanged != rightClickPopUp )
    return;

  RTreeViewNode *selectedItem = rightClickPopUp->getSelectedItem();
  if( selectedItem == NULL )
    return;
  int selectedIdentifier = selectedItem->getNodeIdentifier();

  switch( selectedIdentifier )
  {
  case ENTER_VALUE:   setValue(openModalNumberEntryField(getValue()),        true, false); break;
  case DEFAULT_VALUE: setValue(selectedItem->getNodeText().getDoubleValue(), true, false); break; //?
  default: AutomatableWidget::rPopUpMenuChanged(menuThatHasChanged);
  }
}

void AutomatableSlider::addPopUpMenuItems()
{
  addPopUpEnterValueItem();
  addPopUpDefaultValueItems();
  AutomatableWidget::addPopUpMenuItems();
}

void AutomatableSlider::addPopUpEnterValueItem()
{
  rightClickPopUp->addItem(ENTER_VALUE, "Enter Value");
}

void AutomatableSlider::addPopUpDefaultValueItems()
{
  if( defaultValues.size() > 0 )
  {
    RTreeViewNode *defaultValuesNode = new RTreeViewNode("Default Values");
    for(int i = 0; i < size(defaultValues); i++)
      defaultValuesNode->addChildNode(new RTreeViewNode(String(defaultValues[i]), DEFAULT_VALUE));
    rightClickPopUp->addTreeNodeItem(defaultValuesNode);
  }
}

//=================================================================================================

AutomatableComboBox::AutomatableComboBox()
  : AutomatableWidget(this)
{

}

void AutomatableComboBox::mouseDown(const MouseEvent& e)
{
  if( e.mods.isRightButtonDown() )
  {
    if(!isPopUpOpen())
      openRightClickPopupMenu();
    else
      AutomatableWidget::closePopUp();
  }
  else
    RComboBox::mouseDown(e);
}

void AutomatableComboBox::parameterChanged(Parameter* p)
{
  RWidget::parameterChanged(p);
  // not sure, why that's needed - isn't it supposed to be called anyway, i.e. if we don't override
  // parameterChanged, the RWidget baseclass method would be called? but somehow, it doesn't seem
  // to work
}

//=================================================================================================

AutomatableButton::AutomatableButton(const juce::String& buttonText)
  : RButton(buttonText), AutomatableWidget(this)
{

}

void AutomatableButton::mouseDown(const MouseEvent& e)
{
  if( e.mods.isRightButtonDown() )
  {
    if(!isPopUpOpen())
      openRightClickPopupMenu();
    else
      closePopUp();
  }
  else
    RButton::mouseDown(e);
}

void AutomatableButton::parameterChanged(Parameter* p)
{
  RWidget::parameterChanged(p);
}

//=================================================================================================

void rsModulationDepthSlider::rPopUpMenuChanged(RPopUpMenu* menu)
{
  int id = rightClickPopUp->getSelectedIdentifier();
  switch( id )
  {
  case MOD_DEPTH_MIN: setModDepthMin(openModalNumberEntryField(getModDepthMin())); break;
  case MOD_DEPTH_MAX: setModDepthMax(openModalNumberEntryField(getModDepthMax())); break;
  default: AutomatableSlider::rPopUpMenuChanged(menu);
  }
  juce::String text = menu->getSelectedText();
  typedef ModulationConnection::modModes MM;
  if(text == "Mode: Absolute")       setModMode(MM::ABSOLUTE);
  if(text == "Mode: Relative")       setModMode(MM::RELATIVE);
  if(text == "Mode: Exponential")    setModMode(MM::EXPONENTIAL);
  if(text == "Mode: Multiplicative") setModMode(MM::MULTIPLICATIVE);
}

void rsModulationDepthSlider::addPopUpMenuItems()
{
  AutomatableSlider::addPopUpMenuItems();
  addPopUpMinMaxAndModeItems();
}

void rsModulationDepthSlider::addPopUpMinMaxAndModeItems()
{
  rightClickPopUp->addItem(MOD_DEPTH_MIN, "Mod depth min");
  rightClickPopUp->addItem(MOD_DEPTH_MAX, "Mod depth max");
  typedef ModulationConnection::modModes MM;
  int m = getModMode();
  rightClickPopUp->addItem(MOD_MODE_ABSOLUTE,       "Mode: Absolute",       true, m == MM::ABSOLUTE);
  rightClickPopUp->addItem(MOD_MODE_RELATIVE,       "Mode: Relative",       true, m == MM::RELATIVE);
  rightClickPopUp->addItem(MOD_MODE_EXPONENTIAL,    "Mode: Exponential",    true, m == MM::EXPONENTIAL);
  rightClickPopUp->addItem(MOD_MODE_MULTIPLICATIVE, "Mode: Multiplicative", true, m == MM::MULTIPLICATIVE);
}

//=================================================================================================

rsModulationConnectionWidget::rsModulationConnectionWidget(ModulationConnection* connection,
  rsGarbageCollector* deletor)
  : rsDeletionRequester(deletor)
{
  addChildWidget(depthSlider  = new rsModulationDepthSlider(connection));
  addChildWidget(removeButton = new RClickButtonNotifyOnMouseUp(RButton::CLOSE));
}

void rsModulationConnectionWidget::resized()
{
  int w = getWidth();
  int h = getHeight();
  int buttonWidth  = h;
  int buttonMargin = 2;
  depthSlider->setBounds(0, 0, w-buttonWidth-buttonMargin, h);
  removeButton->setBounds(depthSlider->getRight()+buttonMargin, 0, buttonWidth, h);
}

//=================================================================================================

ModulatableSlider::~ModulatableSlider()
{
  ModulatableParameter* mp = dynamic_cast<ModulatableParameter*> (assignedParameter);
  if(mp)
    mp->deRegisterModulationTargetObserver(this);
}

void ModulatableSlider::modulationsChanged()
{
  repaint();
}

void ModulatableSlider::assignParameter(Parameter* p)
{
  AutomatableSlider::assignParameter(p);
  ModulatableParameter* mp = dynamic_cast<ModulatableParameter*> (p);
  if(mp)
    mp->registerModulationTargetObserver(this);
}

void ModulatableSlider::paint(Graphics& g)
{
  AutomatableSlider::paint(g);
  ModulatableParameter* mp = dynamic_cast<ModulatableParameter*> (assignedParameter);
  if(mp && mp->hasModulation())
    g.fillAll(Colour::fromFloatRGBA(1.f, 0.f, 0.f, 0.125f)); // preliminary
}