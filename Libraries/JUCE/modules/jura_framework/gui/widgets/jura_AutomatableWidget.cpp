AutomatableWidget::AutomatableWidget(RWidget *widgetToWrap)
{
  wrappedWidget = widgetToWrap;
}

AutomatableWidget::~AutomatableWidget()
{
  delete rightClickPopUp;
}

void AutomatableWidget::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  if(menuThatHasChanged != rightClickPopUp)
    return;
  RTreeViewNode *selectedItem = rightClickPopUp->getSelectedItem();
  if(selectedItem == nullptr)
    return;
  int selectedIdentifier = selectedItem->getNodeIdentifier();

  MetaControlledParameter* mcp = getMetaControlledParameter();
  if(mcp != nullptr)
  {
    switch(selectedIdentifier)
    {
    case META_ATTACH: mcp->attachToMetaParameter( 
      (int)wrappedWidget->openModalNumberEntryField() ); break;
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
      int result = (int)wrappedWidget->openModalNumberEntryField();
      result = (int)clip(result, 0, 127);
      ap->assignMidiController(result);
    } break;
    case MIDI_MIN:    ap->setLowerAutomationLimit(ap->getValue());   break;
    case MIDI_MAX:    ap->setUpperAutomationLimit(ap->getValue());   break;
    case MIDI_REVERT: ap->revertToDefaults(false, false, false);     break;
    }
  }
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
}

AutomatableParameter* AutomatableWidget::getAutomatableParameter()
{
  return dynamic_cast<AutomatableParameter*> (wrappedWidget->assignedParameter);
}

MetaControlledParameter* AutomatableWidget::getMetaControlledParameter()
{
  return dynamic_cast<MetaControlledParameter*> (wrappedWidget->assignedParameter);
}

//=================================================================================================

AutomatableSlider::AutomatableSlider() 
  : AutomatableWidget(this)
{

}

void AutomatableSlider::mouseDown(const MouseEvent& e)
{
  if( e.mods.isRightButtonDown() )
    openRightClickPopupMenu();
  else
    RSlider::mouseDown(e);
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
  case ENTER_VALUE:   setValue(openModalNumberEntryField(),                  true, false); break;
  case DEFAULT_VALUE: setValue(selectedItem->getNodeText().getDoubleValue(), true, false); break;
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
    for(int i = 0; i < defaultValues.size(); i++)
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
    openRightClickPopupMenu();
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
    openRightClickPopupMenu();
  else
    RButton::mouseDown(e);
}

void AutomatableButton::parameterChanged(Parameter* p)
{
  RWidget::parameterChanged(p); 
}
