AutomatableWidget::AutomatableWidget()
{

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
  if(selectedItem == NULL)
    return;
  int selectedIdentifier = selectedItem->getNodeIdentifier();

  AutomatableParameter *ap = dynamic_cast<AutomatableParameter*> (assignedParameter);
  if( ap == NULL )
    return;
  if(selectedIdentifier != MIDI_LEARN)
    ap->switchIntoMidiLearnMode(false); // turn off, if currently on and somehting else is selected
  switch( selectedIdentifier )
  {
  case MIDI_LEARN:  ap->switchIntoMidiLearnMode();             break;
  case MIDI_ASSIGN:
  {
    int result = (int) openModalNumberEntryField();
    result = (int) clip(result, 0, 127);
    ap->assignMidiController(result);
  } break;
  case MIDI_MIN:    ap->setLowerAutomationLimit(ap->getValue());   break;
  case MIDI_MAX:    ap->setUpperAutomationLimit(ap->getValue());   break;
  case MIDI_REVERT: ap->revertToDefaults(false, false, false);     break;
  }
}

void AutomatableWidget::updatePopUpMenu()
{
  if(rightClickPopUp == nullptr)
  {
    // popup used the 1st time - we need to create it:
    rightClickPopUp = new RPopUpMenu(this);
    rightClickPopUp->registerPopUpMenuObserver(this);  
    rightClickPopUp->setDismissOnFocusLoss(true);
    addChildWidget(rightClickPopUp, false, false);
  }
  rightClickPopUp->clear();
  addPopUpMenuItems();
}

void AutomatableWidget::addPopUpMenuItems()
{
  addPopUpMidiItems();
}

void AutomatableWidget::addPopUpMidiItems()
{
  if( assignedParameter != NULL )
  {
    AutomatableParameter *ap;
    ap = dynamic_cast<AutomatableParameter*> (assignedParameter);
    if( ap != NULL )
    {
      // prepare some strings for the popup menu:
      int cc = ap->getAssignedMidiController();
      String ccString;
      if( cc > -1 )
        ccString = String("(currently CC") + String(cc) + String(")");
      else
        ccString = String("(currently none)"); 

      int defaultCc = ap->getDefaultMidiController();
      String defaultString;
      if( defaultCc > -1 )
        defaultString = String("CC") + String(defaultCc);
      else
        defaultString = String("none");
      String minString = stringConversionFunction(ap->getLowerAutomationLimit());
      String maxString = stringConversionFunction(ap->getUpperAutomationLimit());

      rightClickPopUp->addItem(MIDI_LEARN,  String("MIDI learn ") + ccString);
      rightClickPopUp->addItem(MIDI_ASSIGN, String("MIDI assign"));
      rightClickPopUp->addItem(MIDI_MIN,    String("use value as lower limit (currently ") + minString + String(")"));
      rightClickPopUp->addItem(MIDI_MAX,    String("use value as upper limit (currently ") + maxString + String(")"));
      rightClickPopUp->addItem(MIDI_REVERT, String("revert MIDI mapping to defaults") );
    }
  }
}

void AutomatableWidget::openRightClickPopupMenu()
{
  updatePopUpMenu();
  int w = jmax(getWidth(), rightClickPopUp->getRequiredWidth(true));
  int h = jmin(200,        rightClickPopUp->getRequiredHeight(true));
  //rightClickPopUp->show(false, RPopUpComponent::BELOW, w, h); // showModally = false
  rightClickPopUp->show(true, RPopUpComponent::BELOW, w, h); // showModally = true
  // If we don't show it modally (1st parameter = true), it will be immediately dismissed
  // after opening (so it appears as if it doesn't open at all). We could avoid it by calling
  // setDismissOnFocusLoss(false) in our constructor, but then it will stay open all the time
  // until we choose some option.
}

double AutomatableWidget::openModalNumberEntryField()
{
  AutomatableParameter *ap = dynamic_cast<AutomatableParameter*> (assignedParameter);
  if( ap == NULL )
    return 0.0;

  RTextEntryField *entryField = new RTextEntryField( String(ap->getValue()) );
  //entryField->setBounds(handleRectangle);
  entryField->setBounds(2, 2, getWidth()-4, getHeight()-4);
  entryField->setColourScheme(getColourScheme());
  addAndMakeVisible(entryField);
  entryField->setPermittedCharacters(String("0123456789.-"));
  entryField->selectAll();

  entryField->runModalLoop(); // should not be used according to doc...
  // entryField->enterModalState(true);  // ...but this doesn't work at all
  // maybe we should keep an RTextEntryField member and register ourselves as observer

  double result = entryField->getText().getDoubleValue();
  removeChildComponent(entryField);
  delete entryField;
  return result;
}
