//#include "rojue_RComboBox.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RComboBox::RComboBox(const String& name) : RTextField(name)
{
  popUpMenu = new RPopUpMenu(this);
  popUpMenu->registerPopUpMenuObserver(this);
  // popUpMenu->setAttachPosition(BELOW); // later
  addChildWidget(popUpMenu, false, false);
  setNoBackgroundAndOutline(false);
  dontOpenPopUpOnNextMouseClick = false;;
}

RComboBox::~RComboBox()
{
  popUpMenu->deRegisterPopUpMenuObserver(this);
  closePopUp();
  delete popUpMenu;
}

//-------------------------------------------------------------------------------------------------
// setup:

void RComboBox::addItem(int itemResultId, const juce::String& itemText, bool isEnabled,
  bool isTicked)
{
  popUpMenu->addItem(itemResultId, itemText, isEnabled, isTicked);
}

void RComboBox::setItemEnabled(int index, bool shouldBeEnabled)
{
  popUpMenu->setItemEnabled(index, shouldBeEnabled);
}

void RComboBox::setItemText(int index, const juce::String& newText)
{
  popUpMenu->setItemText(index, newText);
}

void RComboBox::clear(const bool dontSendChangeMessage)
{
  popUpMenu->clear();
}

void RComboBox::selectItemByIndex(int indexToSelect, bool sendNotification, bool updateParameter)
{
  jassert( indexToSelect >= 0 && indexToSelect < popUpMenu->getNumTopLevelItems() );  // index out of range
  if( indexToSelect < 0 || indexToSelect >= popUpMenu->getNumTopLevelItems() )
    return;

  popUpMenu->selectItemByIndex(indexToSelect, false);

  // old:
  //if( assignedParameter != nullptr )
  //  assignedParameter->setValue((double) indexToSelect, sendNotification, sendNotification);
  // check, if we have a bug here and may run into endless recursive cross-notifications between 
  // the parameter and the combobox - put a breakpoint, look at the stack trace
  // yes: plug in a BreakpointModulator in ToolChain and change the grid setting - endless 
  // recursion
  // possible solutions:
  // -when the parameter notifies the combo-box about a parameter-change, we should recognize that
  //  it is our assigned parameter and not call assignedParameter->setValue here (as it already has
  //  been updated)

  // new:
  if( assignedParameter != nullptr && updateParameter )
    assignedParameter->setValue((double) indexToSelect, true, true);


  setText(getItemText(indexToSelect));

  if( sendNotification == true )
    sendComboBoxChangeNotifications();

  repaintOnMessageThread();
}

void RComboBox::selectItemFromText(const juce::String& textToSelect, bool sendNotification)
{
  // new - todo: check, if this works correctly in all cases:
  popUpMenu->selectItemByText(textToSelect, false); // the popup should not send notifications
  setText(popUpMenu->getSelectedText());

  if(assignedParameter != nullptr)
  {
    rsChoiceParameter* cp = dynamic_cast<rsChoiceParameter*>(assignedParameter);
    if(cp)
      cp->setStringValue(textToSelect, true, true); 
      // can we just call setStringValue for the Parameter baseclass (else-case below) as well and
      // get rid of the dynamic cast and the if? try it! ..that would clean up the code
    else
    {
      for(int i = 0; i < popUpMenu->getNumTopLevelItems(); i++)
        if(getItemText(i) == textToSelect)
          assignedParameter->setValue((double)i, true, true);
    }
  }

  if( sendNotification == true )
    sendComboBoxChangeNotifications();

  /*
  // old - works only for a flat array of options:
  //int numItems = popUpMenu->getNumTopLevelItems(); // for debug
  for(int i=0; i<popUpMenu->getNumTopLevelItems(); i++)
  {
    //juce::String itemText = getItemText(i); // for debug
    if( getItemText(i) == textToSelect )
    {
      selectItemByIndex(i, sendNotification);
      return;
    }
  }
  jassertfalse; // the passed text is not among the items
  // this function works only when the combobox has a flat array of options - with a tree
  // structure, it fails - we need a function getNumLeafNodes instead of getNumItems ...
  // or something like that
  */
}

void RComboBox::setStateFromString(const juce::String &stateString, bool sendChangeMessage)
{
  selectItemFromText(stateString, sendChangeMessage);
}

void RComboBox::registerComboBoxObserver(RComboBoxObserver* const observerToRegister)
{
  comboBoxObservers.addIfNotAlreadyThere(observerToRegister);
}

void RComboBox::deRegisterComboBoxObserver(RComboBoxObserver* const observerToDeRegister)
{
  comboBoxObservers.removeFirstMatchingValue(observerToDeRegister);
}

void RComboBox::assignParameter(Parameter* parameterToAssign)
{
  //RWidget::assignParameter(parameterToAssign);
  if( assignedParameter != NULL )
    assignedParameter->deRegisterParameterObserver(this);
  assignedParameter = NULL;
  clear(true);
  assignedParameter = parameterToAssign;
  if( assignedParameter != NULL )
  {
    assignedParameter->registerParameterObserver(this);
    jassert( assignedParameter->isStringParameter() ); // use Parameter::STRING for the scaling
                                                       // when you attach a combobox to a parameter
    rsChoiceParameter* cp = dynamic_cast<rsChoiceParameter*>(assignedParameter);
    if(cp)
      for(int i = 0; i < assignedParameter->getNumStringValues(); i++)
        addItem(cp->getChoiceEnumValue(i), assignedParameter->getOptionStringAtIndex(i));
    else
      for(int i = 0; i < assignedParameter->getNumStringValues(); i++)
        addItem(i, assignedParameter->getOptionStringAtIndex(i));


    updateWidgetFromAssignedParameter(false);
  }
}

void RComboBox::parameterChanged(Parameter* p)
{
  RTextField::parameterChanged(p);
  selectItemByIndex((int)p->getValue(), true, false);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

const String RComboBox::getItemText(const int index) const
{
  jassert( index >= 0 && index < popUpMenu->getNumTopLevelItems() );  // index out of range
  if( index >= 0 && index < popUpMenu->getNumTopLevelItems() )
    return popUpMenu->getItemByIndex(index)->getNodeText();
  else
    return String();
}

int RComboBox::getSelectedItemIdentifier() const
{
  RTreeViewNode* selectedItem = getSelectedItem();
  if( selectedItem != NULL )
    return selectedItem->getNodeIdentifier();
  else
    return -1;
}

juce::String RComboBox::getSelectedItemText() const
{
  RTreeViewNode* selectedItem = getSelectedItem();
  if( selectedItem != NULL )
    return selectedItem->getNodeText();
  else
    return juce::String("nothing selected");
}

juce::String RComboBox::getStateAsString() const
{
  return getSelectedItemText();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RComboBox::rPopUpDismissedByClickOnOwner(ROwnedPopUpComponent *popUp)
{
  dontOpenPopUpOnNextMouseClick = true;
}

void RComboBox::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  if( menuThatHasChanged == popUpMenu )
  {
    closePopUp();
    //sendComboBoxChangeNotifications();
    selectItemFromText(popUpMenu->getSelectedItem()->getNodeText(), true);
    repaint();
  }
}

void RComboBox::updateWidgetFromAssignedParameter(bool sendNotification)
{
  if(assignedParameter != nullptr)
  {
    rsChoiceParameter* cp = dynamic_cast<rsChoiceParameter*>(assignedParameter);
    if(cp)
      selectItemFromText(assignedParameter->getStringValue(), sendNotification);
    else
      selectItemByIndex((int)assignedParameter->getRawValue(), sendNotification, false);
  }
}

void RComboBox::mouseDown(const MouseEvent& e)
{
  if( isEnabled() && e.eventComponent == this )
  {
    if (assignedParameter != NULL && e.mods.isCommandDown())
    {
      selectItemByIndex((int)assignedParameter->getDefaultValue(), true, true);
      return;
    }

    if( dontOpenPopUpOnNextMouseClick == false )
      openPopUp();
    dontOpenPopUpOnNextMouseClick = false;
  }
}

void RComboBox::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), 2);
  int x = 4;
  int y = getHeight()/2 - font->getFontAscent()/2;
  drawBitmapFontText(g, x, y, getSelectedItemText(), font, getTextColour());

  //grayOutIfDisabled(g);
}

//-------------------------------------------------------------------------------------------------
// others:

void RComboBox::sendComboBoxChangeNotifications()
{
  for(int i=0; i<comboBoxObservers.size(); i++)
    comboBoxObservers[i]->rComboBoxChanged(this);
}

void RComboBox::openPopUp()
{
  int w  = jmin(maxPopUpWidth,  getWidth());
  int h  = jmin(maxPopUpHeight, popUpMenu->getRequiredHeight(true));
  w  = jmax(minPopUpWidth,  w);
  h  = jmax(minPopUpHeight, h);

  popUpMenu->setSize(w, h);
    // maybe introduce a max-height or something - it's strange but if don't set the size here, the popup appears black
  //popUpMenu->showAttachedTo(this, true, RPopUpComponent::BELOW, getWidth(), popUpMenu->getRequiredHeight(), 0, -outlineThickness);
  popUpMenu->show(true, RPopUpComponent::BELOW, w, h, 0, -outlineThickness);
}

void RComboBox::closePopUp()
{
  popUpMenu->dismiss();
}

//=================================================================================================
// class RNamedComboBox:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RNamedComboBox::RNamedComboBox(const String& componentName, const String& comboBoxName)
: RComboBox(componentName)
{
  addAndMakeVisible( nameLabel = new RTextField(comboBoxName) );
  nameLabel->setNoBackgroundAndOutline(true);
  nameLabel->setJustification(Justification::centredLeft);
  nameLabelWidth = font->getTextPixelWidth(comboBoxName, font->getDefaultKerning()) + 8;
  nameLabelPosition = LEFT_TO_BOX;
}

RNamedComboBox::~RNamedComboBox()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void RNamedComboBox::setColourScheme(const WidgetColourScheme& newColourScheme)
{
  RComboBox::setColourScheme(newColourScheme);
  nameLabel->setColourScheme(newColourScheme);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RNamedComboBox::paint(Graphics &g)
{
  int x, y, w, h;
  if( nameLabelPosition == LEFT_TO_BOX )
  {
    x = nameLabel->getRight();
    y = 0;
    w = getWidth()-x;
  }
  else  // above
  {
    x = 0;
    y = getHeight()/2;
    w = getWidth()-x;
  }
  h = nameLabel->getHeight();

  g.setColour(getBackgroundColour());
  g.fillRect(x, y, w, h);
  g.setColour(getOutlineColour());
  g.drawRect(x, y, w, h, 2);

  if( nameLabelPosition == LEFT_TO_BOX )
  {
    x = nameLabel->getRight()+4;
    y = getHeight()/2 - font->getFontAscent()/2;
  }
  else
  {
    x = 4;
    y = 3*getHeight()/4 - font->getFontAscent()/2;
  }

  drawBitmapFontText(g, x, y, getSelectedItemText(), font, getTextColour());

  //if( !isEnabled() )
  //  g.fillAll(Colours::lightgrey.withAlpha(0.75f));
}

void RNamedComboBox::resized()
{
  if( nameLabelPosition == LEFT_TO_BOX )
    nameLabel->setBounds(0, 0, nameLabelWidth, getHeight());
  else if( nameLabelPosition == ABOVE_BOX )
    nameLabel->setBounds(0, 0, nameLabelWidth, getHeight()/2);
  RComboBox::resized();
}

//-------------------------------------------------------------------------------------------------
// others:

void RNamedComboBox::openPopUp()
{
  int x = 0;
  int w = getWidth();
  int h = popUpMenu->getRequiredHeight(true);
  if( nameLabelPosition == LEFT_TO_BOX )
  {
    x += nameLabelWidth+4;
    w -= nameLabelWidth+4;
  }
  popUpMenu->setSize(w, h);

  //popUpMenu->show(false, RPopUpComponent::BELOW, w, h, x, -outlineThickness);
  popUpMenu->show(true, RPopUpComponent::BELOW, w, h, x, -outlineThickness);
    // true for the 1st argument (showModally) causes problems - why? reproduce: in EngineersFilter
    // open/close the same combobox twice ...but with false, the boxes for the color-setup don't
    // respond...
}
