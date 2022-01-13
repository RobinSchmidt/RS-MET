const BitmapFont* RWidget::font = &BitmapFontRoundedBoldA10D0::instance;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RWidget::RWidget(const String& newDescription) 
{
  // set up inherited ParameterObserver:
  setLocalAutomationSwitch(true);
  setIsGuiElement(true);
  notifyPreSmoothing(true);
  notifyPostSmoothing(false);

  // set up appearance stuff:
  noBackgroundAndOutline   = false;
  alphaMultiplier          = 1.f;
  description              = newDescription;
  descriptionField         = NULL;
  assignedParameter        = NULL;
  stringConversionFunction = &valueToString2;

  //setMouseClickGrabsKeyboardFocus(false); 
  // otherwise, the qwerty-keyboard in Fruity Loops (and presumably other hosts as well) won't 
  // respond anymore after clicking on a widget on a plugin-GUI ...but it seems to block 
  // mouseWheelMoved messages as well.. mmm

  // we add ourselves as listener in order to send change messages to ourselves (maybe this could be 
  // better handled with AsyncUpdater):
  //addChangeListener(this);
}

RWidget::~RWidget()
{
  // remove ourselves as listener from the Parameter object, such that it does not try to notify a 
  // nonexistent listener:
  setLocalAutomationSwitch(false);
  if( assignedParameter != NULL )
    assignedParameter->deRegisterParameterObserver(this);

  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void RWidget::assignParameter(Parameter *parameterToAssign)
{
  if(assignedParameter == parameterToAssign)
    return; // nothing to do


  if( assignedParameter != NULL )
    assignedParameter->deRegisterParameterObserver(this);

  assignedParameter = parameterToAssign;
  updateWidgetFromAssignedParameter();

  if( assignedParameter != NULL )
    assignedParameter->registerParameterObserver(this);
}

void RWidget::setColourScheme(const WidgetColourScheme &newColourScheme)
{
  if( noBackgroundAndOutline == true )
  {
    WidgetColourScheme tmpColourScheme = newColourScheme;
    ColourSchemeComponent *colorParent = getOutlyingColourSchemeComponent();
    if( colorParent != NULL )
    {
      int editorAppearance = colorParent->getEditorColourScheme().getAppearance();
      tmpColourScheme.setAppearance(editorAppearance);
    }
    tmpColourScheme.outline    = tmpColourScheme.outline.withAlpha(0.f);
    tmpColourScheme.background = tmpColourScheme.background.withAlpha(0.f);
    colourScheme = tmpColourScheme;
  }
  else
    colourScheme = newColourScheme;

  for(int i=0; i<childWidgets.size(); i++)
    childWidgets[i]->setColourScheme(newColourScheme);

  Component* thisAsComponent = dynamic_cast<Component*> (this);
  if( thisAsComponent != NULL )
    thisAsComponent->repaint();
}

void RWidget::setColourSchemeFromXml(const XmlElement &xml)
{
  WidgetColourScheme tmpColourScheme;
  //tmpColourScheme.setColourSchemeFromXml(xml); // why commented?
  setColourScheme(tmpColourScheme);
  // we don't use colourScheme.setColourSchemeFromXml(xml) directly here such that subclasses 
  // need to override only setColourScheme (and not also setColourSchemeFromXml) when they need 
  // special actions on colour-scheme changes)
}

void RWidget::addChildWidget(RWidget *newChild, bool addAsChildComponent, bool makeVisible)
{
  childWidgets.add(newChild);
  if( addAsChildComponent == true )
    addChildComponent(newChild);
  if( makeVisible == true )
    newChild->setVisible(true);
}

void RWidget::removeChildWidget(RWidget* widgetToRemove, bool removeAsChildComponent, 
  bool deleteObject)
{
  jassert( childWidgets.contains(widgetToRemove) ); // trying to remove a child widget that was not added before?
  //childWidgets.removeValue(widgetToRemove);
  childWidgets.removeFirstMatchingValue(widgetToRemove);
  if( removeAsChildComponent == true )
    removeChildComponent(widgetToRemove);
  if( deleteObject == true )
    delete widgetToRemove;
}

void RWidget::setNoBackgroundAndOutline(bool shouldBeWithout)
{
  noBackgroundAndOutline = shouldBeWithout;
  setColourScheme(colourScheme);  
}

void RWidget::enablementChanged()
{
  if( !isEnabled() )
    alphaMultiplier = 0.25f;
  else
    alphaMultiplier = 1.f;
  repaint();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

ColourSchemeComponent* RWidget::getOutlyingColourSchemeComponent()
{
  Component             *parent      = this;  // ...will soon indeed become a parent
  ColourSchemeComponent *colorParent = NULL;
  while( colorParent == NULL )
  {
    parent = parent->getParentComponent();
    if( parent == NULL )
      return NULL;
    else
    {
      colorParent = dynamic_cast<ColourSchemeComponent*> (parent);
      if( colorParent != NULL )
        return colorParent;
    }
  }
  jassertfalse;  // should never be reached
  return NULL;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RWidget::parameterChanged(Parameter* parameterThatHasChanged)
{
  //triggerAsyncUpdate();
  if( assignedParameter != NULL )
  {
    setLocalAutomationSwitch(false); // to avoid circular notifications and updates
    updateWidgetFromAssignedParameter();
    setLocalAutomationSwitch(true);
  }
}

void RWidget::parameterRangeChanged(Parameter* parameterThatHasChanged)
{
  //triggerAsyncUpdate();
  if( assignedParameter != NULL )
  {
    setLocalAutomationSwitch(false); // to avoid circular notifications and updates
    updateWidgetFromAssignedParameter();
    setLocalAutomationSwitch(true);
  }
}

// was formerly parameterIsGoingToBeDeleted
void RWidget::parameterWillBeDeleted(Parameter* p)
{
  if( assignedParameter == p )
  {
    assignedParameter->deRegisterParameterObserver(this);
    assignedParameter = nullptr;
  }
}
/*
void RWidget::handleAsyncUpdate()
{
  if( assignedParameter != NULL )
  {
    setLocalAutomationSwitch(false); // to avoid circular notifications and updates
    updateWidgetFromAssignedParameter();
    setLocalAutomationSwitch(true);
  }
  //repaint();
}
*/
void RWidget::updateWidgetFromAssignedParameter(bool sendChangeMessage)
{
  // needs to be overriden in the subclasses to - for example - update a slider like this:
  // if( assignedParameter != NULL )
  // {
  //  setValue(assignedParameter->getValue());
  // }
}

void RWidget::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), RWidget::outlineThickness);
  //g.fillAll(Colours::red); test
}

double RWidget::openModalNumberEntryField(double numberToShowInitially)
{
  //if( assignedParameter == NULL )
  //  return 0.0;
  //RTextEntryField *entryField = new RTextEntryField( 
  //  String(assignedParameter->getValue()));

  RTextEntryField *entryField = new RTextEntryField(String(numberToShowInitially));
  entryField->setBounds(2, 2, getWidth()-4, getHeight()-4);
  entryField->setColourScheme(getColourScheme());
  addAndMakeVisible(entryField);
  entryField->setPermittedCharacters(String("0123456789.-"));
  entryField->selectAll();

  entryField->runModalLoop(); // should not be used according to doc...
  // entryField->enterModalState(true);  // ...but this doesn't work at all
  // maybe we should keep an RTextEntryField member and register ourselves as observer

  // see here for Elan's solution:
  // https://github.com/RobinSchmidt/RS-MET/issues/221#issuecomment-427450329

  double result = entryField->getText().getDoubleValue();
  removeChildComponent(entryField);
  delete entryField;
  return result;
}
