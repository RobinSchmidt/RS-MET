//#include "rojue_RWidget.h"
//#include "rojue_RTextField.h"
//#include "../editors/rojue_ColourSchemeComponent.h"
//using namespace rojue;

const BitmapFont* RWidget::font = &BitmapFontRoundedBoldA10D0::instance;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RWidget::RWidget(const String& newDescription) 
{
  localAutomationSwitch    = true;
  isGuiElement             = true;
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
  ParameterObserver::localAutomationSwitch = false;
  if( assignedParameter != NULL )
    assignedParameter->deRegisterParameterObserver(this);

  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void RWidget::assignParameter(Parameter *parameterToAssign)
{
  if( assignedParameter != NULL )
    assignedParameter->deRegisterParameterObserver(this);

  assignedParameter = parameterToAssign;

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
  //tmpColourScheme.setColourSchemeFromXml(xml);
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
  // send changeMessage to ourselves in order to do the update of the slider in the message thread:
  //sendChangeMessage(this);

  triggerAsyncUpdate();
}

void RWidget::parameterRangeChanged(Parameter* parameterThatHasChanged)
{
  triggerAsyncUpdate();
}

void RWidget::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{
  if( assignedParameter == parameterThatWillBeDeleted )
  {
    assignedParameter->deRegisterParameterObserver(this);
    assignedParameter = NULL;
  }
}

void RWidget::handleAsyncUpdate()
{
  if( assignedParameter != NULL )
  {
    localAutomationSwitch = false; // to avoid circular notifications and updates
    updateWidgetFromAssignedParameter();
    localAutomationSwitch = true;
  }
}

/*
void RWidget::changeListenerCallback(void *objectThatHasChanged)
{
  if( objectThatHasChanged == this )
  {
    if( assignedParameter != NULL )
    {
      // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
      // class off to avoid circular notifications and updates:
      localAutomationSwitch = false;
  
      // call the method which updates the widget:
      updateWidgetFromAssignedParameter();
      //setValue(assignedParameter->getValue());
  
      // switch the wantsAutomationNotification flag on again:  
      localAutomationSwitch = true;
    }
  }
}
*/

/*
void RWidget::mouseEnter(const juce::MouseEvent &e)
{
  if( descriptionField != NULL )
    descriptionField->setText(description, false);
}

void RWidget::mouseExit(const MouseEvent &e)
{
  if( descriptionField != NULL )
    descriptionField->setText(String::empty, false);
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

/*
void RWidget::setValueFromString(const juce::String &valueString, bool sendChangeMessage)
{

}
*/

void RWidget::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), RWidget::outlineThickness);
  //g.fillAll(Colours::red); test
}

/*
double RWidget::getValue()
{
  return 1.0;
}
*/