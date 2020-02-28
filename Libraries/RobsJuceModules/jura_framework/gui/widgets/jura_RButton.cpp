//#include "rojue_RButton.h"
//using namespace rojue;

//const BitmapFont* RButton::font = &BitmapFontRoundedBoldA10D0::instance;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RButton::RButton(int newSymbolIndex) : RWidget(String("SymbolButton"))
{
  symbolIndex       = newSymbolIndex;
  text              = String();
  isOn              = false;
  clickTogglesState = true;
}

RButton::RButton(const String& buttonText) : RWidget(buttonText)
{
  symbolIndex       = 0;
  text              = buttonText;
  isOn              = false;
  clickTogglesState = true;
}

RButton::~RButton()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void RButton::setSymbolIndex(int newSymbolIndex)
{
  symbolIndex = newSymbolIndex;
}

void RButton::setButtonText(const String& newText) throw()
{
  if(text != newText)
  {
    text = newText;
    repaint();
  }
}

void RButton::setClickingTogglesState(const bool shouldToggle)
{
  clickTogglesState = shouldToggle;
}

void RButton::setToggleState(const bool shouldBeOn, const bool sendChangeNotification)
{
  if(shouldBeOn != isOn)
  {
    isOn = shouldBeOn;
    repaintOnMessageThread();
    if(sendChangeNotification)
      sendClickMessage();
  }
}

void RButton::setStateFromString(const juce::String &stateString, bool sendChangeMessage)
{
  if( stateString.getIntValue() == 0 )
    setToggleState(false, sendChangeMessage);
  else
    setToggleState(true,  sendChangeMessage);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

juce::String RButton::getStateAsString() const
{
  if( getToggleState() == true )
    return "1";
  else
    return "0";
}

//-------------------------------------------------------------------------------------------------
// others:

void RButton::addRButtonListener (RButtonListener* const newListener) throw()
{
  jassert(newListener != 0);
  jassert(!buttonListeners.contains(newListener)); // trying to add a listener to the list twice!
  if(newListener != 0)
    buttonListeners.add(newListener);
}

void RButton::removeRButtonListener(RButtonListener* const listener) throw()
{
  jassert(buttonListeners.contains(listener)); // trying to remove a listener that isn't on the list!
  buttonListeners.removeFirstMatchingValue(listener);
}

void RButton::sendClickMessage()
{
  for(int i=0; i<buttonListeners.size(); i++)
    buttonListeners[i]->rButtonClicked(this);
}

void RButton::mouseDown(const MouseEvent& e)
{
  clicked();
}

void RButton::enablementChanged()
{
  repaint();
}

void RButton::updateWidgetFromAssignedParameter(bool sendChangeMessage)
{
  if( assignedParameter != NULL )
    setToggleState(assignedParameter->getValue() >= 0.5, sendChangeMessage);
}

void RButton::clicked()
{
  if( clickTogglesState )
    isOn = !isOn;
  if( assignedParameter != NULL )
    assignedParameter->setValue((double) getToggleState(), true, true);
  sendClickMessage();
  repaint();
}

void RButton::paint(Graphics &g)
{
  if(painter != nullptr) { painter->paint(g, this); return; }

  //int x = 0;
  //int y = 0;
  //int w = getWidth();
  //int h = getHeight();

  // draw background and outline (flat 2D looking):
  if( getToggleState() ) // || isDown() )
    g.fillAll(getHandleColour());
  else
    g.fillAll(getBackgroundColour());
  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), 2);

  // draw the text or symbol:
  //int y = getHeight()/2 - 3;
  g.setColour(getTextColour());
  if( symbolIndex <= 0 || symbolIndex > NUM_SYMBOLS )
  {
    int x = getWidth()/2
      - font->getTextPixelWidth(getButtonText(), font->getDefaultKerning()) / 2;
    int y = getHeight()/2 - font->getFontAscent()/2;
    drawBitmapFontText(g, x, y, getButtonText(), font, getTextColour());
  }
  else
    drawSymbol(g);

  // gray out the button if it's disabled:
  //if( !isEnabled() )
  //  g.fillAll(Colours::lightgrey.withAlpha(0.75f));
}

void RButton::drawSymbol(Graphics &g) const
{
  float w = (float) getWidth();
  float h = (float) getHeight();
  float m = 4.f;  // margin in pixels

  if( h >= 20.f )
    m = 5.f;

  float x1, x2, y1, y2;

  switch( symbolIndex )
  {
  case PLUS:
    {
      x1 = 0.5f*getWidth();
      x2 = x1;
      y1 = getHeight()-4.f;
      y2 = 4.f;
      g.drawLine(x1, y1, x2, y2, 2.f);

      x1 = 4.f;
      x2 = getWidth()-4.f;
      y1 = 0.5f*(getHeight());
      y2 = y1;
      g.drawLine(x1, y1, x2, y2, 2.f);
    }
    break;

  case MINUS:
    {
      x1 = 4.f;
      x2 = getWidth()-4.f;
      y1 = 0.5f*(getHeight());
      y2 = y1;
      g.drawLine(x1, y1, x2, y2, 2.f);
    }
    break;

  case ARROW_UP:
    {
      Path path;
      path.addTriangle(m, h-m, 0.5f*w, m, w-m, h-m);
      g.fillPath(path);
    }
    break;

  case ARROW_DOWN:
    {
      Path path;
      path.addTriangle(m, m, w-m, m, 0.5f*w, h-m);
      g.fillPath(path);
    }
    break;

  case ARROW_LEFT:
    {
      Path path;
      path.addTriangle(m, 0.5f*h, w-m, m, w-m, h-m);
      g.fillPath(path);
    }
    break;

  case ARROW_RIGHT:
    {
      Path path;
      path.addTriangle(m, m, m, h-m, w-m, 0.5f*h);
      g.fillPath(path);
    }
    break;

  case PLAY:
    {
      Path path;
      path.addTriangle(12.f, 8.f, getWidth()-12.f, 0.5f*getHeight(), 12.f, getHeight()-8.f);
      g.fillPath(path);
    }
    break;

  case SKIP_FORWARD:
    {
      Path path;
      path.addTriangle(8.f, 4.f, getWidth()-8.f, 0.5f*getHeight(), 8.f, getHeight()-4.f);
      path.addLineSegment(Line<float>(getWidth()-8.f, 4.f, getWidth()-8.f, getHeight()-4.f), 2.f);
      g.fillPath(path);
    }
    break;

  case SKIP_BACK:
    {
      Path path;
      path.addTriangle(getWidth()-8.f, 4.f, 8.f, 0.5f*getHeight(), getWidth()-8.f, getHeight()-4.f);
      path.addLineSegment(Line<float>(8.f, 4.f, 8.f, getHeight()-4.f), 2.f);
      g.fillPath(path);
    }
    break;

  case CLOSE:
    {
      g.drawLine(                  4.f, 4.f, (float)getWidth()-4.f, (float)getHeight()-4.f, 2.f);
      g.drawLine((float)getWidth()-4.f, 4.f,                   4.f, (float)getHeight()-4.f, 2.f);
    }
    break;
  }
}


//=================================================================================================
// class RClickButton:

RClickButton::RClickButton(int newSymbolIndex) : RButton(newSymbolIndex)
{
  clickTogglesState = false;
}

RClickButton::RClickButton(const juce::String& buttonText) : RButton(buttonText)
{
  clickTogglesState = false;
}

void RClickButton::mouseDown(const MouseEvent& e)
{
  isOn = true;
  sendClickMessage();
  repaint();
}

void RClickButton::mouseUp(const MouseEvent& e)
{
  isOn = false;
  //sendClickMessage();
  repaint();
}


//=================================================================================================
// class RClickButtonNotifyOnMouseUp:

RClickButtonNotifyOnMouseUp::RClickButtonNotifyOnMouseUp(int newSymbolIndex)
  : RClickButton(newSymbolIndex)
{

}

RClickButtonNotifyOnMouseUp::RClickButtonNotifyOnMouseUp(const juce::String& buttonText)
  : RClickButton(buttonText)
{

}

void RClickButtonNotifyOnMouseUp::mouseDown(const MouseEvent& e)
{
  isOn = true;
  repaint();
}

void RClickButtonNotifyOnMouseUp::mouseUp(const MouseEvent& e)
{
  if( isOn && contains(Point<int>(e.x, e.y)) )
    sendClickMessage();
  isOn = false;
  repaint();
}


//=================================================================================================
// class RClickButtonWithAutoRepeat:

RClickButtonWithAutoRepeat::RClickButtonWithAutoRepeat(int newSymbolIndex)
  : RClickButton(newSymbolIndex)
{
  initialDelay = 200;
  timeInterval = 50;
}

RClickButtonWithAutoRepeat::RClickButtonWithAutoRepeat(const juce::String& buttonText)
  : RClickButton(buttonText)
{

}

void RClickButtonWithAutoRepeat::mouseDown(const MouseEvent& e)
{
  isOn = true;
  repaint();
  sendClickMessage();
  startTimer(initialDelay);
}

void RClickButtonWithAutoRepeat::mouseUp(const MouseEvent& e)
{
  isOn = false;
  repaint();
  stopTimer();
}

void RClickButtonWithAutoRepeat::timerCallback()
{
  sendClickMessage();
  startTimer(timeInterval);
}


//=================================================================================================
// class RRadioButton:

RRadioButton::RRadioButton(int newSymbolIndex) : RButton(newSymbolIndex)
{
  radioGroupToUse = NULL;
}

RRadioButton::RRadioButton(const juce::String& buttonText) : RButton(buttonText)
{
  radioGroupToUse = NULL;
}

void RRadioButton::addToRadioButtonGroup(RRadioButtonGroup *newGroupToUse)
{
  if( radioGroupToUse != NULL )
    radioGroupToUse->removeButtonFromRadioGroup(this);

  radioGroupToUse = newGroupToUse;

  if( radioGroupToUse != NULL )
    radioGroupToUse->addButtonToRadioGroup(this);
}

void RRadioButton::clicked()
{
  setToggleState(!isOn, true);
}


void RRadioButton::setToggleState(const bool shouldBeOn, const bool sendNotifications)
{
  jassert( radioGroupToUse != NULL ); // forgotten to assign this radio-button to a radio-group?
  if( radioGroupToUse != NULL )
  {
    if( shouldBeOn )
      radioGroupToUse->toggleRadioButtonOn(this, sendNotifications);
  }
}

void RRadioButtonGroup::addButtonToRadioGroup(RRadioButton *buttonToAdd)
{
  radioButtons.addIfNotAlreadyThere(buttonToAdd);
}

void RRadioButtonGroup::removeButtonFromRadioGroup(RRadioButton *buttonToRemove)
{
  radioButtons.removeFirstMatchingValue(buttonToRemove);
}

void RRadioButtonGroup::toggleRadioButtonOn(RRadioButton *buttonToToggleOn, bool sendNotifications)
{
  jassert( radioButtons.contains(buttonToToggleOn) ); // the button passed does not belong into this group
  for(int i=0; i<radioButtons.size(); i++)
  {
    RButton *button = (RButton*) radioButtons[i];     // upcast because we want to invoke the baseclass' method
    if( radioButtons[i] == buttonToToggleOn )
      button->setToggleState(true, sendNotifications);
    else
      button->setToggleState(false, sendNotifications);
  }
}

bool RRadioButtonGroup::isButtonMemberOfGroup(RButton *buttonToCheck)
{
  RRadioButton *radioButton = dynamic_cast<RRadioButton*> (buttonToCheck);
  if(radioButton == nullptr)
    return false;
  else
    return radioButtons.contains(radioButton);
}

//=================================================================================================
// class RHyperlinkButton

RHyperlinkButton::RHyperlinkButton(const String& linkText, const URL& linkURL) : RButton(linkText), url(linkURL)
{
  setMouseCursor(MouseCursor::PointingHandCursor);
}

RHyperlinkButton::~RHyperlinkButton()
{

}

void RHyperlinkButton::setURL (const URL& newURL) throw()
{
  url = newURL;
}

void RHyperlinkButton::clicked()
{
  if(url.isWellFormed())
    url.launchInDefaultBrowser();
}

void RHyperlinkButton::paint(Graphics& g)
{
  drawBitmapFontText(g, 0, 0, getButtonText(), font, getTextColour());
}

//=================================================================================================

void RButtonPainter3D::paint(Graphics& g, RButton *button)
{
  // experimental

  // todo - use w,h variables for width and height

  float threeDeeDepth = 1.0f;
  float x1, x2, y1, y2;
  float thickness = (float) button->outlineThickness;
  if(button->getToggleState())
  {
    g.fillAll(button->getHandleColour());
    g.setColour(button->getHandleColour().brighter(2.f*threeDeeDepth));

    // right line:
    x1 = (float)button->getWidth() - 0.5f * thickness;
    x2 = (float)button->getWidth() - 0.5f * thickness;
    y1 = 0.f;
    y2 = (float)button->getHeight();
    g.drawLine(x1, y1, x2, y2, thickness);

    // bottom line:
    x1 = 0.f;
    x2 = (float)button->getWidth();
    y1 = (float)button->getHeight() - 0.5f * thickness;
    y2 = (float)button->getHeight() - 0.5f * thickness;
    g.drawLine(x1, y1, x2, y2, thickness);

    g.setColour(button->getHandleColour().darker(threeDeeDepth));

    // left line:
    x1 = 0.5f * thickness;
    x2 = 0.5f * thickness;
    y1 = 0.f;
    y2 = (float)button->getHeight();
    g.drawLine(x1, y1, x2, y2, thickness);

    // top line:
    x1 = 0.f;
    x2 = (float)button->getWidth();
    y1 = 0.5f * thickness;
    y2 = 0.5f * thickness;
    g.drawLine(x1, y1, x2, y2, thickness);
  }
  else
  {
    g.fillAll(button->getBackgroundColour());

    g.setColour(button->getHandleColour().darker(threeDeeDepth));

    // right line:
    x1 = (float)button->getWidth() - 0.5f * thickness;
    x2 = (float)button->getWidth() - 0.5f * thickness;
    y1 = 0.f;
    y2 = (float)button->getHeight();
    g.drawLine(x1, y1, x2, y2, thickness);

    // bottom line:
    x1 = 0.f;
    x2 = (float)button->getWidth();
    y1 = (float)button->getHeight() - 0.5f * thickness;
    y2 = (float)button->getHeight() - 0.5f * thickness;
    g.drawLine(x1, y1, x2, y2, thickness);

    g.setColour(button->getHandleColour().brighter(2.f*threeDeeDepth));

    // left line:
    x1 = 0.5f * thickness;
    x2 = 0.5f * thickness;
    y1 = 0.f;
    y2 = (float)button->getHeight();
    g.drawLine(x1, y1, x2, y2, thickness);

    // top line:
    x1 = 0.f;
    x2 = (float)button->getWidth();
    y1 = 0.5f * thickness;
    y2 = 0.5f * thickness;
    g.drawLine(x1, y1, x2, y2, thickness);
  }
}
