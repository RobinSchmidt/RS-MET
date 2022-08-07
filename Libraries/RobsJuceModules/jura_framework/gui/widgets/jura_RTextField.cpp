
//const BitmapFont* RTextField::font = &BitmapFontRoundedBoldA10D0::instance;

RTextField::RTextField(const juce::String& initialText) : justification(Justification::centredLeft)
{ 
  text                   = initialText; 
  justification          = Justification::centredLeft;
  noBackgroundAndOutline = true;
}

void RTextField::setText(const juce::String &newText) 
{ 
  text = newText; 
  repaintOnMessageThread();
  //repaint();
}
   
void RTextField::setJustification(const Justification& newJustification)
{
  justification = newJustification;
  repaint();
}

 void RTextField::setStateFromString(const juce::String &stateString, bool sendChangeMessage)
 {
   setText(stateString);
 }

int RTextField::getTextPixelPositionX() const
{
  //const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  int x = 0;
  int hFlags = justification.getOnlyHorizontalFlags();
  if( justification.testFlags(hFlags & Justification::centredLeft) )
    x = horizontalMargin;
  else if( justification.testFlags(hFlags & Justification::centred) )
    x = getWidth()/2 - font->getTextPixelWidth(text, font->getDefaultKerning())/2;
  else if( justification.testFlags(hFlags & Justification::centredRight) )
    x = getWidth() - font->getTextPixelWidth(text, font->getDefaultKerning()) - horizontalMargin;

  return x;
}

int RTextField::getTextPixelPositionY() const
{
  //const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  return getHeight()/2 - font->getFontAscent()/2;
}

int RTextField::getWidthToFitText() const
{
  //const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  return font->getTextPixelWidth(text, font->getDefaultKerning()) + 2*outlineThickness 
    + 2*horizontalMargin;
}


juce::String RTextField::getStateAsString() const
{
  return getText();
}

void RTextField::assignParameter(Parameter* parameterToAssign)
{
  RWidget::assignParameter(parameterToAssign);
  updateWidgetFromAssignedParameter(false);
}

void RTextField::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( assignedParameter != NULL && parameterThatHasChanged == assignedParameter )
    setText(String(assignedParameter->getValue()));
}

void RTextField::updateWidgetFromAssignedParameter(bool sendChangeMessage)
{
  if( assignedParameter != NULL )
    setText(String(assignedParameter->getValue()));
}

void RTextField::paint(Graphics &g)
{
  //const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  if( noBackgroundAndOutline == false )
    RWidget::paint(g);
  drawBitmapFontText(g, getTextPixelPositionX(), getTextPixelPositionY(), text, font, 
    getTextColour());
}

//=================================================================================================
// class RTextEntryField:

RTextEntryField::RTextEntryField(const juce::String& initialText) : RTextField(initialText)
{
  noBackgroundAndOutline = false;
  caretPosition          = 0;
  selectionStart         = 0;
  selectionEnd           = 0; 
  caretVisible           = false;
  replaceMode            = false;
  textInvalid            = false;
  oldText                = text;
  setWantsKeyboardFocus(true);
}

void RTextEntryField::setPermittedCharacters(const String& newCharacters, bool deleteNonPermittedCharsFromText)
{
  permittedCharacters = newCharacters;
  if( deleteNonPermittedCharsFromText )
  {
    deSelect();
    caretPosition = 0;
    juce::String oldText = text;
    text = text.retainCharacters(permittedCharacters);
    if( text != oldText )
    {
      sendTextChangeNotificationIfChanged();
      repaint();
    }
  }
}

void RTextEntryField::registerTextEntryFieldObserver(RTextEntryFieldObserver *observerToRegister)
{
  textEntryFieldObservers.add(observerToRegister);
}

void RTextEntryField::deRegisterTextEntryFieldObserver(RTextEntryFieldObserver *observerToDeRegister)
{
  textEntryFieldObservers.removeFirstMatchingValue(observerToDeRegister);
}

void RTextEntryField::selectAll()
{
  selectionStart = 0;
  selectionEnd   = text.length();
}
 
void RTextEntryField::deSelect()
{
  selectionStart = 0;
  selectionEnd   = 0;
}

void RTextEntryField::markTextAsInvalid(bool shouldBeMarkedInvalid)
{
  textInvalid = shouldBeMarkedInvalid;
  repaint();
}

Point<int> RTextEntryField::getCaretPixelPosition() const
{
  int x = characterIndexToPixelPosition(caretPosition);
  int y = 2;
  return Point<int>(x, y);
}

int RTextEntryField::characterIndexToPixelPosition(int index) const
{
  //const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  juce::String subString = text.substring(0, caretPosition);
  int x = horizontalMargin + font->getTextPixelWidth(subString, font->getDefaultKerning());
  return x;
}

int RTextEntryField::pixelPositionToCharacterIndex(int pixelX) const
{
  //const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  int i = 0;
  int x = horizontalMargin;
  int w = 0;
  while( x < pixelX && i < text.length() )
  {
    w  = font->getGlyphWidth(text[i]) + font->getDefaultKerning();
    x += w;
    i++;
  }
  if( i > 0 && x-pixelX > w/2 )
    i--; // move one character leftward if we are before the midpoint
  return jlimit(0, text.length(), i);
}

String RTextEntryField::getTextBeforeSelection() const
{
  return text.substring(0, selectionStart);
}

String RTextEntryField::getTextInsideSelection() const
{
  return text.substring(selectionStart, selectionEnd);
}

String RTextEntryField::getTextAfterSelection() const
{
  return text.substring(selectionEnd, text.length());
}

void RTextEntryField::getHeadAndTailString(String &head, String &tail) const
{
  if( isSelectionEmpty() )
  {
    head = text.substring(0, caretPosition);
    tail = text.substring(caretPosition);
  }
  else
  {
    head = getTextBeforeSelection();
    tail = getTextAfterSelection();
  }  
}

void RTextEntryField::setStateFromString(const juce::String &stateString, bool sendChangeMessage)
{
  RTextField::setText(stateString);
  if( sendChangeMessage == true )
    sendTextChangeNotificationIfChanged();
}

void RTextEntryField::inputAttemptWhenModal()
{
  text = oldText;
  if( isCurrentlyModal() )
    exitModalState(0);
}

bool RTextEntryField::keyPressed(const KeyPress &key)
{
  String head, tail;
  if( key.getModifiers().isCommandDown() )
  {
    if( key == KeyPress('c', ModifierKeys::commandModifier, 0) )
    {     
      if( !isSelectionEmpty() )
        SystemClipboard::copyTextToClipboard(getTextInsideSelection());
      return true;
    }
    else if( key == KeyPress('v', ModifierKeys::commandModifier, 0) )
    {
      getHeadAndTailString(head, tail);
      String stringToInsert(SystemClipboard::getTextFromClipboard());
      text          = head + stringToInsert + tail;
      caretPosition = head.length() + stringToInsert.length();
      sendTypingNotification();
      repaint();
    }
  }
  else if( key == KeyPress::insertKey )
  {
    replaceMode = !replaceMode;
  }
  else if( key == KeyPress::returnKey )
  {
    sendTextChangeNotificationIfChanged();
    if( isCurrentlyModal() )
      exitModalState(0);
  }
  else if( key == KeyPress::escapeKey )
  {
    text = oldText;
    if( isCurrentlyModal() )
      exitModalState(0);
  }
  else if( key == KeyPress::backspaceKey )
  {
    getHeadAndTailString(head, tail);
    if( isSelectionEmpty() )
      head = head.substring(0, head.length()-1);
    text = head + tail;
    caretPosition = head.length();
    sendTypingNotification();
    repaint();
  }
  else if( key == KeyPress::deleteKey )
  {
    getHeadAndTailString(head, tail);
    if( isSelectionEmpty() )
      tail = tail.substring(1, tail.length());
    text = head + tail;
    caretPosition = head.length();
    sendTypingNotification();
    repaint();
  }
  else if( key == KeyPress::leftKey )
  {
    caretPosition--;
    caretPosition = jmax(0, caretPosition);
    repaint();
  }
  else if( key == KeyPress::rightKey )
  {
    caretPosition++;
    caretPosition = jmin(caretPosition, text.length());
    repaint();
  }
  else if( key == KeyPress::homeKey )
  {
    caretPosition = 0;
    repaint();
  }
  else if( key == KeyPress::endKey )
  {
    caretPosition = text.length();
    repaint();
  }
  else
  {
    juce_wchar ch = key.getTextCharacter();

    //if(key == KeyPress::numberPadDecimalPoint) ch = '.';    // Allow decimal point be entered on german number pad
    // doesn't work

    if(ch == ',')
      ch = '.';

    if( permittedCharacters != String() && !permittedCharacters.containsChar(ch) )
      return true; // character is not among the permitted ones 
    getHeadAndTailString(head, tail);
    if( isSelectionEmpty() && replaceMode == true )
      tail = tail.substring(1, tail.length());
    text = head + ch + tail;
    caretPosition = head.length() + 1;
    sendTypingNotification();
    repaint();
  }
  deSelect();

  return true; 
  // KeyPress was consumed and will not be passed any further to possibly registered KeyListeners

  // ToDo:
  // -interpret the ',' on the number pad as '.'. JUCE defines
  //  KeyPress::numberPadDecimalPoint
}

void RTextEntryField::focusGained(FocusChangeType cause)
{
  caretVisible = true;
  startTimer(blinkInterval);
  repaint();
}

void RTextEntryField::focusLost(FocusChangeType cause)
{
  caretVisible = false;
  stopTimer();
  sendTextChangeNotificationIfChanged();
  repaint();
  if( isCurrentlyModal() )  
    exitModalState(0);
}

void RTextEntryField::timerCallback()
{
  caretVisible = !caretVisible;
  startTimer(blinkInterval);
  repaint();
}

void RTextEntryField::mouseDown(const MouseEvent &e)
{  
  caretPosition = pixelPositionToCharacterIndex(e.x);
  deSelect();
}

void RTextEntryField::mouseDoubleClick(const MouseEvent &e)
{  
  caretPosition = pixelPositionToCharacterIndex(e.x);
  selectAll();
}

void RTextEntryField::mouseDrag(const MouseEvent &e)
{  
  int x1 = e.getMouseDownX();
  int x2 = x1 + e.getDistanceFromDragStartX();

  if( x1 > x2 )
    RAPT::rsSwap(x1, x2);

  selectionStart = pixelPositionToCharacterIndex(x1);
  selectionEnd   = pixelPositionToCharacterIndex(x2);
}

void RTextEntryField::paint(juce::Graphics &g)
{
  //const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  if( !hasKeyboardFocus(false) )
    RTextField::paint(g);  
  else
  {
    if( isSelectionEmpty() )
    {
      RTextField::paint(g);  
      if( caretVisible )
      {
        Point<int> p = getCaretPixelPosition();
        int caretWidth  = 2;
        int caretHeight = getHeight()-4;
        if( caretPosition < text.length() && replaceMode == true )
          caretWidth = font->getGlyphWidth(text[caretPosition]) + 1;
        g.fillRect(p.getX(), p.getY(), caretWidth, caretHeight);
      }
    }
    else
    {
      if( noBackgroundAndOutline == false )
        RWidget::paint(g);

      int x     = getTextPixelPositionX();
      int y     = getTextPixelPositionY();
      Colour tc = getTextColour();

      x      = drawBitmapFontText(g, x, y, getTextBeforeSelection(), font, tc);
      int xR = x + font->getTextPixelWidth(getTextInsideSelection(), font->getDefaultKerning());
      g.setColour(getHandleColour());
      g.fillRect(x, 2, xR-x, getHeight()-4);
      x      = drawBitmapFontText(g, x, y, getTextInsideSelection(), font, tc);
      x      = drawBitmapFontText(g, x, y, getTextAfterSelection(),  font, tc);
    }
  }

  if( textInvalid )
  {
    g.setColour(Colours::red.withAlpha(0.375f));
    Rectangle<int> r(outlineThickness, outlineThickness, getWidth()-2*outlineThickness, getHeight()-2*outlineThickness);
    g.fillRect(r);
  }
}

void RTextEntryField::sendTextChangeNotificationIfChanged()
{
  if( text != oldText )
  {

    // maybe factor this out into an updateAssignedParameter function:
    if( assignedParameter != NULL )
    {
      if( assignedParameter->isStringParameter() )
        assignedParameter->setStringValue(text, true, true);
      else
        assignedParameter->setValue(text.getDoubleValue(), true, true);
    }


    for(int i=0; i<textEntryFieldObservers.size(); i++)
      textEntryFieldObservers[i]->textChanged(this);
    oldText = text;
  }
}
    
void RTextEntryField::sendTypingNotification()
{
  for(int i=0; i<textEntryFieldObservers.size(); i++)
    textEntryFieldObservers[i]->somethingWasTypedIn(this);
}


//=========================================================================================================================================
// class RLabeledTextEntryField:

RLabeledTextEntryField::RLabeledTextEntryField(const juce::String& labelText, const juce::String& entryFieldText)
{
  setNoBackgroundAndOutline(true);

  labelField = new RTextField(labelText);
  labelField->setNoBackgroundAndOutline(true);
  addChildWidget(labelField);

  entryField = new RTextEntryField(entryFieldText);
  addChildWidget(entryField);

  labelWidth = 60;
}

RLabeledTextEntryField::~RLabeledTextEntryField()
{
  deleteAllChildren();  //necessary?
}

void RLabeledTextEntryField::setDescription(const juce::String &newDescription)
{
  labelField->setDescription(newDescription);
  entryField->setDescription(newDescription);
}

void RLabeledTextEntryField::setLabelWidth(int newWidth)
{
  labelWidth = newWidth;
  resized();
}

void RLabeledTextEntryField::setStateFromString(const juce::String &stateString, bool sendChangeMessage)
{
  entryField->setStateFromString(stateString, sendChangeMessage);
}


void RLabeledTextEntryField::resized()
{
  int w = getWidth();
  int h = getHeight();
  int entryWidth = w-labelWidth;
  labelField->setBounds(0,          0, labelWidth, h);
  entryField->setBounds(labelWidth, 0, entryWidth, h);
}



