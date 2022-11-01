//#include "rosof_AcidSequencerAudioModule.h"
//using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AcidSequencerAudioModule::AcidSequencerAudioModule(CriticalSection *newPlugInLock,
  rosic::AcidSequencer *acidSequencerToWrap) : AudioModule(newPlugInLock)
{
  jassert(acidSequencerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedAcidSequencer = acidSequencerToWrap;
  setModuleTypeName("AcidSequencer");
  createParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void AcidSequencerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedAcidSequencer == NULL )
    return;

  int    index = getIndexOfParameter(parameterThatHasChanged);
  double value = parameterThatHasChanged->getValue();
  switch( index )
  {
  case   0: wrappedAcidSequencer->setStepLength(     value);   break;
  case   1: wrappedAcidSequencer->setMode(      (int)value);   break;
  default:
    {
      // do nothing
    }
  } // end of switch( parameterIndex )

}

void AcidSequencerAudioModule::setStateFromXml(const XmlElement &xmlState,
                                               const juce::String &stateName, bool markAsClean)
{
  if( wrappedAcidSequencer != NULL )
  {
    // set permissibilities...

    for(int p=0; p<wrappedAcidSequencer->getNumPatterns(); p++ )
    {
      AcidPattern *pattern = wrappedAcidSequencer->getPattern(p);
      pattern->clear();
      XmlElement *xmlPattern = xmlState.getChildByName(juce::String("Pattern") + juce::String(p));
      if( xmlPattern != NULL )
      {
        pattern->setStepLength(xmlPattern->getDoubleAttribute("StepLength", 0.5));
        for(int s=0; s<pattern->getMaxNumSteps(); s++)
        {
          AcidNote   *note    = pattern->getNote(s);
          XmlElement *xmlStep = xmlPattern->getChildByName(juce::String("Step") + juce::String(s));
          if( xmlStep != NULL )
          {
            note->gate   = xmlStep->getBoolAttribute("Gate",   false);
            note->accent = xmlStep->getBoolAttribute("Accent", false);
            note->slide  = xmlStep->getBoolAttribute("Slide",  false);
            note->key    = xmlStep->getIntAttribute( "Key",    0);
            note->octave = xmlStep->getIntAttribute( "Octave", 0);
          }
        }
      }
    }
  }

  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
}

XmlElement* AcidSequencerAudioModule::getStateAsXml(const juce::String &stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedAcidSequencer != NULL )
  {
    // retrieve permissibilities...


    for(int p=0; p<wrappedAcidSequencer->getNumPatterns(); p++ )
    {
      AcidPattern *pattern = wrappedAcidSequencer->getPattern(p);
      if( !pattern->isEmpty() )
      {
        XmlElement *xmlPattern = new XmlElement(juce::String("Pattern") + juce::String(p));
        xmlPattern->setAttribute("StepLength", pattern->getStepLength());
        for(int s=0; s<pattern->getMaxNumSteps(); s++)
        {
          AcidNote *note = pattern->getNote(s);
          if( !note->isInDefaultState() )
          {
            XmlElement *xmlStep = new XmlElement(juce::String("Step") + juce::String(s));
            xmlStep->setAttribute("Gate",   note->gate);
            xmlStep->setAttribute("Accent", note->accent);
            xmlStep->setAttribute("Slide",  note->slide);
            xmlStep->setAttribute("Key",    note->key);
            xmlStep->setAttribute("Octave", note->octave);
            xmlPattern->addChildElement(xmlStep);
          }
        }
        xmlState->addChildElement(xmlPattern);
      }
    }
  }

  return xmlState;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void AcidSequencerAudioModule::createParameters()
{
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "StepLength", 0.25, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p);

  // #001:
  // ...this should *not* be automatable!
  p = new AutomatableParameter(lock, juce::String("Mode"), 0.0, 1.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String("Off"));
  p->addStringValue(juce::String("Key"));
  //p->addStringValue(juce::String(T("Host")));
  p->setValue(0.0, false, false);
  addObservedParameter(p);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================


AcidPatternEditor::AcidPatternEditor(rosic::AcidSequencer *sequencerToEdit)
  : jura::ColourSchemeComponent(juce::String("AcidPatternEditor"))
{
  patternToEdit         = nullptr;
  this->sequencerToEdit = sequencerToEdit;
  setDescription("Editor for acid patterns");

  rowHeight     = 12.f;
  columnWidth   = 20.f;
  keyLength     = 48.f;
  topLaneHeight = 16.f;
}

//-------------------------------------------------------------------------------------------------
// setup:

void AcidPatternEditor::setPatternToEdit(rosic::AcidPattern *newPatternToEdit)
{
  patternToEdit = newPatternToEdit;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int AcidPatternEditor::getStepAt(float x)
{
  if( x < keyLength )
    return -1;
  else
  {
    float x2 = x - keyLength;
    return floorInt(x2/columnWidth);
  }
}

int AcidPatternEditor::getNoteAt(float y)
{
  if( y < 4*topLaneHeight )
    return -1;
  else
  {
    float y2 = (float) getHeight() - y;
    return floorInt(y2/rowHeight);
  }
}

int AcidPatternEditor::getVirtualKeyAt(float x, float y)
{
  int key = getNoteAt(y);
  if( key != -1 && x <= keyLength )
    return key;
  else
    return -1;
}

bool AcidPatternEditor::isInGateRow(float y)
{
  if( y < topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInAccentRow(float y)
{
  if( y < 2*topLaneHeight && y >= topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInSlideRow(float y)
{
  if( y < 3*topLaneHeight && y >= 2*topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInOctaveRow(float y)
{
  if( y < 4*topLaneHeight && y >= 3*topLaneHeight )
    return true;
  else
    return false;
}

bool AcidPatternEditor::isInKeyboardColumn(float x)
{
  if( x < keyLength )
    return true;
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void AcidPatternEditor::mouseDown(const MouseEvent &e)
{
  if( patternToEdit == NULL )
    return;

  int step = getStepAt((float) e.x);
  int key  = getNoteAt((float) e.y);
  if( step != -1 && key != -1 )
  {
    if( patternToEdit->getKey(step) == key && patternToEdit->getGate(step) == true )
    {
      patternToEdit->setGate(step, false);
    }
    else
    {
      patternToEdit->setKey( step, key);
      patternToEdit->setGate(step, true);
    }
  }
  else if( isInGateRow((float) e.y) )
  {
    if( step != -1 )
      patternToEdit->setGate(step, !patternToEdit->getGate(step));
  }
  else if( isInAccentRow((float) e.y) )
  {
    if( step != -1 )
      patternToEdit->setAccent(step, !patternToEdit->getAccent(step));
  }
  else if( isInSlideRow((float) e.y) )
  {
    if( step != -1 )
      patternToEdit->setSlide(step, !patternToEdit->getSlide(step));
  }
  else if( isInOctaveRow((float) e.y) )
  {
    if( step != -1 )
    {
      if( e.mods.isLeftButtonDown() )
        patternToEdit->setOctave(step, patternToEdit->getOctave(step)-1);
      else if( e.mods.isRightButtonDown() )
        patternToEdit->setOctave(step, patternToEdit->getOctave(step)+1);
      else if( e.mods.isMiddleButtonDown() )
        patternToEdit->setOctave(step, 0);
    }
  }
  else if( key != -1 && isInKeyboardColumn((float) e.x) )
  {
    if( sequencerToEdit != NULL )
      sequencerToEdit->toggleKeyPermissibility(key);
  }

  repaint();
}

void AcidPatternEditor::paint(juce::Graphics &g)
{
  g.fillAll(getColorBackground());

  // Maybe use ints instead of floats - simplifies code below as well...but maybe there's a reason
  // we use floats, maybe to accomodate for GUI size settings for which things do not add up 
  // nicely?
  float x = 0.f;
  float y = 0.f;
  float w = (float) getWidth();
  float h = (float) getHeight();
  //float s = (float) rowHeight;
  float thickness = 2.f;


  const BitmapFont* font = &BitmapFontRoundedBoldA10D0::instance;
  drawBitmapFontText(g, (int)x+3, (int)y+3, "Gate:", font, getColorText());
  y += topLaneHeight;
  drawBitmapFontText(g, (int)x+3, (int)y+3, "Accent:", font, getColorText());
  y += topLaneHeight;
  drawBitmapFontText(g, (int)x+3, (int)y+3, "Slide:", font, getColorText());
  y += topLaneHeight;
  drawBitmapFontText(g, (int)x+3, (int)y+3, "Octave:", font, getColorText());

  w = keyLength;

  // draw background for white keys and black keys on top of it:
  float keyboardY = 4*topLaneHeight;


  Colour red = Colours::red; // for debugging

  // Draw background for keyboard. The color will eventually become the color of the white keys as
  // all other elements like black keys and key separators are drawn on top:
  //g.setColour(whiteKeyColour);
  //g.setColour(red);
  g.setColour(getColorWhiteKeys());
  g.fillRect(0.f, keyboardY, keyLength, (float)(13*rowHeight));

  // Draw background for sequencer. The color will eventually become the color for the white lanes:
  g.setColour(getColorWhiteLanes());
  g.fillRect(keyLength, keyboardY, (float)getWidth()-keyLength, (float)(13*rowHeight));

  // Draw black keys:
  //g.setColour(red);
  float bkw = 2*keyLength/3;         // width of black the keys
  g.setColour(getColorBlackKeys());
  x = 0;
  //w = 2*keyLength/3;    // use a variable bkw for black key width - is need later again
  h = (float) rowHeight;
  y = keyboardY + 11*rowHeight;
  g.fillRect(x, y, bkw, h); y -= 2*h;
  g.fillRect(x, y, bkw, h); y -= 3*h;
  g.fillRect(x, y, bkw, h); y -= 2*h;
  g.fillRect(x, y, bkw, h); y -= 2*h;
  g.fillRect(x, y, bkw, h);

  // Draw lanes for the black keys:
  x = keyLength;
  w = (float)getWidth()-keyLength;
  y = keyboardY + 11*rowHeight;
  g.setColour(getColorBlackLanes());
  //g.setColour(red);
  g.fillRect(x, y, w, h);  y -= 2*h;
  g.fillRect(x, y, w, h);  y -= 3*h;
  g.fillRect(x, y, w, h);  y -= 2*h;
  g.fillRect(x, y, w, h);  y -= 2*h;
  g.fillRect(x, y, w, h);

  // Draw horizontal lines for the 4 top-lanes:
  thickness = 2.f;
  x = 0.f;
  y = topLaneHeight;
  w = (float) getWidth();
  g.setColour(getColorLines());
  g.drawLine(x, y, w, y, thickness);
  y += topLaneHeight;
  g.drawLine(x, y, w, y, thickness);
  y += topLaneHeight;
  g.drawLine(x, y, w, y, thickness);
  y += topLaneHeight;
  g.drawLine(x, y, w, y, thickness);
  g.drawLine(keyLength, 0, keyLength, (float)getHeight(), thickness); // vertical

  // Draw the lines between the piano-roll rows:
  g.setColour(getColorLines());
  x = keyLength;
  h = (float) rowHeight;
  y = keyboardY + h;
  thickness = 1.f;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);  y += h;
  g.drawLine(x, y, w, y, thickness);
  // maybe use a loop

  // Draw the lines between the piano-roll white keys:
  g.setColour(getColorLines());
  //g.setColour(red);
  x = 0.f;
  w = keyLength;
  h = (float) rowHeight;
  y = keyboardY + 13*rowHeight - (h+h/2);
  float x2 = x+bkw, w2 = x+w;                       // for shortened lines due to black keys
  g.drawLine(x2, y, w2, y, thickness); y -= 2*h;
  g.drawLine(x2, y, w2, y, thickness); y -= h+h/2;
  g.drawLine(x,  y, w,  y, thickness); y -= h+h/2;
  g.drawLine(x2, y, w2, y, thickness); y -= 2*h;
  g.drawLine(x2, y, w2, y, thickness); y -= 2*h;
  g.drawLine(x2, y, w2, y, thickness); y -= h+h/2;
  g.drawLine(x,  y, w,  y, thickness);

  // Draw the vertical lines between the steps:
  int numSteps = rosic::AcidPattern::getMaxNumSteps();
  x = keyLength;
  y = 0.f;
  h = (float) getHeight();
  for(int i=0; i<numSteps; i++)
  {
    x += columnWidth;
    if( (i+1)%4 == 0 )
      thickness = 2.f;
    else
      thickness = 1.f;
    g.drawLine(x, y, x, h, thickness);
  }

  // Draw the pattern data:
  rosic::AcidPattern* ptn = patternToEdit;
  if( ptn != nullptr )
  {
    x            = keyLength;
    w            = columnWidth;
    h            = topLaneHeight;
    float dx     = columnWidth   / 2.f;
    float dy     = topLaneHeight / 2.f;
    int numSteps = ptn->getNumSteps();

    //g.setColour(getColorHandles());
    g.setColour(red);


    for(int i=0; i < ptn->getMaxNumSteps(); i++)
    {
      // Draw the circles that indicate presence of a step's feature (like gate, slide, accent):
      g.setColour(getColorHandles());
      bool slide = ptn->getSlide(i) && ptn->getGate((i+1)%numSteps);
      y = 0;
      if( ptn->getGate(i) == true )  {
        g.fillEllipse(x+5, y+3, w-10, h-6);
        if( slide )
          g.drawLine(x+dx, y+dy, x+3.f*dx, y+dy, 4.f); }
      y += topLaneHeight;
      if( ptn->getAccent(i) == true )
        g.fillEllipse(x+5, y+3, w-10, h-6);
      y += topLaneHeight;
      if( ptn->getSlide(i) == true )
        g.fillEllipse(x+5, y+3, w-10, h-6);

      // Draw the octave shift indicators:
      y += topLaneHeight;
      juce::String octString = valueToStringWithSign0( ptn->getOctave(i) );
      drawBitmapFontText(g, (int)(x+dx), (int)(y+dy), octString, font, getColorText(), 
        -1, Justification::centred);

      // Draw the events in the actual sequencer view, possibly with slide-connectors:
      y = keyboardY + 12*rowHeight;
      if( ptn->getGate(i) == true )
      {
        int key1 = ptn->getKey(i);
        y -= key1 * rowHeight;
        float w2 = (float) ptn->getStepLength();
        if( slide )
        {
          int key2 = ptn->getKey(i+1);
          float x2 = x + columnWidth;
          float y2 = keyboardY + 12*rowHeight - key2 * rowHeight;
          if(key1 > key2)
            g.drawLine(x2, y, x2, y2+rowHeight, 3.f);
          else
            g.drawLine(x2, y+rowHeight, x2, y2, 3.f);
          w2 = 1.f;
        }
        g.fillRect(x, y, w2*columnWidth, rowHeight);
      }
      x += columnWidth;
    }
  }

  // draw enclosing rectangle:
  g.drawRect(0, 0, getWidth(), getHeight(), 2);
}

//=================================================================================================
// class AcidSequencerModuleEditor:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AcidSequencerModuleEditor::AcidSequencerModuleEditor(CriticalSection *newPlugInLock, 
  AcidSequencerAudioModule* newAcidSequencerAudioModule) : AudioModuleEditor(newAcidSequencerAudioModule)
{
  setHeadlineStyle(SUB_HEADLINE);
  setLinkPosition(INVISIBLE);

  // assign the pointer to the rosic::AcidSequencer object to be used as aduio engine:
  jassert(newAcidSequencerAudioModule != NULL ); // you must pass a valid module here
  acidSequencerModuleToEdit = newAcidSequencerAudioModule;

  isTopLevelEditor = false;

  patternEditor = new AcidPatternEditor(acidSequencerModuleToEdit->wrappedAcidSequencer);
  addChildColourSchemeComponent(patternEditor);
  patternEditor->setPatternToEdit(acidSequencerModuleToEdit->wrappedAcidSequencer->getPattern(0));
  patternEditor->setDescriptionField(infoField, true);



  addWidget( modeLabel = new RTextField( juce::String("Mode:")) );
  modeLabel->setDescription("Chooses the sequencer mode");
  modeLabel->setDescriptionField(infoField);

  addWidget( modeBox = new RComboBox(juce::String("ModeComboBox")) );
  modeBox->assignParameter( moduleToEdit->getParameterByName("Mode") );
  modeBox->setDescription(modeLabel->getDescription());
  //modeBox->setNoBackgroundAndOutline(true);
  modeBox->setDescriptionField(infoField);

  addWidget( stepLengthSlider = new RSlider ("StepLengthSlider") );
  stepLengthSlider->assignParameter( acidSequencerModuleToEdit->getParameterByName("StepLength") );
  stepLengthSlider->setDescription(juce::String("Length of the steps in 16th notes"));
  stepLengthSlider->setStringConversionFunction(valueToString2);
  stepLengthSlider->setDescriptionField(infoField);
  stepLengthSlider->addListener(this);

  addWidget( shiftLabel = new RTextField( juce::String("Shift:")) );
  shiftLabel->setDescription("Shift the whole pattern for left or right (circularly)");
  shiftLabel->setDescriptionField(infoField);

  // Helper function to reduce boilerplate for button creation:
  auto addButton = [&](RClickButton** pButton, const String& name, const String& description)
  {
    addWidget( *pButton = new jura::RClickButton(name) );
    if(name == "L") (*pButton)->setSymbolIndex(jura::RButton::buttonSymbols::ARROW_LEFT);
    if(name == "R") (*pButton)->setSymbolIndex(jura::RButton::buttonSymbols::ARROW_RIGHT);

    (*pButton)->setDescription(description);
    (*pButton)->setDescriptionField(infoField);
    (*pButton)->setClickingTogglesState(false);
    (*pButton)->addRButtonListener(this);
  };

  addButton(&shiftLeftButton,         "L", "Shift the whole pattern one postion to the left (circularly)");
  addButton(&shiftRightButton,        "R", "Shift the whole pattern one postion to the right (circularly)");
  addButton(&shiftAccentsLeftButton,  "L", "Shift the accents one postion to the left (circularly)");
  addButton(&shiftAccentsRightButton, "R", "Shift the accents one postion to the right (circularly)");
  addButton(&shiftSlidesLeftButton,   "L", "Shift the slides one postion to the left (circularly)");
  addButton(&shiftSlidesRightButton,  "R", "Shift the slides one postion to the right (circularly)");
  addButton(&shiftNotesLeftButton,    "L", "Shift the notes one postion to the left (circularly)");
  addButton(&shiftNotesRightButton,   "R", "Shift the notes one postion to the right (circularly)");
  addButton(&shiftOctavesLeftButton,  "L", "Shift the octaves one postion to the left (circularly)");
  addButton(&shiftOctavesRightButton, "R", "Shift the octaves one postion to the right (circularly)");
  // reduce boilerplate further by making two functions addShiftLeft/RightButton that takes as 
  // string only "accents", "slides", etc.
  // todo: maybe instead of L and R use the left/right arrow symbols

  addButton(&reverseAllButton,    "Rev", "Reverses the whole pattern");
  addButton(&reverseAccentsButton,"Rev", "Reverses the accents");
  addButton(&reverseSlidesButton, "Rev", "Reverses the slides");
  addButton(&reverseNotesButton,  "Rev", "Reverses the notes");
  addButton(&reverseOctavesButton,"Rev", "Reverses the octaves");

  addButton(&invertAccentsButton, "Inv", "Inverts the accents");
  addButton(&invertSlidesButton,  "Inv", "Inverts the slides");
  addButton(&invertOctavesButton, "Inv", "Inverts the octaves");


  addButton(&swapAccentsSlidesButton, "A2S", "Swaps accents with slides");
  addButton(&xorAccentsSlidesButton,  "AXS", "Xors accents with slides");
  addButton(&xorSlidesAccentsButton,  "SXA", "Xors slides with accents");

  // maybe let new accents be old accenzts xor'ed with slides, smae for slides


  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// setup:



//-------------------------------------------------------------------------------------------------
// callbacks:

void AcidSequencerModuleEditor::rButtonClicked(RButton *b)
{
  if( acidSequencerModuleToEdit==NULL || acidSequencerModuleToEdit->wrappedAcidSequencer==NULL )
    return;

  auto seq = acidSequencerModuleToEdit->wrappedAcidSequencer;

  if(      b == shiftLeftButton         ) seq->circularShiftAll(-1);
  else if( b == shiftRightButton        ) seq->circularShiftAll(+1);
  else if( b == shiftAccentsLeftButton  ) seq->circularShiftAccents(-1);
  else if( b == shiftAccentsRightButton ) seq->circularShiftAccents(+1);
  else if( b == shiftSlidesLeftButton   ) seq->circularShiftSlides(-1);
  else if( b == shiftSlidesRightButton  ) seq->circularShiftSlides(+1);
  else if( b == shiftOctavesLeftButton  ) seq->circularShiftOctaves(-1);
  else if( b == shiftOctavesRightButton ) seq->circularShiftOctaves(+1);
  else if( b == shiftNotesLeftButton    ) seq->circularShiftNotes(-1);
  else if( b == shiftNotesRightButton   ) seq->circularShiftNotes(+1);

  else if( b == reverseAllButton     ) seq->reverseAll();
  else if( b == reverseAccentsButton ) seq->reverseAccents();
  else if( b == reverseSlidesButton  ) seq->reverseSlides();
  else if( b == reverseNotesButton   ) seq->reverseNotes();
  else if( b == reverseOctavesButton ) seq->reverseOctaves();

  else if( b == invertAccentsButton ) seq->invertAccents();
  else if( b == invertSlidesButton  ) seq->invertSlides();
  else if( b == invertOctavesButton ) seq->invertOctaves();

  else if( b == swapAccentsSlidesButton ) seq->swapAccentsWithSlides();
  else if( b == xorAccentsSlidesButton  ) seq->xorAccentsWithSlides();
  else if( b == xorSlidesAccentsButton  ) seq->xorSlidesWithAccents();



  // todo: reverse, swap-halves (swap 1st and 2nd half), exchange, for example, slide for accent,
  // ...all these features are easier with parallel arrays, maybe then, the length of the 
  // individual arrays can be different

  // \todo: randomization stuff....we should give the user the option to define a set of notes,
  // like C, E#, G, G#, Bb with certain probabilities and then a random number generator generates
  // a sequence.  it should also randomly set accents and slides...perhaps the probabilities
  // for on-beat and off-beat accents etc. can be set up



  // BUG: shifting the octaves sometimes produces garbage on the gui

  patternEditor->repaint();
}

void AcidSequencerModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  patternEditor->repaint();
}

void AcidSequencerModuleEditor::rSliderValueChanged(RSlider *rSliderThatHasChanged)
{
  if( rSliderThatHasChanged == stepLengthSlider )
    patternEditor->repaint();
}

void AcidSequencerModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);
}

void AcidSequencerModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  //int h = getHeight();

  y = getPresetSectionBottom()+4;
  x = 4;
  patternEditor->setBounds(x, y, 368, 220);
  x = patternEditor->getRight();
  y = patternEditor->getY();
  setRightKeepLeft(stateWidgetSet, x);
  int h = patternEditor->getTopLaneHeight();
  w = 28; 
  shiftLeftButton->setBounds( x,       y, w, h);
  shiftRightButton->setBounds(x+  w-2, y, w, h);
  reverseAllButton->setBounds(x+2*w-4, y, w, h);
  y += h;
  shiftAccentsLeftButton->setBounds( x,       y, w, h);
  shiftAccentsRightButton->setBounds(x+  w-2, y, w, h);
  reverseAccentsButton->setBounds(   x+2*w-4, y, w, h);
  invertAccentsButton->setBounds(    x+3*w-6, y, w, h); // test - if looks good, use for others, too
  y += h;
  shiftSlidesLeftButton->setBounds( x,       y, w, h);
  shiftSlidesRightButton->setBounds(x+  w-2, y, w, h);
  reverseSlidesButton->setBounds(   x+2*w-4, y, w, h);
  invertSlidesButton->setBounds(    x+3*w-6, y, w, h);
  y += h;
  shiftOctavesLeftButton->setBounds( x,       y, w, h);
  shiftOctavesRightButton->setBounds(x+  w-2, y, w, h);
  reverseOctavesButton->setBounds(   x+2*w-4, y, w, h);
  invertOctavesButton->setBounds(    x+3*w-6, y, w, h);
  y += h;
  shiftNotesLeftButton->setBounds( x,       y, w, h);
  shiftNotesRightButton->setBounds(x+  w-2, y, w, h);
  reverseNotesButton->setBounds(   x+2*w-4, y, w, h);

  x = patternEditor->getRight();
  y = stateWidgetSet->getY();
  modeLabel->setBounds(x, y, 40, 16);
  x = modeLabel->getRight();
  w = getWidth() - x;
  modeBox->setBounds(x, y, w-4, 16);

  y = shiftNotesRightButton->getBottom() + 8;
  x = patternEditor->getRight();
  w = 32;
  swapAccentsSlidesButton->setBounds(x+16,       y, w, 16);
  xorAccentsSlidesButton->setBounds( x+16+  w-2, y, w, 16);
  xorSlidesAccentsButton->setBounds( x+16+2*w-4, y, w, 16);

  y += 24;
  w  = getWidth() - x;
  stepLengthSlider->setBounds(x+4, y+4, w-8, 16);
}


/*
-make step-length available for modulation
-add selector for pattern (maybe a 4x4 array)
-implement copy/paste for patterns
-fix colors (done)
-fix manipulator button positioning
-maybe use triangles instead of L/R letters
-animate the sequencer - highlight the column where we currently are or let a cursor step through. Take other
 animated widgets as reference, such as level-metering widgets


Ideas for sequence manipulations:
-Generally, we want operations that are their own inverses to be able to undo them easily as long 
 as we don't have a proper "Undo" functionality
-Permutations:
 -Swap first and second halves recursively ...well, doing it recursively just results in reversing 
  the array, so that's nothing new. How about non-recursively or with limited recursion depth?
 -De/interleave via viewing the pattern of length 16 as 2 patterns of length 8
-Unary functions: logical not - we call it "Inv" on teh buttons
-Combinations: XOR (done), NXOR, swap accents/slides
 



*/