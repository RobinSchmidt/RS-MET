//-------------------------------------------------------------------------------------------------
// construction/destruction:

BreakpointModulatorEditorMulti::BreakpointModulatorEditorMulti(CriticalSection *newPlugInLock,
  BreakpointModulatorAudioModule* newBreakpointModulatorAudioModule)
 : AudioModuleEditor(newBreakpointModulatorAudioModule)
 , BreakpointModulatorEditor(newPlugInLock, newBreakpointModulatorAudioModule)
{
  // init the pointer to the modulator to be edited to NULL:
  modulatorToEdit = NULL;
  jassert(newBreakpointModulatorAudioModule != NULL);

  editedModulatorIndex = -1;

  // change the headline from the default "Sub-Editor" to "Modulator-Editor":
  //headline->setText(juce::String(T("Modulator")), false);
  setHeadlineText("Modulators");
  setDescription("This is an editor for multi-breakpoint modulation generators");

  // create the breakpoint-editor for multiple modulators, make the inherited single-modulator
  // editor invisible and re-assign the zoomer:
  breakpointEditor->setVisible(false);
  breakpointEditorMulti = new ModulatorCurveEditorMulti("PlotEditor");
  breakpointEditorMulti->addChangeListener(this);
  addPlot(breakpointEditorMulti);
  breakpointZoomer->setCoordinateSystem(breakpointEditorMulti);
  breakpointZoomer->toFront(false);

  // add the inherited widgets for the single modulator as first elements to our arrays:
  globalEditors.add(globalEditor);
  globalEditor->timeScaleByKeySlider->setSliderName("K");
  globalEditor->timeScaleByVelSlider->setSliderName("V");
  globalEditor->depthByKeySlider->setSliderName(    "K");
  globalEditor->depthByVelSlider->setSliderName(    "V");
  globalEditor->setLayout(1);
  globalEditor->editButton->addRButtonListener(this);
  //globalEditor->editButton->setRadioGroupId(1);
  newBreakpointModulatorAudioModule->addStateWatcher(globalEditor->stateWidgetSet);
  globalEditor->stateWidgetSet->addChangeListener(this);

  // customize the inherited stateWidgetSet:
  stateWidgetSet->setVisible(false);

  addModulatorToEdit(newBreakpointModulatorAudioModule, false);
     // this will also set up the widgets according to the state of the modulator
}

/*
BreakpointModulatorEditorMulti::~BreakpointModulatorEditorMulti()
{

}
*/

//-------------------------------------------------------------------------------------------------
// parameter-settings:

void BreakpointModulatorEditorMulti::addModulatorToEdit(
  BreakpointModulatorAudioModule* newModulatorToEdit, bool createWidgets)
{
  modulatorToEdit = newModulatorToEdit->wrappedBreakpointModulator;

  modulatorModules.getLock().enter();
  modulatorModules.addIfNotAlreadyThere(newModulatorToEdit);
  if( createWidgets == true )
    createWidgetsForNewModulator(newModulatorToEdit);
  modulatorModules.getLock().exit();

  breakpointEditorMulti->addModulatorToEdit(modulatorToEdit);
  breakpointZoomer->zoomToAllXY();
}

void BreakpointModulatorEditorMulti::selectModulatorToEdit(int index)
{
  modulatorModules.getLock().enter();
  jassert( index < modulatorModules.size() ); // index out of range
  index %= modulatorModules.size();
  editedModulatorIndex = index;
  modulatorModule      = modulatorModules[index];
  modulatorToEdit      = modulatorModules[index]->wrappedBreakpointModulator;
  breakpointParameterEditor->setModulatorToEdit(modulatorModules[index]);
  modulatorModules.getLock().exit();

  globalEditors.getLock().enter();
  jassert(index < globalEditors.size()); // index out of range

  for(int i=0; i<globalEditors.size(); i++)
    globalEditors[i]->editButton->setToggleState(false, false);
  globalEditors[index]->editButton->setToggleState(true, false);

  breakpointEditorMulti->selectModulatorToEdit(index);


  breakpointParameterEditor->copyColourSettingsFrom(globalEditors[index]);
  breakpointZoomer->copyColourSettingsFrom(globalEditors[index]);
  breakpointEditorMulti->setColourScheme(globalEditors[index]->getPlotColourScheme());

  editorColourScheme = globalEditors[index]->editorColourScheme;
  widgetColourScheme = globalEditors[index]->widgetColourScheme;
  plotColourScheme   = globalEditors[index]->plotColourScheme;
  for(int i=0; i<widgets.size(); i++)
    widgets[i]->setColourScheme(widgetColourScheme);

  // new:
  assignGridAndSnapWidgets(modulatorModule);
  updatePlotAndGridWidgets(modulatorModule, breakpointEditorMulti);

  repaint();

  breakpointZoomer->updateScrollbars();
  //breakpointZoomer->repaint();
  globalEditors.getLock().exit();
}

/*
void BreakpointModulatorEditorMulti::setDescriptionField(RLabel* newDescriptionField)
{
  BreakpointModulatorEditor::setDescriptionField(newDescriptionField);
  breakpointEditorMulti->setDescriptionField(newDescriptionField);
  globalEditors.getLock().enter();
  for(int i= 0; i<globalEditors.size(); i++)
  {
    globalEditors[i]->setDescriptionField(newDescriptionField);
  }
  globalEditors.getLock().exit();
}
*/

void BreakpointModulatorEditorMulti::setModulatorLabel(int index, const juce::String& newLabel)
{
  globalEditors.getLock().enter();
  jassert(index >= 0 && index < globalEditors.size());
  if(index >= 0 && index < globalEditors.size())
  {
    globalEditors[index]->stateWidgetSet->stateLabel->setText(newLabel);
    //globalEditors[index]->setHeadlineText(newLabel);
  }
  globalEditors.getLock().exit();
  /*
  stateWidgetSets.getLock().enter();
  jassert( index < stateWidgetSets.size() ); // index out of range
  if( index  < stateWidgetSets.size() )
    stateWidgetSets[index]->stateLabel->setText(newLabel, false);
  stateWidgetSets.getLock().exit();
  */
}

void BreakpointModulatorEditorMulti::setChildColourScheme(int index,
  const EditorColourScheme& newEditorColourScheme, const WidgetColourScheme& newWidgetColourScheme)
{
  //int dummy = 0;

  /*
  globalEditors.getLock().enter();
  jassert( index < globalEditors.size() ); // index out of range
  if( index  < globalEditors.size() )
  {
    globalEditors[index]->setColourScheme(newEditorColourScheme, true);
    globalEditors[index]->setWidgetColourScheme(newWidgetColourScheme, true);
  }
  globalEditors.getLock().exit();
  */

  // \ todo: call setHue
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void BreakpointModulatorEditorMulti::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( moduleToEdit == NULL )
    return;

  if( buttonThatWasClicked == breakpointParameterEditor->shapeToAllButton )
  {
    if( breakpointParameterEditor->shapeToAllButton->getToggleState() == true )
      breakpointEditorMulti->updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
    moduleToEdit->markStateAsDirty();
    return;
  }
  else if( buttonThatWasClicked == snapXButton || buttonThatWasClicked == snapYButton )
  {
    BreakpointModulatorEditor::rButtonClicked(buttonThatWasClicked);
    //updatePlotAndGridWidgets(modulatorModule, breakpointEditorMulti);
    //breakpointEditorMulti->setSnapToFineGridX(snapXButton->getToggleState());
    //breakpointEditorMulti->setVerticalFineGridVisible(snapXButton->getToggleState());
  }
  /*
  else if( buttonThatWasClicked == snapYButton )
  {
    BreakpointModulatorEditor::rButtonClicked(buttonThatWasClicked);
    //updatePlotAndGridWidgets(modulatorModule, breakpointEditorMulti);
    //breakpointEditorMulti->setSnapToFineGridY(snapYButton->getToggleState());
    //breakpointEditorMulti->setHorizontalFineGridVisible(snapYButton->getToggleState());
  }
  */


  globalEditors.getLock().enter();
  for(int i=0; i<globalEditors.size(); i++)
  {
    if( buttonThatWasClicked == globalEditors[i]->editButton )
    {
      selectModulatorToEdit(i);
      globalEditors.getLock().exit();
      //breakpointEditorMulti->updatePlotImage();
      return;
    }
    else if( buttonThatWasClicked == globalEditors[i]->loopButton )
    {
      // maybe this (ugly) nested if should be in wrapped into a function...
      modulatorModules.getLock().enter();
      if( i < modulatorModules.size() )
      {
        if( modulatorModules[i] != NULL )
        {
          if( modulatorModules[i]->wrappedBreakpointModulator != NULL )
          {
            modulatorModules[i]->wrappedBreakpointModulator->setLoopMode(
              globalEditors[i]->loopButton->getToggleState() );
            modulatorModules[i]->markStateAsDirty();
          }
          else
            jassertfalse;
        }
        else
          jassertfalse;
      }
      else
        jassertfalse;
      modulatorModules.getLock().exit();

      globalEditors.getLock().exit();
      breakpointEditorMulti->updatePlotImage();
      return;
    }
  }
  globalEditors.getLock().exit();


  // let the baseclass handle clicks on the snapX/y buttons:
  //BreakpointModulatorEditor::rButtonClicked(buttonThatWasClicked);
}

void BreakpointModulatorEditorMulti::changeListenerCallback(
  ChangeBroadcaster *objectThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;

  globalEditors.getLock().enter();
  for(int i=0; i<globalEditors.size(); i++)
  {
    if( objectThatHasChanged == globalEditors[i]->stateWidgetSet )
    {
      updateWidgetsAccordingToState(true);
      globalEditors.getLock().exit();
      return;
    }
  }
  globalEditors.getLock().exit();

  /*
  if( objectThatHasChanged == stateWidgetSet )
  {
    // we must check for any one of the stateWidget-sets....

    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
    return;
  }
  */

  if( objectThatHasChanged == breakpointEditorMulti )
  {
    breakpointParameterEditor->selectBreakpoint(
      breakpointEditorMulti->getSelectedBreakpointIndex());
    breakpointZoomer->updateScrollbars();
  }
  else if( objectThatHasChanged == breakpointParameterEditor )
    breakpointEditorMulti->updatePlotCurveData();

  moduleToEdit->markStateAsDirty();
}

void BreakpointModulatorEditorMulti::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;

  if( rComboBoxThatHasChanged == breakpointParameterEditor->shapeComboBox )
    breakpointEditorMulti->updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);

  // !!!NEEDS UPDATE!!!
  else if( rComboBoxThatHasChanged == gridXComboBox )
  {
    int newGridIntervalIndex = gridXComboBox->getSelectedItemIdentifier();
    breakpointEditorMulti->setVerticalFineGrid(gridIntervalFromIndex(newGridIntervalIndex),
      snapXButton->getToggleState());
    breakpointEditorMulti->repaint();
  }
  else if( rComboBoxThatHasChanged == gridYComboBox )
  {
    int newGridIntervalIndex = gridYComboBox->getSelectedItemIdentifier();
    breakpointEditorMulti->setHorizontalFineGrid(gridIntervalFromIndex(newGridIntervalIndex),
      snapYButton->getToggleState());
    breakpointEditorMulti->repaint();
  }
}

void BreakpointModulatorEditorMulti::copyColourSettingsFrom(
  const ColourSchemeComponent *componentToCopyFrom)
{
  AudioModuleEditor::copyColourSettingsFrom(componentToCopyFrom);

  globalEditors.getLock().enter();
  int i;
  breakpointEditorMulti->clearCurveColours();
  ColourAHSL tmpColour = breakpointEditorMulti->getPlotColourScheme().curvesAHSL;
  for(i=0; i<globalEditors.size(); i++)
  {
    float hue = editorColourScheme.getCentralHue() + editorColourScheme.getHueOffset(i);
    globalEditors[i]->setCentralHue(hue);
    breakpointEditorMulti->appendCurveColour(tmpColour.withHue(hue));
  }
  globalEditors.getLock().exit();

  selectModulatorToEdit(editedModulatorIndex);
}

void BreakpointModulatorEditorMulti::paint(Graphics &g)
{
  Editor::paint(g);

  fillRectWithBilinearGradient(g, snapRectangle, editorColourScheme.topLeft,
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);

  Editor::drawHeadline(g);

  g.setColour(editorColourScheme.outline);
  g.drawRect(0, 0, getWidth(), getHeight());
}

void BreakpointModulatorEditorMulti::resized()
{
  presetSectionPosition = INVISIBLE;
  linkPosition          = INVISIBLE;
  AudioModuleEditor::resized();

  int rightSectionWidth = 100;
  int leftSectionWidth  = 282;

  int x = getWidth()-rightSectionWidth;
  int y = 0;
  int w = getWidth()-x;
  int h = getHeight()-y-44;

  // the right section:
  breakpointGroupRectangle.setBounds(x, y, w, h);
  breakpointParameterEditor->setBounds(x-2, y, w+2, h);
  snapRectangle.setBounds(x, breakpointGroupRectangle.getBottom(), w,
    getHeight() - breakpointGroupRectangle.getBottom() );

  x = snapRectangle.getX();
  y = snapRectangle.getY();
  snapXButton->setBounds(x+4, y+4, 32, 16);
  x = snapXButton->getRight();
  gridXComboBox->setBounds(x+4, y+4, snapRectangle.getRight()-x-8, 16);
  x = snapRectangle.getX();
  y = snapXButton->getBottom();
  snapYButton->setBounds(x+4, y+4, 32, 16);
  x = snapXButton->getRight();
  gridYComboBox->setBounds(x+4, y+4, snapRectangle.getRight()-x-8, 16);

  // the middle section (the curve editor):
  x = leftSectionWidth;
  y = 0;
  w = breakpointGroupRectangle.getX()-x+2;
  h = getHeight()-y+2;
  breakpointEditorMulti->setBounds(x, y, w-breakpointZoomer->getZoomerSize(),
    h-breakpointZoomer->getZoomerSize());
  breakpointZoomer->alignWidgetsToCoordinateSystem();

  // the left section:
  globalEditors.getLock().enter();
  int numModulators = globalEditors.size();
  leftSectionRectangles.clear();
  //x = 0;
  x = 2;
  //w = leftSectionWidth+2;
  w = leftSectionWidth-2;
  h = 78;
  //y = getHeight() - h;
  y = getHeight() - h - 2;
  for(int i=numModulators-1; i>=0; i--)
  {
    juce::Rectangle<int> *r = new juce::Rectangle<int>(x, y, w, h);
    leftSectionRectangles.add(r);
    globalEditors[i]->setBounds(*r);
    //y -= h-2;
    y -= h;

  }

  globalEditors.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// others:

void BreakpointModulatorEditorMulti::deSelectBreakpoint()
{
  breakpointEditorMulti->setSelectedBreakpointIndex(-1);
  BreakpointModulatorEditor::deSelectBreakpoint();
}

void BreakpointModulatorEditorMulti::updateWidgetsAccordingToState(bool deSelectBreakpointBefore)
{
  if( modulatorToEdit == NULL )
    return;

  if( deSelectBreakpointBefore == true )
    deSelectBreakpoint();

  globalEditors.getLock().enter();
  for(int m=0; m<globalEditors.size(); m++)
    globalEditors[m]->updateWidgetsAccordingToState();
  globalEditors.getLock().exit();


  // update the plot:
  updatePlotAndGridWidgets(modulatorModule, breakpointEditorMulti);  // new
  breakpointEditorMulti->updateMaximumRange(true);
  breakpointEditorMulti->updateCurveDataForAllPlots(true, true);
  breakpointZoomer->zoomToAllXY();
}

void BreakpointModulatorEditorMulti::updateWidgetsAccordingToState()
{
  updateWidgetsAccordingToState(true);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void BreakpointModulatorEditorMulti::createWidgetsForNewModulator(
  BreakpointModulatorAudioModule *newModulator)
{
  BreakpointModulatorGlobalEditor* tmpEditor =
    new BreakpointModulatorGlobalEditor(lock, newModulator);

  //addWidgetSet(tmpWidgetSet); // this makes the encapsulated widgets child-components and
                                // registers this object as listener to the widgets
  //automatableSliders.addIfNotAlreadyThere(tmpWidgetSet->timeScaleSlider);
  //automatableSliders.addIfNotAlreadyThere(tmpWidgetSet->depthSlider);
  tmpEditor->timeScaleByKeySlider->setSliderName("K");
  tmpEditor->timeScaleByVelSlider->setSliderName("V");
  tmpEditor->depthByKeySlider->setSliderName(    "K");
  tmpEditor->depthByVelSlider->setSliderName(    "V");
  tmpEditor->setLayout(1);
  tmpEditor->loopButton->addRButtonListener(this);
  tmpEditor->editButton->addRButtonListener(this);

  newModulator->addStateWatcher(tmpEditor->stateWidgetSet);
  tmpEditor->stateWidgetSet->addChangeListener(this);
  //tmpEditor->stateWidgetSet->setDescriptionField(infoField);

  //tmpEditor->editButton->setRadioGroupId(1);
  globalEditors.add(tmpEditor);
  //addAndMakeVisible(tmpEditor);
  addChildEditor(tmpEditor);

  /*
  StateLoadSaveWidgetSet* tmpStateWidgetSet = new StateLoadSaveWidgetSet();
  addChildColourSchemeComponent( tmpStateWidgetSet );
  tmpStateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  newModulator->addStateWatcher(tmpStateWidgetSet);
  tmpStateWidgetSet->setDescriptionField(infoField);
  tmpStateWidgetSet->stateLabel->setText(juce::String(T("Preset")), false);
  tmpStateWidgetSet->addChangeListener(this);
  stateWidgetSets.add(tmpStateWidgetSet);
  */
}

void BreakpointModulatorEditorMulti::autoAdjustPlotRangeX()
{
  double maxTime = modulatorToEdit->getEndTime();
  breakpointEditorMulti->setMaximumRangeMaxX(maxTime);
  breakpointEditorMulti->setCurrentRangeMaxX(maxTime);
  breakpointZoomer->zoomToAllX();
}

void BreakpointModulatorEditorMulti::autoAdjustPlotRangeY()
{
  double minLevel = modulatorToEdit->getMinLevel();
  double maxLevel = modulatorToEdit->getMaxLevel();
  breakpointEditorMulti->setMaximumRangeY(minLevel, maxLevel);
  breakpointEditorMulti->setCurrentRangeY(minLevel, maxLevel);
  breakpointZoomer->zoomToAllY();
}
