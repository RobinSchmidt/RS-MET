#include "rosof_ModulatorCurveEditor.h"
using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulatorCurveEditor::ModulatorCurveEditor(const juce::String& name) 
: CurveFamilyPlotOld(name), InteractiveCoordinateSystemOld(name)
{
  setDescription(T("Left-click: insert, right-click: remove, shift-drag: time-shifts subsequent breakpoints"));

  modulatorToEdit            = NULL;
  selectedBreakpoint         = -1;
  breakpointBeingDragged     = -1;
  locatorBeingDragged        = -1;
  //graphColour                = Colours::blueviolet;
  //graphColour                = Colour(0.7f, 1.0f, 0.75f, (uint8) 255);
  //loopIsOn                   = false;
  mostRecentBreakpointIndex  = -1;
  mostRecentAction           = NO_ACTION;

  // set up the plotting options:  
  setAutoReRendering(false);
  setMaximumRange(-0.0, 3.0, 0.0, 1.0);
  setHorizontalCoarseGrid(1.0, false);
  setHorizontalFineGrid(0.25,  false);
  setAxisValuesPositionX(CoordinateSystemOld::ABOVE_AXIS);
  setAxisLabelX(juce::String::empty);
  setVerticalCoarseGrid(1.0,   false);
  setVerticalFineGrid(0.25,    false);
  setAxisValuesPositionY(CoordinateSystemOld::RIGHT_TO_AXIS);
  setAxisLabelY(juce::String::empty);
  setStringConversionForAxisX(&valueToString0);
  setStringConversionForAxisY(&valueToString0);
  setAutoReRendering(true);

  // fill the arrays with the possible snap-values:
  xSnapIntervals.add(0.0);
  xSnapIntervals.add(1.0/128.0);
  xSnapIntervals.add(1.0/96.0);
  xSnapIntervals.add(1.0/64.0);
  xSnapIntervals.add(1.0/48.0);
  xSnapIntervals.add(1.0/32.0);
  xSnapIntervals.add(1.0/24.0);
  xSnapIntervals.add(1.0/16.0);
  xSnapIntervals.add(1.0/12.0);
  xSnapIntervals.add(1.0/10.0);
  xSnapIntervals.add(1.0/8.0);
  xSnapIntervals.add(1.0/6.0);
  xSnapIntervals.add(1.0/4.0);

  ySnapIntervals.add(0.0);
  ySnapIntervals.add(1.0/128.0);
  ySnapIntervals.add(1.0/96.0);
  ySnapIntervals.add(1.0/64.0);
  ySnapIntervals.add(1.0/48.0);
  ySnapIntervals.add(1.0/32.0);
  ySnapIntervals.add(1.0/24.0);
  ySnapIntervals.add(1.0/16.0);
  ySnapIntervals.add(1.0/12.0);
  ySnapIntervals.add(1.0/10.0);
  ySnapIntervals.add(1.0/8.0);
  ySnapIntervals.add(1.0/6.0);
  ySnapIntervals.add(1.0/4.0);

  // initialize the plot data:
  numModulators        = 1;
  numSamplesInPlot     = 0;
  //plotDataX            = NULL;
  plotDataFlatX        = NULL;
  plotDataX            = NULL;  
  plotDataFlatY        = NULL;
  plotDataY            = NULL;   
  editedModulatorIndex = 0;
  setCurveFamilyValues(numSamplesInPlot, numModulators, plotDataX, plotDataY);
  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, false);
}

ModulatorCurveEditor::~ModulatorCurveEditor()
{
  freePlotBufferArrays();
}

//-------------------------------------------------------------------------------------------------
// setup:

void ModulatorCurveEditor::setModulatorToEdit(rosic::BreakpointModulator* newModulatorToEdit)
{
  modulatorToEdit = newModulatorToEdit;
  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
}

/*
void ModulatorCurveEditor::setLoopMode(bool shouldBeLooped)
{ 
  if( modulatorToEdit == NULL )
    return;

  modulatorToEdit->setLoopMode(shouldBeLooped);
  //loopIsOn = shouldBeLooped; 
  repaint();
}
*/

/*
void ModulatorCurveEditor::setSyncMode(bool shouldBeSynced)
{ 
  if( modulatorToEdit == NULL )
    return;

  modulatorToEdit->setSyncMode(shouldBeSynced);
  updateModulationCurve();
  repaint();
}
*/

bool ModulatorCurveEditor::setSelectedBreakpointIndex(int indexToActivate)
{
  if( modulatorToEdit == NULL )
    return false;

  jassert( indexToActivate >= -1 && indexToActivate <= modulatorToEdit->lastBreakpointIndex() );

  if( indexToActivate >= 0 && indexToActivate <= modulatorToEdit->lastBreakpointIndex() )
  {
    selectedBreakpoint = indexToActivate;
    return true;
  }
  else 
  {
    selectedBreakpoint = -1;
    return false;
  }
}

bool ModulatorCurveEditor::setSelectedBreakpointTime(double newTime, bool broadcastMessage)
{
  if( modulatorToEdit == NULL )
    return false;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    bool success = modulatorToEdit->setBreakpointTime(selectedBreakpoint, newTime);
    updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);

    if( broadcastMessage == true ) 
      sendChangeMessage(); 

    // maybe the maximum meaningful plot-range has changed, so we need to
    // inform the zoomer (if there is any)
    return success;
  }
  else
  {
    //jassertfalse
    return false;
  }
}

bool ModulatorCurveEditor::setSelectedBreakpointLevel(double newLevel, bool broadcastMessage)
{
  if( modulatorToEdit == NULL )
    return 0.0;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    bool success = modulatorToEdit->setBreakpointLevel(selectedBreakpoint, newLevel);
    updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);

    if( broadcastMessage == true ) 
      sendChangeMessage(); 

    // maybe the maximum meaningful plot-range has changed, so we need to
    // inform the zoomer (if there is any)
    return success;
  }
  else
  {
    //jassertfalse
    return false;
  }
}

bool ModulatorCurveEditor::setSelectedBreakpointShape(int newShape, bool broadcastMessage)
{
  if( modulatorToEdit == NULL )
    return false;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    bool success = modulatorToEdit->setBreakpointShape(selectedBreakpoint, newShape);

    if( broadcastMessage == true ) 
      sendChangeMessage(); 

    updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
    return success;
  }
  else
  {
    //jassertfalse
    return false;
  }
}

void ModulatorCurveEditor::setAllBreakpointShapes(int newShape, bool broadcastMessage)
{
  if( modulatorToEdit == NULL )
    return;

  for(int p=0; p<=modulatorToEdit->lastBreakpointIndex(); p++)
    modulatorToEdit->setBreakpointShape(p, newShape);
    
  if( broadcastMessage == true )   
    sendChangeMessage(); 

  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
}

bool ModulatorCurveEditor::setSelectedBreakpointShapeAmount(double newShapeAmount, bool broadcastMessage)
{
  if( modulatorToEdit == NULL )
    return false;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    bool success = modulatorToEdit->setBreakpointShapeAmount(
      selectedBreakpoint, newShapeAmount);

    if( broadcastMessage == true ) 
      sendChangeMessage(); 

    updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
    return success;
  }
  else
  {
    //jassertfalse
    return false;
  }
}

void ModulatorCurveEditor::setAllBreakpointShapeAmounts(double newShapeAmount, bool broadcastMessage)
{
  if( modulatorToEdit == NULL )
    return;

  for(int p=0; p<=modulatorToEdit->lastBreakpointIndex(); p++)
    modulatorToEdit->setBreakpointShapeAmount(p, newShapeAmount);

  if( broadcastMessage == true )
    sendChangeMessage(); 

  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
}

void ModulatorCurveEditor::setCurrentRangeX(double newMinX, double newMaxX)
{
  setAutoReRendering(false);
  CurveFamilyPlotOld::setCurrentRangeX(newMinX, newMaxX);
  CurveFamilyPlotOld::updateBackgroundImage();
  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
  setAutoReRendering(true);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int ModulatorCurveEditor::getNumBreakpoints() const
{
  if( modulatorToEdit == NULL )
    return 0;
  else
    return modulatorToEdit->lastBreakpointIndex()+1;
}

/*
CoordinateSystemRange ModulatorCurveEditor::getMaximumMeaningfulRange(
  double relativeMarginLeft, double relativeMarginRight,
  double relativeMarginTop,  double relativeMarginBottom)
{
  if( modulatorToEdit == NULL )
    return CoordinateSystemRange();

  // require at least 1.0% margin:
  jassert( relativeMarginLeft   >= 1.0 );
  jassert( relativeMarginRight  >= 1.0 );
  jassert( relativeMarginTop    >= 1.0 );
  jassert( relativeMarginBottom >= 1.0 );
  if( relativeMarginLeft   < 1.0 )
    relativeMarginLeft   = 1.0;
  if( relativeMarginRight  < 1.0 )
    relativeMarginRight  = 1.0;
  if( relativeMarginTop    < 1.0 )
    relativeMarginTop    = 1.0;
  if( relativeMarginBottom < 1.0 )
    relativeMarginBottom = 1.0;

  CoordinateSystemRange r;

  double startTime = modulatorToEdit->getStartTime();
  double endTime   = modulatorToEdit->getEndTime();
  double width     = endTime-startTime;
  //double minLevel  = modulatorToEdit->getMinLevel();
  //double maxLevel  = modulatorToEdit->getMaxLevel();
  double minLevel = maximumRange.getMinY();
  double maxLevel = maximumRange.getMaxY();
  double height    = maxLevel-minLevel;

  // a somewhat kludgy solution for degenerate modulation curves:
  if( minLevel == maxLevel )
  {
    minLevel -= 0.01;
    maxLevel += 0.01;
  }

  r.setMinX(startTime - 0.01*relativeMarginLeft   * width);
  r.setMaxX(endTime   + 0.01*relativeMarginRight  * width);
  r.setMinY(minLevel  - 0.01*relativeMarginBottom * height);
  r.setMaxY(maxLevel  + 0.01*relativeMarginTop    * height);

  return r;
}
*/

bool ModulatorCurveEditor::getLoopMode()
{ 
  if( modulatorToEdit == NULL )
    return false;
  else
    return (modulatorToEdit->getLoopMode() != 0);
}

/*
bool ModulatorCurveEditor::getSyncMode()
{ 
  if( modulatorToEdit == NULL )
    return false;
  else
    return (modulatorToEdit->isInSyncMode() != 0);
}
*/

double ModulatorCurveEditor::getBreakpointTime(int index)
{
  if( modulatorToEdit == NULL )
    return 0.0;
  else
    return modulatorToEdit->getBreakpointTime(index);
}

double ModulatorCurveEditor::getBreakpointLevel(int index)
{
  if( modulatorToEdit == NULL )
    return 0.0;
  else
    return modulatorToEdit->getBreakpointLevel(index);
}

int ModulatorCurveEditor::getBreakpointShape(int index)
{
  if( modulatorToEdit == NULL )
    return 0;
  else
    return modulatorToEdit->getBreakpointShape(index);
}

double ModulatorCurveEditor::getBreakpointShapeAmount(int index)
{
  if( modulatorToEdit == NULL )
    return 0.0;
  else
    return modulatorToEdit->getBreakpointShapeAmount(index);
}

int ModulatorCurveEditor::getSelectedBreakpointIndex()
{
  return selectedBreakpoint;
}

double ModulatorCurveEditor::getSelectedBreakpointTime()
{
  if( modulatorToEdit == NULL )
    return 0.0;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    return modulatorToEdit->getBreakpointTime(selectedBreakpoint);
  }
  else
  {
    //jassertfalse
    return 0.0;
  }
}

double ModulatorCurveEditor::getSelectedBreakpointMinTime()
{
  if( modulatorToEdit == NULL )
    return 0.0;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    return modulatorToEdit->getBreakpointMinTime(selectedBreakpoint);
  }
  else
  {
    //jassertfalse
    return 0.0;
  }
}

double ModulatorCurveEditor::getSelectedBreakpointMaxTime()
{
  if( modulatorToEdit == NULL )
    return 0.0;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    return modulatorToEdit->getBreakpointMaxTime(selectedBreakpoint);
  }
  else
  {
    //jassertfalse
    return 0.0;
  }
}

double ModulatorCurveEditor::getSelectedBreakpointLevel()
{
  if( modulatorToEdit == NULL )
    return 0.0;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    return modulatorToEdit->getBreakpointLevel(selectedBreakpoint);
  }
  else
  {
    //jassertfalse
    return 0.0;
  }
}

int ModulatorCurveEditor::getSelectedBreakpointShape()
{
  if( modulatorToEdit == NULL )
    return 0;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    return modulatorToEdit->getBreakpointShape(selectedBreakpoint);
  }
  else
  {
    //jassertfalse
    return 0;
  }
}

double ModulatorCurveEditor::getSelectedBreakpointShapeAmount()
{
  if( modulatorToEdit == NULL )
    return 0.0;

  if( selectedBreakpoint >= 0 &&
    selectedBreakpoint <= modulatorToEdit->lastBreakpointIndex() )
  {
    return modulatorToEdit->getBreakpointShapeAmount(selectedBreakpoint);
  }
  else
  {
    //jassertfalse
    return 0.0;
  }
}

int ModulatorCurveEditor::getMostRecentBreakpointIndex()
{
  return mostRecentBreakpointIndex;
}

int ModulatorCurveEditor::getMostRecentBreakpointAction()
{
  return mostRecentAction;
}

int ModulatorCurveEditor::whatIsUnderTheMouseCursor(const MouseEvent &e)
{
  if( modulatorToEdit == NULL )
    return 0;

  // get the position of the event in components coordinates
  mouseX = e.getMouseDownX();
  mouseY = e.getMouseDownY();

  // coordinates (measured in the component coordinate-system) of the point
  // which is currently under examination:
  double pointX;
  double pointY;

  // a margin for the breakpoint-dots
  double marginInPixels = dotRadius;

  // retrieve the scale and offset variables form the Modulator-object
  double scale          = modulatorToEdit->getScaleFactor();
  double offset         = modulatorToEdit->getOffset();

  // loop through the points to identify the point which is being moved (if any)
  int p;
  for(p=0; p<=modulatorToEdit->lastBreakpointIndex(); p++)
  {
    // select a point and transform it to the component-coordinate system:
    pointX = modulatorToEdit->getBreakpointTime(p);
    pointY = modulatorToEdit->getBreakpointLevel(p);
    pointY = scale*pointY + offset;
    transformToComponentsCoordinates(pointX, pointY);

    // measure the distance of the point to the location of the mouse-down event:
    double d = sqrt( (mouseX-pointX)*(mouseX-pointX) + (mouseY-pointY)*(mouseY-pointY)    );

    //  use a larger margin, if the breakpoint is the selected one:
    if( p == selectedBreakpoint )
      marginInPixels = 2*dotRadius;
    else
      marginInPixels = dotRadius;

    // return the index, which indicates that some breakpoint is under the mouse:
    if( d <= marginInPixels )
      return SOME_BREAKPOINT;
  }

  int    p1 = modulatorToEdit->getLoopStartIndex();
  double x1 = modulatorToEdit->getBreakpointTime(p1);
  int    p2 = modulatorToEdit->getLoopEndIndex();
  double x2 = modulatorToEdit->getBreakpointTime(p2);
  double y1 = 0.0; // dummy
  double y2 = 0.0; // dummy
  transformToComponentsCoordinates(x1, y1);
  transformToComponentsCoordinates(x2, y2);

  if( abs(x1+2-mouseX) <= 4.0   &&
      modulatorToEdit->getLoopMode() != rosic::BreakpointModulator::NO_LOOP )
  {
    return LOOP_START_LOCATOR;
  }
  else if( abs(x2+2-mouseX) <= 4.0  &&
           modulatorToEdit->getLoopMode() != rosic::BreakpointModulator::NO_LOOP )
  {
    return LOOP_END_LOCATOR;
  }

  // no object has been identified to be under the mouse cursor:
  return NO_OBJECT;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void ModulatorCurveEditor::parameterChanged(Parameter* parameterThatHasChanged)
{

  // maybe we need to repaint here....

}

void ModulatorCurveEditor::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{

  // invalidate references to any assigned parameter

}

void ModulatorCurveEditor::mouseDown(const MouseEvent &e)
{
  if( modulatorToEdit == NULL )
    return;

  // get the position of the event in components coordinates
  mouseX = e.getMouseDownX();
  mouseY = e.getMouseDownY();

  // coordinates (measured in the component coordinate-system) of the point
  // which is currently under examination:
  double pointX;
  double pointY;

  // add a ControlPoint at the current position or move an existing 
  // Control-Point, when the left button was pressed:
  breakpointBeingDragged = -1;
  double marginInPixels = dotRadius;
  double scale  = modulatorToEdit->getScaleFactor();
  double offset = modulatorToEdit->getOffset();

  if( e.mods.isLeftButtonDown() )
  {
    // loop through the points to identify the point which is being moved (if
    // any), if no existing point is being moved, a new one will be added:
    breakpointBeingDragged = -1;
    int p; // for the breakpoint index
    for(p=0; p<=modulatorToEdit->lastBreakpointIndex(); p++)
    {
      // select a point and transform it to the component-coordinate system:
      pointX = modulatorToEdit->getBreakpointTime(p);
      pointY = modulatorToEdit->getBreakpointLevel(p);
      pointY = scale*pointY + offset;
      transformToComponentsCoordinates(pointX, pointY);

      // measure the distance of the point to the location of the 
      // mouse-down event:
      double d = sqrt( (mouseX-pointX)*(mouseX-pointX) + 
        (mouseY-pointY)*(mouseY-pointY)    );

      // select the point, if the distance is below the margin:
      if( p == selectedBreakpoint )
        marginInPixels = 2*dotRadius;
      else
        marginInPixels = dotRadius;
      if( d <= marginInPixels )
      {
        breakpointBeingDragged    = p;   
        selectedBreakpoint        = breakpointBeingDragged;
        mostRecentBreakpointIndex = breakpointBeingDragged;
        mostRecentAction          = BREAKPOINT_MODIFIED;
      }
    } // end of  for(p=0; p<=modulatorToEdit->lastBreakpointIndex(); p++)

    // check if one of the locators is being dragged (we need to check that only
    // when no point is being dragged):
    if( breakpointBeingDragged == -1 && 
        modulatorToEdit->getLoopMode() != rosic::BreakpointModulator::NO_LOOP)
    {
      if( whatIsUnderTheMouseCursor(e) == LOOP_START_LOCATOR )
      {
        locatorBeingDragged       = LOOP_START_LOCATOR;
        mostRecentAction          = LOOP_START_MODIFIED;
        mostRecentBreakpointIndex = -1;
      }
      else if( whatIsUnderTheMouseCursor(e) == LOOP_END_LOCATOR )
      {
        locatorBeingDragged       = LOOP_END_LOCATOR;
        mostRecentAction          = LOOP_END_MODIFIED;
        mostRecentBreakpointIndex = -1;
      }

      /*
      p = modulatorToEdit->getLoopStartIndex();
      pointX = modulatorToEdit->getBreakpointTime(p);
      pointY = range.getMaxY();
      transformToComponentsCoordinates(pointX, pointY);

      if( abs(pointX+4-mouseX) <= 4.0 &&
        abs(pointY+5-mouseY) <= 4.0 &&
        modulatorToEdit->getLoopMode() != rosic::BreakpointModulator::NO_LOOP )
      {
        locatorBeingDragged       = LOOP_START_LOCATOR;
        mostRecentAction          = LOOP_START_MODIFIED;
        mostRecentBreakpointIndex = -1;
      }

      p = modulatorToEdit->getLoopEndIndex();
      pointX = modulatorToEdit->getBreakpointTime(p);
      pointY = range.getMaxY();
      transformToComponentsCoordinates(pointX, pointY);

      if( abs(pointX-4-mouseX) <= 4.0 &&
        abs(pointY+5-mouseY) <= 4.0 &&
        modulatorToEdit->getLoopMode() != rosic::BreakpointModulator::NO_LOOP )
      {
        locatorBeingDragged       = LOOP_END_LOCATOR;
        mostRecentAction          = LOOP_END_MODIFIED;
        mostRecentBreakpointIndex = -1;
      }
      */
    }

    // insert a new point into the array, if no existing point or locator is
    // being moved:
    if( breakpointBeingDragged == -1 && locatorBeingDragged == -1 )
    {
      double newX = (double) mouseX;
      double newY = (double) mouseY;
      transformFromComponentsCoordinates(newX, newY);

      newY = (newY-offset) / scale;

      // we have to insert a new breakpoint into our modulatorToEdit and keep 
      // track of it:
      breakpointBeingDragged    = modulatorToEdit->insertBreakpoint(newX, newY);

      selectedBreakpoint        = breakpointBeingDragged;
      mostRecentBreakpointIndex = breakpointBeingDragged;
      mostRecentAction          = BREAKPOINT_INSERTED;
    }

  } // end of if( e.mods.isLeftButtonDown() )

  // remove a ControlPoint at the current position, when the right button was 
  // pressed:
  else if( e.mods.isRightButtonDown() )
  {
    bool breakpointWasRemoved = false;

    // loop through the points to identify the point which is being removed (if
    // any) and remove it:
    for(int p=0; p<=modulatorToEdit->lastBreakpointIndex(); p++)
    {
      // select a point and transform it to the component-coordinate system:
      pointX = modulatorToEdit->getBreakpointTime(p);
      pointY = modulatorToEdit->getBreakpointLevel(p);
      pointY = scale*pointY + offset;
      transformToComponentsCoordinates(pointX, pointY);

      // measure the distance of the point to the location of the 
      // mouse-down event:
      double d = sqrt( (mouseX-pointX)*(mouseX-pointX) + 
        (mouseY-pointY)*(mouseY-pointY)    );

      // delete the point, if the distance is below the margin:
      if( p == selectedBreakpoint )
        marginInPixels = 2*dotRadius;
      else
        marginInPixels = dotRadius;
      if( d <= marginInPixels )
      {
        modulatorToEdit->removeBreakpoint(p);
        selectedBreakpoint        = -1;
        mostRecentBreakpointIndex = p;
        mostRecentAction          = BREAKPOINT_REMOVED;
        breakpointWasRemoved      = true;
        sendChangeMessage();

        //updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
        //return; 
        // we must return prematurely to avoid that more than one breakpoint can
        // be removed at a time - that would invalidate the 
        // mostRecentBreakpointIndex
      }

    } // end of for(p=0; p<modulatorToEdit->lastBreakpointIndex(); p++)

    // if a breakpoint was removed by this action, the corresponding flag is true now - otherwise
    // the user did a right-click into an empty area. in this case we fall back to the standard
    // InteactiveCoordinateSystem right-click behaviour (open the right-click context menu):
    if( breakpointWasRemoved == false )
      InteractiveCoordinateSystemOld::mouseDown(e);

  } // end of else if( e.mods.isRightButtonDown() )

  // calculate the new curve (this will not trigger a repaint):
  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);

  // removing breakpoints may change the maximum meaning full range to be displayed, so we change
  // the maximum range of the inherited CoordinateSystem:
  updateMaximumRange();

  // inform our listeners about the change:
  sendChangeMessage();

  // TEST:
  //repaint();
}

void ModulatorCurveEditor::mouseDrag(const MouseEvent &e)
{
  if( modulatorToEdit == NULL )
    return;

  // get the position of the event in components coordinates
  mouseX = e.getMouseDownX() + e.getDistanceFromDragStartX();
  mouseY = e.getMouseDownY() + e.getDistanceFromDragStartY();

  double x, y ;   // for the mapped values
  int    p;

  // drag a point to a new position:
  if( breakpointBeingDragged != -1 &&
    breakpointBeingDragged <= modulatorToEdit->lastBreakpointIndex() )
  {
    x = (double) mouseX;
    y = (double) mouseY;
    transformFromComponentsCoordinates(x,y);

    // take scale and shift into account:
    double scale  = modulatorToEdit->getScaleFactor();
    double offset = modulatorToEdit->getOffset();
    if( scale != 0.0 )
      y = (y - offset) / scale;
    else
      y = 0.0;

    // snap to the grid of the coordinate system, if desired:
    snapToGrid(x,y);

    if( e.mods.isShiftDown() )
      modulatorToEdit->setEditMode(BreakpointModulator::EDIT_WITH_SHIFT);
    else
      modulatorToEdit->setEditMode(BreakpointModulator::EDIT_WITHOUT_SHIFT);

    // pass the new breakpoint to our preview envelope generator:
    modulatorToEdit->modifyBreakpoint(breakpointBeingDragged, x, y);

    mostRecentBreakpointIndex = breakpointBeingDragged;
    mostRecentAction          = BREAKPOINT_MODIFIED;
  }
  // drag the left locator to a new position:
  else if ( locatorBeingDragged == LOOP_START_LOCATOR )
  {
    // find the breakpoint which is closest to the current x-position of the 
    // mouse:
    double minDistance = 1000000.0;
    int    minIndex    = 0;
    for(p=0; p<modulatorToEdit->lastBreakpointIndex(); p++)
    {
      x = modulatorToEdit->getBreakpointTime(p);
      y = 0.f;
      transformToComponentsCoordinates(x,y);
      if( abs(mouseX-x) < minDistance )
      {
        minDistance = abs(mouseX-x);
        minIndex    = p;
      }
    }
    // ...and make this breakpoint to the new loop-start:
    modulatorToEdit->setLoopStartIndex(minIndex);

    mostRecentBreakpointIndex = -1;
    mostRecentAction          = NO_ACTION;
  }
  // drag the right locator to a new position:
  else if ( locatorBeingDragged == LOOP_END_LOCATOR )
  {
    // find the breakpoint which is closest to the current x-position of the 
    // mouse:
    double minDistance = 1000000.0;
    int    minIndex    = 0;
    for(p=0; p<=modulatorToEdit->lastBreakpointIndex(); p++)
    {
      x = modulatorToEdit->getBreakpointTime(p);
      y = 0.f;
      transformToComponentsCoordinates(x,y);
      if( abs(mouseX-x) < minDistance )
      {
        minDistance = abs(mouseX-x);
        minIndex    = p;
      }
    }
    // ...and make this breakpoint to the new loop-start:
    modulatorToEdit->setLoopEndIndex(minIndex);

    mostRecentBreakpointIndex = -1;
    mostRecentAction          = NO_ACTION;
  }

  // moving breakpoints may change the maximum meaning full range to be displayed, so we change
  // the maximum range of the inherited CoordinateSystem:
  //updateMaximumRange();

  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
  sendChangeMessage();

  // TEST:
  //repaint();
}

void ModulatorCurveEditor::mouseMove(const MouseEvent &e)
{
  switch( whatIsUnderTheMouseCursor(e) )
  {
  case NO_OBJECT:       
    currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);          
    break;
  case SOME_BREAKPOINT: 
    currentMouseCursor = MouseCursor(MouseCursor::PointingHandCursor);    
    break;
  case LOOP_START_LOCATOR:    
    currentMouseCursor = MouseCursor(MouseCursor::LeftRightResizeCursor); 
    break;
  case LOOP_END_LOCATOR:    
    currentMouseCursor = MouseCursor(MouseCursor::LeftRightResizeCursor); 
    break;
  }
}

void ModulatorCurveEditor::mouseUp(const MouseEvent &e)
{
  breakpointBeingDragged    = -1;
  locatorBeingDragged       = -1;
  mostRecentBreakpointIndex = -1;
  mostRecentAction          = NO_ACTION;

  // moving breakpoints may change the maximum meaning full range to be displayed, so we change
  // the maximum range of the inherited CoordinateSystem:
  updateMaximumRange();

  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
  sendChangeMessage();
}

void ModulatorCurveEditor::resized()
{
  CurveFamilyPlotOld::resized();

  // (re) allocate memory if necesarry fill the arrays for the magnitude plot
  if( numSamplesInPlot != getWidth()+2 )
  {
    freePlotBufferArrays();
    numSamplesInPlot = getWidth()+2; // one sample overhead at left and right
    allocatePlotBufferArrays();
  }

  setCurveFamilyValues(numSamplesInPlot, numModulators, plotDataX, plotDataY);
  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
}

//-------------------------------------------------------------------------------------------------
// others:

void ModulatorCurveEditor::updatePlotCurveData(int curveIndex, BreakpointModulator* modulator, 
                                               bool updateGUI)
{
  if( modulator == NULL )
    return;

  jassert(curveIndex >= 0 && curveIndex < numModulators);
  if( curveIndex >= numModulators || curveIndex < 0 )
    return;

  // calculate some variables for the plot's samplerate, etc.:
  double visibleRangeX     = getCurrentRangeMaxX()-getCurrentRangeMinX(); // in seconds or beats
  double secondsPerPixel   = visibleRangeX / (double) numSamplesInPlot;
  double pixelsPerSecond   = 1.0 / secondsPerPixel;
  double envStartOffset    = modulator->getBreakpointTime(0);
  int    firstSampleToPlot = jmax(0, (int) floor(getCurrentRangeMinX()*pixelsPerSecond));
  double startTime         = firstSampleToPlot/pixelsPerSecond + envStartOffset;

  // copy the data from the modulator into a temporary object which has mostly the same data 
  // as the passed modulator now, but with some adjustments for plotting:
  BreakpointModulator tmpModulator;
  tmpModulator.copyDataFrom(*modulator);
  tmpModulator.setLoopMode(false);
  tmpModulator.setSyncMode(false);
  tmpModulator.setTimeScale(1.0);
  tmpModulator.setDepth(1.0);
  tmpModulator.setSampleRate(pixelsPerSecond);
  tmpModulator.noteOnAndAdvanceTime(firstSampleToPlot);

  for(int i=0; i<numSamplesInPlot; i++)
  {
    plotDataX[curveIndex][i] = startTime + (i*secondsPerPixel);
    plotDataY[curveIndex][i] = tmpModulator.getSample();
    plotDataY[curveIndex][i] = clip(plotDataY[curveIndex][i], -100.0, +100.0);
  }

  // call CurveFamilyPlot::updatePlotImage() to reflect the new data on the GUI:
  if( updateGUI == true )
    updatePlotImage();
}

void ModulatorCurveEditor::updatePlotCurveData()
{ 
  updatePlotCurveData(editedModulatorIndex, modulatorToEdit, true);
}

void ModulatorCurveEditor::updateMaximumRange(bool alsoUpdateCurrentRange)
{
  if( modulatorToEdit == NULL )
    return;

  double minX    = modulatorToEdit->getStartTime();
  double maxX    = modulatorToEdit->getEndTime();
  double minY    = modulatorToEdit->getMinLevel();
  double maxY    = modulatorToEdit->getMaxLevel();
  minY           = jmin(0.0, minY);
  maxY           = jmax(1.0, maxY);
  double marginX = 0.1 * (maxX-minX);
  double marginY = 0.1 * (maxY-minY);

  if( alsoUpdateCurrentRange == true )
  {
    setAutoReRendering(false);
    setMaximumRange(minX-marginX, maxX+marginX, minY-marginY, maxY+marginY);
    setAutoReRendering(true);
    setCurrentRange(minX-marginX, maxX+marginX, minY-marginY, maxY+marginY);
  }
  else
    setMaximumRange(minX-marginX, maxX+marginX, minY-marginY, maxY+marginY);

}

//-------------------------------------------------------------------------------------------------
// internal functions:

void ModulatorCurveEditor::plotCurveFamily(Graphics &g, Image *targetImage, XmlElement *targetSVG)
{
  // call the paint-method of the CurveFamilyPlot base class, which draws the axes, labels, grids,
  // the curves themselves etc.:
  CurveFamilyPlotOld::plotCurveFamily(g, targetImage, targetSVG);


  // draw the dots and loop locators - \todo: choose colours:
  //Colour dotColour = colourScheme.curves; 
  Colour dotColour = plotColourScheme.getCurveColour(0);
  //if( colourScheme.plotColours.size() > 0 )
  //  dotColour = colourScheme.plotColours[0];


  plotBreakpoints(g, targetImage, modulatorToEdit, dotColour);

  Colour locatorColour = Colours::yellow;
  //if( colourScheme.plotColours.size() > 1 )    
  //  locatorColour = colourScheme.plotColours[1];


  plotLoopLocators(g, targetImage, modulatorToEdit, locatorColour);
}

void ModulatorCurveEditor::plotBreakpoints(Graphics &g, Image *targetImage, 
  BreakpointModulator* modulator, const Colour& dotColour)
{
  if( modulator == NULL )
    return;

  double x, y;
  double scale  = modulator->getScaleFactor();
  double offset = modulator->getOffset();

  g.setColour(dotColour);
  for(int p=0; p<=modulator->lastBreakpointIndex(); p++)
  {
    x = modulator->getBreakpointTime(p);
    //if( modulator->isInSyncMode() )
    //x = beatsToSeconds(x, modulator->getBeatsPerMinute());
    y = modulator->getBreakpointLevel(p);
    y = scale*y + offset;
    transformToImageCoordinates(x, y, targetImage);
    if( p == selectedBreakpoint )
    {
      g.setColour(dotColour.withAlpha(0.5f));
      g.fillEllipse((float) (x-2.5*dotRadius), (float) (y-2.5*dotRadius), 
        (float) (5*dotRadius), (float) (5*dotRadius) );
      g.setColour(dotColour);
      g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius), 
        (float) (2*dotRadius), (float) (2*dotRadius) );
    }
    else
    {
      g.setColour(dotColour);
      g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius), 
        (float) (2*dotRadius), (float) (2*dotRadius) );
    }
  }
}

void ModulatorCurveEditor::plotLoopLocators(Graphics &g, Image *targetImage, 
  BreakpointModulator* modulator, const Colour& locatorColour, bool fullHeight)
{
  if( modulator == NULL )
    return;

  double scale  = modulator->getScaleFactor();
  double offset = modulator->getOffset();

  if( modulator->getLoopMode() == rosic::BreakpointModulator::NO_LOOP )
    return; // nothing to do

  double x1, x2, y1, y2;
  g.setColour(locatorColour);

  int p = modulator->getLoopStartIndex();
  x1 = modulator->getBreakpointTime(p);
  if( fullHeight == true )
  {
    drawLeftLocator(g, (float) x1, InteractiveCoordinateSystemOld::ARROW_AT_TOP, 
      locatorColour, targetImage);
  }
  y1 = modulator->getBreakpointLevel(p);
  CoordinateSystemOld::transformToComponentsCoordinates(x1, y1);
  drawTriangle(g, (float)x1-4.f, (float)y1-6.f, (float)x1-4.f, (float)y1+6.f, (float)x1+6.f, 
               (float)y1, true); 

  p  = modulator->getLoopEndIndex();
  x1 = modulator->getBreakpointTime(p);
  if( fullHeight == true )
  {
    drawRightLocator(g, (float) x1, InteractiveCoordinateSystemOld::ARROW_AT_TOP, 
      locatorColour, targetImage);
  }
  y1 = modulator->getBreakpointLevel(p);
  CoordinateSystemOld::transformToComponentsCoordinates(x1, y1);
  drawTriangle(g, (float)x1-6.f, (float)y1, (float)x1+4.f, (float)y1-6.f, (float)x1+4.f, 
               (float)y1+6.f, true); 



  // draw a preliminary transparent locator of the currently edited modulator, when a locator 
  // is being dragged (we need to check modulator == modulatorToEdit because otherwise the call
  // from subclass ModulatorCurveEditorMulti will draw preliminary locators for all modulators 
  // (which messes with the colour becasue of overlays)):
  g.setColour(Colours::yellow.darker(0.25f).withAlpha(0.5f));
  if( locatorBeingDragged == LOOP_START_LOCATOR && modulator == modulatorToEdit )
  {
    x1 = mouseX;
    y1 = mouseY;
    transformFromComponentsCoordinates(x1, y1);
    drawLeftLocator(g, (float) x1, InteractiveCoordinateSystemOld::ARROW_AT_TOP, 
      locatorColour.withMultipliedAlpha(0.5f));
  }
  else if( locatorBeingDragged == LOOP_END_LOCATOR && modulator == modulatorToEdit )
  {
    x1 = mouseX;
    y1 = mouseY;
    transformFromComponentsCoordinates(x1, y1);
    drawRightLocator(g, (float) x1, InteractiveCoordinateSystemOld::ARROW_AT_TOP, 
      locatorColour.withMultipliedAlpha(0.5f));
  }

  // connect the loop start and end-points with a line (green if start and end-level match, red 
  // otherwise) - also only if the passed modulator is currently currently:
  if( modulator == modulatorToEdit )
  {
    p  = modulator->getLoopStartIndex();
    x1 = modulator->getBreakpointTime(p);
    y1 = modulator->getBreakpointLevel(p);
    y1 = scale*y1 + offset;
    transformToImageCoordinates(x1, y1, targetImage);
    p  = modulator->getLoopEndIndex();
    x2 = modulator->getBreakpointTime(p);
    y2 = modulator->getBreakpointLevel(p);
    y2 = scale*y2 + offset;
    transformToImageCoordinates(x2, y2, targetImage);
    if( abs(y1-y2) < 0.00001 ) // equality check with margin
      g.setColour(matchedLoopConnectorColour);
    else
      g.setColour(unmatchedLoopConnectorColour);
    g.drawLine((float) x1, (float) y1, (float) x2, (float) y2, 1.0);
  }
}

void ModulatorCurveEditor::allocatePlotBufferArrays()
{
  // allocate the memory:
  plotDataFlatX = new double[numModulators*numSamplesInPlot];
  plotDataFlatY = new double[numModulators*numSamplesInPlot];
  plotDataX     = new double*[numModulators];
  plotDataY     = new double*[numModulators];

  // initialize with zeros:
  int i;
  for(i=0; i<numModulators*numSamplesInPlot; i++)
  {
    plotDataFlatX[i] = 0.0;
    plotDataFlatY[i] = 0.0;
  }

  // set up the pointer-pointers:
  for(i=0; i<numModulators; i++)
  {
    plotDataX[i] = &(plotDataFlatX[i*numSamplesInPlot]);
    plotDataY[i] = &(plotDataFlatY[i*numSamplesInPlot]);
  }
}

void ModulatorCurveEditor::freePlotBufferArrays()
{
  deleteAndZero(plotDataX);
  deleteAndZero(plotDataFlatX);
  deleteAndZero(plotDataY);
  deleteAndZero(plotDataFlatY);
}
