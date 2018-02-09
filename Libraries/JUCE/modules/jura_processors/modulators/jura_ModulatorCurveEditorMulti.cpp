//#include "rosof_ModulatorCurveEditorMulti.h"
//using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulatorCurveEditorMulti::ModulatorCurveEditorMulti(const juce::String& name)
: InteractiveCoordinateSystemOld(name), ModulatorCurveEditor(name)
{
  //editedModulatorIndex = -1;
}

ModulatorCurveEditorMulti::~ModulatorCurveEditorMulti()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void ModulatorCurveEditorMulti::addModulatorToEdit(
  rosic::BreakpointModulator* newModulatorToEdit)
{
  modulators.getLock().enter();
  modulators.addIfNotAlreadyThere(newModulatorToEdit);
  numModulators = modulators.size();
  modulators.getLock().exit();
}

void ModulatorCurveEditorMulti::selectModulatorToEdit(int index)
{
  modulators.getLock().enter();
  jassert( index < modulators.size() ); // index out of range
  editedModulatorIndex  = index % modulators.size();
  modulatorToEdit       = modulators[index];
  modulators.getLock().exit();
  setHighlightedCurve(index);
  updateMaximumRange(false);
  updatePlotImage(false);
}

void ModulatorCurveEditorMulti::setCurrentRangeX(double newMinX, double newMaxX)
{
  if( newMinX != getCurrentRangeMinX() || newMaxX != getCurrentRangeMaxX() )
  {
    rsDataPlot::setCurrentRangeX(newMinX, newMaxX);
    updateCurveDataForAllPlots(true, true);
  }
}

void ModulatorCurveEditorMulti::setCurrentRange(rsPlotRange newRange)
{
  if( newRange != getCurrentRange() )
  {
    rsDataPlot::setCurrentRange(newRange);
    updateCurveDataForAllPlots(true, true);
  }
}

void ModulatorCurveEditorMulti::clearCurveColours()
{
  curveColours.clear();
}

void ModulatorCurveEditorMulti::appendCurveColour(const ColourAHSL& colourToAppend)
{
  curveColours.add(colourToAppend);
}

void ModulatorCurveEditorMulti::changeCurveColour(int index, const ColourAHSL& newColour)
{
  jassert(index < curveColours.size());  // index out of range
  if( index >= curveColours.size() )
    return;
  curveColours.set(index, newColour);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

Colour ModulatorCurveEditorMulti::getCurveColour(int index) const
{
  if( curveColours.size() == 0 )
    return ModulatorCurveEditor::getCurveColour(index);
  else
    return curveColours[index % curveColours.size()].getAsJuceColour();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void ModulatorCurveEditorMulti::resized()
{
  rsDataPlot::resized();

  // (re) allocate memory if necesarry fill the arrays for the magnitude plot
  if( numSamplesInPlot != getWidth()+2 )
  {
    freePlotBufferArrays();
    numSamplesInPlot = getWidth()+2; // one sample overhead at left and right
    allocatePlotBufferArrays();
  }
  setCurveFamilyValues(numSamplesInPlot, numModulators, plotDataX, plotDataY);

  updateCurveDataForAllPlots(true, true);
   // this is the only line which is different from the basecalss implementation - maybe the
   // other stuff could be factored out into a function (better style)...
}

void ModulatorCurveEditorMulti::updateCurveDataForAllPlots(bool redrawCurves,
  bool redrawCoordinateSystem)
{
  modulators.getLock().enter();
  for(int i=0; i<modulators.size(); i++)
    updatePlotCurveData(i, modulators[i], false);
  modulators.getLock().exit();
  if( redrawCurves )
    updatePlotImage(redrawCoordinateSystem); // called 2 times on zoom
}

void ModulatorCurveEditorMulti::updateMaximumRange(bool alsoUpdateCurrentRange)
{
  if( modulatorToEdit == NULL )
    return;

  double minX    = modulatorToEdit->getStartTime();
  double maxX    = modulatorToEdit->getEndTime();
  double minY    = modulatorToEdit->getMinLevel();
  double maxY    = modulatorToEdit->getMaxLevel();

  modulators.getLock().enter();
  for(int m=0; m<modulators.size(); m++)
  {
    minX = jmin(minX, modulators[m]->getStartTime());
    maxX = jmax(maxX, modulators[m]->getEndTime());
    minY = jmin(minY, modulators[m]->getMinLevel());
    maxY = jmax(maxY, modulators[m]->getMaxLevel());
  }
  modulators.getLock().exit();

  minY           = jmin(0.0, minY);
  maxY           = jmax(1.0, maxY);
  double marginX = 0.1 * (maxX-minX);
  double marginY = 0.1 * (maxY-minY);

  if( alsoUpdateCurrentRange == true )
  {
    setAutoReRendering(false);
    rsPlotRange range(minX-marginX, maxX+marginX, minY-marginY, maxY+marginY);
    setMaximumRange(range);
    setCurrentRange(range);
    //setMaximumRange(minX-marginX, maxX+marginX, minY-marginY, maxY+marginY);
    //setCurrentRange(minX-marginX, maxX+marginX, minY-marginY, maxY+marginY);
    setAutoReRendering(true);
    //updateCurveDataForAllPlots(true, true);
  }
  else
    setMaximumRange(minX-marginX, maxX+marginX, minY-marginY, maxY+marginY);
}

void ModulatorCurveEditorMulti::plotCurveFamily(Graphics &g, juce::Image *targetImage,
  XmlElement *targetSVG)
{
  rsDataPlot::plotCurveFamily(g, targetImage, targetSVG); // draws axes, labels, grids, curves

  // draw the locators for all modulators but the breakpoint dots only for the focused one:
  modulators.getLock().enter();
  int m;
  for(m=0; m<modulators.size(); m++)
  {
    Colour locatorColour = getCurveColour(m);
    g.setColour(locatorColour);

    //if( colourScheme.plotColours.size() > m )
    //  locatorColour = colourScheme.plotColours[m].withMultipliedAlpha(1.0);

    if( m == editedModulatorIndex )
      plotLoopLocators(g, targetImage, modulators[m], locatorColour, true);
    else
      plotLoopLocators(g, targetImage, modulators[m], locatorColour.withMultipliedAlpha(0.5), false);
  }
  modulators.getLock().exit();

  m = editedModulatorIndex;
  Colour dotColour = getCurveColour(m);

  //if( colourScheme.plotColours.size() > m )
  //  dotColour = colourScheme.plotColours[m];

  plotBreakpoints(g, targetImage, modulatorToEdit, dotColour);

  /*
  Colour locatorColour = Colours::yellow;
  if( colourScheme.plotColours.size() > m )
    locatorColour = colourScheme.plotColours[m].withMultipliedAlpha(0.75);
  plotLoopLocators(g, targetImage, modulators[m], locatorColour);
  */
}
