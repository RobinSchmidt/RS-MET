#include "rojue_InteractiveCoordinateSystem.h"
using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

InteractiveCoordinateSystem::InteractiveCoordinateSystem(const String& name) 
: CoordinateSystem(name)
{
  snapToCoarseGridX            =  false;
  snapToCoarseGridY            =  false;
  snapToFineGridX              =  false;
  snapToFineGridY              =  false;
  mouseX                       =  0;
  mouseY                       =  0;
  dotRadius = 3.f;
  loopLocatorColour            = Colour(0xff0000ff);
  matchedLoopConnectorColour   = Colour(0x0000ff00);
  unmatchedLoopConnectorColour = Colour(0xffff0000);

  setDescription(T("Some coordinate system widget."));
}

InteractiveCoordinateSystem::~InteractiveCoordinateSystem()
{

}

//-------------------------------------------------------------------------------------------------
// callbacks:

void InteractiveCoordinateSystem::mouseDown(const MouseEvent& e)
{
  if( e.mods.isRightButtonDown() )
    openRightClickPopupMenu();
}

/*
void InteractiveCoordinateSystem::mouseEnter(const MouseEvent& e)
{
  CoordinateSystem::mouseEnter(e);
  RWidget::mouseEnter(e);
}

void InteractiveCoordinateSystem::mouseExit(const MouseEvent& e)
{
  CoordinateSystem::mouseExit(e);
  RWidget::mouseExit(e);
}
*/

void InteractiveCoordinateSystem::snapToCoarseGrid(double &x, double &y)
{
  // snap to the (coarse) grid of the coordinate system, if desired:
  double tmp;
  if( snapToCoarseGridX == true )
  {
    tmp = fmod(fabs(x), getVerticalCoarseGridInterval());
    if( tmp < 0.5*getVerticalCoarseGridInterval() )
      x = sign(x) * (fabs(x)-tmp);
    else 
      x = sign(x) * (fabs(x)-tmp+getVerticalCoarseGridInterval());
  }
  if( snapToCoarseGridY == true )
  {
    tmp = fmod(fabs(y), getHorizontalCoarseGridInterval());
    if( tmp < 0.5*getHorizontalCoarseGridInterval() )
      y = sign(y) * (fabs(y)-tmp);
    else 
      y = sign(y) * (fabs(y)-tmp+getHorizontalCoarseGridInterval());
  }
}

void InteractiveCoordinateSystem::snapToFineGrid(double &x, double &y)
{
  // snap to the (fine) grid of the coordinate system, if desired:
  double tmp;
  if( snapToFineGridX == true )
  {
    tmp = fmod(fabs(x), getVerticalFineGridInterval());
    if( tmp < 0.5*getVerticalFineGridInterval() )
      x = sign(x) * (fabs(x)-tmp);
    else 
      x = sign(x) * (fabs(x)-tmp+getVerticalFineGridInterval());
  }
  if( snapToFineGridY == true )
  {
    tmp = fmod(fabs(y), getHorizontalFineGridInterval());
    if( tmp < 0.5*getHorizontalFineGridInterval() )
      y = sign(y) * (fabs(y)-tmp);
    else 
      y = sign(y) * (fabs(y)-tmp+getHorizontalFineGridInterval());
  }
}

void InteractiveCoordinateSystem::snapToGrid(double &x, double &y)
{
  snapToFineGrid(x,y);
  snapToCoarseGrid(x,y);
}

//-----------------------------------------------------------------------------

void InteractiveCoordinateSystem::drawLeftLocator(Graphics &g, float x, int arrowPosition, 
                                                  const Colour &locatorColour, Image *targetImage)
{
  g.setColour(locatorColour);

  double x1, x2, y1, y2;
  x1 = x;
  x2 = x;
  y1 = currentRange.getMinY();
  y2 = currentRange.getMaxY();
  transformToImageCoordinates(x1, y1, targetImage);
  transformToImageCoordinates(x2, y2, targetImage);
  g.drawLine((float) x1, (float) (y1-1), (float) x2, (float) (y2+1), 2.0);
  g.fillRect((float) x1, (float) (y2+1), 4.f, 8.f);
  drawTriangle(g, (float) (x1+4), (float) (y2+1), (float) (x1+4), 
    (float) (y2+9), (float) (x1+4+4), (float) (y2+5), true);
}

void InteractiveCoordinateSystem::drawRightLocator(Graphics &g, float x, int arrowPosition, 
                                                   const Colour &locatorColour, Image *targetImage)
{
  g.setColour(locatorColour);

  double x1, x2, y1, y2;
  x1 = x;
  x2 = x;
  y1 = currentRange.getMinY();
  y2 = currentRange.getMaxY();
  transformToImageCoordinates(x1, y1, targetImage);
  transformToImageCoordinates(x2, y2, targetImage);
  g.drawLine((float) x1, (float) (y1-1), (float) x2, (float) (y2+1), 2.0);
  g.fillRect((float) (x1-4), (float) (y2+1), 4.f, 8.f);
  drawTriangle(g, (float) (x1-4), (float) (y2+1), (float) (x1-4), 
    (float) (y2+9), (float) (x1-4-4), (float) (y2+5), true);
}

void InteractiveCoordinateSystem::drawCurrentPositionLocator(Graphics &g, float x, 
                                                             int arrowPosition, const Colour &locatorColour, Image *targetImage)
{
  g.setColour(locatorColour);

  double x1, x2, y1, y2;
  x1 = x;
  x2 = x;
  y1 = currentRange.getMinY();
  y2 = currentRange.getMaxY();
  transformToImageCoordinates(x1, y1, targetImage);
  transformToImageCoordinates(x2, y2, targetImage);
  g.drawLine((float) x1, (float) (y1-1), (float) x2, (float) (y2+1), 2.0);
}

//-----------------------------------------------------------------------------

XmlElement* InteractiveCoordinateSystem::getStateAsXml(const String& stateName) 
const
{
  // initialize a new XmlElement with all the relevant data from the 
  // CoordinateSystem base-class:
  XmlElement* xmlState = CoordinateSystem::getStateAsXml(stateName);

  // add parameters which are specific to this subclass:
  xmlState->setAttribute(String(T("SnapToCoarseGridX")), snapToCoarseGridX);
  xmlState->setAttribute(String(T("SnapToCoarseGridY")), snapToCoarseGridY);
  xmlState->setAttribute(String(T("SnapToFineGridX")),   snapToFineGridX);
  xmlState->setAttribute(String(T("SnapToFineGridY")),   snapToFineGridY);

  return xmlState;
}

bool InteractiveCoordinateSystem::setStateFromXml(const XmlElement &xmlState)
{
  // restore all the relevant parameters inherited form the CoordinateSystem 
  // base-class:
  bool success = CoordinateSystem::setStateFromXml(xmlState); 

  // restore all the relevant parameters which are specific to this subclass:
  snapToCoarseGridX = xmlState.getBoolAttribute(String(T("SnapToCoarseGridX")), false);
  snapToCoarseGridY = xmlState.getBoolAttribute(String(T("SnapToCoarseGridY")), false);
  snapToFineGridX   = xmlState.getBoolAttribute(String(T("SnapToFineGridX")),   false);
  snapToFineGridY   = xmlState.getBoolAttribute(String(T("SnapToFineGridY")),   false);

  return success; // if everything worked well, this flag is still true
}


void InteractiveCoordinateSystem::openRightClickPopupMenu()
{
  PopupMenu* menu         = NULL;
  PopupMenu* xSnapSubMenu = NULL;
  PopupMenu* ySnapSubMenu = NULL;

  menu = new PopupMenu();

  int index = 1; 
  menu->addItem(index, String(T("Export Image")));
  index++;

  int xSnapIntervalIndicesMin = index;
  if( xSnapIntervals.size() > 0 )
  {
    xSnapSubMenu = new PopupMenu();
    for(int i=0; i<xSnapIntervals.size(); i++)
    {
      xSnapSubMenu->addItem(index, String(xSnapIntervals[i]));
      index++;
    }
    menu->addSubMenu(String(T("Vertical Grid")), *xSnapSubMenu);
  }
  int xSnapIntervalIndicesMax = jmax(xSnapIntervalIndicesMin, index-1);


  int ySnapIntervalIndicesMin = index;
  if( ySnapIntervals.size() > 0 )
  {
    ySnapSubMenu = new PopupMenu();
    for(int i=0; i<ySnapIntervals.size(); i++)
    {
      ySnapSubMenu->addItem(index, String(ySnapIntervals[i]));
      index++;
    }
    menu->addSubMenu(String(T("Horizontal Grid")), *ySnapSubMenu);
  }
  int ySnapIntervalIndicesMax = jmax(ySnapIntervalIndicesMin, index-1);


  // customize the look of the menu and show it:
  /*
  RLookAndFeel tmpLookAndFeel;
  tmpLookAndFeel.setWidgetColours(backgroundColour, outlineColour, handleColour, textColour, 
    specialColour1, specialColour2);
  menu->setLookAndFeel(&tmpLookAndFeel);
  //if( defaultValuesSubMenu != NULL )  
  //  defaultValuesSubMenu->setLookAndFeel(&tmpLookAndFeel);
  */

  int result = menu->show();

  // we have retrieved the result - we don't need the menus anymore:
  if( xSnapSubMenu != NULL )
    delete xSnapSubMenu;
  if( ySnapSubMenu != NULL )
    delete ySnapSubMenu;
  if( menu != NULL )
    delete menu;

  // do some action according to the chosen item:
  if(result == 1)
    openExportDialog(getWidth(), getHeight(), String(T("png")), File::nonexistent);
  else if( result <= xSnapIntervalIndicesMax )
  {
    int xSnapIntervalIndex = result - xSnapIntervalIndicesMin;
    if( xSnapIntervalIndex < xSnapIntervals.size() )
    {
      double snapInterval = xSnapIntervals[xSnapIntervalIndex];
      if( snapInterval > 0.0 )
      {
        setVerticalFineGrid(snapInterval, true);
        setSnapToFineGridX(true);    
        repaint();
      }
      else
      {
        setVerticalFineGridVisible(false);
        setSnapToFineGridX(false);
        repaint();
      }
    }
  }
  else if( result <= ySnapIntervalIndicesMax )
  {
    int ySnapIntervalIndex = result - ySnapIntervalIndicesMin;
    if( ySnapIntervalIndex < ySnapIntervals.size() )
    {
      double snapInterval = ySnapIntervals[ySnapIntervalIndex];
      if( snapInterval > 0.0 )
      {
        setHorizontalFineGrid(snapInterval, true);
        setSnapToFineGridY(true);    
        repaint();
      }
      else
      {
        setHorizontalFineGridVisible(false);
        setSnapToFineGridY(false);
        repaint();
      }
    }
  }
  
}


