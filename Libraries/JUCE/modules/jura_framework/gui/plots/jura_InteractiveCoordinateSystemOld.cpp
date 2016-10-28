//-------------------------------------------------------------------------------------------------
// construction/destruction:

InteractiveCoordinateSystemOld::InteractiveCoordinateSystemOld(const String& name) 
: CoordinateSystemOld(name)
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

  setDescription("Some coordinate system widget.");
}

InteractiveCoordinateSystemOld::~InteractiveCoordinateSystemOld()
{

}

//-------------------------------------------------------------------------------------------------
// callbacks:

void InteractiveCoordinateSystemOld::mouseDown(const MouseEvent& e)
{
  if( e.mods.isRightButtonDown() && showPopUpOnRightClick == true )
    openRightClickPopupMenu();
}

/*
void InteractiveCoordinateSystemOld::mouseEnter(const MouseEvent& e)
{
  CoordinateSystemOld::mouseEnter(e);
  RWidget::mouseEnter(e);
}

void InteractiveCoordinateSystemOld::mouseExit(const MouseEvent& e)
{
  CoordinateSystemOld::mouseExit(e);
  RWidget::mouseExit(e);
}
*/

void InteractiveCoordinateSystemOld::snapToCoarseGrid(double &x, double &y)
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

void InteractiveCoordinateSystemOld::snapToFineGrid(double &x, double &y)
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

void InteractiveCoordinateSystemOld::snapToGrid(double &x, double &y)
{
  snapToFineGrid(x,y);
  snapToCoarseGrid(x,y);
}

//-----------------------------------------------------------------------------

void InteractiveCoordinateSystemOld::drawLeftLocator(Graphics &g, float x, int arrowPosition, 
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

void InteractiveCoordinateSystemOld::drawRightLocator(Graphics &g, float x, int arrowPosition, 
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

void InteractiveCoordinateSystemOld::drawCurrentPositionLocator(Graphics &g, float x, 
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

XmlElement* InteractiveCoordinateSystemOld::getStateAsXml(const String& stateName) 
const
{
  // initialize a new XmlElement with all the relevant data from the 
  // CoordinateSystemOld base-class:
  XmlElement* xmlState = CoordinateSystemOld::getStateAsXml(stateName);

  // add parameters which are specific to this subclass:
  xmlState->setAttribute("SnapToCoarseGridX", snapToCoarseGridX);
  xmlState->setAttribute("SnapToCoarseGridY", snapToCoarseGridY);
  xmlState->setAttribute("SnapToFineGridX",   snapToFineGridX);
  xmlState->setAttribute("SnapToFineGridY",   snapToFineGridY);

  return xmlState;
}

bool InteractiveCoordinateSystemOld::setStateFromXml(const XmlElement &xmlState)
{
  // restore all the relevant parameters inherited form the CoordinateSystemOld 
  // base-class:
  bool success = CoordinateSystemOld::setStateFromXml(xmlState); 

  // restore all the relevant parameters which are specific to this subclass:
  snapToCoarseGridX = xmlState.getBoolAttribute("SnapToCoarseGridX", false);
  snapToCoarseGridY = xmlState.getBoolAttribute("SnapToCoarseGridY", false);
  snapToFineGridX   = xmlState.getBoolAttribute("SnapToFineGridX",   false);
  snapToFineGridY   = xmlState.getBoolAttribute("SnapToFineGridY",   false);

  return success; // if everything worked well, this flag is still true
}


void InteractiveCoordinateSystemOld::openRightClickPopupMenu()
{
  jassertfalse;
  // the code below for opening the context menu is outdated - change it to deal with the new 
  // RPopUpMenu

  /*
  RPopUpMenuOld* menu         = NULL;
  RPopUpMenuOld* xSnapSubMenu = NULL;
  RPopUpMenuOld* ySnapSubMenu = NULL;

  menu = new RPopUpMenuOld();
  menu->setColourScheme(popUpColourScheme);

  int index = 1; 
  menu->addItem(index, String(T("Export Image")));
  index++;

  int xSnapIntervalIndicesMin = index;
  if( xSnapIntervals.size() > 0 )
  {
    xSnapSubMenu = new RPopUpMenuOld();
    xSnapSubMenu->setColourScheme(popUpColourScheme);
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
    ySnapSubMenu = new RPopUpMenuOld();
    ySnapSubMenu->setColourScheme(popUpColourScheme);
    for(int i=0; i<ySnapIntervals.size(); i++)
    {
      ySnapSubMenu->addItem(index, String(ySnapIntervals[i]));
      index++;
    }
    menu->addSubMenu(String(T("Horizontal Grid")), *ySnapSubMenu);
  }
  int ySnapIntervalIndicesMax = jmax(ySnapIntervalIndicesMin, index-1);

  int result = menu->show(0, 130);

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
  */
}


