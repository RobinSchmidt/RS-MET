using namespace RSLib;

//-------------------------------------------------------------------------------------------------
// class rsGraphicsRenderer2D (the abstract baseclass):

// construction/destruction:

rsGraphicsRenderer2D::rsGraphicsRenderer2D()
{

}

rsGraphicsRenderer2D::~rsGraphicsRenderer2D()
{

}

// setup:

void rsGraphicsRenderer2D::setColor(const rsColorRGBA& newColor)
{
  state.color = newColor;
}

void rsGraphicsRenderer2D::setLineThickness(const double newThickness)
{
  state.lineThickness = newThickness;
}

// drawing:

void rsGraphicsRenderer2D::drawFilledEllipseInscribedPolygon(int numVertices, double cx, double cy, 
                                                             double rx, double ry, double a)
{
  rsPolygon2D<double> polygon(numVertices, rsPoint2D<double>(cx, cy), rx, ry, a);  
  drawFilledPolygon(polygon);
}

void rsGraphicsRenderer2D::drawLine(double x1, double y1, double x2, double y2)
{
  double dx = x2 - x1;
  double dy = y2 - y1;
  double d  = rsSqrt(dx*dx + dy*dy);
  double s  = 0.5*state.lineThickness/d;
  dx        = s * dy;
  dy        = s * (x2-x1);  // == s times the old dx (before reassigning it)

  rsPolygon2D<double> p;
  p.reserveVertexMemory(4);
  p.addVertex(rsPoint2D<double>(x1 - dx,  y1 + dy));
  p.addVertex(rsPoint2D<double>(x2 - dx,  y2 + dy));
  p.addVertex(rsPoint2D<double>(x2 + dx,  y2 - dy));
  p.addVertex(rsPoint2D<double>(x1 + dx,  y1 - dy));
  drawFilledPolygon(p);
}

void rsGraphicsRenderer2D::drawOutlinedRectangle(double x, double y, double w, double h)
{
  double t = state.lineThickness;
  double t2 = 0.5*t;
  drawFilledRectangle(x-t2,   y-t2,   t,    h+t2);  // left edge
  drawFilledRectangle(x-t2,   y-t2,   w+t2, t);     // top edge
  drawFilledRectangle(x+w-t2, y-t2,   t,    h+t2);  // right edge
  drawFilledRectangle(x-t2,   y+h-t2, w+t2, t);     // bottom edge
}

void rsGraphicsRenderer2D::drawFilledRectangle(double x, double y, double w, double h)
{
  rsPolygon2D<double> p;
  p.reserveVertexMemory(4);
  p.addVertex(rsPoint2D<double>(x,   y  ));
  p.addVertex(rsPoint2D<double>(x+w, y  ));
  p.addVertex(rsPoint2D<double>(x+w, y+h));
  p.addVertex(rsPoint2D<double>(x,   y+h));
  drawFilledPolygon(p);
}

void rsGraphicsRenderer2D::drawFilledEllipse(double cx, double cy, double rx, double ry)
{
  int n = rsMax((int) ceil(rx), (int) ceil(ry));  // \todo use rsCeilInt
  if( rsIsOdd(n) )
    n += 1;
  n = rsMax(n, 8);                                // at least 8 vertices
  drawFilledEllipseInscribedPolygon(n, cx, cy, rx, ry);

  // \todo maybe an ellipse can be drawn more efficiently using some functions from
  // agg_renderer_primitives.h - look into this, if this turns out to be true, override it
  // rsGraphicsRenderer2DImage
}

void rsGraphicsRenderer2D::drawOutlinedPixelRectangle(int x, int y, int w, int h)
{
  drawFilledRectangle(x,       y,       1.0, h);    // left edge
  drawFilledRectangle(x,       y,       w,   1.0);  // top edge
  drawFilledRectangle(x+w-1.0, y,       1.0, h);    // right edge
  drawFilledRectangle(x,       y+h-1.0, w,   1.0);  // top edge
}

void rsGraphicsRenderer2D::drawFilledPixelRectangle(int x, int y, int w, int h)
{
  drawFilledRectangle(x, y, w, h);
}

//-------------------------------------------------------------------------------------------------
// class rsGraphicsRenderer2DImage:

rsGraphicsRenderer2DImage::rsGraphicsRenderer2DImage(rsImageRegionRGBA &imageRegionToRenderOn)
: imageRegion(imageRegionToRenderOn)
{

}

rsGraphicsRenderer2DImage::~rsGraphicsRenderer2DImage()
{

}

// drawing:

void rsGraphicsRenderer2DImage::fillAll()
{
  imageRegion.fillAll(state.color);
}

void rsGraphicsRenderer2DImage::drawFilledPolygon(const rsPolygon2D<double> &polygon)
{
  rsAGG::drawFilledPolygon(polygon, imageRegion, state.color);
}

void rsGraphicsRenderer2DImage::drawOutlinedPolygon(const rsPolygon2D<double> &polygon)
{
  rsAGG::drawOutlinedPolygon(polygon, imageRegion, state.color, state.lineThickness); 
}

void rsGraphicsRenderer2DImage::drawText(const rsString &text, int x, int y, int width, int height)
{
  const rsPixelFont *fontToUse = rsGlobalFontInstances::getPixelFontRoundedBoldA16D0();
    // \todo shall become reference parameter or state-variable later

  int kerning = fontToUse->getDefaultKerning();

  rsJustification justification;
  //justification.setJustifiedLeft();
  //justification.setJustifiedRight();
  //justification.setHorizontallyCentered();
  //justification.setJustifiedTop();
  //justification.setVerticallyCentered();
  //justification.setJustifiedBottom();
  //justification.setCentered();

  // adjust x, y according to justification (factor out):
  int xText = x;
  int yText = y;
  int textWidth  = fontToUse->getTextPixelWidth(text, kerning);
  int textHeight = fontToUse->getFontHeight();
  if( justification.isHorizontallyCentered() )
    xText += (width - textWidth) / 2;
  else if( justification.isJustifiedRight() )
    xText += (width - textWidth);
  if( justification.isVerticallyCentered() )
    yText += (height - textHeight) / 2;
  else if( justification.isJustifiedBottom() )
    yText += (height - textHeight);

  // draw top-left justified text with the modified x, y (factor out):
  int xMin = x;
  int yMin = y;
  int xMax = x + width;
  int yMax = y + height;
  for(int i = 0; i < text.getLength(); i++)
  {
    rsUint8 c = text.getElement(i);
    const rsImageGray *glyphImage = fontToUse->getGlyphImage(c/*, state.color*/);
    if( glyphImage != NULL )
    {
      // clip glyph when text extends beyond its box:
      int x = rsMax(xText, xMin);
      int y = rsMax(yText, yMin);
      int w = rsMin(glyphImage->getWidth(),  xMax-x);
      int h = rsMin(glyphImage->getHeight(), yMax-y);

      // draw glyph:
      fillRectangleUsingAlphaMask(x, y, w, h, glyphImage, x-xText, y-yText);
      xText += fontToUse->getGlyphWidth(c) + kerning;
    }
  }

  // \todo maybe handle line-break characters by incrementing y
  // \todo factor out a function drawTopLeftJustifiedTextAt
}

void rsGraphicsRenderer2DImage::fillRectangleWithBilinearGradient(int x, int y, int w, int h,
       const rsColorRGBA &topLeftColor, const rsColorRGBA &topRightColor,                                 
       const rsColorRGBA &bottomLeftColor, const rsColorRGBA &bottomRightColor)
{
  rsColorRGBA cL, cR;
  double scaler = 1.0 / (h-1);    
  rsUint8 weight;
  int iStart = rsMax(y, 0);
  int iEnd   = rsMin(y+h-1, imageRegion.getHeight()-1);
  for(int i = iStart; i <= iEnd; i++)
  {
    weight = (rsUint8) (255.f * scaler * (i-y));
    cL.setAsWeightedAverage(topLeftColor,  bottomLeftColor, weight);
    cR.setAsWeightedAverage(topRightColor, bottomRightColor, weight);
    drawHorizontalGradientLine(y+i, x, x+w-1, cL, cR);
  }
}

void rsGraphicsRenderer2DImage::drawOutlinedPixelRectangle(int x, int y, int w, int h)
{
  if( w <= 0 || h <= 0 )
    return;
  int xMax = imageRegion.getWidth()  - 1;
  int yMax = imageRegion.getHeight() - 1;
  if( rsIsInRange(y, 0, yMax) )
    drawHorizontalPixelLine(y, x, x+w-1);      // top line
  if( rsIsInRange(y+h-1, 0, yMax) )
    drawHorizontalPixelLine(y+h-1, x, x+w-1);  // bottom line
  if( rsIsInRange(x, 0, xMax) )
    drawVerticalPixelLine(x, y, y+h-1);        // left line
  if( rsIsInRange(x+w-1, 0, xMax) )
    drawVerticalPixelLine(x+w-1, y, y+h-1);    // right line
}

void rsGraphicsRenderer2DImage::drawFilledPixelRectangle(int x, int y, int w, int h)
{
  clipToValidRange(x, y, w, h);
  for(int iy = y; iy < y+h; iy++)
  {
    for(int ix = x; ix < x+w; ix++)
      imageRegion.setPixelColor(ix, iy, state.color);
  }
}

void rsGraphicsRenderer2DImage::drawHorizontalPixelLine(int y, int x1, int x2)
{
  rsSortAscending(x1, x2);
  for(int x = rsMax(0, x1); x <= rsMin(x2, imageRegion.getWidth()-1); x++)
    imageRegion.setPixelColor(x, y, state.color);
}

void rsGraphicsRenderer2DImage::drawVerticalPixelLine(int x, int y1, int y2)
{
  rsSortAscending(y1, y2);
  for(int y = rsMax(0, y1); y <= rsMin(y2, imageRegion.getHeight()-1); y++)
    imageRegion.setPixelColor(x, y, state.color);
}

void rsGraphicsRenderer2DImage::drawBresenhamLine(int x1, int y1, int x2, int y2)
{
  bool steep = abs(y2 - y1) > abs(x2 - x1);

  if( steep )
  {
    rsSwap(x1, y1);
    rsSwap(x2, y2);
  }
  if( x1 > x2 )
  {
    rsSwap(x1, x2);
    rsSwap(y1, y2);
  }

  int deltaX = x2 - x1;
  int deltaY = abs(y2 - y1);
  int error  = deltaX / 2;
  int ystep;
  int y = y1;

  if( y1 < y2 )
    ystep = 1;
  else
    ystep = -1;

  for(int x=x1; x<=x2; x++)
  {
    if( steep )
      imageRegion.setPixelColor(y, x, state.color);
    else
      imageRegion.setPixelColor(x, y, state.color);

    error = error - deltaY;
    if( error < 0 )
    {
      y = y + ystep;
      error = error + deltaX;
    }
  }
}

void rsGraphicsRenderer2DImage::drawVerticalGradientLine(int x, int yTop, int yBottom,  
                                                         const rsColorRGBA &topColor, 
                                                         const rsColorRGBA &bottomColor)
{
  if( x < 0 || x > imageRegion.getWidth()-1 )
    return;
  int h = yBottom - yTop;
  double scaler = 1.0 / h;  
  rsColorRGBA c;
  rsUint8 weight;
  int iStart = rsMax(yTop, 0);
  int iEnd   = rsMin(yBottom, imageRegion.getHeight()-1);
  for(int i = iStart; i <= iEnd; i++)
  {
    weight = (rsUint8) (255.f * scaler * (i-yTop));
    c.setAsWeightedAverage(topColor, bottomColor, weight);
    imageRegion.setPixelColor(x, i, c);
  }
}

void rsGraphicsRenderer2DImage::drawHorizontalGradientLine(int y, int xLeft, int xRight,  
                                                           const rsColorRGBA &leftColor, 
                                                           const rsColorRGBA &rightColor)
{
  if( y < 0 || y > imageRegion.getHeight()-1 )
    return;
  int w = xRight - xLeft;
  double scaler = 1.0 / w;  
  rsColorRGBA c;
  rsUint8 weight;
  int iStart = rsMax(xLeft, 0);
  int iEnd   = rsMin(xRight, imageRegion.getWidth()-1);
  for(int i = iStart; i <= iEnd; i++)
  {
    weight = (rsUint8) (255.f * scaler * (i-xLeft));
    c.setAsWeightedAverage(leftColor, rightColor, weight);
    imageRegion.setPixelColor(i, y, c);
  }
}

void rsGraphicsRenderer2DImage::fillRectangleUsingAlphaMask(int x, int y, int w, int h,
                                                            const rsImageGray *mask,                                                        
                                                            int xOffset, int yOffset)
{
  // range checking (write a test and factor this out (maybe into 
  // rsImage::restrictImageCopyRange or something))
  if( x < 0 )
  {
    w += x;
    x  = 0;
  }
  if( y < 0 )
  {
    h += y;
    y  = 0;
  }
  if( x+w > imageRegion.getWidth() )
    w = imageRegion.getWidth() - x;
  if( y+h > imageRegion.getHeight() )
    h = imageRegion.getHeight() - y;
  if( w > imageRegion.getWidth()-xOffset )
    w = imageRegion.getWidth() - xOffset;
  if( h > imageRegion.getHeight()-yOffset )
    h = imageRegion.getHeight() - yOffset;

  // OK - range should be safe now

  for(int iy = 0; iy < h; iy++)
  {
    for(int ix = 0; ix < w; ix++)
    {
      rsUint8 alpha          = mask->getPixelColor(ix+xOffset, iy+yOffset);
      rsColorRGBA oldColor = imageRegion.getPixelColor(x+ix, y+iy);
      rsColorRGBA newColor = state.color;
      rsColorRGBA blendColor;
      blendColor.setAsWeightedAverage(oldColor, newColor, alpha);
      imageRegion.setPixelColor(x+ix, y+iy, blendColor);

      // \todo take into account the alpha value of "newColor" - retain alpha of the "oldColor" but
      // mix-in the new color according to its own alpha value and the alpha value from the 
      // alpha-mask
      // write rsColorRGBA.alphaBlendWith(newColor, alpha) and 
      // image.blendPixelWith(x, y, color, alpha)

      // ....maybe it can be streamlined

    }
  }
}

void rsGraphicsRenderer2DImage::clipToValidRange(int &x, int &y, int &w, int &h)
{
  x = rsClipToRange(x, 0, imageRegion.getWidth()-1);
  w = rsClipToRange(w, 0, imageRegion.getWidth()-x);
  y = rsClipToRange(y, 0, imageRegion.getHeight()-1);
  h = rsClipToRange(h, 0, imageRegion.getHeight()-y);
}

//-------------------------------------------------------------------------------------------------
// class rsGraphicsRenderer2DOpenGL:

rsGraphicsRenderer2DOpenGL::rsGraphicsRenderer2DOpenGL()
{

}

rsGraphicsRenderer2DOpenGL::~rsGraphicsRenderer2DOpenGL()
{

}

// setup:

void rsGraphicsRenderer2DOpenGL::setColor(const rsColorRGBA& newColor)
{
  rsGraphicsRenderer2D::setColor(newColor);
  // glSetColor ...or something
}

// drawing:
/*
void rsGraphicsRenderer2DOpenGL::drawLine(rsPoint2D<double> p1, rsPoint2D<double> p2)
{
  // glDrawLine or something
}
*/
