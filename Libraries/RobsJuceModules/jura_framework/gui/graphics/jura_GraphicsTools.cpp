
// some of these functions need to be updated - setPixelData is not available anymore in the
// juce::Image class. maybe, we have to somehow use Image::BimapData::getPixelPointer or something

void copyImage(juce::Image *sourceImage, juce::Image *targetImage)
{
  int w = targetImage->getWidth();
  int h = targetImage->getHeight();

  // preliminary - at some point, we need to update the commentd code below:
  Graphics g(*targetImage);
  g.fillAll(Colours::red);
  g.setColour(Colours::black);
  g.drawText("jura_GraphicsTools.cpp copyImage() needs to be updated", 0, 0, w, h,
    Justification::centred, false);
  return;

  //// the code below doesn't work anymore because the juce::Image class has changed and does not have
  //// the setPixelData function anymore - we need to check out, how this code needs to be updated:
  //if(  sourceImage->getWidth()  != targetImage->getWidth()
  //  || sourceImage->getHeight() != targetImage->getHeight()
  //  || sourceImage->getFormat() != targetImage->getFormat() )
  //{
  //  Graphics g(*targetImage);
  //  g.fillAll(Colours::red);
  //  g.setColour(Colours::black);
  //  g.drawText("Images incompatible in rojue::copyImage", 0, 0, w, h, Justification::centred, false);
  //  return;
  //}
  //targetImage->setPixelData(0, 0, w, h, sourceImage->getSharedImage()->getPixelData(0, 0), sourceImage->getSharedImage()->getLineStride());
  //  // maybe optimize this by retrieving the BitmapData and doing a memcpy - but comparing with "do-nothing", it doesn't seem to be
  //  // much of a difference, so it probably doesn't matter -- hmm - in another try it used more

  // this does not belong to the code that needs to be updated - it was already commented out:
  //targetImage->clear(Rectangle<int>(0,0,w,h), Colours::white); // seems to use more CPU than setPixelData
  //int numPixels = w*h;
  //int byteSize  = targetImage->getSharedImage()->getPixelStride();
  //int dummy = 0;
}

void normalizedDataToImage(float *data, juce::Image &image, const jura::ColorMap& colorMap)
{
  juce::Image::BitmapData bitmap(image, juce::Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);
  uint8 *p8 = bitmap.getPixelPointer(0, 0);
  uint32 *p = reinterpret_cast<uint32*>(p8);
  for(int i = 0; i < bitmap.height * bitmap.width; i++)
  {
    //// for debug:
    //float tmp = data[i];
    jassert(data[i] >= 0);
    jassert(data[i] <= 1);

    p[i] = colorMap.getColorAsUint32(data[i]);
  }
}

void normalizedDataToImage(float *data, juce::Image &image)
{
  juce::Image::BitmapData bitmap(image, juce::Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);

  // indices for RGBA components in target image:
  int ri, gi, bi, ai;
  colorComponentIndices(image, ri, gi, bi, ai);

  uint8 *p = bitmap.getPixelPointer(0, 0);
  for(int i = 0; i < bitmap.height * bitmap.width; i++)
  {
    //// for debug:
    //float tmp = data[i];
    //jassert(data[i] <= 1);

    p[ri] = p[gi] = p[bi] = (uint8)(255 * data[i]);
    p[ai] = 255;  // full opacity
    p += 4;
  }
}

// these 3 are obsolete - superseded by dataToImage - delete soon
void dataMatrixToPixelBrightnessGray(float **data, uint8 *pixels, int width, int height, uint8 gray)
{
  uint8 *p = pixels;
  for(int i = 0; i < height; i++)     // loop over lines
  {
    for(int j = 0; j < width; j++)    // loop over pixels
    {
      // we assume here, that the alpha channel comes last in the byte order of the pixels
      p[0] = p[1] = p[2] = (uint8) (gray * data[i][j]); // data determines gray-value
      p[3] = 255;                                       // set to full opacity ("alpha")
      p   += 4;                                         // jump to next pixel
    }
  }
  // maybe templatize so it can be used for double also
}
void dataMatrixToPixelBrightness(float **data, uint8 *pixels, int width, int height,
  uint8 red, uint8 green, uint8 blue)
{
  uint8 *p = pixels;
  for(int i = 0; i < height; i++)     // loop over lines
  {
    for(int j = 0; j < width; j++)    // loop over pixels
    {
      // we assume here, that the byte order of the pixels is BGRA (seems to be true on PC)
      p[0] = (uint8) (blue  * data[i][j]);
      p[1] = (uint8) (green * data[i][j]);
      p[2] = (uint8) (red   * data[i][j]);
      p[3] = 255; // full opacity ("alpha")
      p   += 4;
    }
  }
  // todo: write a version of this function that uses a colormap
}
void dataMatrixToImage(float **data, juce::Image &image, uint8 red, uint8 green, uint8 blue)
{
  // We assume here that the size of the image (width and height) matches the dimensions of the
  // data matrix)
  juce::Image::BitmapData bitmap(image, juce::Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);
  uint8 *pixelPointer = bitmap.getPixelPointer(0, 0);
  if(red == green && green == blue)
    dataMatrixToPixelBrightnessGray(data, pixelPointer, bitmap.width, bitmap.height, red);
  else
    dataMatrixToPixelBrightness(data, pixelPointer, bitmap.width, bitmap.height, red, green, blue);
}

int drawBitmapFontText(Graphics &g, int x, int y, const String& textToDraw,
  const BitmapFont* fontToUse, const Colour& colourToUse, int kerning, Justification justification)
{
  //jassert(colourToUse != Colours::yellow); // for test only
  if( kerning == -1 )
    kerning = fontToUse->getDefaultKerning();

  int hFlags = justification.getOnlyHorizontalFlags();
  if( justification.testFlags(hFlags & Justification::horizontallyCentred) )
    x -= fontToUse->getTextPixelWidth(textToDraw, kerning)/2;
  else if( justification.testFlags(hFlags & Justification::right) )
    x -= fontToUse->getTextPixelWidth(textToDraw, kerning);

  int vFlags = justification.getOnlyVerticalFlags();
  if( justification.testFlags(vFlags & Justification::verticallyCentred) )
    y -= fontToUse->getFontAscent()/2;
  else if( justification.testFlags(vFlags & Justification::bottom) )
    y -= fontToUse->getFontAscent();

  // Make sure that g's opacity is set up right, otherwise glyph-images may be drawn transparent
  // regardless of their own opacity values (why is that the case?):
  //g.saveState();                // ToDo: try to avoid this - it may be expensive
  g.setColour(colourToUse);
  g.setOpacity(1.f);

  for(int i = 0; i < textToDraw.length(); ++i)
  {
    juce_wchar c = textToDraw[i];
    const Image* image = fontToUse->getGlyphImage(c, colourToUse);
    if( image != nullptr )
    {
      g.drawImageAt(*image, x, y, false);
      x += fontToUse->getGlyphWidth(c) + kerning;
    }
  }

  //g.restoreState();             // Try to avoid, see comment on g.saveState()
  return x;

  // ToDo:
  // -Try to optimize this. This function is very important for the responsiveness of the GUI.
  //  -What about the g.saveState()/g.restoreState() calls? Can we avoid them? Maybe they are 
  //   costly? Seems like the only state-chaning action is g.setColour(), so perhaps, we could 
  //   just retrieve the old color and set it back when done. Write a benchmark and then try
  //   to do that optimization. The saveState creates an object with new and puts it on some 
  //   stack, so yeah - it could plausibly be an expensive operation. But Graphics has no 
  //   getColour() function, so how can we figure out the current color? I actually tried to just
  //   commenting out the save/restore calls with no immediate ill effects, so maybe we don't
  //   actually need them? It's not good to call them in such a low-level function. What if we
  //   comment the setColour call also? ...this messes up the colors. When calling setOpacity()
  //   instead of setColour(), the font seems to be a bit brighter. When calling both, it's also
  //   brighter. Apparently, the opacity is set to a lower value somewhere on a higher level?
  //   OK, let's keep it at full opacity here. Text is a little bit too bright now, though. 
  //   -> Adjust the color elsewhere (in RWidget possibly)

  // -Graphics does not seem to have a member like getCurrentColour - it has getCurrentFont though
  //  can we nonetheless somehow inquire the current color

  // -maybe try to avoid using g.setColor/setOpacity - client code should do this
  //  -> remove the colourToUse parameter...this will ripple through the codebase
}

void colorComponentIndices(juce::Image& image, int &ri, int &gi, int &bi, int &ai)
{
  // preliminary, seems valid on PC - todo: figure these out in a platform specific way from the
  // passed image:
  ai = 3;
  ri = 2;
  gi = 1;
  bi = 0;
}

void dataToImageOpaqueFloat32x4(float *data, juce::Image &image)
{
  juce::Image::BitmapData bitmap(image, juce::Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);

  // indices for RGBA components in target image:
  int ri, gi, bi, ai;
  colorComponentIndices(image, ri, gi, bi, ai);

  uint8 *p = bitmap.getPixelPointer(0, 0);
  for(int i = 0; i < bitmap.height * bitmap.width; i++)
  {
    p[ri] = (uint8)(255 * data[0]);
    p[gi] = (uint8)(255 * data[1]);
    p[bi] = (uint8)(255 * data[2]);
    p[ai] = 255;  // full opacity
    p    += 4;
    data += 4;
  }
}

Colour getMixedColour(const Colour colour1, const Colour colour2, double weight1, double weight2)
{
  float a1 = colour1.getFloatAlpha();
  float r1 = colour1.getFloatRed();
  float g1 = colour1.getFloatGreen();
  float b1 = colour1.getFloatBlue();

  float a2 = colour2.getFloatAlpha();
  float r2 = colour2.getFloatRed();
  float g2 = colour2.getFloatGreen();
  float b2 = colour2.getFloatBlue();

  uint8 a  = (uint8) (255.f * (weight1*a1 + weight2*a2));
  uint8 r  = (uint8) (255.f * (weight1*r1 + weight2*r2));
  uint8 g  = (uint8) (255.f * (weight1*g1 + weight2*g2));
  uint8 b  = (uint8) (255.f * (weight1*b1 + weight2*b2));

  return Colour(r, g, b, a);
}

void fillRectWithBilinearGradientSlow(Graphics &graphics, int x, int y, int w, int h,
  Colour topLeftColour, Colour topRightColour, Colour bottomLeftColour, Colour bottomRightColour)
{
  // require at least 2 pixels width and height:
  if( w <= 1 || h <= 1 )
    return;

  // check, if all four colours are equal, if so just fill the graphics object with that colour:
  if( topLeftColour == topRightColour
    && topRightColour == bottomLeftColour
    && bottomLeftColour == bottomRightColour )
  {
    graphics.fillAll(topLeftColour);
  }

  // ToDo: use the alpha-channels of the colours as well..., conditional compilation for Mac/Win with
  // different byte orderings for the a,r,g,b values

  // allocate an empty image-object:
  Image *background = new Image(juce::Image::RGB, w, h, false);

  // allocate memory for the raw pixel-data (h lines, each line has w colums,
  // each pixel has 3 bytes for blue, green and red (in that order)):
  uint8 *pixelData  = new uint8[h*w*3];

  int baseAddress;
  // address for a blue-value (is followed by addresses for green and red)

  int baseAddressInc = 3*w;
  // increment for the baseAddress per iteration of the inner loop

  // get the colour-components of the 4 edges as float-values:
  // top-left (abbreviated as t_l):
  uint8 t_l_r = topLeftColour.getRed();
  uint8 t_l_g = topLeftColour.getGreen();
  uint8 t_l_b = topLeftColour.getBlue();
  // top-right (abbreviated as t_r):
  uint8 t_r_r = topRightColour.getRed();
  uint8 t_r_g = topRightColour.getGreen();
  uint8 t_r_b = topRightColour.getBlue();
  // bottom-left (abbreviated as b_l):
  uint8 b_l_r = bottomLeftColour.getRed();
  uint8 b_l_g = bottomLeftColour.getGreen();
  uint8 b_l_b = bottomLeftColour.getBlue();
  // bottom-right (abbreviated as b_r):
  uint8 b_r_r = bottomRightColour.getRed();
  uint8 b_r_g = bottomRightColour.getGreen();
  uint8 b_r_b = bottomRightColour.getBlue();

  // declare variables for the top and bottom line colour-components:
  uint8 t_r, t_g, t_b, b_r, b_g, b_b;

  // declare variables for the colour components at the current pixel (in
  // double and int format:
  uint8   r, g, b;

  int p   = 0; // current x-position (in pixels, depends on i)
  int q   = 0; // current y-position (in pixels, depends on j)
  int wm1 = w-1;
  int hm1 = h-1;

  for(int i=0; i<w; i++)
  {
    // colour components of one pixel in the top-line:
    t_r = t_l_r + (p*(t_r_r-t_l_r))/wm1;
    t_g = t_l_g + (p*(t_r_g-t_l_g))/wm1;
    t_b = t_l_b + (p*(t_r_b-t_l_b))/wm1;

    // colour components of one pixel in the bottom-line:
    b_r = b_l_r + (p*(b_r_r-b_l_r))/wm1;
    b_g = b_l_g + (p*(b_r_g-b_l_g))/wm1;
    b_b = b_l_b + (p*(b_r_b-b_l_b))/wm1;

    // draw the current top-line pixel:
    baseAddress = 3*i;
    pixelData[baseAddress+0] = t_b;  // blue comes first,
    pixelData[baseAddress+1] = t_g;  // green comes second,
    pixelData[baseAddress+2] = t_r;  // red comes third in memory

    // draw the current bottom-line pixel:
    baseAddress = 3*i+w*(h-1)*3;
    pixelData[baseAddress+0] = b_b;
    pixelData[baseAddress+1] = b_g;
    pixelData[baseAddress+2] = b_r;

    // increment the x-position:
    p++;

    // fill the column between 'top' and 'bottom':
    baseAddress = 3*i;
    q = 1; // ? -> 1?
    for(int j=1; j<(h-1); j++)
    {
      r = t_r + (q*(b_r-t_r))/hm1;
      g = t_g + (q*(b_g-t_g))/hm1;
      b = t_b + (q*(b_b-t_b))/hm1;

      baseAddress += baseAddressInc;

      pixelData[baseAddress+0] = b;  // the blue-value
      pixelData[baseAddress+1] = g;  // the green-value
      pixelData[baseAddress+2] = r;  // the red-value

      // increment the y-position:
      q++;
    }
  }

  // copy the generated pixel-data into the image-object:
  //background->setPixelData(0, 0, w, h, pixelData, 3*w);
  graphics.drawText("jura_GraphicsTools.cpp fillRectWithBilinearGradientSlow() needs to be updated", 0, 0, w, h, Justification::centred, false);

  // draw the image into the graphics-object:
  graphics.drawImageAt(*background, x, y);

  // free dynamically allocated memory:
  delete   background;
  delete[] pixelData;
}

void fillRectWithBilinearGradient(Graphics &graphics, int x, int y, int w, int h,
  Colour topLeftColour, Colour topRightColour, Colour bottomLeftColour, Colour bottomRightColour)
{
  // ToDo:
  // -use float instead of double

  //// We create a 2x2 pixel image with the given 4 colors and scale it up. The interpolation will
  //// create the gradient:
  //Image img(Image::ARGB, 2, 2, false);
  //img.setPixelAt(0, 0, topLeftColour);
  //img.setPixelAt(0, 1, bottomLeftColour);
  //img.setPixelAt(1, 0, topRightColour);
  //img.setPixelAt(1, 1, bottomRightColour);
  //graphics.drawImage(img, x, y, w, h, 0, 0, 2, 2, Graphics::mediumResamplingQuality);
  //return;

  // require at least 2 pixels width and height:
  if( w <= 1 || h <= 1 )
    return;

  // check, if all four colours are equal, if so just fill the graphics object with that colour:
  if( topLeftColour == topRightColour
    && topRightColour == bottomLeftColour
    && bottomLeftColour == bottomRightColour )
  {
    graphics.fillAll(topLeftColour);
  }

  // ToDo: use the alpha-channels of the colours as well...

  // allocate image, init pointer-variables:
  Image background(juce::Image::RGB, w, h, false);
  juce::Image::BitmapData bitmap(background, juce::Image::BitmapData::writeOnly);

  //jassert(bitmap.pixelStride == 3);
  int pixelStride = bitmap.pixelStride;

  uint8* pixelData = bitmap.getPixelPointer(0, 0);
  int baseAddress;                        // address for a blue (is followed by reen and red)
  //int baseAddressInc = bitmap.lineStride; // increment for  baseAddress per iteration of inner loop
  int lineStride = bitmap.lineStride; // increment for  baseAddress per iteration of inner loop

  // get the colour-components of the 4 edges as float-values:
  // top-left (abbreviated as t_l):
  double t_l_r = topLeftColour.getFloatRed();
  double t_l_g = topLeftColour.getFloatGreen();
  double t_l_b = topLeftColour.getFloatBlue();
  // top-right (abbreviated as t_r):
  double t_r_r = topRightColour.getFloatRed();
  double t_r_g = topRightColour.getFloatGreen();
  double t_r_b = topRightColour.getFloatBlue();
  // bottom-left (abbreviated as b_l):
  double b_l_r = bottomLeftColour.getFloatRed();
  double b_l_g = bottomLeftColour.getFloatGreen();
  double b_l_b = bottomLeftColour.getFloatBlue();
  // bottom-right (abbreviated as b_r):
  double b_r_r = bottomRightColour.getFloatRed();
  double b_r_g = bottomRightColour.getFloatGreen();
  double b_r_b = bottomRightColour.getFloatBlue();

  // declare variables for the top and bottom line colour-components:
  double t_r, t_g, t_b, b_r, b_g, b_b;

  // declare variables for the colour components at the current pixel (in
  // double and int format:
  double r, g, b;
  uint8  r_int, g_int, b_int;

  double p = 0.0; // current relative x-position (from 0.0...1.0, depends on i)
  double q = 0.0; // current relative y-position (from 0.0...1.0, depends on j)
  double pInc = 1.0/(double)(w-1); // increment per iteration for p
  double qInc = 1.0/(double)(h-1); // increment per (inner) iteration for q

  for(int i = 0; i < w; i++)
  {
    // colour components of one pixel in the top-line:
    t_r    = (1.0-p)*t_l_r + p*t_r_r;
    t_g    = (1.0-p)*t_l_g + p*t_r_g;
    t_b    = (1.0-p)*t_l_b + p*t_r_b;

    // colour components of one pixel in the bottom-line:
    b_r    = (1.0-p)*b_l_r + p*b_r_r;
    b_g    = (1.0-p)*b_l_g + p*b_r_g;
    b_b    = (1.0-p)*b_l_b + p*b_r_b;


    // draw the current bottom-line pixel:
    //baseAddress = pixelStride*i + w*(h-1)*pixelStride; 
         // old: 3*i+w*(h-1)*3;
         // maybe needs to be pixelStride*i + lineStride*(h-1)?
    baseAddress = i * pixelStride + (h-1) * lineStride; 
    r_int = (uint8) (255 * b_r);
    g_int = (uint8) (255 * b_g);
    b_int = (uint8) (255 * b_b);
    setPixelRGB(pixelData+baseAddress, r_int, g_int, b_int);

    // draw the current top-line pixel:
    baseAddress = i * pixelStride; // old: 3*i; ...maybe we need lineStride?
    r_int = (uint8) (255 * t_r);
    g_int = (uint8) (255 * t_g);
    b_int = (uint8) (255 * t_b);
    setPixelRGB(pixelData+baseAddress, r_int, g_int, b_int);

    // increment the relative x-position:
    p += pInc;

    // fill the column between 'top' and 'bottom':
    //baseAddress = pixelStride*i; // old: 3*i;
    q = 0.0;
    for(int j = 1; j < (h-1); j++)
    {
      r = (1.0-q)*t_r + q*b_r;
      g = (1.0-q)*t_g + q*b_g;
      b = (1.0-q)*t_b + q*b_b;

      r_int = (uint8) (255 * r);
      g_int = (uint8) (255 * g);
      b_int = (uint8) (255 * b);

      baseAddress += lineStride;
      setPixelRGB(pixelData+baseAddress, r_int, g_int, b_int);
      q += qInc; // increment the relative y-position:
    }
  }

  // set the alpha channel to 255, if the image has one:
  if(pixelStride == 4)
    for(int i = 0; i < h; i++)
      for(int j = 0; j < w; j++)
        pixelData[lineStride*i + 4*j + 3] = 255;

  graphics.drawImageAt(background, x, y);  // draw the image into the graphics-object

  //juce::String debugInfo = "w=" + String(w) + ", h=" + String(h);
  //graphics.drawText(debugInfo, 0, 0, w, h, Justification::centred, false);

  // Notes:
  // -When fillig a big rectangle, artifacts can be seen: areas of same color separated by 
  //  hyperbolas. Maybe using floating point arithemtic could avodi this? Try it!
}

void fillRectWithBilinearGradient(Graphics &graphics, Rectangle<int> r,  Colour topLeftColour,
  Colour topRightColour, Colour bottomLeftColour, Colour bottomRightColour)
{
  fillRectWithBilinearGradient(graphics,
    r.getX(), r.getY(), r.getWidth(), r.getHeight(),
    topLeftColour, topRightColour, bottomLeftColour, bottomRightColour);
}

void fillRectWithBilinearGrayScaleGradient(Graphics &graphics, int x, int y, int w, int h,
  float topLeftWhite, float topRightWhite, float bottomLeftWhite, float bottomRightWhite)
{
  // allocate an empty image-object:
  Image *background = new Image(juce::Image::RGB, w, h, false);

  // allocate memory for the raw pixel-data (h lines, each line has w colums,
  // each pixel has 3 bytes for blue, green and red (in that order)):
  uint8 *pixelData  = new uint8[h*w*3];

  int    baseAddress;
  // address for a blue-value (is followed by addresses for green and red)

  int    baseAddressInc = 3*w;
  // increment for the baseAddress per iteration of the inner loop

  // get the white values of the 4 edges as float-values:
  double t_l_w = topLeftWhite;
  double t_r_w = topRightWhite;
  double b_l_w = bottomLeftWhite;
  double b_r_w = bottomRightWhite;


  // declare variables for the top and bottom line white-values:
  double t_w, b_w;

  // declare variables for the colour components at the current pixel (in
  // double and int format:
  double white;
  uint8  white_int;

  double p = 0.0; // current relative x-position (from 0.0...1.0, depends on i)
  double q = 0.0; // current relative y-position (from 0.0...1.0, depends on j)
  double pInc = 1.0/(double)(w-1); // increment per iteration for p
  double qInc = 1.0/(double)(h-1); // increment per (inner) iteration for q

  for(int i=0; i<w; i++)
  {
    // white-value of one pixel in the top-line:
    t_w    = (1.0-p)*t_l_w + p*t_r_w;

    // white-value of one pixel in the bottom-line:
    b_w    = (1.0-p)*b_l_w + p*b_r_w;

    // draw the current top-line pixel:
    baseAddress = 3*i;
    white_int = (uint8) (255 * t_w);
    pixelData[baseAddress+0] = white_int;  // blue comes first,
    pixelData[baseAddress+1] = white_int;  // green comes second,
    pixelData[baseAddress+2] = white_int;  // red comes third in memory

    // draw the current bottom-line pixel:
    baseAddress = 3*i+w*(h-1)*3;
    white_int = (uint8) (255 * b_w);
    pixelData[baseAddress+0] = white_int;
    pixelData[baseAddress+1] = white_int;
    pixelData[baseAddress+2] = white_int;

    // increment the relative x-position:
    p += pInc;

    // fill the column between 'top' and 'bottom':
    baseAddress = 3*i;
    q = 0.0;
    for(int j=1; j<(h-1); j++)
    {
      white     = (1.0-q)*t_w + q*b_w;
      white_int = (uint8) (255 * white);

      baseAddress += baseAddressInc;

      pixelData[baseAddress+0] = white_int;  // the blue-value
      pixelData[baseAddress+1] = white_int;  // the green-value
      pixelData[baseAddress+2] = white_int;  // the red-value

      q += qInc; // increment the relative y-position:
    }
  }

  // copy the generated pixel-data into the image-object:
  //background->setPixelData(0, 0, w, h, pixelData, 3*w);
  graphics.drawText("jura_GraphicsTools.cpp fillRectWithBilinearGrayScaleGradient() needs to be updated", 0, 0, w, h, Justification::centred, false);

  // draw the image into the graphics-object:
  graphics.drawImageAt(*background, x, y);

  // free dynamically allocated memory:
  delete   background;
  delete[] pixelData;
}

void fillRectWithDefaultBackground(Graphics &g, int x, int y, int w, int h)
{
  fillRectWithBilinearGradient(g, x, y, w, h,
    Colours::white, Colour(190,200,225), Colour(190,200,225), Colours::white);
}

void fillRectWithDefaultBackground(Graphics &g, Rectangle<int> r)
{
  fillRectWithDefaultBackground(g, r.getX(), r.getY(), r.getWidth(), r.getHeight());
}

void drawTriangle(Graphics &g, float x1, float y1, float x2, float y2, float x3, float y3,
  bool fill)
{
  // make a path from the 3 given points:
  Path	trianglePath;
  trianglePath.startNewSubPath(x1, y1);
  trianglePath.lineTo(x2, y2);
  trianglePath.lineTo(x3, y3);
  trianglePath.lineTo(x1, y1);

  // draw or fill the path:
  if( fill == true )
    g.fillPath(trianglePath);
  else
    g.strokePath(trianglePath, PathStrokeType(1.0));
}

void drawBlockDiagramPlus(Graphics &g, float x, float y, float w, float h, float thickness)
{
  g.drawEllipse(x, y, w, h, thickness);
  g.drawLine(x, y+0.5f*h, x+w, y+0.5f*h, thickness);
  g.drawLine(x+0.5f*w, y, x+0.5f*w, y+h, thickness);
}

void fitLineToRectangle(double &x1, double &y1, double &x2, double &y2, double xMin, double yMin,
  double xMax, double yMax)
{
  // catch some special cases:
  if( x1 == x2 ) // vertical line, infinite slope
  {
    y1 = yMin;
    y2 = yMax;
    return;
  }
  if( y1 == y2 ) // horizontal line, zero slope
  {
    x1 = xMin;
    x2 = xMax;
    return;
  }

  // calculate slope 'a' and constant term 'b' for the line equation y = a*x + b:
  double a = (y2-y1)/(x2-x1);
  double b = y1-a*x1;

  // calculate x- and y-values at the rectangle's boundaries:
  double yAtMinX = a*xMin+b;
  double yAtMaxX = a*xMax+b;
  double xAtMinY = (yMin-b)/a;
  double xAtMaxY = (yMax-b)/a;

  if( yAtMinX > yMin && yAtMinX < yMax )
    // line intercepts left boundary
  {
    x1 = xMin;
    y1 = a*xMin+b;
    if( xAtMaxY > xMin && xAtMaxY < xMax )
      // line intercepts left and top boundary (chops off the top-left corner)
    {
      x2 = (yMax-b)/a;
      y2 = yMax;
    }
    else if( xAtMinY > xMin && xAtMinY < xMax )
      // line intercepts left and bottom boundary (chops off the bottom-left corner)
    {
      x2 = (yMin-b)/a;
      y2 = yMin;
    }
    else
      // line intercepts right boundary (divides the rectangle into top and bottom)
    {
      x2 = xMax;
      y2 = a*xMax+b;
    }
  }
  else if( yAtMaxX > yMin && yAtMaxX < yMax )
    // line intercepts right boundary
  {
    x2 = xMax;
    y2 = a*xMax+b;
    if( xAtMaxY > xMin && xAtMaxY < xMax )
      // line intercepts right and top boundary (chops off top right corner)
    {
      x1 = (yMax-b)/a;
      y1 = yMax;
    }
    else if( xAtMinY > xMin && xAtMinY < xMax )
      // line intercepts right and bottom boundary (chops off bottom right corner)
    {
      x1 = (yMin-b)/a;
      y1 = yMin;
    }
  }
  else if( xAtMaxY > xMin && xAtMaxY < xMax && xAtMinY > xMin && xAtMinY < xMax )
    // line intercepts top and bottom boundary (divides the rectangle into left and right)
  {
    x1 = (yMin-b)/a;
    y1 = yMin;
    x2 = (yMax-b)/a;
    y2 = yMax;
  }

}

void clipLineToRectangle(double &x1, double &y1, double &x2, double &y2, double xMin, double yMin,
  double xMax, double yMax)
{
  bool   firstIsInside  = false;
  bool   secondIsInside = false;
  double xIn, yIn, xOut, yOut;

  if( x1 >= xMin && x1 <= xMax && y1 >= yMin && y1 <= yMax )
    firstIsInside = true;
  if( x2 >= xMin && x2 <= xMax && y2 >= yMin && y2 <= yMax )
    secondIsInside = true;

  if( firstIsInside && secondIsInside )
    return;
  if( firstIsInside == false && secondIsInside == false )
  {
    fitLineToRectangle(x1, y1, x2, y2, xMin, yMin, xMax, yMax);
    return;
  }

  if( firstIsInside && secondIsInside == false )
  {
    xIn  = x1;
    yIn  = y1;
    xOut = x2;
    yOut = y2;
  }
  else if( secondIsInside && firstIsInside == false )
  {
    xIn  = x2;
    yIn  = y2;
    xOut = x1;
    yOut = y1;
  }

  // we have the case that one point is inside and the other one outside the rectangle and we have
  // already determined which is which - now we must find out the clipped outside value

  if( x1 == x2 ) // vertical line, infinite slope
  {
    if( yOut > yMax )
    {
      yOut = yMax;
      if( firstIsInside && secondIsInside == false )
      {
        y2 = yOut;
        return;
      }
      else if( secondIsInside && firstIsInside == false )
      {
        y1 = yOut;
        return;
      }
    }
    else if( yOut < yMin )
    {
      yOut = yMin;
      if( firstIsInside && secondIsInside == false )
      {
        y2 = yOut;
        return;
      }
      else if( secondIsInside && firstIsInside == false )
      {
        y1 = yOut;
        return;
      }
    }

  }
  if( y1 == y2 ) // horizontal line, zero slope
  {
    // something to do? i dont thibk so, but...

  }

  // calculate y-values at the rectangle's left and right boundaries:
  double a = (yOut-yIn)/(xOut-xIn);
  double b = yIn-a*xIn;
  double yAtMinX = a*xMin+b;
  double yAtMaxX = a*xMax+b;
  double xAtMinY = (yMin-b)/a;
  double xAtMaxY = (yMax-b)/a;

  if( xOut >= xMax ) // line intercepts right, top or bottom boundary
  {
    if( yAtMaxX > yMax )      // line intercepts top boundary
    {
      xOut = xAtMaxY;
      yOut = yMax;
    }
    else if( yAtMaxX < yMin ) // line intercepts bottom boundary
    {
      xOut = xAtMinY;
      yOut = yMin;
    }
    else                      // line intercepts right boundary
    {
      xOut = xMax;
      yOut = yAtMaxX;
    }
  }
  else if( xOut <= xMin ) // line intercepts left, top or bottom boundary
  {
    if( yAtMinX > yMax )      // line intercepts top boundary
    {
      xOut = xAtMaxY;
      yOut = yMax;
    }
    else if( yAtMinX < yMin ) // line intercepts bottom boundary
    {
      xOut = xAtMinY;
      yOut = yMin;
    }
    else                      // line intercepts left boundary
    {
      xOut = xMin;
      yOut = yAtMinX;
    }
  }
  else if( yOut > yMax ) // line intercepts top boundary
  {
    xOut = xAtMaxY;
    yOut = yMax;
  }
  else if( yOut > yMax ) // line intercepts bottom boundary
  {
    xOut = xAtMinY;
    yOut = yMin;
  }


  if( firstIsInside && secondIsInside == false )
  {
    x2 = xOut;
    y2 = yOut;
  }
  else if( secondIsInside && firstIsInside == false )
  {
    x1 = xOut;
    y1 = yOut;
  }
}

//=================================================================================================
// coordinate system drawing:

void setupCoordinateMapper(RAPT::rsCoordinateMapper2D<double>& mapper, double w, double h)
{
  //mapper.setOutputRange(0.5, w-0.5, h-0.5, 0.5); // inf/nan wehn w or h <= 1
  mapper.setOutputRange(0.5, jmax(w-0.5, 1.0), jmax(h-0.5, 1.0), 0.5);
}

void setupCoordinateMapper(RAPT::rsCoordinateMapper2D<double>& mapper, const Component* cmp)
{
  int w = cmp->getWidth();
  int h = cmp->getHeight();
  setupCoordinateMapper(mapper, w, h);
}

void setupCoordinateMapper(RAPT::rsCoordinateMapper2D<double>& mapper, const Image* img)
{
  int w = img->getWidth();
  int h = img->getHeight();
  setupCoordinateMapper(mapper, w, h);
}

void setupCoordinateMapper(RAPT::rsCoordinateMapper2D<double>& mapper, const XmlElement* svg)
{
  double w = svg->getDoubleAttribute("width", 0);
  double h = svg->getDoubleAttribute("height", 0);
  jassert(w != 0 && h != 0); // svg must have width and height attributes
  setupCoordinateMapper(mapper, w, h);
}

void drawBitmapText(Graphics &g, const String &text, double x, double y, double w, double h,
  BitmapFont const* font, Justification justification, Colour color)
{
  int hFlags = justification.getOnlyHorizontalFlags();
  if(justification.testFlags(hFlags & Justification::horizontallyCentred))
    x = (x+x+w)/2.0; // average between x and x+w
  else if(justification.testFlags(hFlags & Justification::right))
    x = x+w;

  int vFlags = justification.getOnlyVerticalFlags();
  if(justification.testFlags(vFlags & Justification::verticallyCentred))
    y = (y+y+h)/2.0;
  else if(justification.testFlags(vFlags & Justification::bottom))
    y = y+h;

  int xInt = roundToInt(x);
  int yInt = roundToInt(y);
  drawBitmapFontText(g, xInt, yInt, text, font, color, -1, justification);
}

void initGridDrawing(const RAPT::rsCoordinateMapper<double>& mapper, double spacing,
  double& start, double& delta)
{
  // todo: for log-scaling, the current formulas lead to frequencies that are powers of two when
  // the freq-axis is log scaled with octave spacing, like ...,512,1024,2048,... but we want
  // ...500,1000,2000,...
  // i think, we need to define a grid anchor, like a = 1000, then compute the smallest integer k,
  // such that a * b^k <= f_min, our start is then a * b^k for that found k

  if(mapper.isLogScaled()) {
    double k = ceil(RAPT::rsLogB(mapper.getInMin(), spacing)) + 1;
    start = pow(spacing, k);
    start = mapper.map(start);
    delta = mapper.map(spacing) - mapper.map(1); }
  else {
    start = spacing * floor(mapper.getInMin() / spacing + 0.5); // use ceil
    start = mapper.map(start);
    delta = mapper.map(spacing) - mapper.map(0); }
}

void drawHorizontalGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness)
{
  float xL = (float) mapper.mapX(mapper.getInMinX());
  float xR = (float) mapper.mapX(mapper.getInMaxX());
  double y, dy; 
  initGridDrawing(mapper.mapperY, spacing, y, dy);
  while(y > mapper.getOutMaxY()) {                     // > because pixel y-axis goes down
    g.drawLine(xL, float(y), xR, float(y), thickness); // juce doc says, it's better to use fillRect
    y += dy; }
}

void drawAxisValuesY(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, double xPos, juce::String (*yToString) (double y), Colour textColor)
{
  float  x = (float) mapper.mapX(xPos); // in pixels
  double xt = x+4;                      // text-position
  Justification just(Justification::centredLeft);
  if(x > 10) { xt = x-25; just = Justification::centredRight; }
  double y, dy;
  initGridDrawing(mapper.mapperY, spacing, y, dy);
  while(y > mapper.getOutMaxY()) {                     // loop over the tics
    g.drawLine(x-4.f, float(y), x+4.f, float(y), 1.f); // juce doc says, it's better to use fillRect
    drawBitmapText(g, yToString(mapper.unmapY(y)), xt, y-10, 20, 20, &normalFont7px, just, 
      textColor);
    y += dy; }
}

void drawVerticalGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness)
{
  float yB = (float) mapper.mapY(mapper.getInMinY());
  float yT = (float) mapper.mapY(mapper.getInMaxY());
  double x, dx; 
  initGridDrawing(mapper.mapperX, spacing, x, dx);
  while(x < mapper.getOutMaxX()) {
    g.drawLine(float(x), yB, float(x), yT, thickness);
    x += dx; }
}

void drawAxisValuesX(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, double yPos, juce::String (*xToString) (double x), Colour textColor)
{
  float  y = (float) mapper.mapY(yPos); // in pixels
  double x, dx;
  initGridDrawing(mapper.mapperX, spacing, x, dx);
  double yt = y+4; 
  Justification just(Justification::centred);
  if(y > mapper.getOutMinY() - 10) 
    yt = y-20;
  while(x < mapper.getOutMaxX()) {
    g.drawLine(float(x), y-4.f, float(x), y+4.f, 1.f);
    drawBitmapText(g, xToString(mapper.unmapX(x)), x-32, yt, 64, 16, &normalFont7px, just, 
      textColor);
    x += dx; }
}

double getMaxRadius(const RAPT::rsCoordinateMapper2D<double>& mapper)
{
  double x = jmax(fabs(mapper.getInMinX()), fabs(mapper.getInMaxX()));
  double y = jmax(fabs(mapper.getInMinY()), fabs(mapper.getInMaxY()));
  return sqrt(x*x + y*y);
}

void drawRadialGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness)
{
  int i = 1; 
  double radius = spacing; 
  double maxRadius = getMaxRadius(mapper);
  while(radius <= maxRadius) {
    float xL = (float) mapper.mapX(-radius); float xR = (float) mapper.mapX(+radius);
    float yB = (float) mapper.mapY(-radius); float yT = (float) mapper.mapY(+radius);
    g.drawEllipse(xL, yT, xR-xL, yB-yT, thickness); // circle may deform into ellipse
    i++; radius = spacing * (double) i; }
}

void angularLineEndPoints(double angle, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double& startX, double& startY, double& endX, double& endY)
{
  endX = cos(angle);
  endY = sin(angle);
  startX = -endX;
  startY = -endY;
  fitLineToRectangle(startX, startY, endX, endY, 
    mapper.getInMinX(), mapper.getInMinY(), mapper.getInMaxX(), mapper.getInMaxY() );
  startX = mapper.mapX(startX);
  endX   = mapper.mapX(endX);
  startY = mapper.mapY(startY);
  endY   = mapper.mapY(endY);
}

void drawAngularGrid(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness)
{
  spacing = spacing*(PI/180.0); // degree-to-radiant
  double angle = 0.0; int i = 0;
  double startX, endX, startY, endY;
  while(angle <= PI) {
    angularLineEndPoints(angle, mapper, startX, startY, endX, endY);
    g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, thickness);
    i++; angle = spacing * (double) i; }
}

void pixelCoordsForAxisX(const RAPT::rsCoordinateMapper2D<double>& mapper, float& xs, float& xe, 
  float& y, float& yt) 
{
  y  = (float) mapper.mapY(y); // y is input in model- and output in pixel coords
  xs = (float) mapper.getOutMinX();
  xe = (float) mapper.getOutMaxX();
  if(y > mapper.getOutMinY() - 10)  
    yt = y - 20; // label text above x-axis
  else
    yt = y + 4;  // label text below x-axis
}

void drawAxisX(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, double yPos, 
  const juce::String& label, Colour labelColor)
{
  float y = (float) yPos;
  float xs, xe, yt; // start/end x, y for text
  pixelCoordsForAxisX(mapper, xs, xe, y, yt);
  g.drawArrow(Line<float>(xs, y, xe, y), 1.0, 6.0, 6.0);
  static const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  drawBitmapText(g, label, xe-100, yt, 96, 16, font, Justification::centredRight, labelColor);
}

void pixelCoordsForAxisY(const RAPT::rsCoordinateMapper2D<double>& mapper, float& ys, float& ye, 
  float& x, float& xt) 
{
  x  = (float) mapper.mapX(x);
  ys = (float) mapper.getOutMinY();
  ye = (float) mapper.getOutMaxY();
  if(x < 10)  
    xt = x + 8;  // label text right to y-axis
  else
    xt = x - 8;  // label text left to y-axis
}

void drawAxisY(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& mapper, double xPos, 
  const juce::String& label, Colour labelColor)
{
  float x = (float) xPos;
  float ys, ye, xt;
  pixelCoordsForAxisY(mapper, ys, ye, x, xt);
  g.drawArrow(Line<float>(x, ys, x, ye), 1.0, 6.0, 6.0);

  // label positioning may not yet be ideal in all circumstances:
  int width = 100;
  Justification just(Justification::centredRight);
  if(x < width+8)
    just = Justification::centredLeft;
  else
    xt -= width;

  static const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  drawBitmapText(g, label, xt, ye, width, 16, font, just, labelColor);
}

//=================================================================================================
// coordinate system drawing for svg export:

void addLineToSvgDrawing(XmlElement* svg, float x1, float y1, float x2, float y2,
  float thickness, Colour color, bool withArrowHead)
{
  XmlElement* line = new XmlElement(String("line"));
  line->setAttribute(String("x1"), x1);
  line->setAttribute(String("y1"), y1);
  if( withArrowHead == true && y1 == y2 )
    line->setAttribute(String("x2"), x2-8);
  else
    line->setAttribute(String("x2"), x2);

  if( withArrowHead == true && x1 == x2 )
    line->setAttribute(String("y2"), y2+8);
  else
    line->setAttribute(String("y2"), y2);

  line->setAttribute(String("style"), String("stroke-width: ") + String(thickness) + 
    String("; stroke: #") + color.toString().substring(2) + String(";") );
  svg->addChildElement(line);

  if( withArrowHead == true )
  {
    XmlElement* triangle = new XmlElement(String("path"));

    if( y1 == y2 ) // this is a horizontal rightward arrow 
    {
      triangle->setAttribute(String("d"), 
        String("M ")   + String(x2-8) + String(" ") + String(y2-4) + 
        String(", L ") + String(x2-8) + String(" ") + String(y2+4) + 
        String(", L ") + String(x2)   + String(" ") + String(y2)   + 
        String(", Z") );
      triangle->setAttribute(String("style"), String("stroke: none, fill: #") 
        + color.toString().substring(2) + String(";") );
    }
    else if( x1 == x2 ) // this is an upward verzical rightward arrow 
    {
      triangle->setAttribute(String("d"), 
        String("M ")   + String(x2-4) + String(" ") + String(y2+8) + 
        String(", L ") + String(x2+4) + String(" ") + String(y2+8) + 
        String(", L ") + String(x2)   + String(" ") + String(y2)   + 
        String(", Z") );
      triangle->setAttribute(String("style"), String("stroke: none, fill: #") 
        + color.toString().substring(2) + String(";") );
    }
    svg->addChildElement(triangle);
  }
}

void addTextToSvgDrawing(XmlElement* svg, juce::String theText, float x, float y,
  Justification justification, Colour color)
{
  XmlElement* textContainer = new XmlElement(String("text"));
  XmlElement* text = XmlElement::createTextElement(theText);

  String jString = String();
  if( justification.getFlags() == Justification::centredLeft )
    jString = String("start");
  else if( justification.getFlags() == Justification::centred )
    jString = String("middle");
  else if( justification.getFlags() == Justification::centredRight )
    jString = String("end");

  textContainer->setAttribute(String("x"), x);
  textContainer->setAttribute(String("y"), y);
  textContainer->setAttribute(String("style"), String("font-family: sans-serif;") +  
    String(" font-size: 12px;") + String(" stroke: none;") + String(" fill: black;") +
    String(" text-anchor: ") + jString + String(";") );
  textContainer->addChildElement(text);
  svg->addChildElement(textContainer);
}

void drawHorizontalGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness, Colour colour)
{
  float xL = (float) mapper.mapX(mapper.getInMinX());
  float xR = (float) mapper.mapX(mapper.getInMaxX());
  double y, dy; 
  initGridDrawing(mapper.mapperY, spacing, y, dy);
  String gridPathDataString;

  while(y > mapper.getOutMaxY()) {
    gridPathDataString += String("M ") + String(xL) + String(" ") + String(y) + String(" ");
    gridPathDataString += String("L ") + String(xR) + String(" ") + String(y) + String(" ");
    y += dy; }

  XmlElement* gridPath = new XmlElement(String("path"));
  gridPath->setAttribute(String("d"), gridPathDataString);
  gridPath->setAttribute(String("style"), String("stroke-width: ") + String(thickness) + 
    String("; stroke: #") + colour.toString().substring(2) + String(";") );
  svg->addChildElement(gridPath);
}

void drawAxisValuesY(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, double xPos, juce::String (*yToString) (double y), Colour color)
{
  float  x = (float) mapper.mapX(xPos); // in pixels
  double xt = x+8;                      // text-position...why not float?
  Justification just(Justification::centredLeft);
  if(x > 10) { xt = x-8; just = Justification::centredRight; }
  double y, dy;
  initGridDrawing(mapper.mapperY, spacing, y, dy);
  while(y > mapper.getOutMaxY()) {
    addLineToSvgDrawing(svg, x-4.f, float(y), x+4.f, float(y), 1.0, color, false);
    addTextToSvgDrawing(svg, yToString(mapper.unmapY(y)), float(xt), float(y)+4.f, just, color);
    y += dy; }
}

void drawVerticalGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness, Colour colour)
{
  float yB = (float) mapper.mapY(mapper.getInMinY());
  float yT = (float) mapper.mapY(mapper.getInMaxY());
  double x, dx; 
  initGridDrawing(mapper.mapperX, spacing, x, dx);
  String gridPathDataString;

  while(x < mapper.getOutMaxX()) {
    gridPathDataString += String("M ") + String(x) + String(" ") + String(yB) + String(" ");
    gridPathDataString += String("L ") + String(x) + String(" ") + String(yT) + String(" ");
    x += dx; }

  // factor out (duplicated from drawHorizontalGrid):
  XmlElement* gridPath = new XmlElement(String("path"));
  gridPath->setAttribute(String("d"), gridPathDataString);
  gridPath->setAttribute(String("style"), String("stroke-width: ") + String(thickness) + 
    String("; stroke: #") + colour.toString().substring(2) + String(";") );
  svg->addChildElement(gridPath);
}

void drawAxisValuesX(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, double yPos, juce::String (*xToString) (double x), Colour color)
{
  float y = (float) mapper.mapY(yPos); // in pixels
  double x, dx;
  initGridDrawing(mapper.mapperX, spacing, x, dx);
  double yt = y+12; 
  Justification just(Justification::centred);
  if(y > mapper.getOutMinY() - 10) 
    yt = y-28;
  while(x < mapper.getOutMaxX()) {
    addLineToSvgDrawing(svg, float(x), y-4.f, float(x), y+4.f, 1.f, color, false);
    addTextToSvgDrawing(svg, xToString(mapper.unmapX(x)), float(x), float(yt), just, color);
    x += dx; }
}

void drawRadialGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness, Colour color)
{
  int i = 1;
  double maxRadius = getMaxRadius(mapper); double radius = spacing;
  float xC = (float) mapper.mapX(0); float yC = (float) mapper.mapY(0);
  while(radius <= maxRadius) {
    float xL = (float) mapper.mapX(-radius); float xR = (float) mapper.mapX(+radius);
    float yB = (float) mapper.mapY(-radius); float yT = (float) mapper.mapY(+radius);
    XmlElement* ellipse = new XmlElement("ellipse"); // circle may deform into ellipse
    ellipse->setAttribute("cx", xC);
    ellipse->setAttribute("cy", yC);
    ellipse->setAttribute("rx", 0.5f*(xR-xL));
    ellipse->setAttribute("ry", 0.5f*(yB-yT));
    ellipse->setAttribute("style", "stroke-width: " + String(thickness) + 
      "; stroke: #" + color.toString().substring(2) + ";" + "fill: none;" );
    svg->addChildElement(ellipse);
    i++; radius = spacing * (double) i; }
}

void drawAngularGrid(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double spacing, float thickness, Colour color)
{
  spacing = spacing*(PI/180.0); // degree-to-radiant
  double angle = 0.0; int i = 0;
  double startX, endX, startY, endY;
  String pathString;
  while(angle <= PI) 
  {
    angularLineEndPoints(angle, mapper, startX, startY, endX, endY);
    pathString += String("M ") + String(startX) + String(" ") + String(startY) + String(" ");
    pathString += String("L ") + String(endX) + String(" ")   + String(endY) + String(" ");
    i++; angle = spacing * (double) i; 
  }
  XmlElement* gridPath = new XmlElement(String("path"));
  gridPath->setAttribute(String("d"), pathString);
  gridPath->setAttribute(String("style"), String("stroke-width: ") + String(thickness) 
    + String("; stroke: #") + color.toString().substring(2) + String(";") );
  svg->addChildElement(gridPath);
}

void drawAxisX(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double yPos, const juce::String& label, Colour color)
{
  float y = (float) yPos;
  float xs, xe, yt; // start/end x, y for text
  pixelCoordsForAxisX(mapper, xs, xe, y, yt);
  addLineToSvgDrawing(svg, xs, y, xe, y, 2.0, color, true);
  addTextToSvgDrawing(svg, label, xe-4, yt+8, Justification::centredRight, color);
}

void drawAxisY(XmlElement* svg, const RAPT::rsCoordinateMapper2D<double>& mapper,
  double xPos, const juce::String& label, Colour color)
{
  float x = (float) xPos;
  float ys, ye, xt;
  pixelCoordsForAxisY(mapper, ys, ye, x, xt);
  addLineToSvgDrawing(svg, x, ys, x, ye, 2.0, color, true);

  addTextToSvgDrawing(svg, label, xt-8, ye+12, Justification::centredRight, color); 
    // still wrongly positioned
}