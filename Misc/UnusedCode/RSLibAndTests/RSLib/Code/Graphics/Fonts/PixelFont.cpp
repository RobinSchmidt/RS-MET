using namespace RSLib;

// construction/destruction:

rsPixelFont::rsPixelFont() 
{
  ascent         = 0;
  descent        = 0;
  defaultKerning = 1;

  for(int g=0; g<numGlyphs; g++)
  {
    glyphImages[g] = NULL;
    glyphWidths[g]  = 0;
  }

  createGlyphImages();
}

rsPixelFont::~rsPixelFont()
{
  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphImages[g] != NULL )
    {
      delete glyphImages[g];
      glyphImages[g] = NULL;
    }
  }
}

// inquiry:

int rsPixelFont::getTextPixelWidth(const rsString &text, int kerning) const
{
  int   result = 0;
  rsUint8 c; 
  for(int i = 0; i < text.getLength(); i++)
  {
    c       = text.getElement(i);
    result += getGlyphWidth(c) + kerning;
  }
  return result - kerning; // we added one too much in the loop
}

int rsPixelFont::glyphIndexToX(const rsString &text, int index, int kerning) const
{
  rsAssert( index < text.getLength() ); // index is out of range
  if( index >= text.getLength() )
    return 0;

  int   result = 0;
  rsUint8 c; 
  for(int i = 0; i < index; i++)
  {
    c       = text.getElement(i);
    result += getGlyphWidth(c) + kerning;
  }
  return result; 
}

int rsPixelFont::xToGlyphIndex(const rsString &text, int xToFindGlyphIndexAt, int kerning) const
{
  int x     = 0;
  int index = 0;
  rsUint8 c; 
  while( x < xToFindGlyphIndexAt && index < text.getLength() )
  {
    c  = text.getElement(index);
    x += getGlyphWidth(c) + kerning;
    index++;
  }
  return index;
}

// internal functions:

void rsPixelFont::createGlyphImages()
{

}
