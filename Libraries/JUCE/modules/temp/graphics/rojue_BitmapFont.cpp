#include "rojue_BitmapFont.h"
using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

//BitmapFont* BitmapFont::instance = NULL;

BitmapFont::BitmapFont() 
{
  ascent         = 0;
  descent        = 0;
  defaultKerning = 1;

  for(int g=0; g<numGlyphs; g++)
  {
    glyphBitmaps[g] = NULL;
    glyphWidths[g]  = 0;
  }

  createGlyphBitmaps();
}

BitmapFont::~BitmapFont()
{
  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphBitmaps[g] != NULL )
    {
      delete glyphBitmaps[g];
      glyphBitmaps[g] = NULL;
    }
  }
}



//-------------------------------------------------------------------------------------------------
// setup:






//-------------------------------------------------------------------------------------------------
// inquiry:

int BitmapFont::getTextPixelWidth(const juce::String &text, int kerning) const
{
  int   result = 0;
  tchar c; 
  for(int i=0; i<text.length(); i++)
  {
    c       = text[i];
    result += getGlyphWidth(c) + kerning;
  }
  return result - kerning; // we added one too much in the loop
}

int BitmapFont::glyphIndexToX(const juce::String &text, int index, int kerning) const
{
  jassert( index < text.length() ); // index is out of range
  if( index >= text.length() )
    return 0;

  int   result = 0;
  tchar c; 
  for(int i=0; i<index; i++)
  {
    c       = text[i];
    result += getGlyphWidth(c) + kerning;
  }
  return result; 
}

int BitmapFont::xToGlyphIndex(const String &text, int xToFindGlyphIndexAt, int kerning) const
{
  int x     = 0;
  int index = 0;
  tchar c; 
  while( x < xToFindGlyphIndexAt && index < text.length() )
  {
    c  = text[index];
    x += getGlyphWidth(c) + kerning;
    index++;
  }
  return index;
}


//-------------------------------------------------------------------------------------------------
// internal functions:

void BitmapFont::createGlyphBitmaps()
{

}
