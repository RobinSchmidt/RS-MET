#ifndef rojue_BitmapFont_h
#define rojue_BitmapFont_h

#include "rojue_ColourizableBitmap.h"

namespace rojue
{

  /**

  This class 

  */

  class BitmapFont //: public Font
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    //BitmapFont* getInstance();

    /** Constructor. */
    BitmapFont();   

    /** Destructor. */
    virtual ~BitmapFont(); 

    //---------------------------------------------------------------------------------------------
    // setup:


    //---------------------------------------------------------------------------------------------
    // inquiry:

    virtual int getFontHeight()     const { return ascent+descent; }
    virtual int getFontAscent()     const { return ascent; }
    virtual int getFontDescent()    const { return descent; }
    virtual int getDefaultKerning() const { return defaultKerning; }

    /** Returns the width of the text (in pixels) when the text is rendered with this font at a 
    given kerning. */
    virtual int getTextPixelWidth(const juce::String& text, int kerning) const;

    /** Returns the x coordinate (in pixels) of the (left side of) a glyph in some text when the 
    text is rendered with this font at a given kerning. */
    virtual int glyphIndexToX(const juce::String& text, int index, int kerning) const;

    virtual int xToGlyphIndex(const juce::String& text, int xToFindGlyphIndexAt, int kerning) const;


    /** Returns a pointer to an image that represents the glyph at a specific Colour. The pointer 
    points to a member-variable of this class, so don't delete it. */
    virtual const Image* getGlyphImage(const tchar charOrAsciiCode, const Colour& colour) const
    {
      if( charOrAsciiCode >= numGlyphs )
      {
        jassertfalse;  // no glyph available for this character
        return NULL; 
      }

      if( glyphBitmaps[charOrAsciiCode] != NULL )
      {
        glyphBitmaps[charOrAsciiCode]->setColour(colour);
        return glyphBitmaps[charOrAsciiCode];
      }
      else
      {
        //jassertfalse;  // no glyph available for this character
        return NULL;
      }
    }

    /** Returns the width of one of the glyphs in pixels. */
    virtual const int getGlyphWidth(const tchar charOrAsciiCode) const
    { 
      if( charOrAsciiCode < numGlyphs )
        return glyphWidths[charOrAsciiCode]; 
      else
        return 0;
    }

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    virtual void createGlyphBitmaps();

    int                 ascent, descent, defaultKerning;
    static const int    numGlyphs = 256;
    ColourizableBitmap* glyphBitmaps[numGlyphs];
    int                 glyphWidths[numGlyphs];


  };

}

#endif  
