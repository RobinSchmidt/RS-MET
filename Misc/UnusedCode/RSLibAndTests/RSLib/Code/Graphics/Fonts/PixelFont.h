#ifndef RS_PIXELFONT_H
#define RS_PIXELFONT_H

namespace RSLib
{

  /**

  This class ...

  \todo: maybe have a true "bitmap" font class that really stores one bit per pixel (currently 
  there's a 'char' per pixel - on the other hand, with the char, we could use the in-between 
  values to have "gray" pixels in the letter - maybe use this to make the font look better)
  -> reduces dll size
  ->but maybe someday we'll have vector-fonts anyway (and render pixel-fonts from them on startup)
  consolidate all embedded pixel-fonts into one pair of h/cpp files EmbeddedPixelFonts.h/cpp

  */

  class RSLib_API rsPixelFont // : public rsFont 
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsPixelFont();

    /** Destructor. */
    virtual ~rsPixelFont();


    /** \name Inquiry */

    virtual int getFontHeight()     const { return ascent+descent; }
    virtual int getFontAscent()     const { return ascent; }
    virtual int getFontDescent()    const { return descent; }
    virtual int getDefaultKerning() const { return defaultKerning; }

    /** Returns the width of the text (in pixels) when the text is rendered with this font at a
    given kerning. */
    virtual int getTextPixelWidth(const rsString &text, int kerning) const;

    /** Returns the x coordinate (in pixels) of the (left side of) a glyph in some text when the
    text is rendered with this font at a given kerning. */
    virtual int glyphIndexToX(const rsString &text, int index, int kerning) const;

    /** .... */
    virtual int xToGlyphIndex(const rsString &text, int xToFindGlyphIndexAt,
      int kerning) const;


    /** Returns a pointer to an image that represents the glyph at a specific Colour. The pointer
    points to a member-variable of this class, so don't delete it. */
    virtual const rsImageGray* getGlyphImage(const rsUint8 charOrAsciiCode
      /*, const rsColorRGBA &colour // remnant from rojue */) const
    {
      if( glyphImages[charOrAsciiCode] != NULL )
      {
        //glyphImages[charOrAsciiCode]->setColour(colour);
        return glyphImages[charOrAsciiCode];
      }
      else
      {
        //rassertfalse;  // no glyph available for this character
        return NULL;
      }
    }

    /** Returns the width of one of the glyphs in pixels. */
    virtual int getGlyphWidth(const rsUint8 charOrAsciiCode) const
    {
      return glyphWidths[charOrAsciiCode];
    }

  protected:

    /** \name Internal Functions */

    virtual void createGlyphImages();


    /** \name Data */

    int              ascent, descent, defaultKerning;
    static const int numGlyphs = 256;
    rsImageGray*     glyphImages[numGlyphs];
    int              glyphWidths[numGlyphs];

  };

}

#endif
