#ifndef jura_BitmapFont_h
#define jura_BitmapFont_h

/**

This is a baseclass for representing the bitmap fonts that are used in the GUI components. I'm not
really happy with JUCE's font rendering - it looks blurry, especially for bright-on-dark text, so
I'm defining bitmap fonts which have the advantage of always being ultra-sharp. The disadvantage is
that they do not allow for scaling, transforming, etc. The actual fonts are subclasses of this
class which define a set of related fonts (different sizes, weights, etc.) that are used on the
GUI components. Of each such font subclass, there's a single global object.

*/

class JUCE_API BitmapFont //: public Font
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
  virtual const Image* getGlyphImage(const char charOrAsciiCode, const Colour& colour) const
  {
    int i = charOrAsciiCode; // using the char directly gives a warning in gcc/win
    if(i >= numGlyphs)
    {
      jassertfalse;  // no glyph available for this character
      return NULL;
    }

    if(glyphBitmaps[i] != NULL)
    {
      glyphBitmaps[i]->setColour(colour);
      return glyphBitmaps[i];
    }
    else
    {
      //jassertfalse;  // no glyph available for this character
      return NULL;
    }
  }

  /** Returns the width of one of the glyphs in pixels. */
  //virtual const int getGlyphWidth(const char charOrAsciiCode) const
  virtual const int getGlyphWidth(const juce_wchar charOrAsciiCode) const
  {
    if(charOrAsciiCode < numGlyphs)
      return glyphWidths[charOrAsciiCode];
    else
      return 0;
  }

  juce_UseDebuggingNewOperator;

protected:

  virtual void createGlyphBitmaps();

  int                 ascent, descent, defaultKerning;
  static const int    numGlyphs = 256;
  ColourizableBitmap* glyphBitmaps[numGlyphs];
  int                 glyphWidths[numGlyphs];

};

//=================================================================================================
// the actual fonts - in the cpp file, the glyph bitmaps are created (maybe this can be done
// better using binary resources at some point):


/** This is the smallest bitmap font. It has an ascent of 7 pixels and no descent below the
baseline. It is used mainly in plots for the axis labels, ticks, etc. */
class BitmapFontRoundedA7D0 : public BitmapFont
{
public:
  BitmapFontRoundedA7D0();
protected:
  virtual void createGlyphBitmaps();
  juce_UseDebuggingNewOperator;
};
//static const BitmapFontRoundedA7D0 normalFont7px;

/** This is a medium sized bitmap font. It has an ascent of 9 pixels and no descent below the
baseline. It is currently not used anywhere (i think).
Maybe it was once used in the alpha versions of Liberty? ...better to not delete it, just in
case...
*/
class BitmapFontRoundedBoldA9D0 : public BitmapFont
{
public:
  BitmapFontRoundedBoldA9D0();
protected:
  virtual void createGlyphBitmaps();
  juce_UseDebuggingNewOperator;
};
//static const BitmapFontRoundedBoldA9D0 normalFont9px;

/** This is a bold font with an acsent of 10 pixels and no descent. It is used in most GUI widgets
such as sliders and comboboxes. */
class BitmapFontRoundedBoldA10D0 : public BitmapFont
{
public:

  BitmapFontRoundedBoldA10D0();

  static const BitmapFontRoundedBoldA10D0 instance;
  // we have one static member of this class as instance variable - this is the global object
  // which we can access from all widgets, editors, etc. - that seems to work - we should use that
  // strategy for the other fonts, too - doing i this way makes them available in other modules,
  // for example in the jura_audioprocessors module

protected:
  virtual void createGlyphBitmaps();
  juce_UseDebuggingNewOperator;
};
//static const BitmapFontRoundedBoldA10D0 boldFont10px;

/** This is a large font that is used mainly for headlines of GUI editors. */
class BitmapFontRoundedBoldA16D0 : public BitmapFont
{
public:
  BitmapFontRoundedBoldA16D0();
  static const BitmapFontRoundedBoldA10D0 instance;
protected:
  virtual void createGlyphBitmaps();
  juce_UseDebuggingNewOperator;
};
//static const BitmapFontRoundedBoldA16D0 boldFont16px;
// this is still defined in the .cpp file - but we should use the "instance" member instead 
// everywhere - such that it can be used from dependent library modules, too - then delete
// the definition in the .cpp file

#endif
