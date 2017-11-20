#ifndef RS_PIXELFONTROUNDEDBOLDA16D0_H
#define RS_PIXELFONTROUNDEDBOLDA16D0_H

namespace RSLib
{

  /**

  This class represents a pixel font with 16 pixels ascent and 0 pixels descent.

  */

  class RSLib_API rsPixelFontRoundedBoldA16D0 : public rsPixelFont
  {

  public:

    rsPixelFontRoundedBoldA16D0();

  protected:

    virtual void createGlyphImages();

  };

}

#endif
