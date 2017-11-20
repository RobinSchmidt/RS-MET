#ifndef RS_PIXELFONTROUNDEDBOLDA9D0_H
#define RS_PIXELFONTROUNDEDBOLDA9D0_H

namespace RSLib
{

  /**

  This class represents a pixel font with 9 pixels ascent and 0 pixels descent.

  */

  class RSLib_API rsPixelFontRoundedBoldA9D0 : public rsPixelFont
  {

  public:

    rsPixelFontRoundedBoldA9D0();

  protected:

    virtual void createGlyphImages();

  };

}

#endif
