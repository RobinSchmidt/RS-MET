#ifndef RS_PIXELFONTROUNDEDBOLDA10D0_H
#define RS_PIXELFONTROUNDEDBOLDA10D0_H

namespace RSLib
{

  /**

  This class represents a pixel font with 10 pixels ascent and 0 pixels descent.

  */

  class RSLib_API rsPixelFontRoundedBoldA10D0 : public rsPixelFont
  {

  public:

    rsPixelFontRoundedBoldA10D0();

  protected:

    virtual void createGlyphImages();

  };

}

#endif
