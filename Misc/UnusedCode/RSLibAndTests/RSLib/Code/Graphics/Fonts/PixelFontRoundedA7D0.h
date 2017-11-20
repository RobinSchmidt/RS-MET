#ifndef RS_PIXELFONTROUNDEDA7D0_H
#define RS_PIXELFONTROUNDEDA7D0_H

namespace RSLib
{

  /**

  This class represents a pixels font with 7 pixels ascent and 0 pixels descent.

  */

  class RSLib_API rsPixelFontRoundedA7D0 : public rsPixelFont
  {

  public:

    rsPixelFontRoundedA7D0();

  protected:

    virtual void createGlyphImages();

  };

}

#endif
