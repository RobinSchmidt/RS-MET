#ifndef rojue_BitmapFontBlockyA10D3_h
#define rojue_BitmapFontBlockyA10D3_h

#include "rojue_BitmapFont.h"

namespace rojue
{

  /**

  This class 

  */

  class BitmapFontBlockyA10D3 : public BitmapFont
  {

  public:

    BitmapFontBlockyA10D3::BitmapFontBlockyA10D3();



    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    virtual void createGlyphBitmaps();

  };

}

#endif  
