#ifndef rojue_GlobalFontInstances_h
#define rojue_GlobalFontInstances_h

#include "rojue_BitmapFont.h"

#include "rojue_BitmapFontRoundedA7D0.h"        // for plot axes
#include "rojue_BitmapFontRoundedBoldA9D0.h"        // for ...
#include "rojue_BitmapFontRoundedBoldA10D0.h"   // for widget texts
#include "rojue_BitmapFontRoundedBoldA16D0.h"   // for headlines

namespace rojue
{

  /**

  This file declares some font objects for global access. The singleton pattern approach turned out 
  to be unwieldy because of undefined destruction of the objects.

  */

  static const BitmapFontRoundedA7D0      normalFont7px;
  static const BitmapFontRoundedBoldA9D0      normalFont9px;
  static const BitmapFontRoundedBoldA10D0 boldFont10px;
  static const BitmapFontRoundedBoldA16D0 boldFont16px;  

}

#endif  
