#ifndef RS_GLOBALFONTINSTANCES_H
#define RS_GLOBALFONTINSTANCES_H

namespace RSLib
{

  /**

  This class provides global access to some objects of subclasses of rsPixelFont objects to be 
  used for text rendering.

  \todo refactor this (but how?) - it contains lots of ugly code duplication
  ...maybe have a single access function static rsPixelFont* getPixelFont(); 
  with some parameters for height and boldness ...this may also facilitate a later switch to
  runtime time rendered pixel-fonts from vector-fonts ...we'll see

  */
    
  class RSLib_API rsGlobalFontInstances
  {

  public:

    /** Returns a pointer to a valid rsPixelFontRoundedA7D0 object that can be used to render text.
    The function will either create a new object (thereby internally storing it) and return a 
    pointer to it or - if such a font object was already created by a previous call - will return 
    the same pointer to the internally stored object again. However, client code should not make 
    any assumptions about the validity of the pointer later on because it will be valid only up to 
    the next call to cleanUpFonts, so it's probably not a good idea to keep the pointer around for 
    a long time - you can use it only as long as you can assure, that there are no interfering
    calls to cleanUpFonts. */
    static rsPixelFont* getPixelFontRoundedA7D0();
    static rsPixelFont* getPixelFontRoundedBoldA10D0();
    static rsPixelFont* getPixelFontRoundedBoldA16D0();

    /** Deletes all instances of fonts. */
    static void cleanUpFonts();

    /*
    \todo maybe make cleanUpFonts private and add RSLib::cleanUp() as friend and require that 
    cleanUpFonts should only be called from the global RSLib::cleanUp() on program shutdown
    ...but no: what if two programs use the same dll, one program suts down (calling cleanUp) but 
    the other one still needs the library? ...seems better to just require client code to not
    assume the pointer to be valid for a long time...
    */

  private:
  
    static rsPixelFontRoundedA7D0      *normalFont7px;
    static rsPixelFontRoundedBoldA10D0 *boldFont10px;
    static rsPixelFontRoundedBoldA16D0 *boldFont16px; 

    //friend void RSLib::cleanUp();

  };


  /**

  This file declares some font objects for global access. The singleton pattern approach turned out 
  to be unwieldy because of undefined destruction of the objects.

  */

  /*
  // ..."already defined in RSFC.hpp" error:
  rsPixelFontRoundedA7D0      *normalFont7px;
  rsPixelFontRoundedBoldA9D0  *normalFont9px;
  rsPixelFontRoundedBoldA10D0 *boldFont10px;
  rsPixelFontRoundedBoldA16D0 *boldFont16px; 
  */

  /*
  // "unresolved external symbol" error (even when defined in the .cpp file):
  extern rsPixelFontRoundedA7D0      *normalFont7px;
  extern rsPixelFontRoundedBoldA9D0  *normalFont9px;
  extern rsPixelFontRoundedBoldA10D0 *boldFont10px;
  extern rsPixelFontRoundedBoldA16D0 *boldFont16px; 
  */
  
  /*
  // ..."already defined in RSFC.hpp" error:
  RSLib_API rsPixelFontRoundedA7D0      *normalFont7px;
  RSLib_API rsPixelFontRoundedBoldA9D0  *normalFont9px;
  RSLib_API rsPixelFontRoundedBoldA10D0 *boldFont10px;
  RSLib_API rsPixelFontRoundedBoldA16D0 *boldFont16px; 
  */
  
  /*
  // pointers are NULL when used, despite being intialized:
  static rsPixelFontRoundedA7D0      *normalFont7px;
  static rsPixelFontRoundedBoldA9D0  *normalFont9px;
  static rsPixelFontRoundedBoldA10D0 *boldFont10px;
  static rsPixelFontRoundedBoldA16D0 *boldFont16px; 
  */

  /*
  static const rsPixelFontRoundedA7D0      normalFont7px;
  static const rsPixelFontRoundedBoldA9D0  normalFont9px;
  static const rsPixelFontRoundedBoldA10D0 boldFont10px;
  static const rsPixelFontRoundedBoldA16D0 boldFont16px;  
  */

}

#endif  
