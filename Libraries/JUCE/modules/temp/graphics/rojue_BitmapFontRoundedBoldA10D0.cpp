#include "rojue_BitmapFontRoundedBoldA10D0.h"
using namespace rojue;

BitmapFontRoundedBoldA10D0::BitmapFontRoundedBoldA10D0()
{
  ascent  = 10;
  descent = 0;

  createGlyphBitmaps();
}

//-------------------------------------------------------------------------------------------------
// create the glyphs:



void BitmapFontRoundedBoldA10D0::createGlyphBitmaps()
{
  unsigned char _ = 0;
  unsigned char X = 255;

  /*
  unsigned char glyph_dummy[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_
  };
  glyphBitmaps[60] = new ColourizableBitmap(6, 10, glyph_glyph_dummy);
  */


  unsigned char glyph_space[4*10] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
  };
  glyphBitmaps[32] = new ColourizableBitmap(4, 10, glyph_space);

  unsigned char glyph_exclamationMark[2*10] = 
  {
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    _,_,
    _,_,
    X,X,
    X,X
  };
  glyphBitmaps[33] = new ColourizableBitmap(2, 10, glyph_exclamationMark);

  unsigned char glyph_quotes[6*10] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_
  };
  glyphBitmaps[34] = new ColourizableBitmap(6, 10, glyph_quotes);

  unsigned char glyph_sharp[8*10] = 
  {
    _,X,X,_,_,X,X,_,
    _,X,X,_,_,X,X,_,
    X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,
    _,X,X,_,_,X,X,_,
    _,X,X,_,_,X,X,_,
    X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,
    _,X,X,_,_,X,X,_,
    _,X,X,_,_,X,X,_,
  };
  glyphBitmaps[35] = new ColourizableBitmap(8, 10, glyph_sharp);

  unsigned char glyph_dollar[6*10] = 
  {
    _,_,X,X,_,_,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    _,_,X,X,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    _,_,X,X,_,_
  };
  glyphBitmaps[36] = new ColourizableBitmap(6, 10, glyph_dollar);

  unsigned char glyph_percent[6*10] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphBitmaps[37] = new ColourizableBitmap(6, 10, glyph_percent);

  unsigned char glyph_ampersand[8*10] = 
  {
    _,X,X,X,_,_,_,_,
    X,X,X,X,X,_,_,_,
    X,X,_,X,X,_,_,_,
    X,X,_,X,X,_,_,_,
    X,X,X,X,_,_,X,X,
    _,X,X,X,X,X,X,X,
    X,X,_,X,X,X,_,_,
    X,X,_,X,X,X,_,_,
    X,X,X,X,X,X,X,X,
    _,X,X,X,_,_,X,X
  };
  glyphBitmaps[38] = new ColourizableBitmap(8, 10, glyph_ampersand);

  unsigned char glyph_singlequote[2*10] = 
  {
    X,X,
    X,X,
    X,X,
    X,X,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_
  };
  glyphBitmaps[39] = new ColourizableBitmap(2, 10, glyph_singlequote);

  unsigned char glyph_openingBrace[3*10] = 
  {
    _,_,X,
    _,X,X,
    _,X,X,
    X,X,_,
    X,X,_,
    X,X,_,
    X,X,_,
    _,X,X,
    _,X,X,
    _,_,X
  };
  glyphBitmaps[40] = new ColourizableBitmap(3, 10, glyph_openingBrace);

  unsigned char glyph_closingBrace[3*10] = 
  {
    X,_,_,
    X,X,_,
    X,X,_,
    _,X,X,
    _,X,X,
    _,X,X,
    _,X,X,
    X,X,_,
    X,X,_,
    X,_,_
  };
  glyphBitmaps[41] = new ColourizableBitmap(3, 10, glyph_closingBrace);

  unsigned char glyph_multiply[8*10] = 
  {
    _,_,_,_,_,_,_,_,
    _,_,_,X,X,_,_,_,
    _,_,_,X,X,_,_,_,
    X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,_,_,
    _,X,X,X,X,X,X,_,
    _,X,X,_,_,X,X,_,
    _,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_
  };
  glyphBitmaps[42] = new ColourizableBitmap(8, 10, glyph_multiply);

  unsigned char glyph_plus[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
  };
  glyphBitmaps[43] = new ColourizableBitmap(6, 10, glyph_plus);


  unsigned char glyph_comma[2*10] = 
  {
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    X,X,
    X,X,
    X,X,
    X,X
  };
  glyphBitmaps[44] = new ColourizableBitmap(2, 10, glyph_comma);

  unsigned char glyph_minus[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
  };
  glyphBitmaps[45] = new ColourizableBitmap(6, 10, glyph_minus);

  unsigned char glyph_dot[2*10] = 
  {
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    X,X,
    X,X
  };
  glyphBitmaps[46] = new ColourizableBitmap(2, 10, glyph_dot);

  unsigned char glyph_slash[6*10] = 
  {
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphBitmaps[47] = new ColourizableBitmap(6, 10, glyph_slash);

  unsigned char glyph_0[6*10] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphBitmaps[48] = new ColourizableBitmap(6, 10, glyph_0);

  unsigned char glyph_1[4*10] = 
  {
    _,_,X,X,
    _,X,X,X,
    X,X,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X
  };
  glyphBitmaps[49] = new ColourizableBitmap(4, 10, glyph_1);

  unsigned char glyph_2[6*10] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,X,X,X,X,
    _,X,X,X,X,_,
    X,X,X,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X
  };
  glyphBitmaps[50] = new ColourizableBitmap(6, 10, glyph_2);

  unsigned char glyph_3[6*10] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,X,X,X,X,X,
    _,X,X,X,X,_,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[51] = new ColourizableBitmap(6, 10, glyph_3);


  unsigned char glyph_4[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X
  };
  glyphBitmaps[52] = new ColourizableBitmap(6, 10, glyph_4);

  /*
  unsigned char glyph_5[6*10] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_
  };
  glyphBitmaps[53] = new ColourizableBitmap(6, 10, glyph_5);
  */

  unsigned char glyph_5[6*10] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[53] = new ColourizableBitmap(6, 10, glyph_5);

  unsigned char glyph_6[6*10] = 
  {
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphBitmaps[54] = new ColourizableBitmap(6, 10, glyph_6);

  unsigned char glyph_7[6*10] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_
  };
  glyphBitmaps[55] = new ColourizableBitmap(6, 10, glyph_7);

  unsigned char glyph_8[6*10] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    _,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphBitmaps[56] = new ColourizableBitmap(6, 10, glyph_8);

  unsigned char glyph_9[6*10] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[57] = new ColourizableBitmap(6, 10, glyph_9);

  unsigned char glyph_colon[2*10] = 
  {
    _,_,
    _,_,
    X,X,
    X,X,
    _,_,
    _,_,
    X,X,
    X,X,
    _,_,
    _,_
  };
  glyphBitmaps[58] = new ColourizableBitmap(2, 10, glyph_colon);


  unsigned char glyph_semicolon[2*10] = 
  {
    _,_,
    _,_,
    _,_,
    X,X,
    X,X,
    _,_,
    X,X,
    X,X,
    X,X,
    X,X
  };
  glyphBitmaps[59] = new ColourizableBitmap(2, 10, glyph_semicolon);

  unsigned char glyph_openingAngleBrace[6*10] = 
  {
    _,_,_,_,X,X,
    _,_,_,X,X,X,
    _,_,X,X,X,_,
    _,X,X,X,_,_,
    X,X,X,_,_,_,
    X,X,X,_,_,_,
    _,X,X,X,_,_,
    _,_,X,X,X,_,
    _,_,_,X,X,X,
    _,_,_,_,X,X
  };
  glyphBitmaps[60] = new ColourizableBitmap(6, 10, glyph_openingAngleBrace);

  unsigned char glyph_equalsSign[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,_,_,_,_,
    _,_,_,_,_,_
  };
  glyphBitmaps[61] = new ColourizableBitmap(6, 10, glyph_equalsSign);

  unsigned char glyph_closingAngleBrace[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,X,_,_,_,
    _,X,X,X,_,_,
    _,_,X,X,X,_,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,X,X,X,_,
    _,X,X,X,_,_,
    X,X,X,_,_,_,
    X,X,_,_,_,_
  };
  glyphBitmaps[62] = new ColourizableBitmap(6, 10, glyph_closingAngleBrace);

  unsigned char glyph_questionMark[6*10] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphBitmaps[63] = new ColourizableBitmap(6, 10, glyph_questionMark);


  unsigned char glyph_at[10*10] = 
  {
    _,X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,X,
    X,X,_,_,_,_,_,_,X,X,
    X,X,_,X,X,X,X,_,X,X,
    X,X,_,X,_,_,X,_,X,X,
    X,X,_,X,_,_,X,_,X,X,
    X,X,_,X,X,X,X,X,X,_,
    X,X,_,_,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,X
  };
  glyphBitmaps[64] = new ColourizableBitmap(10, 10, glyph_at);

  unsigned char glyph_openingBracket[4*10] = 
  {
    X,X,X,X,
    X,X,X,X,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,X,X,
    X,X,X,X
  };
  glyphBitmaps[91] = new ColourizableBitmap(4, 10, glyph_openingBracket);

  unsigned char glyph_backslash[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,_,_,X,X,
    _,_,_,_,X,X
  };
  glyphBitmaps[92] = new ColourizableBitmap(6, 10, glyph_backslash);

  unsigned char glyph_closingBracket[4*10] = 
  {
    X,X,X,X,
    X,X,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    X,X,X,X,
    X,X,X,X
  };
  glyphBitmaps[93] = new ColourizableBitmap(4, 10, glyph_closingBracket);


  unsigned char glyph_power[7*10] = 
  {
    _,_,_,X,_,_,_,
    _,_,X,X,X,_,_,
    _,X,X,X,X,X,_,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_
  };
  glyphBitmaps[94] = new ColourizableBitmap(7, 10, glyph_power);

  unsigned char glyph_underscore[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X
  };
  glyphBitmaps[95] = new ColourizableBitmap(6, 10, glyph_underscore);

  unsigned char glyph_96[3*10] = 
  {
    X,X,_,
    X,X,X,
    X,X,X,
    _,X,X,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_
  };
  glyphBitmaps[96] = new ColourizableBitmap(3, 10, glyph_96);

  unsigned char glyph_openingCurlyBrace[4*10] = 
  {
    _,X,X,X,
    _,X,X,X,
    _,X,X,_,
    _,X,X,_,
    X,X,X,_,
    X,X,X,_,
    _,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,X,X,X
  };
  glyphBitmaps[123] = new ColourizableBitmap(4, 10, glyph_openingCurlyBrace);

  unsigned char glyph_verticalLine[2*10] = 
  {
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X
  };
  glyphBitmaps[124] = new ColourizableBitmap(2, 10, glyph_verticalLine);

  unsigned char glyph_closingCurlyBrace[4*10] = 
  {
    X,X,X,_,
    X,X,X,_,
    _,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,X,X,X,
    _,X,X,_,
    _,X,X,_,
    X,X,X,_,
    X,X,X,_
  };
  glyphBitmaps[125] = new ColourizableBitmap(4, 10, glyph_closingCurlyBrace);

  unsigned char glyph_tilde[7*10] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,X,X,X,_,X,X,
    X,X,X,X,X,X,X,
    X,X,_,X,X,X,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_
  };
  glyphBitmaps[126] = new ColourizableBitmap(7, 10, glyph_tilde);

  unsigned char glyph_A[6*10] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphBitmaps[65] = new ColourizableBitmap(6, 10, glyph_A);

  unsigned char glyph_a[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
  };
  glyphBitmaps[97] = new ColourizableBitmap(6, 10, glyph_a);

  unsigned char glyph_B[6*10] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[66] = new ColourizableBitmap(6, 10, glyph_B);

  unsigned char glyph_b[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[98] = new ColourizableBitmap(6, 10, glyph_b);

  unsigned char glyph_C[6*10] = 
  {
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    _,X,X,X,X,X
  };
  glyphBitmaps[67] = new ColourizableBitmap(6, 10, glyph_C);

  unsigned char glyph_c[5*10] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,X,X,X,X,
    X,X,X,X,X,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,X,X,X,
    _,X,X,X,X
  };
  glyphBitmaps[99] = new ColourizableBitmap(5, 10, glyph_c);

  unsigned char glyph_D[6*10] = 
  {
    X,X,X,X,_,_,
    X,X,X,X,X,_,
    X,X,_,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_
  };
  glyphBitmaps[68] = new ColourizableBitmap(6, 10, glyph_D);

  unsigned char glyph_d[6*10] = 
  {
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X
  };
  glyphBitmaps[100] = new ColourizableBitmap(6, 10, glyph_d);

  unsigned char glyph_E[6*10] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X
  };
  glyphBitmaps[69] = new ColourizableBitmap(6, 10, glyph_E);

  unsigned char glyph_e[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    _,X,X,X,X,X
  };
  glyphBitmaps[101] = new ColourizableBitmap(6, 10, glyph_e);

  unsigned char glyph_F[6*10] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
  };
  glyphBitmaps[70] = new ColourizableBitmap(6, 10, glyph_F);

  unsigned char glyph_f[5*10] = 
  {
    _,_,X,X,X,
    _,X,X,X,X,
    _,X,X,_,_,
    _,X,X,_,_,
    X,X,X,X,_,  
    X,X,X,X,_,
    _,X,X,_,_,
    _,X,X,_,_,
    _,X,X,_,_,
    _,X,X,_,_
  };
  glyphBitmaps[102] = new ColourizableBitmap(5, 10, glyph_f);

  unsigned char glyph_G[6*10] = 
  {
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,X,X,X,
    X,X,_,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_,
  };
  glyphBitmaps[71] = new ColourizableBitmap(6, 10, glyph_G);

  unsigned char glyph_g[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
  };
  glyphBitmaps[103] = new ColourizableBitmap(6, 10, glyph_g);

  unsigned char glyph_H[6*10] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
  };
  glyphBitmaps[72] = new ColourizableBitmap(6, 10, glyph_H);

  unsigned char glyph_h[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
  };
  glyphBitmaps[104] = new ColourizableBitmap(6, 10, glyph_h);

  unsigned char glyph_I[2*10] = 
  {
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X
  };
  glyphBitmaps[73] = new ColourizableBitmap(2, 10, glyph_I);

  unsigned char glyph_i[2*10] = 
  {
    X,X,
    X,X,
    _,_,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X
  };
  glyphBitmaps[105] = new ColourizableBitmap(2, 10, glyph_i);

  unsigned char glyph_J[6*10] = 
  {
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_,
  };
  glyphBitmaps[74] = new ColourizableBitmap(6, 10, glyph_J);

  unsigned char glyph_j[4*10] = 
  {
    _,_,X,X,
    _,_,X,X,
    _,_,_,_,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    X,X,X,X,
    X,X,X,_
  };
  glyphBitmaps[106] = new ColourizableBitmap(4, 10, glyph_j);

  unsigned char glyph_K[6*10] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_,
    X,X,X,X,X,_,
    X,X,_,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphBitmaps[75] = new ColourizableBitmap(6, 10, glyph_K);

  unsigned char glyph_k[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_,
    X,X,X,X,X,_,
    X,X,_,X,X,X,
    X,X,_,_,X,X,
  };
  glyphBitmaps[107] = new ColourizableBitmap(6, 10, glyph_k);

  unsigned char glyph_L[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X
  };
  glyphBitmaps[76] = new ColourizableBitmap(6, 10, glyph_L);

  unsigned char glyph_l[3*10] = 
  {
    X,X,_,
    X,X,_,
    X,X,_,
    X,X,_,
    X,X,_,
    X,X,_,
    X,X,_,
    X,X,_,
    X,X,X,
    X,X,X
  };
  glyphBitmaps[108] = new ColourizableBitmap(3, 10, glyph_l);

  unsigned char glyph_M[9*10] = 
  {
    X,X,_,_,_,_,_,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,_,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,_,X,X,X,_,X,X,
    X,X,_,_,X,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X
  };
  glyphBitmaps[77] = new ColourizableBitmap(9, 10, glyph_M);

  unsigned char glyph_m[8*10] = 
  {
    _,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,
    X,X,X,_,X,X,X,_,
    X,X,X,X,X,X,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X
  };
  glyphBitmaps[109] = new ColourizableBitmap(8, 10, glyph_m);

  unsigned char glyph_N[7*10] = 
  {
    X,X,_,_,_,X,X,
    X,X,X,_,_,X,X,
    X,X,X,X,_,X,X,
    X,X,X,X,X,X,X,
    X,X,_,X,X,X,X,
    X,X,_,_,X,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X
  };
  glyphBitmaps[78] = new ColourizableBitmap(7, 10, glyph_N);

  unsigned char glyph_n[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphBitmaps[110] = new ColourizableBitmap(6, 10, glyph_n);

  unsigned char glyph_O[6*10] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphBitmaps[79] = new ColourizableBitmap(6, 10, glyph_O);

  unsigned char glyph_o[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphBitmaps[111] = new ColourizableBitmap(6, 10, glyph_o);

  unsigned char glyph_P[6*10] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphBitmaps[80] = new ColourizableBitmap(6, 10, glyph_P);

  unsigned char glyph_p[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphBitmaps[112] = new ColourizableBitmap(6, 10, glyph_p);

  unsigned char glyph_Q[7*10] = 
  {
    _,X,X,X,X,_,_,
    X,X,X,X,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,X,X,X,_,
    X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,
  };
  glyphBitmaps[81] = new ColourizableBitmap(7, 10, glyph_Q);

  unsigned char glyph_q[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X
  };
  glyphBitmaps[113] = new ColourizableBitmap(6, 10, glyph_q);

  unsigned char glyph_R[6*10] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_,
    X,X,X,X,X,_,
    X,X,_,X,X,X,
    X,X,_,_,X,X
  };
  glyphBitmaps[82] = new ColourizableBitmap(6, 10, glyph_R);

  unsigned char glyph_r[5*10] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,_,X,X,
    X,X,X,X,X,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_
  };
  glyphBitmaps[114] = new ColourizableBitmap(5, 10, glyph_r);

  unsigned char glyph_S[6*10] = 
  {
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[83] = new ColourizableBitmap(6, 10, glyph_S);

  unsigned char glyph_s[5*10] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,X,X,X,X,
    X,X,X,X,X,
    X,X,_,_,_,
    X,X,X,X,_,
    _,X,X,X,X,
    _,_,_,X,X,
    X,X,X,X,X,
    X,X,X,X,_
  };
  glyphBitmaps[115] = new ColourizableBitmap(5, 10, glyph_s);

  unsigned char glyph_T[6*10] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_
  };
  glyphBitmaps[84] = new ColourizableBitmap(6, 10, glyph_T);

  unsigned char glyph_t[4*10] = 
  {
    _,X,X,_,
    _,X,X,_,
    X,X,X,X,
    X,X,X,X,
    _,X,X,_,
    _,X,X,_,
    _,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,_,X,X
  };
  glyphBitmaps[116] = new ColourizableBitmap(4, 10, glyph_t);

  unsigned char glyph_U[6*10] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphBitmaps[85] = new ColourizableBitmap(6, 10, glyph_U);

  unsigned char glyph_u[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphBitmaps[117] = new ColourizableBitmap(6, 10, glyph_u);

  unsigned char glyph_V[7*10] = 
  {
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,X,_,X,X,X,
    _,X,X,X,X,X,_,
    _,_,X,X,X,_,_,
    _,_,_,X,_,_,_
  };
  glyphBitmaps[86] = new ColourizableBitmap(7, 10, glyph_V);

  unsigned char glyph_v[7*10] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,X,_,X,X,X,
    _,X,X,X,X,X,_,
    _,_,X,X,X,_,_,
    _,_,_,X,_,_,_
  };
  glyphBitmaps[118] = new ColourizableBitmap(7, 10, glyph_v);

  unsigned char glyph_W[9*10] = 
  {
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,X,_,_,X,X,
    X,X,_,X,X,X,_,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,_,_,_,_,_,X,X
  };
  glyphBitmaps[87] = new ColourizableBitmap(9, 10, glyph_W);

  /*
  unsigned char glyph_w[9*10] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,X,_,_,X,X,
    X,X,_,X,X,X,_,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,_,_,_,_,_,X,X
  };
  glyphBitmaps[119] = new ColourizableBitmap(9, 10, glyph_w);
  */
  unsigned char glyph_w[7*10] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,X,_,X,X,
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X
  };
  glyphBitmaps[119] = new ColourizableBitmap(7, 10, glyph_w);

  unsigned char glyph_X[7*10] = 
  {
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,X,_,X,X,X,
    _,X,X,X,X,X,_,
    _,_,X,X,X,_,_,
    _,X,X,X,X,X,_,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X
  };
  glyphBitmaps[88] = new ColourizableBitmap(7, 10, glyph_X);

  unsigned char glyph_x[7*10] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,X,_,_,_,X,X,
    X,X,X,_,X,X,X,
    _,X,X,X,X,X,_,
    _,_,X,X,X,_,_,
    _,_,X,X,X,_,_,
    _,X,X,X,X,X,_,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X
  };
  glyphBitmaps[120] = new ColourizableBitmap(7, 10, glyph_x);

  unsigned char glyph_Y[6*10] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[89] = new ColourizableBitmap(6, 10, glyph_Y);

  unsigned char glyph_y[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphBitmaps[121] = new ColourizableBitmap(6, 10, glyph_y);

  /*
  unsigned char glyph_Y[6*10] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_
  };
  glyphBitmaps[89] = new ColourizableBitmap(6, 10, glyph_Y);

  unsigned char glyph_y[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_
  };
  glyphBitmaps[121] = new ColourizableBitmap(6, 10, glyph_y);
  */

  unsigned char glyph_Z[7*10] = 
  {
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
    _,_,_,_,_,X,X,
    _,_,_,_,X,X,X,
    _,_,_,X,X,X,_,
    _,_,X,X,X,_,_,
    _,X,X,X,_,_,_,
    X,X,X,_,_,_,_,
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
  };
  glyphBitmaps[90] = new ColourizableBitmap(7, 10, glyph_Z);

  unsigned char glyph_z[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,_,X,X,X,
    _,_,X,X,X,_,
    _,X,X,X,_,_,
    X,X,X,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X
  };
  glyphBitmaps[122] = new ColourizableBitmap(6, 10, glyph_z);


  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphBitmaps[g] != NULL )
      glyphWidths[g] = glyphBitmaps[g]->getWidth();
  }
}
