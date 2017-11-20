#include "rojue_BitmapFontRoundedA7D0.h"
using namespace rojue;

BitmapFontRoundedA7D0::BitmapFontRoundedA7D0()
{
  ascent  = 7;
  descent = 0;

  createGlyphBitmaps();
}

//-------------------------------------------------------------------------------------------------
// create the glyphs:



void BitmapFontRoundedA7D0::createGlyphBitmaps()
{
  unsigned char _ = 0;
  unsigned char X = 255;

  /*
  unsigned char glyph_dummy[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphBitmaps[60] = new ColourizableBitmap(5, 7, glyph_glyph_dummy);
  */


  unsigned char glyph_space[2*7] = 
  {
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_
  };
  glyphBitmaps[32] = new ColourizableBitmap(2, 7, glyph_space);

  unsigned char glyph_exclamationMark[1*7] = 
  {
    X,
    X,
    X,
    X,
    _,
    _,
    X
  };
  glyphBitmaps[33] = new ColourizableBitmap(1, 7, glyph_exclamationMark);

  unsigned char glyph_quotes[3*7] = 
  {
    X,_,X,
    X,_,X,
    X,_,X,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
  };
  glyphBitmaps[34] = new ColourizableBitmap(3, 7, glyph_quotes);

  unsigned char glyph_sharp[5*7] = 
  {
    _,X,_,X,_,
    _,X,_,X,_,
    X,X,X,X,X,
    _,X,_,X,_,
    X,X,X,X,X,
    _,X,_,X,_,
    _,X,_,X,_
  };
  glyphBitmaps[35] = new ColourizableBitmap(5, 7, glyph_sharp);

  unsigned char glyph_dollar[3*7] = 
  {
    _,X,_,
    X,X,X,
    X,_,_,
    _,X,_,
    _,_,X,
    X,X,X,
    _,X,_,
  };
  glyphBitmaps[36] = new ColourizableBitmap(3, 7, glyph_dollar);

  unsigned char glyph_percent[3*7] = 
  {
    X,_,X,
    _,_,X,
    _,X,X,
    _,X,_,
    X,X,_,
    X,_,_,
    X,_,X,
  };
  glyphBitmaps[37] = new ColourizableBitmap(3, 7, glyph_percent);

  unsigned char glyph_ampersand[5*7] = 
  {
    _,X,X,_,_,
    X,_,_,X,_,
    X,_,X,_,_,
    _,X,X,_,X,
    X,_,X,X,_,
    X,_,_,X,_,
    _,X,X,_,X
  };
  glyphBitmaps[38] = new ColourizableBitmap(5, 7, glyph_ampersand);

  unsigned char glyph_singlequote[1*7] = 
  {
    X,
    X,
    X,
    _,
    _,
    _,
    _,
  };
  glyphBitmaps[39] = new ColourizableBitmap(1, 7, glyph_singlequote);

  unsigned char glyph_openingBrace[2*7] = 
  {
    _,X,
    X,X,
    X,_,
    X,_,
    X,_,
    X,X,
    _,X,
  };
  glyphBitmaps[40] = new ColourizableBitmap(2, 7, glyph_openingBrace);

  unsigned char glyph_closingBrace[2*7] = 
  {
    X,_,
    X,X,
    _,X,
    _,X,
    _,X,
    X,X,
    X,_,
  };
  glyphBitmaps[41] = new ColourizableBitmap(2, 7, glyph_closingBrace);

  unsigned char glyph_multiply[5*7] = 
  {
    _,_,_,_,_,
    _,_,X,_,_,
    X,X,X,X,X,
    _,_,X,_,_,
    _,X,_,X,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphBitmaps[42] = new ColourizableBitmap(5, 7, glyph_multiply);

  unsigned char glyph_plus[5*7] = 
  {
    _,_,_,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    X,X,X,X,X,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,_,_,_
  };
  glyphBitmaps[43] = new ColourizableBitmap(5, 7, glyph_plus);

  unsigned char glyph_comma[1*7] = 
  {
    _,
    _,
    _,
    _,
    X,
    X,
    X,
  };
  glyphBitmaps[44] = new ColourizableBitmap(1, 7, glyph_comma);

  unsigned char glyph_minus[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_
  };
  glyphBitmaps[45] = new ColourizableBitmap(4, 7, glyph_minus);

  unsigned char glyph_dot[1*7] = 
  {
    _,
    _,
    _,
    _,
    _,
    _,
    X
  };
  glyphBitmaps[46] = new ColourizableBitmap(1, 7, glyph_dot);

  unsigned char glyph_slash[3*7] = 
  {
    _,_,X,
    _,_,X,
    _,X,X,
    _,X,_,
    X,X,_,
    X,_,_,
    X,_,_,
  };
  glyphBitmaps[47] = new ColourizableBitmap(3, 7, glyph_slash);

  unsigned char glyph_0[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[48] = new ColourizableBitmap(4, 7, glyph_0);

  unsigned char glyph_1[2*7] = 
  {
    _,X,
    X,X,
    _,X,
    _,X,
    _,X,
    _,X,
    _,X
  };
  glyphBitmaps[49] = new ColourizableBitmap(2, 7, glyph_1);

  unsigned char glyph_2[4*70] = 
  {
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    _,_,X,_,
    _,X,_,_,
    X,_,_,_,
    X,X,X,X
  };
  glyphBitmaps[50] = new ColourizableBitmap(4, 7, glyph_2);

  unsigned char glyph_3[4*7] = 
  {
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphBitmaps[51] = new ColourizableBitmap(4, 7, glyph_3);


  unsigned char glyph_4[4*7] = 
  {
    X,_,_,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,X,
    _,_,_,X,
    _,_,_,X,
    _,_,_,X
  };
  glyphBitmaps[52] = new ColourizableBitmap(4, 7, glyph_4);


  unsigned char glyph_5[4*7] = 
  {
    X,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphBitmaps[53] = new ColourizableBitmap(4, 7, glyph_5);

  unsigned char glyph_6[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[54] = new ColourizableBitmap(4, 7, glyph_6);

  unsigned char glyph_7[4*7] = 
  {
    X,X,X,X,
    _,_,_,X,
    _,_,X,_,
    _,_,X,_,
    _,X,_,_,
    _,X,_,_,
    X,_,_,_
  };
  glyphBitmaps[55] = new ColourizableBitmap(4, 7, glyph_7);

  unsigned char glyph_8[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[56] = new ColourizableBitmap(4, 7, glyph_8);

  unsigned char glyph_9[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,X,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphBitmaps[57] = new ColourizableBitmap(4, 7, glyph_9);

  unsigned char glyph_colon[1*7] = 
  {
    _,
    X,
    X,
    _,
    X,
    X,
    _,
  };
  glyphBitmaps[58] = new ColourizableBitmap(1, 7, glyph_colon);


  unsigned char glyph_semicolon[1*7] = 
  {
    _,
    X,
    X,
    _,
    X,
    X,
    X,
  };
  glyphBitmaps[59] = new ColourizableBitmap(1, 7, glyph_semicolon);

  unsigned char glyph_openingAngleBrace[4*7] = 
  {
    _,_,_,X,
    _,_,X,_,
    _,X,_,_,
    X,_,_,_,
    _,X,_,_,
    _,_,X,_,
    _,_,_,X,
  };
  glyphBitmaps[60] = new ColourizableBitmap(4, 7, glyph_openingAngleBrace);

  unsigned char glyph_equalsSign[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    _,_,_,_,
    X,X,X,X,
    _,_,_,_,
    _,_,_,_,
  };
  glyphBitmaps[61] = new ColourizableBitmap(4, 7, glyph_equalsSign);

  unsigned char glyph_closingAngleBrace[4*7] = 
  {
    X,_,_,_,
    _,X,_,_,
    _,_,X,_,
    _,_,_,X,
    _,_,X,_,
    _,X,_,_,
    X,_,_,_,
  };
  glyphBitmaps[62] = new ColourizableBitmap(4, 7, glyph_closingAngleBrace);

  unsigned char glyph_questionMark[3*7] = 
  {
    X,X,X,
    _,_,X,
    _,_,X,
    X,X,X,
    X,_,_,
    _,_,_,
    X,_,_,
  };
  glyphBitmaps[63] = new ColourizableBitmap(3, 7, glyph_questionMark);


  unsigned char glyph_at[5*7] = 
  {
    _,X,X,X,X,
    X,_,_,_,X,
    X,_,X,X,X,
    X,_,X,_,X,
    X,_,X,X,X,
    X,_,_,_,_,
    _,X,X,X,X
  };
  glyphBitmaps[64] = new ColourizableBitmap(5, 7, glyph_at);

  unsigned char glyph_openingBracket[2*7] = 
  {
    X,X,
    X,_,
    X,_,
    X,_,
    X,_,
    X,_,
    X,X,
  };
  glyphBitmaps[91] = new ColourizableBitmap(2, 7, glyph_openingBracket);

  unsigned char glyph_backslash[3*7] = 
  {
    X,_,_,
    X,_,_,
    X,X,_,
    _,X,_,
    _,X,X,
    _,_,X,
    _,_,X,
  };
  glyphBitmaps[92] = new ColourizableBitmap(3, 7, glyph_backslash);

  unsigned char glyph_closingBracket[2*7] = 
  {
    X,X,
    _,X,
    _,X,
    _,X,
    _,X,
    _,X,
    X,X,
  };
  glyphBitmaps[93] = new ColourizableBitmap(2, 7, glyph_closingBracket);


  unsigned char glyph_power[5*7] = 
  {
    _,_,X,_,_,
    _,X,_,X,_,
    X,_,_,_,X,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphBitmaps[94] = new ColourizableBitmap(5, 7, glyph_power);

  unsigned char glyph_underscore[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,X,X,X
  };
  glyphBitmaps[95] = new ColourizableBitmap(5, 7, glyph_underscore);

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

  unsigned char glyph_openingCurlyBrace[3*7] = 
  {
    _,X,X,
    _,X,_,
    _,X,_,
    X,X,_,
    _,X,_,
    _,X,_,
    _,X,X
  };
  glyphBitmaps[123] = new ColourizableBitmap(3, 7, glyph_openingCurlyBrace);

  unsigned char glyph_verticalLine[1*7] = 
  {
    X,
    X,
    X,
    X,
    X,
    X,
    X
  };
  glyphBitmaps[124] = new ColourizableBitmap(1, 7, glyph_verticalLine);

  unsigned char glyph_closingCurlyBrace[3*7] = 
  {
    X,X,_,
    _,X,_,
    _,X,_,
    _,X,X,
    _,X,_,
    _,X,_,
    X,X,_
  };
  glyphBitmaps[125] = new ColourizableBitmap(3, 7, glyph_closingCurlyBrace);

  unsigned char glyph_tilde[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,X,_,_,_,
    X,_,X,_,X,
    _,_,_,X,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphBitmaps[126] = new ColourizableBitmap(5, 7, glyph_tilde);

  unsigned char glyph_A[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphBitmaps[65] = new ColourizableBitmap(4, 7, glyph_A);

  unsigned char glyph_a[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    _,_,_,X,
    _,X,X,X,
    X,_,_,X,
    X,X,X,X
  };
  glyphBitmaps[97] = new ColourizableBitmap(4, 7, glyph_a);

  unsigned char glyph_B[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_
  };
  glyphBitmaps[66] = new ColourizableBitmap(4, 7, glyph_B);

  unsigned char glyph_b[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_
  };
  glyphBitmaps[98] = new ColourizableBitmap(4, 7, glyph_b);

  unsigned char glyph_C[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    _,X,X,X
  };
  glyphBitmaps[67] = new ColourizableBitmap(4, 7, glyph_C);

  unsigned char glyph_c[3*7] = 
  {
    _,_,_,
    _,_,_,
    _,X,X,
    X,_,_,
    X,_,_,
    X,_,_,
    _,X,X
  };
  glyphBitmaps[99] = new ColourizableBitmap(3, 7, glyph_c);

  unsigned char glyph_D[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_
  };
  glyphBitmaps[68] = new ColourizableBitmap(4, 7, glyph_D);

  unsigned char glyph_d[4*7] = 
  {
    _,_,_,X,
    _,_,_,X,
    _,X,X,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,X
  };
  glyphBitmaps[100] = new ColourizableBitmap(4, 7, glyph_d);

  unsigned char glyph_E[4*7] = 
  {
    X,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,_,
    X,_,_,_,
    X,X,X,X
  };
  glyphBitmaps[69] = new ColourizableBitmap(4, 7, glyph_E);

  unsigned char glyph_e[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,X,X,_,
    X,_,_,X,
    X,X,X,X,
    X,_,_,_,
    _,X,X,X
  };
  glyphBitmaps[101] = new ColourizableBitmap(4, 7, glyph_e);

  unsigned char glyph_F[4*7] = 
  {
    X,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_
  };
  glyphBitmaps[70] = new ColourizableBitmap(4, 7, glyph_F);

  unsigned char glyph_f[3*7] = 
  {
    _,X,X,
    _,X,_,
    _,X,_,
    X,X,X,
    _,X,_,
    _,X,_,
    _,X,_
  };
  glyphBitmaps[102] = new ColourizableBitmap(3, 7, glyph_f);

  unsigned char glyph_G[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,_,X,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[71] = new ColourizableBitmap(4, 7, glyph_G);

  unsigned char glyph_g[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,X,X,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    X,X,X,X
  };
  glyphBitmaps[103] = new ColourizableBitmap(4, 7, glyph_g);

  unsigned char glyph_H[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphBitmaps[72] = new ColourizableBitmap(4, 7, glyph_H);

  unsigned char glyph_h[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphBitmaps[104] = new ColourizableBitmap(4, 7, glyph_h);

  unsigned char glyph_I[1*7] = 
  {
    X,
    X,
    X,
    X,
    X,
    X,
    X
  };
  glyphBitmaps[73] = new ColourizableBitmap(1, 7, glyph_I);

  unsigned char glyph_i[1*7] = 
  {
    _,
    X,
    _,
    X,
    X,
    X,
    X
  };
  glyphBitmaps[105] = new ColourizableBitmap(1, 7, glyph_i);

  unsigned char glyph_J[4*7] = 
  {
    _,_,_,X,
    _,_,_,X,
    _,_,_,X,
    _,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[74] = new ColourizableBitmap(4, 7, glyph_J);

  unsigned char glyph_j[2*7] = 
  {
    _,X,
    _,_,
    _,X,
    _,X,
    _,X,
    _,X,
    X,X
  };
  glyphBitmaps[106] = new ColourizableBitmap(2, 7, glyph_j);

  unsigned char glyph_K[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,X,_,
    X,X,_,_,
    X,_,X,_,
    X,_,_,X,
    X,_,_,X
  };
  glyphBitmaps[75] = new ColourizableBitmap(4, 7, glyph_K);

  unsigned char glyph_k[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,_,_,X,
    X,_,X,_,
    X,X,_,_,
    X,_,X,_,
    X,_,_,X
  };
  glyphBitmaps[107] = new ColourizableBitmap(4, 7, glyph_k);

  unsigned char glyph_L[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,X,X,X
  };
  glyphBitmaps[76] = new ColourizableBitmap(4, 7, glyph_L);

  unsigned char glyph_l[2*7] = 
  {
    X,_,
    X,_,
    X,_,
    X,_,
    X,_,
    X,_,
    X,X
  };
  glyphBitmaps[108] = new ColourizableBitmap(2, 7, glyph_l);

  unsigned char glyph_M[5*7] = 
  {
    X,_,_,_,X,
    X,X,_,X,X,
    X,_,X,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X
  };
  glyphBitmaps[77] = new ColourizableBitmap(5, 7, glyph_M);

  unsigned char glyph_m[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,X,X,_,
    X,_,X,_,X,
    X,_,X,_,X,
    X,_,X,_,X,
    X,_,X,_,X
  };
  glyphBitmaps[109] = new ColourizableBitmap(5, 7, glyph_m);

  unsigned char glyph_N[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    X,X,_,_,X,
    X,_,X,_,X,
    X,_,_,X,X,
    X,_,_,_,X,
    X,_,_,_,X
  };
  glyphBitmaps[78] = new ColourizableBitmap(5, 7, glyph_N);

  unsigned char glyph_n[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphBitmaps[110] = new ColourizableBitmap(4, 7, glyph_n);

  unsigned char glyph_O[5*7] = 
  {
    _,X,X,X,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,X,X,_
  };
  glyphBitmaps[79] = new ColourizableBitmap(5, 7, glyph_O);

  unsigned char glyph_o[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[111] = new ColourizableBitmap(4, 7, glyph_o);

  unsigned char glyph_P[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_
  };
  glyphBitmaps[80] = new ColourizableBitmap(4, 7, glyph_P);

  unsigned char glyph_p[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    X,_,_,X,
    X,X,X,X,
    X,_,_,_,
    X,_,_,_
  };
  glyphBitmaps[112] = new ColourizableBitmap(4, 7, glyph_p);

  unsigned char glyph_Q[5*7] = 
  {
    _,X,X,X,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,X,_,X,
    X,_,_,X,_,
    _,X,X,_,X
  };
  glyphBitmaps[81] = new ColourizableBitmap(5, 7, glyph_Q);

  unsigned char glyph_q[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    _,_,_,X
  };
  glyphBitmaps[113] = new ColourizableBitmap(4, 7, glyph_q);

  unsigned char glyph_R[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_,
    X,X,_,_,
    X,_,X,_,
    X,_,_,X
  };
  glyphBitmaps[82] = new ColourizableBitmap(4, 7, glyph_R);

  unsigned char glyph_r[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,_,X,X,
    X,X,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_
  };
  glyphBitmaps[114] = new ColourizableBitmap(4, 7, glyph_r);

  unsigned char glyph_S[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    _,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphBitmaps[83] = new ColourizableBitmap(4, 7, glyph_S);

  unsigned char glyph_s[3*7] = 
  {
    _,_,_,
    _,_,_,
    _,X,X,
    X,_,_,
    _,X,_,
    _,_,X,
    X,X,_
  };
  glyphBitmaps[115] = new ColourizableBitmap(3, 7, glyph_s);

  unsigned char glyph_T[5*7] = 
  {
    X,X,X,X,X,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_
  };
  glyphBitmaps[84] = new ColourizableBitmap(5, 7, glyph_T);

  unsigned char glyph_t[3*7] = 
  {
    _,X,_,
    _,X,_,
    X,X,X,
    _,X,_,
    _,X,_,
    _,X,_,
    _,X,X
  };
  glyphBitmaps[116] = new ColourizableBitmap(3, 7, glyph_t);

  unsigned char glyph_U[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[85] = new ColourizableBitmap(4, 7, glyph_U);

  unsigned char glyph_u[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphBitmaps[117] = new ColourizableBitmap(4, 7, glyph_u);

  unsigned char glyph_V[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_
  };
  glyphBitmaps[86] = new ColourizableBitmap(5, 7, glyph_V);

  unsigned char glyph_v[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_
  };
  glyphBitmaps[118] = new ColourizableBitmap(5, 7, glyph_v);

  unsigned char glyph_W[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,X,_,X,
    X,X,_,X,X,
    X,_,_,_,X
  };
  glyphBitmaps[87] = new ColourizableBitmap(5, 7, glyph_W);

  unsigned char glyph_w[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,X,_,X,
    X,X,_,X,X,
    X,_,_,_,X
  };
  glyphBitmaps[119] = new ColourizableBitmap(5, 7, glyph_w);

  unsigned char glyph_X[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_,
    _,X,_,X,_,
    X,_,_,_,X,
    X,_,_,_,X
  };
  glyphBitmaps[88] = new ColourizableBitmap(5, 7, glyph_X);

  unsigned char glyph_x[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_,
    _,X,_,X,_,
    X,_,_,_,X
  };
  glyphBitmaps[120] = new ColourizableBitmap(5, 7, glyph_x);

  unsigned char glyph_Y[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    _,_,_,X,
    X,X,X,X
  };
  glyphBitmaps[89] = new ColourizableBitmap(4, 7, glyph_Y);

  unsigned char glyph_y[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    X,X,X,X
  };
  glyphBitmaps[121] = new ColourizableBitmap(4, 7, glyph_y);


  unsigned char glyph_Z[5*7] = 
  {
    X,X,X,X,X,
    _,_,_,_,X,
    _,_,_,X,_,
    _,_,X,_,_,
    _,X,_,_,_,
    X,_,_,_,_,
    X,X,X,X,X
  };
  glyphBitmaps[90] = new ColourizableBitmap(5, 7, glyph_Z);

  unsigned char glyph_z[3*7] = 
  {
    _,_,_,
    _,_,_,
    X,X,X,
    _,_,X,
    _,X,_,
    X,_,_,
    X,X,X
  };
  glyphBitmaps[122] = new ColourizableBitmap(3, 7, glyph_z);


  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphBitmaps[g] != NULL )
      glyphWidths[g] = glyphBitmaps[g]->getWidth();
  }
}
