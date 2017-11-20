#include "rojue_BitmapFontRoundedBoldA16D0.h"
using namespace rojue;

BitmapFontRoundedBoldA16D0::BitmapFontRoundedBoldA16D0()
{
  ascent  = 16;
  descent = 0;

  createGlyphBitmaps();
}

//-------------------------------------------------------------------------------------------------
// create the glyphs:



void BitmapFontRoundedBoldA16D0::createGlyphBitmaps()
{
  unsigned char _ = 0;
  unsigned char X = 255;

  /*
  unsigned char glyph_A[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_
  };
  */

  unsigned char glyph_32[3*16] = 
  {
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_
  };
  glyphBitmaps[32] = new ColourizableBitmap(3, 16, glyph_32);

  unsigned char glyph040[4*16] = 
  {
    _,_,X,X,
    _,X,X,X,
    _,X,X,X,
    X,X,X,_,
    X,X,X,_,
    X,X,X,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,X,_,
    X,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,_,X,X
  };
  glyphBitmaps[40] = new ColourizableBitmap(4, 16, glyph040);

  unsigned char glyph041[4*16] = 
  {
    X,X,_,_,
    X,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,X,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,X,X,X,
    _,X,X,X,
    _,X,X,_,
    X,X,X,_,
    X,X,_,_
  };
  glyphBitmaps[41] = new ColourizableBitmap(4, 16, glyph041);


  unsigned char glyph045[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_
  };
  glyphBitmaps[45] = new ColourizableBitmap(9, 16, glyph045);

  unsigned char glyph_A[9*16] = 
  {
    _,_,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    _,X,X,X,X,X,X,X,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X
  };
  glyphBitmaps[65] = new ColourizableBitmap(9, 16, glyph_A);

  unsigned char glyph_a[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X
  };
  glyphBitmaps[97] = new ColourizableBitmap(9, 16, glyph_a);

  unsigned char glyph_B[9*16] = 
  {
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_
  };
  glyphBitmaps[66] = new ColourizableBitmap(9, 16, glyph_B);

  unsigned char glyph_b[9*16] = 
  {
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_
  };
  glyphBitmaps[98] = new ColourizableBitmap(9, 16, glyph_b);

  unsigned char glyph_C[9*16] = 
  {
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,_,_,_,_,_,
    _,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X
  };
  glyphBitmaps[67] = new ColourizableBitmap(9, 16, glyph_C);

  unsigned char glyph_c[8*16] = 
  {
    _,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,
    _,_,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,
    X,X,X,X,_,_,_,_,
    X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,
    X,X,X,X,_,_,_,_,
    _,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X
  };
  glyphBitmaps[99] = new ColourizableBitmap(8, 16, glyph_c);

  unsigned char glyph_D[9*16] = 
  {
    X,X,X,X,X,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,_,X,X,X,X,X,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,_,_,_,_
  };
  glyphBitmaps[68] = new ColourizableBitmap(9, 16, glyph_D);

  unsigned char glyph_d[9*16] = 
  {
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X
  };
  glyphBitmaps[100] = new ColourizableBitmap(9, 16, glyph_d);

  unsigned char glyph_E[9*16] = 
  {
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X
  };
  glyphBitmaps[69] = new ColourizableBitmap(9, 16, glyph_E);

  unsigned char glyph_e[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X
  };
  glyphBitmaps[101] = new ColourizableBitmap(9, 16, glyph_e);

  unsigned char glyph_F[9*16] = 
  {
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
  };
  glyphBitmaps[70] = new ColourizableBitmap(9, 16, glyph_F);

  unsigned char glyph_f[7*16] = 
  {
    _,_,_,X,X,X,X,
    _,_,X,X,X,X,X,
    _,X,X,X,X,X,X,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    X,X,X,X,X,X,_,
    X,X,X,X,X,X,_,
    X,X,X,X,X,X,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_
  };
  glyphBitmaps[102] = new ColourizableBitmap(7, 16, glyph_f);

  unsigned char glyph_G[9*16] = 
  {
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,X,X,X,X,X,
    X,X,X,_,X,X,X,X,X,
    X,X,X,_,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X
  };
  glyphBitmaps[71] = new ColourizableBitmap(9, 16, glyph_G);

  unsigned char glyph_g[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_
  };
  glyphBitmaps[103] = new ColourizableBitmap(9, 16, glyph_g);

  unsigned char glyph_H[9*16] = 
  {
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X
  };
  glyphBitmaps[72] = new ColourizableBitmap(9, 16, glyph_H);

  unsigned char glyph_h[9*16] = 
  {
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X
  };
  glyphBitmaps[104] = new ColourizableBitmap(9, 16, glyph_h);

  unsigned char glyph_I[3*116] = 
  {
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X
  };
  glyphBitmaps[73] = new ColourizableBitmap(3, 16, glyph_I);

  unsigned char glyph_i[3*16] = 
  {
    X,X,X,
    X,X,X,
    X,X,X,
    _,_,_,
    _,_,_,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X,
    X,X,X
  };
  glyphBitmaps[105] = new ColourizableBitmap(3, 16, glyph_i);

  unsigned char glyph_J[9*16] = 
  {
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_
  };
  glyphBitmaps[74] = new ColourizableBitmap(9, 16, glyph_J);

  unsigned char glyph_j[6*16] = 
  {
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    _,_,_,X,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_
  };
  glyphBitmaps[106] = new ColourizableBitmap(6, 16, glyph_j);

  unsigned char glyph_K[9*16] = 
  {
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,_,_,_,
    X,X,X,X,X,_,_,_,_,
    X,X,X,X,X,_,_,_,_,
    X,X,X,X,X,X,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,_,X,X,X,X,_,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X
  };
  glyphBitmaps[75] = new ColourizableBitmap(9, 16, glyph_K);

  unsigned char glyph_k[9*16] = 
  {
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,_,X,X,X,X,X,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X
  };
  glyphBitmaps[107] = new ColourizableBitmap(9, 16, glyph_k);

  unsigned char glyph_L[9*16] = 
  {
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X
  };
  glyphBitmaps[76] = new ColourizableBitmap(9, 16, glyph_L);

  unsigned char glyph_l[5*16] = 
  {
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,_,_,
    X,X,X,X,X,
    X,X,X,X,X,
    _,X,X,X,X
  };
  glyphBitmaps[108] = new ColourizableBitmap(5, 16, glyph_l);

  unsigned char glyph_M[11*16] = 
  {
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,X,X,_,X,X,X,X,X,
    X,X,X,X,X,_,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,X,X,
    X,X,X,_,X,X,X,_,X,X,X,
    X,X,X,_,X,X,X,_,X,X,X,
    X,X,X,_,_,X,_,_,X,X,X,
    X,X,X,_,_,X,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X
  };
  glyphBitmaps[77] = new ColourizableBitmap(11, 16, glyph_M);

  unsigned char glyph_m[13*16] = 
  {
    _,_,_,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,_,_,
    X,X,X,_,X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,_,X,X,X
  };
  glyphBitmaps[109] = new ColourizableBitmap(13, 16, glyph_m);

  unsigned char glyph_N[10*16] = 
  {
    X,X,X,_,_,_,_,X,X,X,
    X,X,X,X,_,_,_,X,X,X,
    X,X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,_,_,X,X,X,
    X,X,X,X,X,_,_,X,X,X,
    X,X,X,X,X,X,_,X,X,X,
    X,X,X,_,X,X,_,X,X,X,
    X,X,X,_,X,X,X,X,X,X,
    X,X,X,_,_,X,X,X,X,X,
    X,X,X,_,_,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,X,
    X,X,X,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,X,X,X
  };
  glyphBitmaps[78] = new ColourizableBitmap(10, 16, glyph_N);

  unsigned char glyph_n[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,_,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X
  };
  glyphBitmaps[110] = new ColourizableBitmap(9, 16, glyph_n);

  unsigned char glyph_O[9*16] = 
  {
    _,_,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    _,X,X,X,X,X,X,X,_,
    X,X,X,X,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,_,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_
  };
  glyphBitmaps[79] = new ColourizableBitmap(9, 16, glyph_O);

  unsigned char glyph_o[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_
  };
  glyphBitmaps[111] = new ColourizableBitmap(9, 16, glyph_o);

  unsigned char glyph_P[9*16] = 
  {
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_
  };
  glyphBitmaps[80] = new ColourizableBitmap(9, 16, glyph_P);

  unsigned char glyph_p[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_
  };
  glyphBitmaps[112] = new ColourizableBitmap(9, 16, glyph_p);

  unsigned char glyph_Q[11*16] = 
  {
    _,_,X,X,X,X,X,_,_,_,_,
    _,X,X,X,X,X,X,X,_,_,_,
    _,X,X,X,X,X,X,X,_,_,_,
    X,X,X,X,_,X,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,_,_,_,X,X,X,_,_,
    X,X,X,X,_,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X,X,X
  };
  glyphBitmaps[81] = new ColourizableBitmap(11, 16, glyph_Q);

  unsigned char glyph_q[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X
  };
  glyphBitmaps[113] = new ColourizableBitmap(9, 16, glyph_q);

  unsigned char glyph_R[9*16] = 
  {
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,X,X,X,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,_,_,_,_,
    X,X,X,X,X,X,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,_,X,X,X,X,_,
    X,X,X,_,_,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X
  };
  glyphBitmaps[82] = new ColourizableBitmap(9, 16, glyph_R);

  unsigned char glyph_r[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,_,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_
  };
  glyphBitmaps[114] = new ColourizableBitmap(9, 16, glyph_r);

  unsigned char glyph_S[9*16] = 
  {
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,X,X,
    _,_,_,_,_,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_
  };
  glyphBitmaps[83] = new ColourizableBitmap(9, 16, glyph_S);

  unsigned char glyph_s[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_
  };
  glyphBitmaps[115] = new ColourizableBitmap(9, 16, glyph_s);

  unsigned char glyph_T[9*16] = 
  {
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_
  };
  glyphBitmaps[84] = new ColourizableBitmap(9, 16, glyph_T);

  unsigned char glyph_t[7*16] = 
  {
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,_,_,_,
    _,X,X,X,X,X,X,
    _,_,X,X,X,X,X,
    _,_,_,X,X,X,X
  };
  glyphBitmaps[116] = new ColourizableBitmap(7, 16, glyph_t);

  unsigned char glyph_U[9*16] = 
  {
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_
  };
  glyphBitmaps[85] = new ColourizableBitmap(9, 16, glyph_U);

  unsigned char glyph_u[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_
  };
  glyphBitmaps[117] = new ColourizableBitmap(9, 16, glyph_u);

  unsigned char glyph_V[9*16] = 
  {
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,_,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,_,X,_,_,_,_
  };
  glyphBitmaps[86] = new ColourizableBitmap(9, 16, glyph_V);

  unsigned char glyph_v[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,_,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,_,X,_,_,_,_
  };
  glyphBitmaps[118] = new ColourizableBitmap(9, 16, glyph_v);

  unsigned char glyph_W[11*16] = 
  {
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,X,_,_,X,X,X,
    X,X,X,_,_,X,_,_,X,X,X,
    X,X,X,_,X,X,X,_,X,X,X,
    X,X,X,_,X,X,X,_,X,X,X,
    X,X,X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,_,X,X,X,X,X,
    X,X,X,X,X,_,X,X,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X
  };
  glyphBitmaps[87] = new ColourizableBitmap(11, 16, glyph_W);

  unsigned char glyph_w[11*16] = 
  {
    _,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,X,_,_,X,X,X,
    X,X,X,_,_,X,_,_,X,X,X,
    X,X,X,_,X,X,X,_,X,X,X,
    X,X,X,_,X,X,X,_,X,X,X,
    X,X,X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,_,X,X,X,X,X,
    X,X,X,X,X,_,X,X,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X
  };
  glyphBitmaps[119] = new ColourizableBitmap(11, 16, glyph_w);


  unsigned char glyph_X[11*16] = 
  {
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    _,X,X,X,X,_,X,X,X,X,_,
    _,_,X,X,X,X,X,X,X,_,_,
    _,_,_,X,X,X,X,X,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,X,X,X,X,X,_,_,_,
    _,_,X,X,X,X,X,X,X,_,_,
    _,X,X,X,X,_,X,X,X,X,_,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X
  };
  glyphBitmaps[88] = new ColourizableBitmap(11, 16, glyph_X);

  unsigned char glyph_x[11*16] = 
  {
    _,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    _,X,X,X,X,_,X,X,X,X,_,
    _,_,X,X,X,X,X,X,X,_,_,
    _,_,_,X,X,X,X,X,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,X,X,X,X,X,_,_,_,
    _,_,X,X,X,X,X,X,X,_,_,
    _,X,X,X,X,_,X,X,X,X,_,
    X,X,X,X,_,_,_,X,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X
  };
  glyphBitmaps[120] = new ColourizableBitmap(11, 16, glyph_x);


  unsigned char glyph_Y[11*16] = 
  {
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    _,X,X,X,X,_,X,X,X,X,_,
    _,_,X,X,X,X,X,X,X,_,_,
    _,_,_,X,X,X,X,X,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_
  };
  glyphBitmaps[89] = new ColourizableBitmap(11, 16, glyph_Y);

  unsigned char glyph_y[11*16] = 
  {
    _,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,_,_,_,_,_,X,X,X,
    X,X,X,X,_,_,_,X,X,X,X,
    _,X,X,X,X,_,X,X,X,X,_,
    _,_,X,X,X,X,X,X,X,_,_,
    _,_,_,X,X,X,X,X,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_,
    _,_,_,_,X,X,X,_,_,_,_
  };
  glyphBitmaps[121] = new ColourizableBitmap(11, 16, glyph_y);

  unsigned char glyph_Z[9*16] = 
  {
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,X,X,X,X,
    _,_,_,_,X,X,X,X,_,
    _,_,_,X,X,X,X,_,_,
    _,_,X,X,X,X,_,_,_,
    _,X,X,X,X,_,_,_,_,
    X,X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X
  };
  glyphBitmaps[90] = new ColourizableBitmap(9, 16, glyph_Z);


  unsigned char glyph_z[9*16] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,X,X,X,X,
    _,_,_,_,X,X,X,X,_,
    _,_,_,X,X,X,X,_,_,
    _,_,X,X,X,X,_,_,_,
    _,X,X,X,X,_,_,_,_,
    X,X,X,X,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X
  };
  glyphBitmaps[122] = new ColourizableBitmap(9, 16, glyph_z);


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

  unsigned char glyph_1[7*16] = 
  {
    _,_,_,_,X,X,X,
    _,_,_,X,X,X,X,
    _,_,X,X,X,X,X,
    _,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,_,X,X,X
  };
  glyphBitmaps[49] = new ColourizableBitmap(7, 16, glyph_1);

  unsigned char glyph_2[9*16] = 
  {
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,X,X,X,X,X,
    _,_,_,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_,
    _,X,X,X,X,_,_,_,_,
    X,X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X
  };
  glyphBitmaps[50] = new ColourizableBitmap(9, 16, glyph_2);

  unsigned char glyph_3[9*16] = 
  {
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,X,X,X,X,X,_,
    _,_,_,X,X,X,X,_,_,
    _,_,_,X,X,X,X,X,_,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_
  };
  glyphBitmaps[51] = new ColourizableBitmap(9, 16, glyph_3);


  unsigned char glyph_4[9*16] = 
  {
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X
  };
  glyphBitmaps[52] = new ColourizableBitmap(9, 16, glyph_4);

  unsigned char glyph_5[9*16] = 
  {
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,X,X,X,X,
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,_,_,
    X,X,X,X,X,X,_,_,_
  };
  glyphBitmaps[53] = new ColourizableBitmap(9, 16, glyph_5);

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
    _,_,_,X,X,X,
    _,_,_,X,X,_,
    _,_,X,X,X,_,
    _,_,X,X,_,_,
    _,X,X,X,_,_,
    _,X,X,_,_,_,
    X,X,X,_,_,_
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

  unsigned char glyph_slash[7*10] = 
  {
    _,_,_,_,_,X,X,
    _,_,_,_,X,X,_,
    _,_,_,_,X,X,_,
    _,_,_,X,X,_,_,
    _,_,_,X,X,_,_,
    _,_,X,X,_,_,_,
    _,_,X,X,_,_,_,
    _,X,X,_,_,_,_,
    _,X,X,_,_,_,_,
    X,X,_,_,_,_,_
  };
  glyphBitmaps[47] = new ColourizableBitmap(7, 10, glyph_slash);


  unsigned char glyph_percent[7*10] = 
  {
    X,X,X,_,_,X,X,
    X,_,X,_,X,X,_,
    X,X,X,_,X,X,_,
    _,_,_,X,X,_,_,
    _,_,_,X,X,_,_,
    _,_,X,X,_,_,_,
    _,_,X,X,_,_,_,
    _,X,X,_,X,X,X,
    _,X,X,_,X,_,X,
    X,X,_,_,X,X,X
  };
  glyphBitmaps[37] = new ColourizableBitmap(7, 10, glyph_percent);


  unsigned char glyph_space[2*10] = 
  {
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
  };
  glyphBitmaps[60] = new ColourizableBitmap(2, 10, glyph_space);


  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphBitmaps[g] != NULL )
      glyphWidths[g] = glyphBitmaps[g]->getWidth();
  }
}
