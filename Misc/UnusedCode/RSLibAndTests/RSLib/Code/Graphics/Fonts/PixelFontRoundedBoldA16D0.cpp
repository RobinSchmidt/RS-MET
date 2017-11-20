using namespace RSLib;

rsPixelFontRoundedBoldA16D0::rsPixelFontRoundedBoldA16D0()
{
  ascent  = 16;
  descent = 0;

  createGlyphImages();
}

// create the glyphs:

void rsPixelFontRoundedBoldA16D0::createGlyphImages()
{
  rsUint8 _ = 0;
  rsUint8 X = 255;

  /*
  rsUint8 glyph_A[9*16] = 
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

  rsUint8 glyph_32[3*16] = 
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
  glyphImages[32] = new rsImageGray(3, 16, glyph_32);

  rsUint8 glyph040[4*16] = 
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
  glyphImages[40] = new rsImageGray(4, 16, glyph040);

  rsUint8 glyph041[4*16] = 
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
  glyphImages[41] = new rsImageGray(4, 16, glyph041);


  rsUint8 glyph045[9*16] = 
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
  glyphImages[45] = new rsImageGray(9, 16, glyph045);

  rsUint8 glyph_A[9*16] = 
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
  glyphImages[65] = new rsImageGray(9, 16, glyph_A);

  rsUint8 glyph_a[9*16] = 
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
  glyphImages[97] = new rsImageGray(9, 16, glyph_a);

  rsUint8 glyph_B[9*16] = 
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
  glyphImages[66] = new rsImageGray(9, 16, glyph_B);

  rsUint8 glyph_b[9*16] = 
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
  glyphImages[98] = new rsImageGray(9, 16, glyph_b);

  rsUint8 glyph_C[9*16] = 
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
  glyphImages[67] = new rsImageGray(9, 16, glyph_C);

  rsUint8 glyph_c[8*16] = 
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
  glyphImages[99] = new rsImageGray(8, 16, glyph_c);

  rsUint8 glyph_D[9*16] = 
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
  glyphImages[68] = new rsImageGray(9, 16, glyph_D);

  rsUint8 glyph_d[9*16] = 
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
  glyphImages[100] = new rsImageGray(9, 16, glyph_d);

  rsUint8 glyph_E[9*16] = 
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
  glyphImages[69] = new rsImageGray(9, 16, glyph_E);

  rsUint8 glyph_e[9*16] = 
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
  glyphImages[101] = new rsImageGray(9, 16, glyph_e);

  rsUint8 glyph_F[9*16] = 
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
  glyphImages[70] = new rsImageGray(9, 16, glyph_F);

  rsUint8 glyph_f[7*16] = 
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
  glyphImages[102] = new rsImageGray(7, 16, glyph_f);

  rsUint8 glyph_G[9*16] = 
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
  glyphImages[71] = new rsImageGray(9, 16, glyph_G);

  rsUint8 glyph_g[9*16] = 
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
  glyphImages[103] = new rsImageGray(9, 16, glyph_g);

  rsUint8 glyph_H[9*16] = 
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
  glyphImages[72] = new rsImageGray(9, 16, glyph_H);

  rsUint8 glyph_h[9*16] = 
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
  glyphImages[104] = new rsImageGray(9, 16, glyph_h);

  rsUint8 glyph_I[3*116] = 
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
  glyphImages[73] = new rsImageGray(3, 16, glyph_I);

  rsUint8 glyph_i[3*16] = 
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
  glyphImages[105] = new rsImageGray(3, 16, glyph_i);

  rsUint8 glyph_J[9*16] = 
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
  glyphImages[74] = new rsImageGray(9, 16, glyph_J);

  rsUint8 glyph_j[6*16] = 
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
  glyphImages[106] = new rsImageGray(6, 16, glyph_j);

  rsUint8 glyph_K[9*16] = 
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
  glyphImages[75] = new rsImageGray(9, 16, glyph_K);

  rsUint8 glyph_k[9*16] = 
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
  glyphImages[107] = new rsImageGray(9, 16, glyph_k);

  rsUint8 glyph_L[9*16] = 
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
  glyphImages[76] = new rsImageGray(9, 16, glyph_L);

  rsUint8 glyph_l[5*16] = 
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
  glyphImages[108] = new rsImageGray(5, 16, glyph_l);

  rsUint8 glyph_M[11*16] = 
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
  glyphImages[77] = new rsImageGray(11, 16, glyph_M);

  rsUint8 glyph_m[13*16] = 
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
  glyphImages[109] = new rsImageGray(13, 16, glyph_m);

  rsUint8 glyph_N[10*16] = 
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
  glyphImages[78] = new rsImageGray(10, 16, glyph_N);

  rsUint8 glyph_n[9*16] = 
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
  glyphImages[110] = new rsImageGray(9, 16, glyph_n);

  rsUint8 glyph_O[9*16] = 
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
  glyphImages[79] = new rsImageGray(9, 16, glyph_O);

  rsUint8 glyph_o[9*16] = 
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
  glyphImages[111] = new rsImageGray(9, 16, glyph_o);

  rsUint8 glyph_P[9*16] = 
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
  glyphImages[80] = new rsImageGray(9, 16, glyph_P);

  rsUint8 glyph_p[9*16] = 
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
  glyphImages[112] = new rsImageGray(9, 16, glyph_p);

  rsUint8 glyph_Q[11*16] = 
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
  glyphImages[81] = new rsImageGray(11, 16, glyph_Q);

  rsUint8 glyph_q[9*16] = 
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
  glyphImages[113] = new rsImageGray(9, 16, glyph_q);

  rsUint8 glyph_R[9*16] = 
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
  glyphImages[82] = new rsImageGray(9, 16, glyph_R);

  rsUint8 glyph_r[9*16] = 
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
  glyphImages[114] = new rsImageGray(9, 16, glyph_r);

  rsUint8 glyph_S[9*16] = 
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
  glyphImages[83] = new rsImageGray(9, 16, glyph_S);

  rsUint8 glyph_s[9*16] = 
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
  glyphImages[115] = new rsImageGray(9, 16, glyph_s);

  rsUint8 glyph_T[9*16] = 
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
  glyphImages[84] = new rsImageGray(9, 16, glyph_T);

  rsUint8 glyph_t[7*16] = 
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
  glyphImages[116] = new rsImageGray(7, 16, glyph_t);

  rsUint8 glyph_U[9*16] = 
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
  glyphImages[85] = new rsImageGray(9, 16, glyph_U);

  rsUint8 glyph_u[9*16] = 
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
  glyphImages[117] = new rsImageGray(9, 16, glyph_u);

  rsUint8 glyph_V[9*16] = 
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
  glyphImages[86] = new rsImageGray(9, 16, glyph_V);

  rsUint8 glyph_v[9*16] = 
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
  glyphImages[118] = new rsImageGray(9, 16, glyph_v);

  rsUint8 glyph_W[11*16] = 
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
  glyphImages[87] = new rsImageGray(11, 16, glyph_W);

  rsUint8 glyph_w[11*16] = 
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
  glyphImages[119] = new rsImageGray(11, 16, glyph_w);


  rsUint8 glyph_X[11*16] = 
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
  glyphImages[88] = new rsImageGray(11, 16, glyph_X);

  rsUint8 glyph_x[11*16] = 
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
  glyphImages[120] = new rsImageGray(11, 16, glyph_x);


  rsUint8 glyph_Y[11*16] = 
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
  glyphImages[89] = new rsImageGray(11, 16, glyph_Y);

  rsUint8 glyph_y[11*16] = 
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
  glyphImages[121] = new rsImageGray(11, 16, glyph_y);

  rsUint8 glyph_Z[9*16] = 
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
  glyphImages[90] = new rsImageGray(9, 16, glyph_Z);


  rsUint8 glyph_z[9*16] = 
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
  glyphImages[122] = new rsImageGray(9, 16, glyph_z);


  rsUint8 glyph_0[6*10] = 
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
  glyphImages[48] = new rsImageGray(6, 10, glyph_0);

  rsUint8 glyph_1[7*16] = 
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
  glyphImages[49] = new rsImageGray(7, 16, glyph_1);

  rsUint8 glyph_2[9*16] = 
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
  glyphImages[50] = new rsImageGray(9, 16, glyph_2);

  rsUint8 glyph_3[9*16] = 
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
  glyphImages[51] = new rsImageGray(9, 16, glyph_3);


  rsUint8 glyph_4[9*16] = 
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
  glyphImages[52] = new rsImageGray(9, 16, glyph_4);

  rsUint8 glyph_5[9*16] = 
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
  glyphImages[53] = new rsImageGray(9, 16, glyph_5);

  rsUint8 glyph_6[6*10] = 
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
  glyphImages[54] = new rsImageGray(6, 10, glyph_6);


  rsUint8 glyph_7[6*10] = 
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
  glyphImages[55] = new rsImageGray(6, 10, glyph_7);

  rsUint8 glyph_8[6*10] = 
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
  glyphImages[56] = new rsImageGray(6, 10, glyph_8);

  rsUint8 glyph_9[6*10] = 
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
  glyphImages[57] = new rsImageGray(6, 10, glyph_9);



  rsUint8 glyph_sharp[8*10] = 
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
  glyphImages[35] = new rsImageGray(8, 10, glyph_sharp);


  rsUint8 glyph_dot[2*10] = 
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
  glyphImages[46] = new rsImageGray(2, 10, glyph_dot);

  rsUint8 glyph_slash[7*10] = 
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
  glyphImages[47] = new rsImageGray(7, 10, glyph_slash);


  rsUint8 glyph_percent[7*10] = 
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
  glyphImages[37] = new rsImageGray(7, 10, glyph_percent);


  rsUint8 glyph_space[2*10] = 
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
  glyphImages[60] = new rsImageGray(2, 10, glyph_space);


  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphImages[g] != NULL )
      glyphWidths[g] = glyphImages[g]->getWidth();
  }
}
