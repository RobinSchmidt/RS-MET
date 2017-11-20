#include "rojue_BitmapFontRoundedTopLeftA10D0.h"
using namespace rojue;

BitmapFontRoundedTopLeftA10D0::BitmapFontRoundedTopLeftA10D0()
{
  ascent  = 10;
  descent = 0;

  createGlyphBitmaps();
}

//-------------------------------------------------------------------------------------------------
// create the glyphs:



void BitmapFontRoundedTopLeftA10D0::createGlyphBitmaps()
{
  float _ = 0.f;
  float X = 0.85f;
  float x = 0.6f; 

  float glyph_A[6*10] = 
  {
    _,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x
  };
  glyphBitmaps[65] = new ColourizableBitmap(6, 10, glyph_A);

  float glyph_a[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,_,
    x,x,x,x,X,x,
    _,_,_,_,X,x,
    _,X,X,X,X,x,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,x,
  };
  glyphBitmaps[97] = new ColourizableBitmap(6, 10, glyph_a);

  float glyph_B[6*10] = 
  {
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    X,x,x,x,X,_,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[66] = new ColourizableBitmap(6, 10, glyph_B);

  float glyph_b[6*10] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[98] = new ColourizableBitmap(6, 10, glyph_b);

  float glyph_C[6*10] = 
  {
    _,X,X,X,X,X,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,X,
    _,x,x,x,x,x
  };
  glyphBitmaps[67] = new ColourizableBitmap(6, 10, glyph_C);

  float glyph_c[5*10] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,X,X,X,X,
    X,x,x,x,x,
    X,x,_,_,_,
    X,x,_,_,_,
    X,x,_,_,_,
    X,x,_,_,_,
    X,X,X,X,X,
    _,x,x,x,x
  };
  glyphBitmaps[99] = new ColourizableBitmap(5, 10, glyph_c);

  float glyph_D[6*10] = 
  {
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,X,X,x,
    X,X,X,X,x,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[68] = new ColourizableBitmap(6, 10, glyph_D);

  float glyph_d[6*10] = 
  {
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,X,X,X,X,x,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,x
  };
  glyphBitmaps[100] = new ColourizableBitmap(6, 10, glyph_d);

  float glyph_E[6*10] = 
  {
    X,X,X,X,X,X,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,_,
    X,x,x,x,x,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,X,
    x,x,x,x,x,x
  };
  glyphBitmaps[69] = new ColourizableBitmap(6, 10, glyph_E);

  float glyph_e[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,X,X,X,X,X,
    _,x,x,x,x,x
  };
  glyphBitmaps[101] = new ColourizableBitmap(6, 10, glyph_e);

  float glyph_F[6*10] = 
  {
    X,X,X,X,X,X,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,_,
    X,x,x,x,x,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
  };
  glyphBitmaps[70] = new ColourizableBitmap(6, 10, glyph_F);

  float glyph_f[5*10] = 
  {
    _,_,X,X,X,
    _,X,x,x,x,
    _,X,x,_,_,
    _,X,x,_,_,
    X,X,X,X,_,  
    X,X,x,x,_,
    _,X,x,_,_,
    _,X,x,_,_,
    _,X,x,_,_,
    _,X,x,_,_
  };
  glyphBitmaps[102] = new ColourizableBitmap(5, 10, glyph_f);

  float glyph_G[6*10] = 
  {
    _,X,X,X,X,X,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,X,X,x,
    X,x,_,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,x,
  };
  glyphBitmaps[71] = new ColourizableBitmap(6, 10, glyph_G);

  float glyph_g[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,x,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,X,x,
    _,_,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,x,
  };
  glyphBitmaps[103] = new ColourizableBitmap(6, 10, glyph_g);

  float glyph_H[6*10] = 
  {
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
  };
  glyphBitmaps[72] = new ColourizableBitmap(6, 10, glyph_H);

  float glyph_h[6*10] = 
  {
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
  };
  glyphBitmaps[104] = new ColourizableBitmap(6, 10, glyph_h);

  float glyph_I[2*10] = 
  {
    X,x,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x
  };
  glyphBitmaps[73] = new ColourizableBitmap(2, 10, glyph_I);

  float glyph_i[2*10] = 
  {
    X,x,
    X,x,
    _,_,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x,
    X,x
  };
  glyphBitmaps[105] = new ColourizableBitmap(2, 10, glyph_i);

  float glyph_J[6*10] = 
  {
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,_,
  };
  glyphBitmaps[74] = new ColourizableBitmap(6, 10, glyph_J);

  float glyph_j[4*10] = 
  {
    _,_,_,_,
    _,_,X,x,
    _,_,X,x,
    _,_,_,_,
    _,_,X,x,
    _,_,X,x,
    _,_,X,x,
    _,_,X,x,
    X,X,X,x,
    x,x,x,_
  };
  glyphBitmaps[106] = new ColourizableBitmap(4, 10, glyph_j);

  float glyph_K[6*10] = 
  {
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,X,X,x,
    X,X,X,X,x,_,
    X,X,X,X,_,_,
    X,x,x,X,x,_,
    X,x,_,X,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x
  };
  glyphBitmaps[75] = new ColourizableBitmap(6, 10, glyph_K);

  float glyph_k[6*10] = 
  {
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,X,X,x,
    X,X,X,X,x,_,
    X,X,X,x,_,_,
    X,x,X,x,x,_,
    X,x,_,X,X,x,
    X,x,_,_,X,x,
  };
  glyphBitmaps[107] = new ColourizableBitmap(6, 10, glyph_k);

  float glyph_L[6*10] = 
  {
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,X,
    x,x,x,x,x,x
  };
  glyphBitmaps[76] = new ColourizableBitmap(6, 10, glyph_L);

  float glyph_l[3*10] = 
  {
    X,x,_,
    X,x,_,
    X,x,_,
    X,x,_,
    X,x,_,
    X,x,_,
    X,x,_,
    X,x,_,
    X,X,X,
    x,x,x
  };
  glyphBitmaps[108] = new ColourizableBitmap(3, 10, glyph_l);

  float glyph_M[9*10] = 
  {
    X,x,_,_,_,_,_,x,x,
    X,X,x,_,_,_,x,X,x,
    X,x,X,x,_,x,X,X,x,
    X,x,x,X,x,X,x,X,x,
    X,x,_,x,X,x,_,X,x,
    X,x,_,_,x,_,_,X,x,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,_,_,_,X,x
  };
  glyphBitmaps[77] = new ColourizableBitmap(9, 10, glyph_M);

  float glyph_m[8*10] = 
  {
    _,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,
    X,X,X,_,X,X,X,_,
    X,x,x,X,x,x,X,x,
    X,x,_,X,x,_,X,x,
    X,x,_,X,x,_,X,x,
    X,x,_,X,x,_,X,x,
    X,x,_,X,x,_,X,x,
    X,x,_,X,x,_,X,x,
    X,x,_,X,x,_,X,x
  };
  glyphBitmaps[109] = new ColourizableBitmap(8, 10, glyph_m);

  float glyph_N[7*10] = 
  {
    X,x,_,_,_,X,x,
    X,X,x,_,_,X,x,
    X,x,X,x,_,X,x,
    X,x,x,X,x,X,x,
    X,x,_,x,X,X,x,
    X,x,_,_,x,X,x,
    X,x,_,_,_,X,x,
    X,x,_,_,_,X,x,
    X,x,_,_,_,X,x,
    X,x,_,_,_,X,x
  };
  glyphBitmaps[78] = new ColourizableBitmap(7, 10, glyph_N);

  float glyph_n[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x
  };
  glyphBitmaps[110] = new ColourizableBitmap(6, 10, glyph_n);

  float glyph_O[6*10] = 
  {
    _,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,_
  };
  glyphBitmaps[79] = new ColourizableBitmap(6, 10, glyph_O);

  float glyph_o[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,_
  };
  glyphBitmaps[111] = new ColourizableBitmap(6, 10, glyph_o);

  float glyph_P[6*10] = 
  {
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    X,x,x,x,x,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_
  };
  glyphBitmaps[80] = new ColourizableBitmap(6, 10, glyph_P);

  float glyph_p[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    X,x,x,x,x,_,
    X,x,_,_,_,_,
    X,x,_,_,_,_
  };
  glyphBitmaps[112] = new ColourizableBitmap(6, 10, glyph_p);

  float glyph_Q[7*10] = 
  {
    _,X,X,X,X,_,_,
    X,x,x,x,X,x,_,
    X,x,_,_,X,x,_,
    X,x,_,_,X,x,_,
    X,x,_,_,X,x,_,
    X,x,_,_,X,x,_,
    X,x,_,_,X,x,_,
    X,x,_,X,X,x,_,
    X,X,X,X,X,X,X,
    _,x,x,x,x,x,x,
  };
  glyphBitmaps[81] = new ColourizableBitmap(7, 10, glyph_Q);

  float glyph_q[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,x,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x
  };
  glyphBitmaps[113] = new ColourizableBitmap(6, 10, glyph_q);

  float glyph_R[6*10] = 
  {
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    X,X,X,x,x,_,
    X,X,X,x,_,_,
    X,x,x,X,x,_,
    X,x,_,x,X,x,
    X,x,_,_,x,X
  };
  glyphBitmaps[82] = new ColourizableBitmap(6, 10, glyph_R);

  float glyph_r[5*10] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,_,X,X,
    X,x,X,x,x,
    X,x,_,_,_,
    X,x,_,_,_,
    X,x,_,_,_,
    X,x,_,_,_,
    X,x,_,_,_,
    X,x,_,_,_
  };
  glyphBitmaps[114] = new ColourizableBitmap(5, 10, glyph_r);

  float glyph_S[6*10] = 
  {
    _,X,X,X,X,X,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,_,
    _,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[83] = new ColourizableBitmap(6, 10, glyph_S);

  float glyph_s[5*10] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,X,X,X,X,
    X,x,x,x,x,
    X,x,_,_,_,
    X,X,X,X,_,
    _,x,x,X,x,
    _,_,_,X,x,
    X,X,X,X,x,
    x,x,x,x,_
  };
  glyphBitmaps[115] = new ColourizableBitmap(5, 10, glyph_s);

  float glyph_T[6*10] = 
  {
    X,X,X,X,X,X,
    x,x,X,x,x,x,
    _,_,X,x,_,_,
    _,_,X,x,_,_,
    _,_,X,x,_,_,
    _,_,X,x,_,_,
    _,_,X,x,_,_,
    _,_,X,x,_,_,
    _,_,X,x,_,_,
    _,_,X,x,_,_
  };
  glyphBitmaps[84] = new ColourizableBitmap(6, 10, glyph_T);

  float glyph_t[4*10] = 
  {
    _,X,x,_,
    _,X,x,_,
    X,X,X,x,
    x,X,x,x,
    _,X,x,_,
    _,X,x,_,
    _,X,x,_,
    _,X,x,_,
    _,X,X,x,
    _,x,x,x
  };
  glyphBitmaps[116] = new ColourizableBitmap(4, 10, glyph_t);

  float glyph_U[6*10] = 
  {
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,x
  };
  glyphBitmaps[85] = new ColourizableBitmap(6, 10, glyph_U);

  float glyph_u[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,x
  };
  glyphBitmaps[117] = new ColourizableBitmap(6, 10, glyph_u);

  float glyph_V[7*10] = 
  {
    X,x,_,_,_,X,x,
    X,x,_,_,_,X,x,
    X,x,_,_,_,X,x,
    X,x,_,_,_,X,x,
    X,X,_,_,_,X,x,
    _,X,x,_,X,X,_,
    _,X,X,x,X,x,_,
    _,_,X,X,X,_,_,
    _,_,x,X,x,_,_,
    _,_,_,x,_,_,_
  };
  glyphBitmaps[86] = new ColourizableBitmap(7, 10, glyph_V);

  float glyph_v[7*10] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,x,_,_,_,X,x,
    X,x,_,_,_,X,x,
    X,X,_,_,_,X,x,
    _,X,x,_,X,X,_,
    _,X,X,x,X,x,_,
    _,_,X,X,X,_,_,
    _,_,x,X,x,_,_,
    _,_,_,x,_,_,_
  };
  glyphBitmaps[118] = new ColourizableBitmap(7, 10, glyph_v);

  float glyph_W[9*10] = 
  {
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,x,_,_,X,x,
    X,x,_,x,X,x,_,X,x,
    X,x,x,X,x,X,x,X,x,
    X,x,X,x,_,x,X,X,x,
    X,X,x,_,_,_,x,X,x,
    X,x,_,_,_,_,_,X,x
  };
  glyphBitmaps[87] = new ColourizableBitmap(9, 10, glyph_W);

  float glyph_w[9*10] = 
  {
    _,_,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,_,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,_,_,_,X,x,
    X,x,_,_,x,_,_,X,x,
    X,x,_,x,X,x,_,X,x,
    X,x,x,X,x,X,x,X,x,
    X,x,X,x,_,x,X,X,x,
    X,X,x,_,_,_,x,X,x,
    X,x,_,_,_,_,_,X,x
  };
  glyphBitmaps[119] = new ColourizableBitmap(9, 10, glyph_w);


  float glyph_X[7*10] = 
  {
    X,X,_,_,_,X,X,
    _,X,x,_,X,X,_,
    _,X,x,_,X,x,_,
    _,X,X,x,X,x,_,
    _,_,X,x,X,_,_,
    _,_,X,X,X,_,_,
    _,X,X,x,X,x,_,
    _,X,x,_,X,x,_,
    X,X,x,_,X,X,x,
    X,x,_,_,_,X,x
  };
  glyphBitmaps[88] = new ColourizableBitmap(7, 10, glyph_X);

  float glyph_x[7*10] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,x,_,_,_,x,X,
    x,X,x,_,x,X,x,
    _,x,X,x,X,x,_,
    _,_,x,X,x,_,_,
    _,_,X,x,X,_,_,
    _,X,x,_,x,X,_,
    X,x,x,_,x,x,X,
    X,x,_,_,_,x,X
  };
  glyphBitmaps[120] = new ColourizableBitmap(7, 10, glyph_x);


  float glyph_Y[6*10] = 
  {
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[89] = new ColourizableBitmap(6, 10, glyph_Y);

  float glyph_y[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,X,x,
    _,_,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[121] = new ColourizableBitmap(6, 10, glyph_y);

  float glyph_Z[7*10] = 
  {
    X,X,X,X,X,X,x,
    x,x,x,x,x,X,x,
    _,_,_,_,_,X,x,
    _,_,_,_,X,x,x,
    _,_,_,X,x,x,_,
    _,_,X,x,x,_,_,
    _,X,x,x,_,_,_,
    X,x,x,_,_,_,_,
    X,X,X,X,X,X,X,
    x,x,x,x,x,x,x,
  };
  glyphBitmaps[90] = new ColourizableBitmap(7, 10, glyph_Z);


  float glyph_z[6*10] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,x,
    x,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,X,x,x,
    _,_,X,x,x,_,
    _,X,x,x,_,_,
    X,X,X,X,X,X,
    x,x,x,x,x,x
  };
  glyphBitmaps[122] = new ColourizableBitmap(6, 10, glyph_z);


  float glyph_0[6*10] = 
  {
    _,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,_
  };
  glyphBitmaps[48] = new ColourizableBitmap(6, 10, glyph_0);

  float glyph_1[4*10] = 
  {
    _,_,X,x,
    _,X,X,x,
    x,x,X,x,
    _,_,X,x,
    _,_,X,x,
    _,_,X,x,
    _,_,X,x,
    _,_,X,x,
    _,_,X,x,
    _,_,X,x
  };
  glyphBitmaps[49] = new ColourizableBitmap(4, 10, glyph_1);

  float glyph_2[6*10] = 
  {
    X,X,X,X,X,_,
    x,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,X,X,X,x,
    _,X,x,x,x,_,
    X,x,x,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,X,
    x,x,x,x,x,x
  };
  glyphBitmaps[50] = new ColourizableBitmap(6, 10, glyph_2);

  float glyph_3[6*10] = 
  {
    X,X,X,X,X,_,
    x,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,X,X,X,X,x,
    _,x,x,x,X,_,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[51] = new ColourizableBitmap(6, 10, glyph_3);


  float glyph_4[6*10] = 
  {
    X,x,_,_,_,_,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x
  };
  glyphBitmaps[52] = new ColourizableBitmap(6, 10, glyph_4);

  float glyph_5[6*10] = 
  {
    X,X,X,X,X,X,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,X,X,X,X,_,
    x,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    _,_,_,X,X,x,
    X,X,X,X,x,_,
    x,x,x,x,_,_
  };
  glyphBitmaps[53] = new ColourizableBitmap(6, 10, glyph_5);

  float glyph_6[6*10] = 
  {
    _,X,X,X,X,X,
    X,x,x,x,x,x,
    X,x,_,_,_,_,
    X,x,_,_,_,_,
    X,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,_
  };
  glyphBitmaps[54] = new ColourizableBitmap(6, 10, glyph_6);


  float glyph_7[6*10] = 
  {
    X,X,X,X,X,x,
    x,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,X,X,x,
    _,_,_,X,x,_,
    _,_,X,X,x,_,
    _,_,X,x,_,_,
    _,X,X,x,_,_,
    _,X,x,_,_,_,
    X,X,x,_,_,_
  };
  glyphBitmaps[55] = new ColourizableBitmap(6, 10, glyph_7);

  float glyph_8[6*10] = 
  {
    _,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    _,X,X,X,X,x,
    X,x,x,x,X,_,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,x,_
  };
  glyphBitmaps[56] = new ColourizableBitmap(6, 10, glyph_8);

  float glyph_9[6*10] = 
  {
    _,X,X,X,X,_,
    X,x,x,x,X,x,
    X,x,_,_,X,x,
    X,x,_,_,X,x,
    X,X,X,X,X,x,
    _,x,x,x,X,x,
    _,_,_,_,X,x,
    _,_,_,_,X,x,
    X,X,X,X,X,x,
    x,x,x,x,x,_
  };
  glyphBitmaps[57] = new ColourizableBitmap(6, 10, glyph_9);



  float glyph_sharp[8*10] = 
  {
    _,X,x,_,_,X,x,_,
    _,X,x,_,_,X,x,_,
    X,X,X,X,X,X,X,X,
    x,X,x,x,x,X,x,x,
    _,X,x,_,_,X,x,_,
    _,X,x,_,_,X,x,_,
    X,X,X,X,X,X,X,X,
    x,X,x,x,x,X,x,x,
    _,X,x,_,_,X,x,_,
    _,X,x,_,_,X,x,_,
  };
  glyphBitmaps[35] = new ColourizableBitmap(8, 10, glyph_sharp);


  float glyph_dot[2*10] = 
  {
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    X,x,
    X,x
  };
  glyphBitmaps[46] = new ColourizableBitmap(2, 10, glyph_dot);

  float glyph_slash[7*10] = 
  {
    _,_,_,_,_,X,x,
    _,_,_,_,X,X,_,
    _,_,_,_,X,x,_,
    _,_,_,X,X,_,_,
    _,_,_,X,x,_,_,
    _,_,X,X,_,_,_,
    _,_,X,x,_,_,_,
    _,X,X,_,_,_,_,
    _,X,x,_,_,_,_,
    X,X,_,_,_,_,_
  };
  glyphBitmaps[47] = new ColourizableBitmap(7, 10, glyph_slash);


  float glyph_percent[7*10] = 
  {
    X,X,X,_,_,X,x,
    X,_,X,_,X,X,_,
    X,X,X,_,X,x,_,
    _,_,_,X,X,_,_,
    _,_,_,X,x,_,_,
    _,_,X,X,_,_,_,
    _,_,X,x,_,_,_,
    _,X,X,_,X,X,X,
    _,X,x,_,X,_,X,
    X,X,_,_,X,X,X
  };
  glyphBitmaps[37] = new ColourizableBitmap(7, 10, glyph_percent);


  float glyph_space[2*10] = 
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
