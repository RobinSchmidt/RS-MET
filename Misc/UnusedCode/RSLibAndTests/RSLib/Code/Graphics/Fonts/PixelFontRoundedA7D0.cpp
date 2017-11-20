using namespace RSLib;

rsPixelFontRoundedA7D0::rsPixelFontRoundedA7D0()
{
  ascent  = 7;
  descent = 0;

  createGlyphImages();
}

// create the glyphs:

void rsPixelFontRoundedA7D0::createGlyphImages()
{
  rsUint8 _ = 0;
  rsUint8 X = 255;

  /*
  rsUint8 glyph_dummy[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphImages[60] = new rsImageGray(5, 7, glyph_glyph_dummy);
  */


  rsUint8 glyph_space[2*7] = 
  {
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_
  };
  glyphImages[32] = new rsImageGray(2, 7, glyph_space);

  rsUint8 glyph_exclamationMark[1*7] = 
  {
    X,
    X,
    X,
    X,
    _,
    _,
    X
  };
  glyphImages[33] = new rsImageGray(1, 7, glyph_exclamationMark);

  rsUint8 glyph_quotes[3*7] = 
  {
    X,_,X,
    X,_,X,
    X,_,X,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
  };
  glyphImages[34] = new rsImageGray(3, 7, glyph_quotes);

  rsUint8 glyph_sharp[5*7] = 
  {
    _,X,_,X,_,
    _,X,_,X,_,
    X,X,X,X,X,
    _,X,_,X,_,
    X,X,X,X,X,
    _,X,_,X,_,
    _,X,_,X,_
  };
  glyphImages[35] = new rsImageGray(5, 7, glyph_sharp);

  rsUint8 glyph_dollar[3*7] = 
  {
    _,X,_,
    X,X,X,
    X,_,_,
    _,X,_,
    _,_,X,
    X,X,X,
    _,X,_,
  };
  glyphImages[36] = new rsImageGray(3, 7, glyph_dollar);

  rsUint8 glyph_percent[3*7] = 
  {
    X,_,X,
    _,_,X,
    _,X,X,
    _,X,_,
    X,X,_,
    X,_,_,
    X,_,X,
  };
  glyphImages[37] = new rsImageGray(3, 7, glyph_percent);

  rsUint8 glyph_ampersand[5*7] = 
  {
    _,X,X,_,_,
    X,_,_,X,_,
    X,_,X,_,_,
    _,X,X,_,X,
    X,_,X,X,_,
    X,_,_,X,_,
    _,X,X,_,X
  };
  glyphImages[38] = new rsImageGray(5, 7, glyph_ampersand);

  rsUint8 glyph_singlequote[1*7] = 
  {
    X,
    X,
    X,
    _,
    _,
    _,
    _,
  };
  glyphImages[39] = new rsImageGray(1, 7, glyph_singlequote);

  rsUint8 glyph_openingBrace[2*7] = 
  {
    _,X,
    X,X,
    X,_,
    X,_,
    X,_,
    X,X,
    _,X,
  };
  glyphImages[40] = new rsImageGray(2, 7, glyph_openingBrace);

  rsUint8 glyph_closingBrace[2*7] = 
  {
    X,_,
    X,X,
    _,X,
    _,X,
    _,X,
    X,X,
    X,_,
  };
  glyphImages[41] = new rsImageGray(2, 7, glyph_closingBrace);

  rsUint8 glyph_multiply[5*7] = 
  {
    _,_,_,_,_,
    _,_,X,_,_,
    X,X,X,X,X,
    _,_,X,_,_,
    _,X,_,X,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphImages[42] = new rsImageGray(5, 7, glyph_multiply);

  rsUint8 glyph_plus[5*7] = 
  {
    _,_,_,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    X,X,X,X,X,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,_,_,_
  };
  glyphImages[43] = new rsImageGray(5, 7, glyph_plus);

  rsUint8 glyph_comma[1*7] = 
  {
    _,
    _,
    _,
    _,
    X,
    X,
    X,
  };
  glyphImages[44] = new rsImageGray(1, 7, glyph_comma);

  rsUint8 glyph_minus[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    _,_,_,_,
    _,_,_,_,
    _,_,_,_
  };
  glyphImages[45] = new rsImageGray(4, 7, glyph_minus);

  rsUint8 glyph_dot[1*7] = 
  {
    _,
    _,
    _,
    _,
    _,
    _,
    X
  };
  glyphImages[46] = new rsImageGray(1, 7, glyph_dot);

  rsUint8 glyph_slash[3*7] = 
  {
    _,_,X,
    _,_,X,
    _,X,X,
    _,X,_,
    X,X,_,
    X,_,_,
    X,_,_,
  };
  glyphImages[47] = new rsImageGray(3, 7, glyph_slash);

  rsUint8 glyph_0[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[48] = new rsImageGray(4, 7, glyph_0);

  rsUint8 glyph_1[2*7] = 
  {
    _,X,
    X,X,
    _,X,
    _,X,
    _,X,
    _,X,
    _,X
  };
  glyphImages[49] = new rsImageGray(2, 7, glyph_1);

  rsUint8 glyph_2[4*70] = 
  {
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    _,_,X,_,
    _,X,_,_,
    X,_,_,_,
    X,X,X,X
  };
  glyphImages[50] = new rsImageGray(4, 7, glyph_2);

  rsUint8 glyph_3[4*7] = 
  {
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphImages[51] = new rsImageGray(4, 7, glyph_3);


  rsUint8 glyph_4[4*7] = 
  {
    X,_,_,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,X,
    _,_,_,X,
    _,_,_,X,
    _,_,_,X
  };
  glyphImages[52] = new rsImageGray(4, 7, glyph_4);


  rsUint8 glyph_5[4*7] = 
  {
    X,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphImages[53] = new rsImageGray(4, 7, glyph_5);

  rsUint8 glyph_6[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[54] = new rsImageGray(4, 7, glyph_6);

  rsUint8 glyph_7[4*7] = 
  {
    X,X,X,X,
    _,_,_,X,
    _,_,X,_,
    _,_,X,_,
    _,X,_,_,
    _,X,_,_,
    X,_,_,_
  };
  glyphImages[55] = new rsImageGray(4, 7, glyph_7);

  rsUint8 glyph_8[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[56] = new rsImageGray(4, 7, glyph_8);

  rsUint8 glyph_9[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    _,X,X,X,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphImages[57] = new rsImageGray(4, 7, glyph_9);

  rsUint8 glyph_colon[1*7] = 
  {
    _,
    X,
    X,
    _,
    X,
    X,
    _,
  };
  glyphImages[58] = new rsImageGray(1, 7, glyph_colon);


  rsUint8 glyph_semicolon[1*7] = 
  {
    _,
    X,
    X,
    _,
    X,
    X,
    X,
  };
  glyphImages[59] = new rsImageGray(1, 7, glyph_semicolon);

  rsUint8 glyph_openingAngleBrace[4*7] = 
  {
    _,_,_,X,
    _,_,X,_,
    _,X,_,_,
    X,_,_,_,
    _,X,_,_,
    _,_,X,_,
    _,_,_,X,
  };
  glyphImages[60] = new rsImageGray(4, 7, glyph_openingAngleBrace);

  rsUint8 glyph_equalsSign[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    _,_,_,_,
    X,X,X,X,
    _,_,_,_,
    _,_,_,_,
  };
  glyphImages[61] = new rsImageGray(4, 7, glyph_equalsSign);

  rsUint8 glyph_closingAngleBrace[4*7] = 
  {
    X,_,_,_,
    _,X,_,_,
    _,_,X,_,
    _,_,_,X,
    _,_,X,_,
    _,X,_,_,
    X,_,_,_,
  };
  glyphImages[62] = new rsImageGray(4, 7, glyph_closingAngleBrace);

  rsUint8 glyph_questionMark[3*7] = 
  {
    X,X,X,
    _,_,X,
    _,_,X,
    X,X,X,
    X,_,_,
    _,_,_,
    X,_,_,
  };
  glyphImages[63] = new rsImageGray(3, 7, glyph_questionMark);


  rsUint8 glyph_at[5*7] = 
  {
    _,X,X,X,X,
    X,_,_,_,X,
    X,_,X,X,X,
    X,_,X,_,X,
    X,_,X,X,X,
    X,_,_,_,_,
    _,X,X,X,X
  };
  glyphImages[64] = new rsImageGray(5, 7, glyph_at);

  rsUint8 glyph_openingBracket[2*7] = 
  {
    X,X,
    X,_,
    X,_,
    X,_,
    X,_,
    X,_,
    X,X,
  };
  glyphImages[91] = new rsImageGray(2, 7, glyph_openingBracket);

  rsUint8 glyph_backslash[3*7] = 
  {
    X,_,_,
    X,_,_,
    X,X,_,
    _,X,_,
    _,X,X,
    _,_,X,
    _,_,X,
  };
  glyphImages[92] = new rsImageGray(3, 7, glyph_backslash);

  rsUint8 glyph_closingBracket[2*7] = 
  {
    X,X,
    _,X,
    _,X,
    _,X,
    _,X,
    _,X,
    X,X,
  };
  glyphImages[93] = new rsImageGray(2, 7, glyph_closingBracket);


  rsUint8 glyph_power[5*7] = 
  {
    _,_,X,_,_,
    _,X,_,X,_,
    X,_,_,_,X,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphImages[94] = new rsImageGray(5, 7, glyph_power);

  rsUint8 glyph_underscore[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,X,X,X
  };
  glyphImages[95] = new rsImageGray(5, 7, glyph_underscore);

  rsUint8 glyph_96[3*10] = 
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
  glyphImages[96] = new rsImageGray(3, 10, glyph_96);

  rsUint8 glyph_openingCurlyBrace[3*7] = 
  {
    _,X,X,
    _,X,_,
    _,X,_,
    X,X,_,
    _,X,_,
    _,X,_,
    _,X,X
  };
  glyphImages[123] = new rsImageGray(3, 7, glyph_openingCurlyBrace);

  rsUint8 glyph_verticalLine[1*7] = 
  {
    X,
    X,
    X,
    X,
    X,
    X,
    X
  };
  glyphImages[124] = new rsImageGray(1, 7, glyph_verticalLine);

  rsUint8 glyph_closingCurlyBrace[3*7] = 
  {
    X,X,_,
    _,X,_,
    _,X,_,
    _,X,X,
    _,X,_,
    _,X,_,
    X,X,_
  };
  glyphImages[125] = new rsImageGray(3, 7, glyph_closingCurlyBrace);

  rsUint8 glyph_tilde[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,X,_,_,_,
    X,_,X,_,X,
    _,_,_,X,_,
    _,_,_,_,_,
    _,_,_,_,_
  };
  glyphImages[126] = new rsImageGray(5, 7, glyph_tilde);

  rsUint8 glyph_A[4*7] = 
  {
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphImages[65] = new rsImageGray(4, 7, glyph_A);

  rsUint8 glyph_a[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    _,_,_,X,
    _,X,X,X,
    X,_,_,X,
    X,X,X,X
  };
  glyphImages[97] = new rsImageGray(4, 7, glyph_a);

  rsUint8 glyph_B[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_
  };
  glyphImages[66] = new rsImageGray(4, 7, glyph_B);

  rsUint8 glyph_b[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_
  };
  glyphImages[98] = new rsImageGray(4, 7, glyph_b);

  rsUint8 glyph_C[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    _,X,X,X
  };
  glyphImages[67] = new rsImageGray(4, 7, glyph_C);

  rsUint8 glyph_c[3*7] = 
  {
    _,_,_,
    _,_,_,
    _,X,X,
    X,_,_,
    X,_,_,
    X,_,_,
    _,X,X
  };
  glyphImages[99] = new rsImageGray(3, 7, glyph_c);

  rsUint8 glyph_D[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_
  };
  glyphImages[68] = new rsImageGray(4, 7, glyph_D);

  rsUint8 glyph_d[4*7] = 
  {
    _,_,_,X,
    _,_,_,X,
    _,X,X,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,X
  };
  glyphImages[100] = new rsImageGray(4, 7, glyph_d);

  rsUint8 glyph_E[4*7] = 
  {
    X,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,_,
    X,_,_,_,
    X,X,X,X
  };
  glyphImages[69] = new rsImageGray(4, 7, glyph_E);

  rsUint8 glyph_e[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,X,X,_,
    X,_,_,X,
    X,X,X,X,
    X,_,_,_,
    _,X,X,X
  };
  glyphImages[101] = new rsImageGray(4, 7, glyph_e);

  rsUint8 glyph_F[4*7] = 
  {
    X,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_
  };
  glyphImages[70] = new rsImageGray(4, 7, glyph_F);

  rsUint8 glyph_f[3*7] = 
  {
    _,X,X,
    _,X,_,
    _,X,_,
    X,X,X,
    _,X,_,
    _,X,_,
    _,X,_
  };
  glyphImages[102] = new rsImageGray(3, 7, glyph_f);

  rsUint8 glyph_G[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    X,_,X,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[71] = new rsImageGray(4, 7, glyph_G);

  rsUint8 glyph_g[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,X,X,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    X,X,X,X
  };
  glyphImages[103] = new rsImageGray(4, 7, glyph_g);

  rsUint8 glyph_H[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphImages[72] = new rsImageGray(4, 7, glyph_H);

  rsUint8 glyph_h[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphImages[104] = new rsImageGray(4, 7, glyph_h);

  rsUint8 glyph_I[1*7] = 
  {
    X,
    X,
    X,
    X,
    X,
    X,
    X
  };
  glyphImages[73] = new rsImageGray(1, 7, glyph_I);

  rsUint8 glyph_i[1*7] = 
  {
    _,
    X,
    _,
    X,
    X,
    X,
    X
  };
  glyphImages[105] = new rsImageGray(1, 7, glyph_i);

  rsUint8 glyph_J[4*7] = 
  {
    _,_,_,X,
    _,_,_,X,
    _,_,_,X,
    _,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[74] = new rsImageGray(4, 7, glyph_J);

  rsUint8 glyph_j[2*7] = 
  {
    _,X,
    _,_,
    _,X,
    _,X,
    _,X,
    _,X,
    X,X
  };
  glyphImages[106] = new rsImageGray(2, 7, glyph_j);

  rsUint8 glyph_K[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,X,_,
    X,X,_,_,
    X,_,X,_,
    X,_,_,X,
    X,_,_,X
  };
  glyphImages[75] = new rsImageGray(4, 7, glyph_K);

  rsUint8 glyph_k[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,_,_,X,
    X,_,X,_,
    X,X,_,_,
    X,_,X,_,
    X,_,_,X
  };
  glyphImages[107] = new rsImageGray(4, 7, glyph_k);

  rsUint8 glyph_L[4*7] = 
  {
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_,
    X,X,X,X
  };
  glyphImages[76] = new rsImageGray(4, 7, glyph_L);

  rsUint8 glyph_l[2*7] = 
  {
    X,_,
    X,_,
    X,_,
    X,_,
    X,_,
    X,_,
    X,X
  };
  glyphImages[108] = new rsImageGray(2, 7, glyph_l);

  rsUint8 glyph_M[5*7] = 
  {
    X,_,_,_,X,
    X,X,_,X,X,
    X,_,X,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X
  };
  glyphImages[77] = new rsImageGray(5, 7, glyph_M);

  rsUint8 glyph_m[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,X,X,_,
    X,_,X,_,X,
    X,_,X,_,X,
    X,_,X,_,X,
    X,_,X,_,X
  };
  glyphImages[109] = new rsImageGray(5, 7, glyph_m);

  rsUint8 glyph_N[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    X,X,_,_,X,
    X,_,X,_,X,
    X,_,_,X,X,
    X,_,_,_,X,
    X,_,_,_,X
  };
  glyphImages[78] = new rsImageGray(5, 7, glyph_N);

  rsUint8 glyph_n[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X
  };
  glyphImages[110] = new rsImageGray(4, 7, glyph_n);

  rsUint8 glyph_O[5*7] = 
  {
    _,X,X,X,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,X,X,_
  };
  glyphImages[79] = new rsImageGray(5, 7, glyph_O);

  rsUint8 glyph_o[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    _,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[111] = new rsImageGray(4, 7, glyph_o);

  rsUint8 glyph_P[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_
  };
  glyphImages[80] = new rsImageGray(4, 7, glyph_P);

  rsUint8 glyph_p[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    X,_,_,X,
    X,X,X,X,
    X,_,_,_,
    X,_,_,_
  };
  glyphImages[112] = new rsImageGray(4, 7, glyph_p);

  rsUint8 glyph_Q[5*7] = 
  {
    _,X,X,X,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,X,_,X,
    X,_,_,X,_,
    _,X,X,_,X
  };
  glyphImages[81] = new rsImageGray(5, 7, glyph_Q);

  rsUint8 glyph_q[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,X,X,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    _,_,_,X
  };
  glyphImages[113] = new rsImageGray(4, 7, glyph_q);

  rsUint8 glyph_R[4*7] = 
  {
    X,X,X,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,_,
    X,X,_,_,
    X,_,X,_,
    X,_,_,X
  };
  glyphImages[82] = new rsImageGray(4, 7, glyph_R);

  rsUint8 glyph_r[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,_,X,X,
    X,X,_,_,
    X,_,_,_,
    X,_,_,_,
    X,_,_,_
  };
  glyphImages[114] = new rsImageGray(4, 7, glyph_r);

  rsUint8 glyph_S[4*7] = 
  {
    _,X,X,X,
    X,_,_,_,
    X,_,_,_,
    _,X,X,_,
    _,_,_,X,
    _,_,_,X,
    X,X,X,_
  };
  glyphImages[83] = new rsImageGray(4, 7, glyph_S);

  rsUint8 glyph_s[3*7] = 
  {
    _,_,_,
    _,_,_,
    _,X,X,
    X,_,_,
    _,X,_,
    _,_,X,
    X,X,_
  };
  glyphImages[115] = new rsImageGray(3, 7, glyph_s);

  rsUint8 glyph_T[5*7] = 
  {
    X,X,X,X,X,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_,
    _,_,X,_,_
  };
  glyphImages[84] = new rsImageGray(5, 7, glyph_T);

  rsUint8 glyph_t[3*7] = 
  {
    _,X,_,
    _,X,_,
    X,X,X,
    _,X,_,
    _,X,_,
    _,X,_,
    _,X,X
  };
  glyphImages[116] = new rsImageGray(3, 7, glyph_t);

  rsUint8 glyph_U[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[85] = new rsImageGray(4, 7, glyph_U);

  rsUint8 glyph_u[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    _,X,X,_
  };
  glyphImages[117] = new rsImageGray(4, 7, glyph_u);

  rsUint8 glyph_V[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_
  };
  glyphImages[86] = new rsImageGray(5, 7, glyph_V);

  rsUint8 glyph_v[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_
  };
  glyphImages[118] = new rsImageGray(5, 7, glyph_v);

  rsUint8 glyph_W[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,X,_,X,
    X,X,_,X,X,
    X,_,_,_,X
  };
  glyphImages[87] = new rsImageGray(5, 7, glyph_W);

  rsUint8 glyph_w[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,_,_,_,X,
    X,_,_,_,X,
    X,_,X,_,X,
    X,X,_,X,X,
    X,_,_,_,X
  };
  glyphImages[119] = new rsImageGray(5, 7, glyph_w);

  rsUint8 glyph_X[5*7] = 
  {
    X,_,_,_,X,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_,
    _,X,_,X,_,
    X,_,_,_,X,
    X,_,_,_,X
  };
  glyphImages[88] = new rsImageGray(5, 7, glyph_X);

  rsUint8 glyph_x[5*7] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,_,_,_,X,
    _,X,_,X,_,
    _,_,X,_,_,
    _,X,_,X,_,
    X,_,_,_,X
  };
  glyphImages[120] = new rsImageGray(5, 7, glyph_x);

  rsUint8 glyph_Y[4*7] = 
  {
    X,_,_,X,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    _,_,_,X,
    X,X,X,X
  };
  glyphImages[89] = new rsImageGray(4, 7, glyph_Y);

  rsUint8 glyph_y[4*7] = 
  {
    _,_,_,_,
    _,_,_,_,
    X,_,_,X,
    X,_,_,X,
    X,X,X,X,
    _,_,_,X,
    X,X,X,X
  };
  glyphImages[121] = new rsImageGray(4, 7, glyph_y);


  rsUint8 glyph_Z[5*7] = 
  {
    X,X,X,X,X,
    _,_,_,_,X,
    _,_,_,X,_,
    _,_,X,_,_,
    _,X,_,_,_,
    X,_,_,_,_,
    X,X,X,X,X
  };
  glyphImages[90] = new rsImageGray(5, 7, glyph_Z);

  rsUint8 glyph_z[3*7] = 
  {
    _,_,_,
    _,_,_,
    X,X,X,
    _,_,X,
    _,X,_,
    X,_,_,
    X,X,X
  };
  glyphImages[122] = new rsImageGray(3, 7, glyph_z);


  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphImages[g] != NULL )
      glyphWidths[g] = glyphImages[g]->getWidth();
  }
}
