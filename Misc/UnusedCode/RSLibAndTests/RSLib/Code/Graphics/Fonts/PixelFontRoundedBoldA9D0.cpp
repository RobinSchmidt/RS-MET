using namespace RSLib;

rsPixelFontRoundedBoldA9D0::rsPixelFontRoundedBoldA9D0()
{
  ascent  = 9;
  descent = 0;

  createGlyphImages();
}

// create the glyphs:

void rsPixelFontRoundedBoldA9D0::createGlyphImages()
{
  rsUint8 _ = 0;
  rsUint8 X = 255;

  /*
  rsUint8 glyph_dummy[6*10] = 
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
  glyphImages[60] = new rsImageGray(6, 10, glyph_glyph_dummy);
  */


  rsUint8 glyph_space[4*9] = 
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
  };
  glyphImages[32] = new rsImageGray(4, 9, glyph_space);

  rsUint8 glyph_exclamationMark[2*9] = 
  {
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
  glyphImages[33] = new rsImageGray(2, 9, glyph_exclamationMark);

  rsUint8 glyph_quotes[6*9] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,_,_,_,_,_
  };
  glyphImages[34] = new rsImageGray(6, 9, glyph_quotes);

  rsUint8 glyph_sharp[8*9] = 
  {
    _,X,X,_,_,X,X,_,
    _,X,X,_,_,X,X,_,
    X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,
    _,X,X,_,_,X,X,_,
    X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,
    _,X,X,_,_,X,X,_,
    _,X,X,_,_,X,X,_,
  };
  glyphImages[35] = new rsImageGray(8, 9, glyph_sharp);

  rsUint8 glyph_dollar[6*9] = 
  {
    _,_,X,X,_,_,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    _,_,X,X,_,_
  };
  glyphImages[36] = new rsImageGray(6, 9, glyph_dollar);

  rsUint8 glyph_percent[6*9] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,X,X,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphImages[37] = new rsImageGray(6, 9, glyph_percent);

  rsUint8 glyph_ampersand[8*10] = 
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
  glyphImages[38] = new rsImageGray(8, 10, glyph_ampersand);

  rsUint8 glyph_singlequote[2*9] = 
  {
    X,X,
    X,X,
    X,X,
    X,X,
    _,_,
    _,_,
    _,_,
    _,_,
    _,_
  };
  glyphImages[39] = new rsImageGray(2, 9, glyph_singlequote);

  rsUint8 glyph_openingBrace[3*9] = 
  {
    _,_,X,
    _,X,X,
    _,X,X,
    X,X,_,
    X,X,_,
    X,X,_,
    _,X,X,
    _,X,X,
    _,_,X
  };
  glyphImages[40] = new rsImageGray(3, 9, glyph_openingBrace);

  rsUint8 glyph_closingBrace[3*9] = 
  {
    X,_,_,
    X,X,_,
    X,X,_,
    _,X,X,
    _,X,X,
    _,X,X,
    X,X,_,
    X,X,_,
    X,_,_
  };
  glyphImages[41] = new rsImageGray(3, 9, glyph_closingBrace);

  rsUint8 glyph_multiply[8*9] = 
  {
    _,_,_,_,_,_,_,_,
    _,_,_,X,X,_,_,_,
    _,_,_,X,X,_,_,_,
    X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,
    _,_,X,X,X,X,_,_,
    _,X,X,X,X,X,X,_,
    _,X,X,_,_,X,X,_,
    _,_,_,_,_,_,_,_
  };
  glyphImages[42] = new rsImageGray(8, 9, glyph_multiply);

  rsUint8 glyph_plus[6*9] = 
  {
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
  glyphImages[43] = new rsImageGray(6, 9, glyph_plus);


  rsUint8 glyph_comma[2*9] = 
  {
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
  glyphImages[44] = new rsImageGray(2, 9, glyph_comma);

  rsUint8 glyph_minus[6*9] = 
  {
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
  glyphImages[45] = new rsImageGray(6, 9, glyph_minus);

  rsUint8 glyph_dot[2*9] = 
  {
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
  glyphImages[46] = new rsImageGray(2, 9, glyph_dot);

  rsUint8 glyph_slash[6*9] = 
  {
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,X,X,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphImages[47] = new rsImageGray(6, 9, glyph_slash);

  rsUint8 glyph_0[6*9] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphImages[48] = new rsImageGray(6, 9, glyph_0);

  rsUint8 glyph_1[4*9] = 
  {
    _,_,X,X,
    _,X,X,X,
    X,X,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X
  };
  glyphImages[49] = new rsImageGray(4, 9, glyph_1);

  rsUint8 glyph_2[6*9] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,X,X,X,X,
    _,X,X,X,X,_,
    X,X,X,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X
  };
  glyphImages[50] = new rsImageGray(6, 9, glyph_2);

  rsUint8 glyph_3[6*9] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,X,X,X,X,X,
    _,X,X,X,X,_,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphImages[51] = new rsImageGray(6, 9, glyph_3);


  rsUint8 glyph_4[6*9] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X
  };
  glyphImages[52] = new rsImageGray(6, 9, glyph_4);

  /*
  rsUint8 glyph_5[6*10] = 
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
  glyphImages[53] = new rsImageGray(6, 10, glyph_5);
  */

  rsUint8 glyph_5[6*9] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphImages[53] = new rsImageGray(6, 9, glyph_5);

  rsUint8 glyph_6[6*9] = 
  {
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphImages[54] = new rsImageGray(6, 9, glyph_6);

  rsUint8 glyph_7[6*9] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_
  };
  glyphImages[55] = new rsImageGray(6, 9, glyph_7);

  rsUint8 glyph_8[6*9] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    _,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphImages[56] = new rsImageGray(6, 9, glyph_8);

  rsUint8 glyph_9[6*9] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphImages[57] = new rsImageGray(6, 9, glyph_9);

  rsUint8 glyph_colon[2*9] = 
  {
    _,_,
    _,_,
    X,X,
    X,X,
    _,_,
    _,_,
    X,X,
    X,X,
    _,_
  };
  glyphImages[58] = new rsImageGray(2, 9, glyph_colon);


  rsUint8 glyph_semicolon[2*9] = 
  {
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
  glyphImages[59] = new rsImageGray(2, 9, glyph_semicolon);

  rsUint8 glyph_openingAngleBrace[6*9] = 
  {
    _,_,_,_,X,X,
    _,_,_,X,X,X,
    _,_,X,X,X,_,
    _,X,X,X,_,_,
    X,X,X,_,_,_,
    _,X,X,X,_,_,
    _,_,X,X,X,_,
    _,_,_,X,X,X,
    _,_,_,_,X,X
  };
  glyphImages[60] = new rsImageGray(6, 9, glyph_openingAngleBrace);

  rsUint8 glyph_equalsSign[6*9] = 
  {
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
  glyphImages[61] = new rsImageGray(6, 9, glyph_equalsSign);

  rsUint8 glyph_closingAngleBrace[6*9] = 
  {
    X,X,_,_,_,_,
    X,X,X,_,_,_,
    _,X,X,X,_,_,
    _,_,X,X,X,_,
    _,_,_,X,X,X,
    _,_,X,X,X,_,
    _,X,X,X,_,_,
    X,X,X,_,_,_,
    X,X,_,_,_,_
  };
  glyphImages[62] = new rsImageGray(6, 9, glyph_closingAngleBrace);

  rsUint8 glyph_questionMark[6*9] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphImages[63] = new rsImageGray(6, 9, glyph_questionMark);


  rsUint8 glyph_at[10*9] = 
  {
    _,X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,X,
    X,X,_,_,_,_,_,_,X,X,
    X,X,_,X,X,X,X,_,X,X,
    X,X,_,X,_,_,X,_,X,X,
    X,X,_,X,X,X,X,X,X,_,
    X,X,_,_,_,_,_,_,_,_,
    X,X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,X,X
  };
  glyphImages[64] = new rsImageGray(10, 9, glyph_at);

  rsUint8 glyph_openingBracket[4*9] = 
  {
    X,X,X,X,
    X,X,X,X,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,_,_,
    X,X,X,X,
    X,X,X,X
  };
  glyphImages[91] = new rsImageGray(4, 9, glyph_openingBracket);

  rsUint8 glyph_backslash[6*9] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    _,X,X,_,_,_,
    _,X,X,_,_,_,
    _,_,X,X,_,_,
    _,_,_,X,X,_,
    _,_,_,X,X,_,
    _,_,_,_,X,X,
    _,_,_,_,X,X
  };
  glyphImages[92] = new rsImageGray(6, 9, glyph_backslash);

  rsUint8 glyph_closingBracket[4*9] = 
  {
    X,X,X,X,
    X,X,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    X,X,X,X,
    X,X,X,X
  };
  glyphImages[93] = new rsImageGray(4, 9, glyph_closingBracket);


  rsUint8 glyph_power[7*9] = 
  {
    _,_,_,X,_,_,_,
    _,_,X,X,X,_,_,
    _,X,X,X,X,X,_,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_
  };
  glyphImages[94] = new rsImageGray(7, 9, glyph_power);

  rsUint8 glyph_underscore[6*9] = 
  {
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
  glyphImages[95] = new rsImageGray(6, 9, glyph_underscore);

  rsUint8 glyph_96[3*9] = 
  {
    X,X,_,
    X,X,X,
    X,X,X,
    _,X,X,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_,
    _,_,_
  };
  glyphImages[96] = new rsImageGray(3, 9, glyph_96);

  rsUint8 glyph_openingCurlyBrace[4*9] = 
  {
    _,X,X,X,
    _,X,X,X,
    _,X,X,_,
    X,X,X,_,
    X,X,X,_,
    _,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,X,X,X
  };
  glyphImages[123] = new rsImageGray(4, 9, glyph_openingCurlyBrace);

  rsUint8 glyph_verticalLine[2*9] = 
  {
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
  glyphImages[124] = new rsImageGray(2, 9, glyph_verticalLine);

  rsUint8 glyph_closingCurlyBrace[4*9] = 
  {
    X,X,X,_,
    X,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,X,X,X,
    _,X,X,_,
    _,X,X,_,
    X,X,X,_,
    X,X,X,_
  };
  glyphImages[125] = new rsImageGray(4, 9, glyph_closingCurlyBrace);

  rsUint8 glyph_tilde[7*9] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,X,X,X,_,X,X,
    X,X,X,X,X,X,X,
    X,X,_,X,X,X,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_
  };
  glyphImages[126] = new rsImageGray(7, 9, glyph_tilde);

  rsUint8 glyph_A[6*9] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphImages[65] = new rsImageGray(6, 9, glyph_A);

  rsUint8 glyph_a[6*9] = 
  {
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
  glyphImages[97] = new rsImageGray(6, 9, glyph_a);

  rsUint8 glyph_B[6*10] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphImages[66] = new rsImageGray(6, 9, glyph_B);

  rsUint8 glyph_b[6*9] = 
  {
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphImages[98] = new rsImageGray(6, 9, glyph_b);

  rsUint8 glyph_C[5*9] = 
  {
    _,X,X,X,X,
    X,X,X,X,X,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,X,X,X,
    _,X,X,X,X
  };
  glyphImages[67] = new rsImageGray(5, 9, glyph_C);

  rsUint8 glyph_c[5*9] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    _,X,X,X,X,
    X,X,X,X,X,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,X,X,X,
    _,X,X,X,X
  };
  glyphImages[99] = new rsImageGray(5, 9, glyph_c);

  rsUint8 glyph_D[6*9] = 
  {
    X,X,X,X,_,_,
    X,X,X,X,X,_,
    X,X,_,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_
  };
  glyphImages[68] = new rsImageGray(6, 9, glyph_D);

  rsUint8 glyph_d[6*9] = 
  {
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X
  };
  glyphImages[100] = new rsImageGray(6, 9, glyph_d);

  rsUint8 glyph_E[6*9] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,X,X,X,X,
    X,X,X,X,X,X
  };
  glyphImages[69] = new rsImageGray(6, 9, glyph_E);

  rsUint8 glyph_e[6*9] = 
  {
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
  glyphImages[101] = new rsImageGray(6, 9, glyph_e);

  rsUint8 glyph_F[6*9] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
  };
  glyphImages[70] = new rsImageGray(6, 9, glyph_F);

  rsUint8 glyph_f[5*9] = 
  {
    _,_,X,X,X,
    _,X,X,X,X,
    _,X,X,_,_,
    _,X,X,_,_,
    X,X,X,X,_,  
    X,X,X,X,_,
    _,X,X,_,_,
    _,X,X,_,_,
    _,X,X,_,_
  };
  glyphImages[102] = new rsImageGray(5, 9, glyph_f);

  rsUint8 glyph_G[6*9] = 
  {
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,_,X,X,X,
    X,X,_,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_,
  };
  glyphImages[71] = new rsImageGray(6, 9, glyph_G);

  rsUint8 glyph_g[6*9] = 
  {
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
  glyphImages[103] = new rsImageGray(6, 9, glyph_g);

  rsUint8 glyph_H[6*9] = 
  {
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
  glyphImages[72] = new rsImageGray(6, 9, glyph_H);

  rsUint8 glyph_h[6*9] = 
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
  };
  glyphImages[104] = new rsImageGray(6, 9, glyph_h);

  rsUint8 glyph_I[2*9] = 
  {
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
  glyphImages[73] = new rsImageGray(2, 9, glyph_I);

  rsUint8 glyph_i[2*9] = 
  {
    X,X,
    X,X,
    _,_,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X,
    X,X
  };
  glyphImages[105] = new rsImageGray(2, 9, glyph_i);

  rsUint8 glyph_J[6*9] = 
  {
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
  glyphImages[74] = new rsImageGray(6, 9, glyph_J);

  rsUint8 glyph_j[4*9] = 
  {
    _,_,X,X,
    _,_,X,X,
    _,_,_,_,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    _,_,X,X,
    X,X,X,X,
    X,X,X,_
  };
  glyphImages[106] = new rsImageGray(4, 9, glyph_j);

  rsUint8 glyph_K[6*9] = 
  {
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_,
    X,X,X,X,X,_,
    X,X,_,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphImages[75] = new rsImageGray(6, 9, glyph_K);

  rsUint8 glyph_k[6*9] = 
  {
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
  glyphImages[107] = new rsImageGray(6, 9, glyph_k);

  rsUint8 glyph_L[6*9] = 
  {
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
  glyphImages[76] = new rsImageGray(6, 9, glyph_L);

  rsUint8 glyph_l[3*9] = 
  {
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
  glyphImages[108] = new rsImageGray(3, 9, glyph_l);

  rsUint8 glyph_M[9*9] = 
  {
    X,X,_,_,_,_,_,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,_,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,_,X,X,X,_,X,X,
    X,X,_,_,X,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X,
    X,X,_,_,_,_,_,X,X
  };
  glyphImages[77] = new rsImageGray(9, 9, glyph_M);

  rsUint8 glyph_m[8*9] = 
  {
    _,_,_,_,_,_,_,_,
    _,_,_,_,_,_,_,_,
    X,X,X,_,X,X,X,_,
    X,X,X,X,X,X,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X,
    X,X,_,X,X,_,X,X
  };
  glyphImages[109] = new rsImageGray(8, 9, glyph_m);

  rsUint8 glyph_N[7*9] = 
  {
    X,X,_,_,_,X,X,
    X,X,X,_,_,X,X,
    X,X,X,X,_,X,X,
    X,X,X,X,X,X,X,
    X,X,_,X,X,X,X,
    X,X,_,_,X,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X
  };
  glyphImages[78] = new rsImageGray(7, 9, glyph_N);

  rsUint8 glyph_n[6*9] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X
  };
  glyphImages[110] = new rsImageGray(6, 9, glyph_n);

  rsUint8 glyph_O[6*9] = 
  {
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphImages[79] = new rsImageGray(6, 9, glyph_O);

  rsUint8 glyph_o[6*9] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphImages[111] = new rsImageGray(6, 9, glyph_o);

  rsUint8 glyph_P[6*9] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphImages[80] = new rsImageGray(6, 9, glyph_P);

  rsUint8 glyph_p[6*9] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,_,_,_,_,
    X,X,_,_,_,_
  };
  glyphImages[112] = new rsImageGray(6, 9, glyph_p);

  rsUint8 glyph_Q[7*9] = 
  {
    _,X,X,X,X,_,_,
    X,X,X,X,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,_,X,X,_,
    X,X,_,X,X,X,_,
    X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,
  };
  glyphImages[81] = new rsImageGray(7, 9, glyph_Q);

  rsUint8 glyph_q[6*9] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X
  };
  glyphImages[113] = new rsImageGray(6, 9, glyph_q);

  rsUint8 glyph_R[6*9] = 
  {
    X,X,X,X,X,_,
    X,X,X,X,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_,
    X,X,X,X,_,_,
    X,X,X,X,X,_,
    X,X,_,X,X,X,
    X,X,_,_,X,X
  };
  glyphImages[82] = new rsImageGray(6, 9, glyph_R);

  rsUint8 glyph_r[5*9] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,_,X,X,
    X,X,X,X,X,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_,
    X,X,_,_,_
  };
  glyphImages[114] = new rsImageGray(5, 9, glyph_r);

  rsUint8 glyph_S[6*9] = 
  {
    _,X,X,X,X,X,
    X,X,X,X,X,X,
    X,X,_,_,_,_,
    X,X,X,X,X,_,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphImages[83] = new rsImageGray(6, 9, glyph_S);

  rsUint8 glyph_s[5*9] = 
  {
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
  glyphImages[115] = new rsImageGray(5, 9, glyph_s);

  rsUint8 glyph_T[6*9] = 
  {
    X,X,X,X,X,X,
    X,X,X,X,X,X,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_,
    _,_,X,X,_,_
  };
  glyphImages[84] = new rsImageGray(6, 9, glyph_T);

  rsUint8 glyph_t[4*9] = 
  {
    _,X,X,_,
    _,X,X,_,
    X,X,X,X,
    X,X,X,X,
    _,X,X,_,
    _,X,X,_,
    _,X,X,_,
    _,X,X,X,
    _,_,X,X
  };
  glyphImages[116] = new rsImageGray(4, 9, glyph_t);

  rsUint8 glyph_U[6*9] = 
  {
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
  glyphImages[85] = new rsImageGray(6, 9, glyph_U);

  rsUint8 glyph_u[6*9] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,_
  };
  glyphImages[117] = new rsImageGray(6, 9, glyph_u);

  rsUint8 glyph_V[7*9] = 
  {
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
  glyphImages[86] = new rsImageGray(7, 9, glyph_V);

  rsUint8 glyph_v[7*9] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,X,_,X,X,X,
    _,X,X,X,X,X,_,
    _,_,X,X,X,_,_,
    _,_,_,X,_,_,_
  };
  glyphImages[118] = new rsImageGray(7, 9, glyph_v);

  rsUint8 glyph_W[9*9] = 
  {
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
  glyphImages[87] = new rsImageGray(9, 9, glyph_W);

  /*
  rsUint8 glyph_w[9*10] = 
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
  glyphImages[119] = new rsImageGray(9, 10, glyph_w);
  */
  rsUint8 glyph_w[7*9] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,_,X,_,X,X,
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X
  };
  glyphImages[119] = new rsImageGray(7, 9, glyph_w);

  rsUint8 glyph_X[7*9] = 
  {
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X,
    X,X,X,_,X,X,X,
    _,X,X,X,X,X,_,
    _,_,X,X,X,_,_,
    _,X,X,X,X,X,_,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X,
    X,X,_,_,_,X,X
  };
  glyphImages[88] = new rsImageGray(7, 9, glyph_X);

  rsUint8 glyph_x[7*9] = 
  {
    _,_,_,_,_,_,_,
    _,_,_,_,_,_,_,
    X,X,_,_,_,X,X,
    X,X,X,_,X,X,X,
    _,X,X,X,X,X,_,
    _,_,X,X,X,_,_,
    _,X,X,X,X,X,_,
    X,X,X,_,X,X,X,
    X,X,_,_,_,X,X
  };
  glyphImages[120] = new rsImageGray(7, 9, glyph_x);

  rsUint8 glyph_Y[6*9] = 
  {
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
  glyphImages[89] = new rsImageGray(6, 9, glyph_Y);

  rsUint8 glyph_y[6*9] = 
  {
    _,_,_,_,_,_,
    _,_,_,_,_,_,
    X,X,_,_,X,X,
    X,X,_,_,X,X,
    X,X,X,X,X,X,
    _,X,X,X,X,X,
    _,_,_,_,X,X,
    X,X,X,X,X,X,
    X,X,X,X,X,_
  };
  glyphImages[121] = new rsImageGray(6, 9, glyph_y);

  /*
  rsUint8 glyph_Y[6*10] = 
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
  glyphImages[89] = new rsImageGray(6, 10, glyph_Y);

  rsUint8 glyph_y[6*10] = 
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
  glyphImages[121] = new rsImageGray(6, 10, glyph_y);
  */

  rsUint8 glyph_Z[7*9] = 
  {
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
    _,_,_,_,X,X,X,
    _,_,_,X,X,X,_,
    _,_,X,X,X,_,_,
    _,X,X,X,_,_,_,
    X,X,X,_,_,_,_,
    X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,
  };
  glyphImages[90] = new rsImageGray(7, 9, glyph_Z);

  rsUint8 glyph_z[5*9] = 
  {
    _,_,_,_,_,
    _,_,_,_,_,
    X,X,X,X,X,
    X,X,X,X,X,
    _,_,X,X,X,
    _,X,X,X,_,
    X,X,X,_,_,
    X,X,X,X,X,
    X,X,X,X,X
  };
  glyphImages[122] = new rsImageGray(5, 9, glyph_z);


  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphImages[g] != NULL )
      glyphWidths[g] = glyphImages[g]->getWidth();
  }
}
