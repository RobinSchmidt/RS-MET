
//-------------------------------------------------------------------------------------------------
// construction/destruction:

//BitmapFont* BitmapFont::instance = NULL;

BitmapFont::BitmapFont() 
{
  ascent         = 0;
  descent        = 0;
  defaultKerning = 1;

  for(int g=0; g<numGlyphs; g++)
  {
    glyphBitmaps[g] = NULL;
    glyphWidths[g]  = 0;
  }

  createGlyphBitmaps();
}

BitmapFont::~BitmapFont()
{
  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphBitmaps[g] != NULL )
    {
      delete glyphBitmaps[g];
      glyphBitmaps[g] = NULL;
    }
  }
}

//-------------------------------------------------------------------------------------------------
// setup:

//-------------------------------------------------------------------------------------------------
// inquiry:

int BitmapFont::getTextPixelWidth(const juce::String &text, int kerning) const
{
  int  result = 0;
  juce_wchar c; 
  for(int i=0; i<text.length(); i++)
  {
    c       = text[i];
    result += getGlyphWidth(c) + kerning;
  }
  return result - kerning; // we added one too much in the loop
}

int BitmapFont::glyphIndexToX(const juce::String &text, int index, int kerning) const
{
  jassert( index < text.length() ); // index is out of range
  if( index >= text.length() )
    return 0;

  int   result = 0;
  juce_wchar c; 
  for(int i=0; i<index; i++)
  {
    c       = text[i];
    result += getGlyphWidth(c) + kerning;
  }
  return result; 
}

int BitmapFont::xToGlyphIndex(const String &text, int xToFindGlyphIndexAt, int kerning) const
{
  int x     = 0;
  int index = 0;
  juce_wchar c; 
  while( x < xToFindGlyphIndexAt && index < text.length() )
  {
    c  = text[index];
    x += getGlyphWidth(c) + kerning;
    index++;
  }
  return index;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void BitmapFont::createGlyphBitmaps()
{

}

//=================================================================================================

const BitmapFontRoundedA7D0 BitmapFontRoundedA7D0::instance; 

BitmapFontRoundedA7D0::BitmapFontRoundedA7D0()
{
  ascent  = 7;
  descent = 0;

  createGlyphBitmaps();
}

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

static const BitmapFontRoundedA7D0 normalFont7px;  // the global object

//=================================================================================================

const BitmapFontRoundedBoldA9D0 BitmapFontRoundedBoldA9D0::instance; 

BitmapFontRoundedBoldA9D0::BitmapFontRoundedBoldA9D0()
{
  ascent  = 9;
  descent = 0;

  createGlyphBitmaps();
}

void BitmapFontRoundedBoldA9D0::createGlyphBitmaps()
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


  unsigned char glyph_space[4*9] = 
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
  glyphBitmaps[32] = new ColourizableBitmap(4, 9, glyph_space);

  unsigned char glyph_exclamationMark[2*9] = 
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
  glyphBitmaps[33] = new ColourizableBitmap(2, 9, glyph_exclamationMark);

  unsigned char glyph_quotes[6*9] = 
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
  glyphBitmaps[34] = new ColourizableBitmap(6, 9, glyph_quotes);

  unsigned char glyph_sharp[8*9] = 
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
  glyphBitmaps[35] = new ColourizableBitmap(8, 9, glyph_sharp);

  unsigned char glyph_dollar[6*9] = 
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
  glyphBitmaps[36] = new ColourizableBitmap(6, 9, glyph_dollar);

  unsigned char glyph_percent[6*9] = 
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
  glyphBitmaps[37] = new ColourizableBitmap(6, 9, glyph_percent);

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

  unsigned char glyph_singlequote[2*9] = 
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
  glyphBitmaps[39] = new ColourizableBitmap(2, 9, glyph_singlequote);

  unsigned char glyph_openingBrace[3*9] = 
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
  glyphBitmaps[40] = new ColourizableBitmap(3, 9, glyph_openingBrace);

  unsigned char glyph_closingBrace[3*9] = 
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
  glyphBitmaps[41] = new ColourizableBitmap(3, 9, glyph_closingBrace);

  unsigned char glyph_multiply[8*9] = 
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
  glyphBitmaps[42] = new ColourizableBitmap(8, 9, glyph_multiply);

  unsigned char glyph_plus[6*9] = 
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
  glyphBitmaps[43] = new ColourizableBitmap(6, 9, glyph_plus);


  unsigned char glyph_comma[2*9] = 
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
  glyphBitmaps[44] = new ColourizableBitmap(2, 9, glyph_comma);

  unsigned char glyph_minus[6*9] = 
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
  glyphBitmaps[45] = new ColourizableBitmap(6, 9, glyph_minus);

  unsigned char glyph_dot[2*9] = 
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
  glyphBitmaps[46] = new ColourizableBitmap(2, 9, glyph_dot);

  unsigned char glyph_slash[6*9] = 
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
  glyphBitmaps[47] = new ColourizableBitmap(6, 9, glyph_slash);

  unsigned char glyph_0[6*9] = 
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
  glyphBitmaps[48] = new ColourizableBitmap(6, 9, glyph_0);

  unsigned char glyph_1[4*9] = 
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
  glyphBitmaps[49] = new ColourizableBitmap(4, 9, glyph_1);

  unsigned char glyph_2[6*9] = 
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
  glyphBitmaps[50] = new ColourizableBitmap(6, 9, glyph_2);

  unsigned char glyph_3[6*9] = 
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
  glyphBitmaps[51] = new ColourizableBitmap(6, 9, glyph_3);


  unsigned char glyph_4[6*9] = 
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
  glyphBitmaps[52] = new ColourizableBitmap(6, 9, glyph_4);

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

  unsigned char glyph_5[6*9] = 
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
  glyphBitmaps[53] = new ColourizableBitmap(6, 9, glyph_5);

  unsigned char glyph_6[6*9] = 
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
  glyphBitmaps[54] = new ColourizableBitmap(6, 9, glyph_6);

  unsigned char glyph_7[6*9] = 
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
  glyphBitmaps[55] = new ColourizableBitmap(6, 9, glyph_7);

  unsigned char glyph_8[6*9] = 
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
  glyphBitmaps[56] = new ColourizableBitmap(6, 9, glyph_8);

  unsigned char glyph_9[6*9] = 
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
  glyphBitmaps[57] = new ColourizableBitmap(6, 9, glyph_9);

  unsigned char glyph_colon[2*9] = 
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
  glyphBitmaps[58] = new ColourizableBitmap(2, 9, glyph_colon);


  unsigned char glyph_semicolon[2*9] = 
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
  glyphBitmaps[59] = new ColourizableBitmap(2, 9, glyph_semicolon);

  unsigned char glyph_openingAngleBrace[6*9] = 
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
  glyphBitmaps[60] = new ColourizableBitmap(6, 9, glyph_openingAngleBrace);

  unsigned char glyph_equalsSign[6*9] = 
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
  glyphBitmaps[61] = new ColourizableBitmap(6, 9, glyph_equalsSign);

  unsigned char glyph_closingAngleBrace[6*9] = 
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
  glyphBitmaps[62] = new ColourizableBitmap(6, 9, glyph_closingAngleBrace);

  unsigned char glyph_questionMark[6*9] = 
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
  glyphBitmaps[63] = new ColourizableBitmap(6, 9, glyph_questionMark);


  unsigned char glyph_at[10*9] = 
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
  glyphBitmaps[64] = new ColourizableBitmap(10, 9, glyph_at);

  unsigned char glyph_openingBracket[4*9] = 
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
  glyphBitmaps[91] = new ColourizableBitmap(4, 9, glyph_openingBracket);

  unsigned char glyph_backslash[6*9] = 
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
  glyphBitmaps[92] = new ColourizableBitmap(6, 9, glyph_backslash);

  unsigned char glyph_closingBracket[4*9] = 
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
  glyphBitmaps[93] = new ColourizableBitmap(4, 9, glyph_closingBracket);


  unsigned char glyph_power[7*9] = 
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
  glyphBitmaps[94] = new ColourizableBitmap(7, 9, glyph_power);

  unsigned char glyph_underscore[6*9] = 
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
  glyphBitmaps[95] = new ColourizableBitmap(6, 9, glyph_underscore);

  unsigned char glyph_96[3*9] = 
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
  glyphBitmaps[96] = new ColourizableBitmap(3, 9, glyph_96);

  unsigned char glyph_openingCurlyBrace[4*9] = 
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
  glyphBitmaps[123] = new ColourizableBitmap(4, 9, glyph_openingCurlyBrace);

  unsigned char glyph_verticalLine[2*9] = 
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
  glyphBitmaps[124] = new ColourizableBitmap(2, 9, glyph_verticalLine);

  unsigned char glyph_closingCurlyBrace[4*9] = 
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
  glyphBitmaps[125] = new ColourizableBitmap(4, 9, glyph_closingCurlyBrace);

  unsigned char glyph_tilde[7*9] = 
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
  glyphBitmaps[126] = new ColourizableBitmap(7, 9, glyph_tilde);

  unsigned char glyph_A[6*9] = 
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
  glyphBitmaps[65] = new ColourizableBitmap(6, 9, glyph_A);

  unsigned char glyph_a[6*9] = 
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
  glyphBitmaps[97] = new ColourizableBitmap(6, 9, glyph_a);

  unsigned char glyph_B[6*10] = 
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
  glyphBitmaps[66] = new ColourizableBitmap(6, 9, glyph_B);

  unsigned char glyph_b[6*9] = 
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
  glyphBitmaps[98] = new ColourizableBitmap(6, 9, glyph_b);

  unsigned char glyph_C[5*9] = 
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
  glyphBitmaps[67] = new ColourizableBitmap(5, 9, glyph_C);

  unsigned char glyph_c[5*9] = 
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
  glyphBitmaps[99] = new ColourizableBitmap(5, 9, glyph_c);

  unsigned char glyph_D[6*9] = 
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
  glyphBitmaps[68] = new ColourizableBitmap(6, 9, glyph_D);

  unsigned char glyph_d[6*9] = 
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
  glyphBitmaps[100] = new ColourizableBitmap(6, 9, glyph_d);

  unsigned char glyph_E[6*9] = 
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
  glyphBitmaps[69] = new ColourizableBitmap(6, 9, glyph_E);

  unsigned char glyph_e[6*9] = 
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
  glyphBitmaps[101] = new ColourizableBitmap(6, 9, glyph_e);

  unsigned char glyph_F[6*9] = 
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
  glyphBitmaps[70] = new ColourizableBitmap(6, 9, glyph_F);

  unsigned char glyph_f[5*9] = 
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
  glyphBitmaps[102] = new ColourizableBitmap(5, 9, glyph_f);

  unsigned char glyph_G[6*9] = 
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
  glyphBitmaps[71] = new ColourizableBitmap(6, 9, glyph_G);

  unsigned char glyph_g[6*9] = 
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
  glyphBitmaps[103] = new ColourizableBitmap(6, 9, glyph_g);

  unsigned char glyph_H[6*9] = 
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
  glyphBitmaps[72] = new ColourizableBitmap(6, 9, glyph_H);

  unsigned char glyph_h[6*9] = 
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
  glyphBitmaps[104] = new ColourizableBitmap(6, 9, glyph_h);

  unsigned char glyph_I[2*9] = 
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
  glyphBitmaps[73] = new ColourizableBitmap(2, 9, glyph_I);

  unsigned char glyph_i[2*9] = 
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
  glyphBitmaps[105] = new ColourizableBitmap(2, 9, glyph_i);

  unsigned char glyph_J[6*9] = 
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
  glyphBitmaps[74] = new ColourizableBitmap(6, 9, glyph_J);

  unsigned char glyph_j[4*9] = 
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
  glyphBitmaps[106] = new ColourizableBitmap(4, 9, glyph_j);

  unsigned char glyph_K[6*9] = 
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
  glyphBitmaps[75] = new ColourizableBitmap(6, 9, glyph_K);

  unsigned char glyph_k[6*9] = 
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
  glyphBitmaps[107] = new ColourizableBitmap(6, 9, glyph_k);

  unsigned char glyph_L[6*9] = 
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
  glyphBitmaps[76] = new ColourizableBitmap(6, 9, glyph_L);

  unsigned char glyph_l[3*9] = 
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
  glyphBitmaps[108] = new ColourizableBitmap(3, 9, glyph_l);

  unsigned char glyph_M[9*9] = 
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
  glyphBitmaps[77] = new ColourizableBitmap(9, 9, glyph_M);

  unsigned char glyph_m[8*9] = 
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
  glyphBitmaps[109] = new ColourizableBitmap(8, 9, glyph_m);

  unsigned char glyph_N[7*9] = 
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
  glyphBitmaps[78] = new ColourizableBitmap(7, 9, glyph_N);

  unsigned char glyph_n[6*9] = 
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
  glyphBitmaps[110] = new ColourizableBitmap(6, 9, glyph_n);

  unsigned char glyph_O[6*9] = 
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
  glyphBitmaps[79] = new ColourizableBitmap(6, 9, glyph_O);

  unsigned char glyph_o[6*9] = 
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
  glyphBitmaps[111] = new ColourizableBitmap(6, 9, glyph_o);

  unsigned char glyph_P[6*9] = 
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
  glyphBitmaps[80] = new ColourizableBitmap(6, 9, glyph_P);

  unsigned char glyph_p[6*9] = 
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
  glyphBitmaps[112] = new ColourizableBitmap(6, 9, glyph_p);

  unsigned char glyph_Q[7*9] = 
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
  glyphBitmaps[81] = new ColourizableBitmap(7, 9, glyph_Q);

  unsigned char glyph_q[6*9] = 
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
  glyphBitmaps[113] = new ColourizableBitmap(6, 9, glyph_q);

  unsigned char glyph_R[6*9] = 
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
  glyphBitmaps[82] = new ColourizableBitmap(6, 9, glyph_R);

  unsigned char glyph_r[5*9] = 
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
  glyphBitmaps[114] = new ColourizableBitmap(5, 9, glyph_r);

  unsigned char glyph_S[6*9] = 
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
  glyphBitmaps[83] = new ColourizableBitmap(6, 9, glyph_S);

  unsigned char glyph_s[5*9] = 
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
  glyphBitmaps[115] = new ColourizableBitmap(5, 9, glyph_s);

  unsigned char glyph_T[6*9] = 
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
  glyphBitmaps[84] = new ColourizableBitmap(6, 9, glyph_T);

  unsigned char glyph_t[4*9] = 
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
  glyphBitmaps[116] = new ColourizableBitmap(4, 9, glyph_t);

  unsigned char glyph_U[6*9] = 
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
  glyphBitmaps[85] = new ColourizableBitmap(6, 9, glyph_U);

  unsigned char glyph_u[6*9] = 
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
  glyphBitmaps[117] = new ColourizableBitmap(6, 9, glyph_u);

  unsigned char glyph_V[7*9] = 
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
  glyphBitmaps[86] = new ColourizableBitmap(7, 9, glyph_V);

  unsigned char glyph_v[7*9] = 
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
  glyphBitmaps[118] = new ColourizableBitmap(7, 9, glyph_v);

  unsigned char glyph_W[9*9] = 
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
  glyphBitmaps[87] = new ColourizableBitmap(9, 9, glyph_W);

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
  unsigned char glyph_w[7*9] = 
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
  glyphBitmaps[119] = new ColourizableBitmap(7, 9, glyph_w);

  unsigned char glyph_X[7*9] = 
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
  glyphBitmaps[88] = new ColourizableBitmap(7, 9, glyph_X);

  unsigned char glyph_x[7*9] = 
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
  glyphBitmaps[120] = new ColourizableBitmap(7, 9, glyph_x);

  unsigned char glyph_Y[6*9] = 
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
  glyphBitmaps[89] = new ColourizableBitmap(6, 9, glyph_Y);

  unsigned char glyph_y[6*9] = 
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
  glyphBitmaps[121] = new ColourizableBitmap(6, 9, glyph_y);

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

  unsigned char glyph_Z[7*9] = 
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
  glyphBitmaps[90] = new ColourizableBitmap(7, 9, glyph_Z);

  unsigned char glyph_z[5*9] = 
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
  glyphBitmaps[122] = new ColourizableBitmap(5, 9, glyph_z);


  for(int g=0; g<numGlyphs; g++)
  {
    if( glyphBitmaps[g] != NULL )
      glyphWidths[g] = glyphBitmaps[g]->getWidth();
  }
}

//static const BitmapFontRoundedA9D0 normalFont9px;  // the global object

//=================================================================================================

const BitmapFontRoundedBoldA10D0 BitmapFontRoundedBoldA10D0::instance; 
// creates the static global instance

BitmapFontRoundedBoldA10D0::BitmapFontRoundedBoldA10D0()
{
  ascent  = 10;
  descent = 0;
  createGlyphBitmaps();
}

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

//static const BitmapFontRoundedBoldA10D0 boldFont10px; // the global object

//=================================================================================================

const BitmapFontRoundedBoldA16D0 BitmapFontRoundedBoldA16D0::instance; 

BitmapFontRoundedBoldA16D0::BitmapFontRoundedBoldA16D0()
{
  ascent  = 16;
  descent = 0;
  createGlyphBitmaps();
}

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


  // Here, the lowercase m is wider than the uppercase M. In the smaller font, it's the other 
  // way around. That's inconsistent! Maybe narrow the lowercase m to 11x16

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



  // These Ys are inconsistent with those of the smaller font:

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



  // Needs to be bigger:

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
    X,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,_,_,_
  };
  glyphBitmaps[53] = new ColourizableBitmap(9, 16, glyph_5);
  /*
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
    X,X,X,X,X,X,X,_,_,   // This line is different 
    X,X,X,X,X,X,_,_,_
  };
  glyphBitmaps[53] = new ColourizableBitmap(9, 16, glyph_5);
  */

  unsigned char glyph_6[9*16] = 
  {
    _,_,X,X,X,X,X,X,_,
    _,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
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
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_
  };
  glyphBitmaps[54] = new ColourizableBitmap(9, 16, glyph_6);


  unsigned char glyph_7[9*16] = 
  {
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,_,_,_,_,_,X,X,X,
    _,_,_,_,_,X,X,X,_,
    _,_,_,_,_,X,X,X,_,
    _,_,_,_,X,X,X,_,_,
    _,_,_,_,X,X,X,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,_,X,X,X,_,_,_,
    _,_,X,X,X,_,_,_,_,
    _,_,X,X,X,_,_,_,_,
    _,X,X,X,_,_,_,_,_,
    _,X,X,X,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
    X,X,X,_,_,_,_,_,_,
  };
  glyphBitmaps[55] = new ColourizableBitmap(9, 16, glyph_7);


  unsigned char glyph_8[9*16] = 
  {
    _,_,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    X,X,X,X,X,X,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_,
    _,X,X,X,X,X,X,X,_,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,_,_,_,X,X,X,
    X,X,X,X,X,X,X,X,X,
    _,X,X,X,X,X,X,X,_,
    _,_,X,X,X,X,X,_,_
  };
  glyphBitmaps[56] = new ColourizableBitmap(9, 16, glyph_8);




  // These are still under construction:


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

static const BitmapFontRoundedBoldA16D0 boldFont16px; // the global object