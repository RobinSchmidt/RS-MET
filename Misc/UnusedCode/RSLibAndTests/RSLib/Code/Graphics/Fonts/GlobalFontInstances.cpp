using namespace RSLib;

rsPixelFontRoundedA7D0*      rsGlobalFontInstances::normalFont7px = NULL;
rsPixelFontRoundedBoldA10D0* rsGlobalFontInstances::boldFont10px  = NULL;
rsPixelFontRoundedBoldA16D0* rsGlobalFontInstances::boldFont16px  = NULL; 

rsPixelFont* rsGlobalFontInstances::getPixelFontRoundedA7D0()
{
  if( normalFont7px == NULL )
    normalFont7px = new rsPixelFontRoundedA7D0;
  return normalFont7px;
}

rsPixelFont* rsGlobalFontInstances::getPixelFontRoundedBoldA10D0()
{
  if( boldFont10px == NULL )
    boldFont10px = new rsPixelFontRoundedBoldA10D0;
  return boldFont10px;
}

rsPixelFont* rsGlobalFontInstances::getPixelFontRoundedBoldA16D0()
{
  if( boldFont16px == NULL )
    boldFont16px = new rsPixelFontRoundedBoldA16D0;
  return boldFont16px;
}

void rsGlobalFontInstances::cleanUpFonts()
{
  delete normalFont7px; normalFont7px = NULL;
  delete boldFont10px;  boldFont10px  = NULL;
  delete boldFont16px;  boldFont16px  = NULL;
}

/*
rsPixelFontRoundedA7D0      *normalFont7px = NULL;
rsPixelFontRoundedBoldA9D0  *normalFont9px = NULL;
rsPixelFontRoundedBoldA10D0 *boldFont10px  = NULL;
rsPixelFontRoundedBoldA16D0 *boldFont16px  = NULL; 
*/
/*
static const rsPixelFontRoundedA7D0      normalFont7px;
static const rsPixelFontRoundedBoldA10D0 boldFont10px;
static const rsPixelFontRoundedBoldA16D0 boldFont16px;
*/