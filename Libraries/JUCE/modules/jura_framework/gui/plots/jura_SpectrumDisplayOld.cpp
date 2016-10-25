
SpectrumDisplayOld::SpectrumDisplayOld(const String& name) 
: CurveFamilyPlotOld(name)
{
  // initialize the appearance settings:
  setAutoReRendering(false);

  setHorizontalCoarseGrid(6.0,           true);
  setHorizontalFineGrid(  1.0,           false);
  setVerticalCoarseGrid(  2.0,           true);
  setVerticalFineGrid(pow(2.0, 1.0/3.0), false);
  setMaximumRange(15.625, 32000.0, -96.0, 6.0);
  setCurrentRange(15.625, 32000.0, -96.0, 6.0);
  setAxisPositionX(CoordinateSystemOld::BOTTOM);
  setAxisLabelX(String::empty);
  useLogarithmicScaleX(true, 2.0);
  setAxisPositionY(CoordinateSystemOld::LEFT);
  setAxisValuesPositionY(LEFT);
  setAxisLabelY(String::empty);

  //setStringConversionForInfoLineX(hertzToStringWithUnitTotal5);  
  //setStringConversionForInfoLineX(rojue::frequencyToNoteString);  
  setStringConversionForInfoLineX(frequencyInHzAndAsNote); 
  setStringConversionForInfoLineY(decibelsToStringWithUnit2);

  setAutoReRendering(true);
}

SpectrumDisplayOld::~SpectrumDisplayOld(void)
{
  deleteAllChildren();
}

void SpectrumDisplayOld::useLogarithmicScaleX(bool shouldBeLogScaledX, double newLogBase)
{
  setAutoReRendering(false);

  if( shouldBeLogScaledX == true )
  {
    setMaximumRange(15.625, 32000.0, -96.0, 6.0);
    setCurrentRange(15.625, 32000.0, -96.0, 6.0);
    setVerticalCoarseGrid(2.0, true);
    setVerticalFineGrid(pow(2.0, 1.0/3.0), false);
    CurveFamilyPlotOld::useLogarithmicScaleX(true);
    setAxisPositionY(CoordinateSystemOld::LEFT);
    setStringConversionForAxisX(&valueToString0);
  }
  else
  {
    setMaximumRange(0.0, 21000.0, -96.0, 6.0);
    setCurrentRange(0.0, 21000.0, -96.0, 6.0);
    setVerticalCoarseGrid(1000.0, true);
    setVerticalFineGrid(100.0, false);
    CurveFamilyPlotOld::useLogarithmicScaleX(false);
    setAxisPositionY(CoordinateSystemOld::ZERO);
    setStringConversionForAxisX(&kiloToString0);
  }

  updateBackgroundImage();
  setAutoReRendering(true);
}

void SpectrumDisplayOld::setSpectrum(int newNumBins, double *newBinFrequencies, 
  double *newBinMagnitudes)
{
  CurveFamilyPlotOld::setCurveValues(newNumBins, newBinFrequencies, newBinMagnitudes);
}

void SpectrumDisplayOld::setSpectra(int newNumBins, int newNumSpectra, double* newBinFrequencies, 
  double** newBinMagnitudes)
{
  CurveFamilyPlotOld::setFunctionFamilyValues(newNumBins, newNumSpectra, newBinFrequencies, 
    newBinMagnitudes);
}

/*
void SpectrumDisplayOld::mouseMove(const MouseEvent &e)
{
  setDescription( getCoordinateStringAtPixelPosition(e.x, e.y) );
}
*/

void SpectrumDisplayOld::getDisplayedFrequencies(double *frequencies, int numBins)
{
  double x, y;
  for(int i=0; i<numBins; i++)
  {
    x = (double) i * getWidth() / (double) numBins;
    y = 1.0; // only a dummy
    transformFromComponentsCoordinates(x, y);

    frequencies[i] = x;
  }
}

/*
String SpectrumDisplayOld::getCoordinateStringAtPixelPosition(int x, int y)
{
  double frequency = (double) x;
  double level     = (double) y;
  transformFromComponentsCoordinates(frequency, level);
  String fString = hertzToStringWithUnitTotal5(frequency);
  String aString = decibelsToStringWithUnit2(level);
  return fString + String(T(", ")) + aString;
}
*/