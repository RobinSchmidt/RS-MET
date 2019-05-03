
rsSpectrumPlot::rsSpectrumPlot(const String& name) 
: rsDataPlot(name)
{
  // initialize the appearance settings:
  setAutoReRendering(false);

  setHorizontalCoarseGrid(6.0,           true);
  setHorizontalFineGrid(  1.0,           false);
  setVerticalCoarseGrid(  2.0,           true);
  setVerticalFineGrid(pow(2.0, 1.0/3.0), false);
  setMaximumRange(15.625, 32000.0, -96.0, 6.0);
  setCurrentRange(15.625, 32000.0, -96.0, 6.0);
  setAxisPositionX(rsPlotSettings::BOTTOM);
  setAxisLabelX(String::empty);
  useLogarithmicScaleX(true, 2.0);
  setAxisPositionY(rsPlotSettings::LEFT);
  setAxisValuesPositionY(rsPlotSettings::LEFT);
  setAxisLabelY(String::empty);

  //setStringConversionForInfoLineX(hertzToStringWithUnitTotal5);  
  //setStringConversionForInfoLineX(rojue::frequencyToNoteString);  
  setStringConversionForInfoLineX(frequencyInHzAndAsNote); 
  setStringConversionForInfoLineY(decibelsToStringWithUnit2);

  setAutoReRendering(true);
}

rsSpectrumPlot::~rsSpectrumPlot(void)
{
  deleteAllChildren();
}

void rsSpectrumPlot::useLogarithmicScaleX(bool shouldBeLogScaledX, double newLogBase)
{
  setAutoReRendering(false);

  if( shouldBeLogScaledX == true )
  {
    setMaximumRange(15.625, 32000.0, -96.0, 6.0);
    setCurrentRange(15.625, 32000.0, -96.0, 6.0);
    setVerticalCoarseGrid(2.0, true);
    setVerticalFineGrid(pow(2.0, 1.0/3.0), false);
    rsDataPlot::useLogarithmicScaleX(true);
    setAxisPositionY(rsPlotSettings::LEFT);
    setStringConversionForAxisX(&valueToString0);
  }
  else
  {
    setMaximumRange(0.0, 21000.0, -96.0, 6.0);
    setCurrentRange(0.0, 21000.0, -96.0, 6.0);
    setVerticalCoarseGrid(1000.0, true);
    setVerticalFineGrid(100.0, false);
    rsDataPlot::useLogarithmicScaleX(false);
    setAxisPositionY(rsPlotSettings::ZERO);
    setStringConversionForAxisX(&kiloToString0);
  }

  updateBackgroundImage();
  setAutoReRendering(true);
}

void rsSpectrumPlot::setSpectrum(int newNumBins, double *newBinFrequencies, 
  double *newBinMagnitudes)
{
  rsDataPlot::setCurveValues(newNumBins, newBinFrequencies, newBinMagnitudes);
}

void rsSpectrumPlot::setSpectra(int newNumBins, int newNumSpectra, double* newBinFrequencies, 
  double** newBinMagnitudes)
{
  rsDataPlot::setFunctionFamilyValues(newNumBins, newNumSpectra, newBinFrequencies, 
    newBinMagnitudes);
}

/*
void rsSpectrumPlot::mouseMove(const MouseEvent &e)
{
  setDescription( getCoordinateStringAtPixelPosition(e.x, e.y) );
}
*/

void rsSpectrumPlot::getDisplayedFrequencies(double *frequencies, int numBins)
{
  double x, y;
  for(int i=0; i<numBins; i++)
  {
    x = (double) i * getWidth() / (double) numBins;
    y = 1.0; // only a dummy
    fromPixelCoordinates(x, y);

    frequencies[i] = x;
  }
}

/*
String rsSpectrumPlot::getCoordinateStringAtPixelPosition(int x, int y)
{
  double frequency = (double) x;
  double level     = (double) y;
  fromPixelCoordinates(frequency, level);
  String fString = hertzToStringWithUnitTotal5(frequency);
  String aString = decibelsToStringWithUnit2(level);
  return fString + String(T(", ")) + aString;
}
*/