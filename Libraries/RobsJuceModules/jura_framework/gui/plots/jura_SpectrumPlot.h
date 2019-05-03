#ifndef jura_SpectrumPlot_h
#define jura_SpectrumPlot_h

/** This class is for plotting spectra and frequency responses. */

class rsSpectrumPlot : virtual public rsDataPlot
{

public:

  /** Constructor. */
  rsSpectrumPlot(const juce::String& name = juce::String("SpectrumDisplay"));

  /** Destructor. */
  virtual ~rsSpectrumPlot();

  /** Overrides the inherited function from CoordinateSystem in order to adjust the 
  grid-spacing. */
  virtual void useLogarithmicScaleX(bool shouldBeLogScaledX, double newLogBase = 2.0);

  /** Accepts values for a single spectrum.  */
  virtual void setSpectrum(int newNumBins, double* newBinFrequencies, double* newBinMagnitudes);

  /** Accepts values for a family of spectra.  */
  virtual void setSpectra(int newNumBins, int newNumSpectra, double* newBinFrequencies, 
    double** newBinMagnitudes);

  /** Overrides mouseMove. \todo: call getCoordinateStringAtPixelPosition() already in baseclass, 
  get rid of the popup. */
  //virtual void mouseMove(const MouseEvent& e);

  /** Fills the 'numBins' long array with the frequencies which are ccurrently visible on the 
  display. This can be used to create frequency responses which match with the display. */
  virtual void getDisplayedFrequencies(double* frequencies, int numBins);

  /** Returns a string that represents the coordinates at pixel (x,y). */
  //virtual juce::String getCoordinateStringAtPixelPosition(int x, int y);

  juce_UseDebuggingNewOperator;
};

#endif  
