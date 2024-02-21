#ifndef jura_DataPlot_h
#define jura_DataPlot_h

/** A subclass for plots that are based on data, i.e. arrays of x,y values. The data will NOT be 
copied into internal buffers, an object of class "rsDataPlot" will only hold pointers to the data. 
If the outlying class for some reason deletes the array in which the data is stored, it should call
invalidatePointers() such that the "rsDataPlot" object does not try to aceess the array anymore. */

class JUCE_API rsDataPlot : virtual public rsPlot
{

public:

  rsDataPlot(const juce::String& name = juce::String("rsDataPlot") );
  virtual ~rsDataPlot();

  virtual void paint(Graphics &g) override;
  virtual void resized() override;

  //-----------------------------------------------------------------------------------------------
  // \name Data setup:

  /** Accepts the new function or curve values to be drawn.  */
  virtual void setCurveFamilyValues(int newNumValues, int newNumCurves,
    double**  newFamilyValuesX, double** newFamilyValuesY);

  /** Accepts the new function or curve values to be drawn. The newValuesX and newValuesY are
  assumed to arrays of length newNumValues and this class will treat this function-call as if
  you would have called curveFamilyValues with a number of curves of one. */
  virtual void setCurveValues(int newNumValues, double* newValuesX, double* newValuesY);

  /** Accepts values for a function family.  */
  virtual void setFunctionFamilyValues(int newNumValues, int newNumFunctions,
    double* newValuesX, double** newFamilyValuesY);
  // I think, the difference between a function family and a curve family is that the former uses
  // shared x-axis values for all the graphs - document this better!

  /** Sets the number of curves to be drawn. Can be used to draw not all but only the first few
  curves from the familyValuesY-array. But the number passed here should never be larger than
  the actual number of arrays residing in the memory. */
  virtual void setNumCurves(int newNumCurves) { numCurves = newNumCurves; }

  /** Invalidates the pointers to the data. */
  virtual void invalidatePointers();

  /** Chooses one of the curves to be highlighted, causing the other curves to be drawn more
  faintly (using transparency). If -1 is passed here, all curves will be drawn fully opaque. */
  virtual void setHighlightedCurve(int newHighlightedCurve);


  // \todo: insert setup functions here later that control the appearance of the plot (colors, 
  // line-widths, etc.)


  //-----------------------------------------------------------------------------------------------
  // \name Plotting:

  /** Returns the drawing as SVG compliant XmlElement. The caller must take care to delete the
  pointer to the XmlElement when it's not needed anymore. */
  XmlElement* getPlotAsSVG(int width, int height) override;

  /** Renders the plot to an image object of given width and height. The caller must take
  care to delete the image when it's not needed anymore. */
  juce::Image* getPlotAsImage(int width, int height) override;

  /** Overrides the inherited method from the rsPlot base-class. 
  ...ooohh...no - this is wrong there is such baseclass method...we need to update the comment */
  virtual void updatePlotImage(bool redrawBackground = false);

  /** Overrides the inherited method from the CoordinateSystem base-class. */
  void updateBackgroundImage() override;

  /** Returns the colour for one of the curves to be drawn. */
  virtual Colour getCurveColour(int index) const { return plotColourScheme.getCurveColour(index); }


protected:

  /** Plots the curve family either on a graphics drawing canvas or an image. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL,
    XmlElement *targetSVG = NULL);


  //OwnedArray<Colour> curveColours;	  // array which holds the colurs for the graphs
  //Colour  graphColour;

  int numCurves = 0;        // number of curves in memory
  int numValues = 0;        // number of values per curve
  int numCurvesToDraw = 0;  // the number of curves which should be drawn (should be <= numCurves)

  double** familyValuesX = nullptr;
  double** familyValuesY = nullptr;
  // Pointers to two-dimensional arrays. The first index indicates the curve, second index 
  // indicates the particular x- or y-value.
  // ToDo: use float instead of double

  double* valuesX1 = nullptr;        // pointer to the first array of x-values
  double* valuesY1 = nullptr;        // pointer to the first array of y-values
  // Why do we need these? I think, they may be there for convenience for the common case to draw
  // only one curve or multiple curves with shared x-values for a function family.
  
  //double**  decimatedValuesX;
  //double*** decimatedFamilyValuesY;
  // similar to the arrays above but with recursively decimated peak-data - 
  // the first index represents the decimation level - maybe to this in subclass

  bool fillAreaUnderFunction = false;   // good for spectra - not yet used
  bool isFunctionFamily = false;        // indicates that x-axis is shared - rename to xAxisShared
  int  highlightedCurve = -1;           // -1 means: non
  juce::Image*  plotImage = nullptr;

  juce_UseDebuggingNewOperator;
};

#endif
