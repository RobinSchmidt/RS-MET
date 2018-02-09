#ifndef jura_CurveFamilyPlotOld_h
#define jura_CurveFamilyPlotOld_h

class JUCE_API CurveFamilyPlotOld : virtual public CoordinateSystemOld
{

public:

  CurveFamilyPlotOld(const juce::String& name = juce::String("CurveFamilyPlotOld") );
  virtual ~CurveFamilyPlotOld();

  virtual void paint(Graphics &g) override;
  virtual void resized() override;

  //-----------------------------------------------------------------------------------------------
  // data setup:

  /** Accepts the new function or curve values to be drawn. The data will NOT be copied into
  internal buffers, an object of class "CurveFamilyPlotOld" will only hold pointers to the data. If
  the outlying class for some reason deletes the array in which the data is stored, it should call
  invalidatePointers() - such that the "CurveFamilyPlotOld"-object does not try to aceess the
  array anymore. */
  virtual void setCurveFamilyValues(int newNumValues, int newNumCurves,
    double**  newFamilyValuesX, double** newFamilyValuesY);

  /** Accepts the new function or curve values to be drawn. The newValuesX and newValuesY are
  assumed to arrays of length newNumValues and this class will treat this function-call as if
  you would have called curveFamilyValues with a number of curves of one. */
  virtual void setCurveValues(int newNumValues, double* newValuesX, double* newValuesY);

  /** Accepts values for a fucntion family.  */
  virtual void setFunctionFamilyValues(int newNumValues, int newNumFunctions,
    double* newValuesX, double** newFamilyValuesY);

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
  // plotting:

  /** Returns the drawing as SVG compliant XmlElement. The caller must take care to delete the
  pointer to the XmlElement when it's not needed anymore. */
  XmlElement* getPlotAsSVG(int width, int height) override;

  /** Renders the plot to an image object of given width and height. The caller must take
  care to delete the image when it's not needed anymore. */
  juce::Image* getPlotAsImage(int width, int height) override;

  /** Overrides the inherited method from the CoordinateSystem base-class. */
  virtual void updatePlotImage(bool redrawBackground = false);

  /** Overrides the inherited method from the CoordinateSystem base-class. */
  void updateBackgroundImage() override;

  /** Returns the colour for one of the curves to be drawn. */
  virtual Colour getCurveColour(int index) const { return plotColourScheme.getCurveColour(index); }


protected:

  /** Plots the curve family either on a graphics drawing canvas or an image. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL,
    XmlElement *targetSVG = NULL);

  /** Plots the curve with index 'index'. */
  virtual void plotCurve(Graphics &g, juce::Image *targetImage, XmlElement *targetSV, int index);


  //OwnedArray<Colour> curveColours;	  // array which holds the colurs for the graphs
  //Colour  graphColour;

  int numCurves;	      // number of curves in memory
  int numValues;       	// number of values per curve
  int numCurvesToDraw;  // the number of curves which should be drawn (should be <= numCurves)

  double** familyValuesX;
  double** familyValuesY;
  // pointer to a two-dimensional array. first index indicates the curve, second index indicates 
  // the particular y-value

  double*  valuesX1;	      // pointer to the first array of x-values
  double*  valuesY1;	      // pointer to the first array of x-values

  //double**  decimatedValuesX;
  //double*** decimatedFamilyValuesY;
  // similar to the arrays above but with recursively decimated peak-data - 
  // the first index represents the decimation level

  bool	  fillAreaUnderFunction;   // good for spectra
  bool    isFunctionFamily;
  int     highlightedCurve;
  juce::Image*  plotImage;

  juce_UseDebuggingNewOperator;
};

#endif
