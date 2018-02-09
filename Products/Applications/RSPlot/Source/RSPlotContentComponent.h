/**

This program can visualize mathematical functions based on an expression or imported data in the 
.dat format. In order to work properly, the dat-files must adhere to the following conventions:

-arrays are separated by lines-feeds
-values are separated by whitespaces (1 or 2 are allowed)
-for curve-families the interpretation of the data is: 
 1. line: x-values for the 1st curve
 2. line: y-values for the 1st curve
 3. line: x-values for the 2nd curve
 4. line: y-values for the 2nd curve 
 ...and so on - this implies that the number of lines must be even!

 The data import has been written in such a way as to easily deal with data generated from 
 MatLab/Octave saved as ASCII double, i.e. with the instruction
 save 'ganzonk.dat' x1 y1 x2 y2 x3 y3 -ASCII -DOUBLE

*/

#ifndef __JUCE_RSPLOTCONTENTCOMPONENT_JUCEHEADER__
#define __JUCE_RSPLOTCONTENTCOMPONENT_JUCEHEADER__

#include "RSPlotAxesSetup.h"
#include "RSPlotDataSetup.h"
#include "RSPlotFileSetup.h"

class RSPlotContentComponent : public Component, public StateFileManager, 
  public ImageFileManager, 
  public RButtonListener, public RTextEntryFieldObserver,
  public ComboBoxListener, public LabelListener, public SliderListener // get rid
{  

public:

  enum plotterComponents
  {
    CURVE_FAMILY_PLOT = 0,
    SURFACE_PLOT
  };

  enum visualizationModes
  {
    FUNCTION_ON_PLANE = 1,
    FUNCTION_FAMILY_ON_PLANE,
    CURVE_ON_PLANE,
    CURVE_FAMILY_ON_PLANE,
    SCATTER_PLOT_ON_PLANE,
    CURVE_IN_SPACE,
    SURFACE_IN_SPACE,
    CLOUD_IN_SPACE,
    COMPLEX_MAPPING
  };

  RSPlotContentComponent(const String &newEditorName);  
  ///< Constructor.

  virtual ~RSPlotContentComponent();                             
  ///< Destructor.

  //virtual void buttonClicked(Button *buttonThatWasClicked) override;
  // // get rid - use RButtons everywhere...

  virtual void rButtonClicked(RButton *buttonThatWasClicked) override;
  /**< Implements the purely virtual buttonClicked()-method of the 
       ButtonListener base-class. */

  virtual void comboBoxChanged(ComboBox *comboBoxThatHasChanged);
  /**< Implements the purely virtual comboBoxChanged()-method of the 
       ComboBoxListener base-class. */

  virtual void labelTextChanged(Label *labelThatHasChanged); 
  /**< Implements the purely virtual labelTextChanged()-method of the 
       LablelListener base-class. */

  virtual void textChanged(RTextEntryField *field) override;

  virtual void sliderValueChanged(Slider *sliderThatHasChanged); 
  /**< Implements the purely virtual sliderValueChanged()-method of the 
       SliderListener base-class. */

  virtual void setBounds(int x, int y, int width, int height);
  /**< Overrides the setBounds()-method of the Component base-class in order
       to arrange the widgets according to the size. */

 /** Overrides the getStateAsXml()-method of the PresetFileManager base-class. */
 virtual XmlElement* getStateAsXml() const;

 /** Overrides the setStateFromXml()-method of the PresetFileManager base-class. */
 virtual bool setStateFromXml(const XmlElement& xmlState);


 // this is a bit preliminary:
 virtual XmlElement* getStateAsXml(const String& stateName, bool markAsClean)
 { return getStateAsXml(); }

 virtual void setStateFromXml(const XmlElement& xmlState, const String& stateName, 
                              bool markAsClean)
 { setStateFromXml(xmlState); }



 //=============================================================================================
 juce_UseDebuggingNewOperator;

protected:

  /** This function will be invoked to calculate the actual data from the expression strings. */
  void calculateData();

  /** Updates the result fields of the calculator. */
  void updateCalculatorResult();

  /** Fills all the x- and y-arrays with function family data. */
  void fillArraysWithFamilyData();

  /** Fills one of the x-arrays with data according to the passed expression-string. */
  void fillArrayWithDataX(const String &xString, int curveIndex);

  /** Fills one of the y-arrays with data according to the passed expression-string. */
  void fillArrayWithDataY(const String &yString, int curveIndex);


  // the components for the actual curve- and function-plots:
  rsDataPlot*          curveFamilyPlot;
  rsPlotZoomer*   zoomer2D;    
  //SurfacePlot*              surfacePlot;
  //CoordinateSystem3DZoomer* zoomer3D;                

  // the different setups:
  RSPlotAxesSetup* axesSetup;
  RSPlotDataSetup* dataSetup;
  RSPlotFileSetup* fileSetup;

  // the tabber to choose which of the setups should be shown:
  TabbedComponent* setupTabber; 

  int currentPlotterComponent;

  static const int maxNumCurves = 9;
  static const int maxNumValues = 1024;

  // arrays in which we store the actual data:
  double xValues[maxNumCurves][maxNumValues];
  double yValues[maxNumCurves][maxNumValues];

  // we also need arrays pointers to double - each array element is a pointer to the first value
  // of the respective data-array:
  double** xFamilyPointer;
  double** yFamilyPointer;

  void openDataLoadingDialog();
  void loadDataFromFile(const File &dataFileToLoad);
  void parseDataString(const String& dataString);
  String dataFileName;

  rosic::ExpressionEvaluator evaluator;
   // embedded rosic::ExpressionEvaluator object will perform the expression evaluation.

};

#endif  