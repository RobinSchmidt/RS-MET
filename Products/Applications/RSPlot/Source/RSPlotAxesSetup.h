#ifndef __JUCE_RSPLOTAXESSETUP_JUCEHEADER__
#define __JUCE_RSPLOTAXESSETUP_JUCEHEADER__

// old:
////#include "../../../../third_party_code/juce/juce.h"
//#include "../../../rojue/includesForRojue.h"

//#include "../../../rojue/components/widgets/rojue_RButton.h"
//#include "../../../rojue/components/widgets/rojue_RTextField.h"
//using namespace rojue;

#include "../JuceLibraryCode/JuceHeader.h"
using namespace juce;
using namespace jura;

class RSPlotAxesSetup	: public Component
{

public:

  RSPlotAxesSetup();
  /**< Constructor. */

  virtual ~RSPlotAxesSetup();
  /**< Destructor. */


  virtual void paint(Graphics& g) override;
  virtual void resized() override;


 //=============================================================================================
 juce_UseDebuggingNewOperator;


  // the widgets are public in order to let the outlying 
  // RSPlotContentComponent listen to them and access them

  //-----------------------------------------------------------------------------------------------
  // widgets for controlling the axes:

  Label*            xAxisSetupLabel;
  Label*            xAnnotationLabel;
  Label*            xAnnotationEditLabel;
  Label*            xMinLabel;
  Label*            xMinEditLabel;
  Label*            xMaxLabel;
  Label*            xMaxEditLabel;
  Label*            xDigitsLabel;
  ComboBox*         xDigitsComboBox;
  Label*            xPosLabel;
  ComboBox*         xPosComboBox;
  RButton*          xLogScaleButton;

  Label*            yAxisSetupLabel;
  Label*            yAnnotationLabel;
  Label*            yAnnotationEditLabel;
  Label*            yMinLabel;
  Label*            yMinEditLabel;
  Label*            yMaxLabel;
  Label*            yMaxEditLabel;
  Label*            yDigitsLabel;
  ComboBox*         yDigitsComboBox;
  Label*            yPosLabel;
  ComboBox*         yPosComboBox;
  RButton*          yLogScaleButton;

  /*
  Label*            zAxisSetupLabel;
  Label*            zAnnotationLabel;
  Label*            zAnnotationEditLabel;
  Label*            zMinLabel;
  Label*            zMinEditLabel;
  Label*            zMaxLabel;
  Label*            zMaxEditLabel;
  ToggleButton*     zLogScaleButton;
  ComboBox*         zAxisComboBox;

  Label*            scaleSetupLabel;
  Label*            scaleLabelX;
  Slider*           scaleSliderX;
  Label*            scaleLabelY;
  Slider*           scaleSliderY;
  Label*            scaleLabelZ;
  Slider*           scaleSliderZ;
  */

  //-----------------------------------------------------------------------------------------------
  // widgets for controlling the grids in 2D coordinate-systems:

  Label*            gridSetupLabel;

  Label*            horizontalGridSetupLabel;
  RButton*          horizontalCoarseGridButton;
  Label*            horizontalCoarseGridIntervalLabel;
  RButton*          horizontalFineGridButton;
  Label*            horizontalFineGridIntervalLabel;

  Label*            verticalGridSetupLabel;
  RButton*          verticalCoarseGridButton;
  Label*            verticalCoarseGridIntervalLabel;
  RButton*          verticalFineGridButton;
  Label*            verticalFineGridIntervalLabel;

  Label*            radialGridSetupLabel;
  RButton*          radialCoarseGridButton;
  Label*            radialCoarseGridIntervalLabel;
  RButton*          radialFineGridButton;
  Label*            radialFineGridIntervalLabel;

  Label*            angularGridSetupLabel;
  RButton*          angularCoarseGridButton;
  Label*            angularCoarseGridIntervalLabel;
  RButton*          angularFineGridButton;
  Label*            angularFineGridIntervalLabel;

  //-----------------------------------------------------------------------------------------------
  // widgets for controlling the 3D parameters:

  /*
  Label*            rotationSetupLabel;
  Label*            rotationLabelX;
  Slider*           rotationSliderX;
  Label*            rotationLabelY;
  Slider*           rotationSliderY;
  Label*            rotationLabelZ;
  Slider*           rotationSliderZ;
  */
};

#endif	
