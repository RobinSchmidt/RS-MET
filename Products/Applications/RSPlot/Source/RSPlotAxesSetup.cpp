#include "RSPlotAxesSetup.h"

RSPlotAxesSetup::RSPlotAxesSetup()
{
  //-------------------------------------------------------------------------------------
  // x-axis setup widgets:

  xAxisSetupLabel = new Label( String("xAxisSetupLabel"), String(("x-axis setup:")) );
  xAxisSetupLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( xAxisSetupLabel );

  xAnnotationLabel = new Label( String(("xAnnotationLabel")), String(("Label:")) );
  xAnnotationLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( xAnnotationLabel );

  xAnnotationEditLabel = new Label( String(("xAnnotationEditLabel")), String(("x")) );
  xAnnotationEditLabel->setJustificationType(Justification::centredLeft);
  xAnnotationEditLabel->setColour(Label::outlineColourId, Colours::black);
  xAnnotationEditLabel->setEditable(true, true);
  addAndMakeVisible( xAnnotationEditLabel );

  xMinLabel = new Label( String(("MinXLabel")), String(("Min:")) );
  xMinLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( xMinLabel );

  xMinEditLabel = new Label( String(("MinXEditLabel")), String(("-4.5")) );
  xMinEditLabel->setJustificationType(Justification::centredLeft);
  xMinEditLabel->setColour(Label::outlineColourId, Colours::black);
  xMinEditLabel->setEditable(true, true);
  addAndMakeVisible( xMinEditLabel );

  xMaxLabel = new Label( String(("MaxXLabel")), String(("Max:")) );
  xMaxLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( xMaxLabel );

  xMaxEditLabel = new Label( String(("MaxXEditLabel")), String(("4.5")) );
  xMaxEditLabel->setJustificationType(Justification::centredLeft);
  xMaxEditLabel->setColour(Label::outlineColourId, Colours::black);
  xMaxEditLabel->setEditable(true, true);
  addAndMakeVisible( xMaxEditLabel );

  xDigitsLabel = new Label( String(("DigitsXLabel")), String(("Digits:")) );
  xDigitsLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( xDigitsLabel );

  xDigitsComboBox = new ComboBox(String(("xDigitsComboBox")));
  xDigitsComboBox->addItem(String(("0")),       1);
  xDigitsComboBox->addItem(String(("1")),       2);
  xDigitsComboBox->addItem(String(("2")),       3);
  xDigitsComboBox->setSelectedId(1);
  addAndMakeVisible( xDigitsComboBox );

  xPosLabel = new Label( String(("xPosLabel")), String(("Pos:")) );
  xPosLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( xPosLabel );

  xPosComboBox = new ComboBox(String(("xPosComboBox")));
  xPosComboBox->addItem(String(("Hide")),       1);
  xPosComboBox->addItem(String(("Zero")),       2);
  xPosComboBox->addItem(String(("Top")),        3);
  xPosComboBox->addItem(String(("Bottom")),     4);
  xPosComboBox->setSelectedId(2);
  addAndMakeVisible( xPosComboBox );

  xLogScaleButton = new RButton(String("Log"));
  xLogScaleButton->setClickingTogglesState(true);
  xLogScaleButton->setToggleState(false, true);
  addAndMakeVisible(xLogScaleButton);

  //-------------------------------------------------------------------------------------
  // y-axis setup widgets:

  yAxisSetupLabel = new Label( String(("yAxisSetupLabel")), String(("y-axis setup:")) );
  yAxisSetupLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( yAxisSetupLabel );

  yAnnotationLabel = new Label( String(("yAnnotationLabel")), String(("Label:")) );
  yAnnotationLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( yAnnotationLabel );

  yAnnotationEditLabel = new Label( String(("yAnnotationEditLabel")), String(("y")) );
  yAnnotationEditLabel->setJustificationType(Justification::centredLeft);
  yAnnotationEditLabel->setColour(Label::outlineColourId, Colours::black);
  yAnnotationEditLabel->setEditable(true, true);
  addAndMakeVisible( yAnnotationEditLabel );

  yMinLabel = new Label( String(("MinYLabel")), String(("Min:")) );
  yMinLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( yMinLabel );

  yMinEditLabel = new Label( String(("MinYEditLabel")), String(("-4.5")) );
  yMinEditLabel->setJustificationType(Justification::centredLeft);
  yMinEditLabel->setColour(Label::outlineColourId, Colours::black);
  yMinEditLabel->setEditable(true, true);
  addAndMakeVisible( yMinEditLabel );

  yMaxLabel = new Label( String(("MaxYLabel")), String(("Max:")) );
  yMaxLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( yMaxLabel );

  yMaxEditLabel = new Label( String(("MaxYEditLabel")), String(("4.5")) );
  yMaxEditLabel->setJustificationType(Justification::centredLeft);
  yMaxEditLabel->setColour(Label::outlineColourId, Colours::black);
  yMaxEditLabel->setEditable(true, true);
  addAndMakeVisible( yMaxEditLabel );

  yDigitsLabel = new Label( String(("DigitsXLabel")), String(("Digits:")) );
  yDigitsLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( yDigitsLabel );

  yDigitsComboBox = new ComboBox(String(("yDigitsComboBox")));
  yDigitsComboBox->addItem(String(("0")),       1);
  yDigitsComboBox->addItem(String(("1")),       2);
  yDigitsComboBox->addItem(String(("2")),       3);
  yDigitsComboBox->setSelectedId(1);
  addAndMakeVisible( yDigitsComboBox );

  yPosLabel = new Label( String(("yPosLabel")), String(("Pos:")) );
  yPosLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( yPosLabel );

  yPosComboBox = new ComboBox(String(("yPosComboBox")));
  yPosComboBox->addItem(String(("Hide")),       1);
  yPosComboBox->addItem(String(("Zero")),       2);
  yPosComboBox->addItem(String(("Left")),       3);
  yPosComboBox->addItem(String(("Right")),      4);
  yPosComboBox->setSelectedId(2);
  addAndMakeVisible( yPosComboBox );

  yLogScaleButton = new RButton(String("Log"));
  yLogScaleButton->setClickingTogglesState(true);
  yLogScaleButton->setToggleState(false, true);
  addAndMakeVisible(yLogScaleButton);

  //-------------------------------------------------------------------------------------
  // z-axis setup widgets:

  /*
  zAxisSetupLabel = new Label( String(("zAxisSetupLabel")), String(("z-axis setup:")) );
  zAxisSetupLabel->setJustificationType(Justification::centred);
  addChildComponent( zAxisSetupLabel );

  zAnnotationLabel = new Label( String(("zAnnotationLabel")), String(("annotation:")) );
  zAnnotationLabel->setJustificationType(Justification::centredLeft);
  addChildComponent( zAnnotationLabel );

  zAnnotationEditLabel = new Label( String(("zAnnotationEditLabel")), String(("z")) );
  zAnnotationEditLabel->setJustificationType(Justification::centredLeft);
  zAnnotationEditLabel->setColour(Label::outlineColourId, Colours::black);
  zAnnotationEditLabel->setEditable(true, true);
  addChildComponent( zAnnotationEditLabel );

  minZLabel = new Label( String(("MinZLabel")), String(("min:")) );
  minZLabel->setJustificationType(Justification::centredLeft);
  addChildComponent( minZLabel );

  minZEditLabel = new Label( String(("MinZEditLabel")), String(("-4.5")) );
  minZEditLabel->setJustificationType(Justification::centredLeft);
  minZEditLabel->setColour(Label::outlineColourId, Colours::black);
  minZEditLabel->setEditable(true, true);
  addChildComponent( minZEditLabel );

  maxZLabel = new Label( String(("MaxZLabel")), String(("max:")) );
  maxZLabel->setJustificationType(Justification::centredLeft);
  addChildComponent( maxZLabel );

  maxZEditLabel = new Label( String(("MaxZEditLabel")), String(("4.5")) );
  maxZEditLabel->setJustificationType(Justification::centredLeft);
  maxZEditLabel->setColour(Label::outlineColourId, Colours::black);
  maxZEditLabel->setEditable(true, true);
  addChildComponent( maxZEditLabel );

  logScaleZButton = new ToggleButton(String("log"));
  logScaleZButton->setClickingTogglesState(true);
  logScaleZButton->setToggleState(false, true);
  addChildComponent(logScaleZButton);

  zAxisComboBox = new ComboBox(String(("zAxisComboBox")));
  zAxisComboBox->addItem(String(("no axis")),     1);
  zAxisComboBox->addItem(String(("at zero")),     2);
  zAxisComboBox->addItem(String(("at left")),     3);
  zAxisComboBox->addItem(String(("at right")),    4);
  zAxisComboBox->setSelectedId(2);
  addChildComponent( zAxisComboBox );
  */

  //-----------------------------------------------------------------------------------------------
  // grid setup widgets:

  gridSetupLabel = new Label(String(("gridSetupLabel")), 
    String(("grid setup:")) );
  gridSetupLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( gridSetupLabel );

  //-------------------------------------------------------------------------------------
  // horizontal grid setup widgets:

  horizontalGridSetupLabel = new Label(String(("horizontalGridSetupLabel")), 
    String(("Horizontal Grid:")) );
  horizontalGridSetupLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( horizontalGridSetupLabel );

  horizontalCoarseGridButton = new RButton(String("Coarse:"));
  horizontalCoarseGridButton->setClickingTogglesState(true);
  horizontalCoarseGridButton->setToggleState(false, true);
  addAndMakeVisible(horizontalCoarseGridButton);

  horizontalCoarseGridIntervalLabel = 
    new Label( String(("horizontalCoarseGridIntervalLabel")), String(("1.0")) );
  horizontalCoarseGridIntervalLabel->setJustificationType(Justification::centredLeft);
  horizontalCoarseGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  horizontalCoarseGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( horizontalCoarseGridIntervalLabel );

  horizontalFineGridButton = new RButton(String("Fine:"));
  horizontalFineGridButton->setClickingTogglesState(true);
  horizontalFineGridButton->setToggleState(false, true);
  addAndMakeVisible(horizontalFineGridButton);

  horizontalFineGridIntervalLabel = 
    new Label( String(("horizontalFineGridIntervalLabel")), String(("0.1")) );
  horizontalFineGridIntervalLabel->setJustificationType(Justification::centredLeft);
  horizontalFineGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  horizontalFineGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( horizontalFineGridIntervalLabel );

  //-------------------------------------------------------------------------------------
  // vertical grid setup widgets:

  verticalGridSetupLabel = new Label(String(("verticalGridSetupLabel")), 
    String(("Vertical Grid:")) );
  verticalGridSetupLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( verticalGridSetupLabel );

  verticalCoarseGridButton = new RButton(String("Coarse:"));
  verticalCoarseGridButton->setClickingTogglesState(true);
  verticalCoarseGridButton->setToggleState(false, true);
  addAndMakeVisible(verticalCoarseGridButton);

  verticalCoarseGridIntervalLabel = 
    new Label( String(("verticalCoarseGridIntervalLabel")), String(("1.0")) );
  verticalCoarseGridIntervalLabel->setJustificationType(Justification::centredLeft);
  verticalCoarseGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  verticalCoarseGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( verticalCoarseGridIntervalLabel );

  verticalFineGridButton = new RButton(String("Fine:"));
  verticalFineGridButton->setClickingTogglesState(true);
  verticalFineGridButton->setToggleState(false, true);
  addAndMakeVisible(verticalFineGridButton);

  verticalFineGridIntervalLabel = 
    new Label( String(("verticalFineGridIntervalLabel")), String(("0.1")) );
  verticalFineGridIntervalLabel->setJustificationType(Justification::centredLeft);
  verticalFineGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  verticalFineGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( verticalFineGridIntervalLabel );

  //-------------------------------------------------------------------------------------
  // radial grid setup widgets:

  radialGridSetupLabel = new Label(String(("radialGridSetupLabel")), 
    String(("Radial Grid:")) );
  radialGridSetupLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( radialGridSetupLabel );

  radialCoarseGridButton = new RButton(String("Coarse:"));
  radialCoarseGridButton->setClickingTogglesState(true);
  radialCoarseGridButton->setToggleState(false, true);
  addAndMakeVisible(radialCoarseGridButton);

  radialCoarseGridIntervalLabel = 
    new Label( String(("radialCoarseGridIntervalLabel")), String(("1.0")) );
  radialCoarseGridIntervalLabel->setJustificationType(Justification::centredLeft);
  radialCoarseGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  radialCoarseGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( radialCoarseGridIntervalLabel );

  radialFineGridButton = new RButton(String("Fine:"));
  radialFineGridButton->setClickingTogglesState(true);
  radialFineGridButton->setToggleState(false, true);
  addAndMakeVisible(radialFineGridButton);

  radialFineGridIntervalLabel = 
    new Label( String(("radialFineGridIntervalLabel")), String(("0.1")) );
  radialFineGridIntervalLabel->setJustificationType(Justification::centredLeft);
  radialFineGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  radialFineGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( radialFineGridIntervalLabel );

  //-------------------------------------------------------------------------------------
  // angular grid setup widgets:

  angularGridSetupLabel = new Label(String(("angularGridSetupLabel")), 
    String(("Angular Grid:")) );
  angularGridSetupLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( angularGridSetupLabel );

  angularCoarseGridButton = new RButton(String("Coarse:"));
  angularCoarseGridButton->setClickingTogglesState(true);
  angularCoarseGridButton->setToggleState(false, true);
  addAndMakeVisible(angularCoarseGridButton);

  angularCoarseGridIntervalLabel = 
    new Label( String(("angularCoarseGridIntervalLabel")), String(("15")) );
  angularCoarseGridIntervalLabel->setJustificationType(Justification::centredLeft);
  angularCoarseGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  angularCoarseGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( angularCoarseGridIntervalLabel );

  angularFineGridButton = new RButton(String("Fine:"));
  angularFineGridButton->setClickingTogglesState(true);
  angularFineGridButton->setToggleState(false, true);
  addAndMakeVisible(angularFineGridButton);

  angularFineGridIntervalLabel = 
    new Label( String(("angularFineGridIntervalLabel")), String(("5")) );
  angularFineGridIntervalLabel->setJustificationType(Justification::centredLeft);
  angularFineGridIntervalLabel->setColour(Label::outlineColourId, Colours::black);
  angularFineGridIntervalLabel->setEditable(true, true);
  addAndMakeVisible( angularFineGridIntervalLabel );

  //-----------------------------------------------------------------------------------------------
  // rotation setup widgets:

  /*
  rotationSetupLabel = new Label( String(("rotationSetupLabel")), String(("rotation setup:")) );
  rotationSetupLabel->setJustificationType(Justification::centred);
  addChildComponent( rotationSetupLabel );


  rotationLabelX = new Label( String(("rotationLabelX")), String(("x:")) );
  rotationLabelX->setJustificationType(Justification::centred);
  addChildComponent( rotationLabelX );

  rotationSliderX = new Slider( String(("rotationSliderX")) );
  rotationSliderX->setRange(-90.0, 90.0, 1.0);
  rotationSliderX->setValue(35.0);
  rotationSliderX->setSliderStyle(Slider::LinearBar);
  addChildComponent( rotationSliderX );

  rotationLabelY = new Label( String(("rotationLabelY")), String(("y:")) );
  rotationLabelY->setJustificationType(Justification::centred);
  addChildComponent( rotationLabelY );

  rotationSliderY = new Slider( String(("rotationSliderY")) );
  rotationSliderY->setRange(-90.0, 90.0, 1.0);
  rotationSliderY->setValue(-45.0);
  rotationSliderY->setSliderStyle(Slider::LinearBar);
  addChildComponent( rotationSliderY );

  rotationLabelZ = new Label( String(("rotationLabelZ")), String(("z:")) );
  rotationLabelZ->setJustificationType(Justification::centred);
  addChildComponent( rotationLabelZ );

  rotationSliderZ = new Slider( String(("rotationSliderZ")) );
  rotationSliderZ->setRange(-90.0, 90.0, 1.0);
  rotationSliderZ->setValue(0.0);
  rotationSliderZ->setSliderStyle(Slider::LinearBar);
  addChildComponent( rotationSliderZ );

  */
}

RSPlotAxesSetup::~RSPlotAxesSetup()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// appearance:

void RSPlotAxesSetup::paint(Graphics& g)
{
  g.fillAll(Colours::black);
}

void RSPlotAxesSetup::resized()
{
  Component::resized();

  int xL, xR, yT, yB; // left, right, top and bottom for the child-component
  int w, h;           // width and height for the child component

  xL = 0;
  xR = getWidth()/4;
  w  = xR-xL;
  yT = 0;
  yB = getHeight();
  h  = yB-yT;

  xAxisSetupLabel->setBounds(xL, yT, w, 16);
  yT += 16;
  xAnnotationLabel->setBounds(xL+4, yT, 56-8, 16);
  xAnnotationEditLabel->setBounds(xL+56, yT, w-56-4, 16);
  yT += 16;
  xMinLabel->setBounds(xL+4, yT, 56-8, 16);
  xMinEditLabel->setBounds(xL+56, yT, w-56-4, 16);
  yT += 16;
  xMaxLabel->setBounds(xL+4, yT, 56-8, 16);
  xMaxEditLabel->setBounds(xL+56, yT, w-56-4, 16);
  yT += 16;
  xDigitsLabel->setBounds(xL+4, yT, 56-8, 20);
  xDigitsComboBox->setBounds(xL+56, yT, w-56-4, 20);
  yT += 20;
  xPosLabel->setBounds(xL+4, yT, 56-8, 20);
  xPosComboBox->setBounds(xL+56, yT, w-56-4, 20);
  yT += 20;
  xLogScaleButton->setBounds(xL+w-40-4, yT, 40, 20);

  yT = 0;
  xL += w;

  yAxisSetupLabel->setBounds(xL, yT, w, 16);
  yT += 16;
  yAnnotationLabel->setBounds(xL+4, yT, 56-8, 16);
  yAnnotationEditLabel->setBounds(xL+56, yT, w-56-4, 16);
  yT += 16;
  yMinLabel->setBounds(xL+4, yT, 56-8, 16);
  yMinEditLabel->setBounds(xL+56, yT, w-56-4, 16);
  yT += 16;
  yMaxLabel->setBounds(xL+4, yT, 56-8, 16);
  yMaxEditLabel->setBounds(xL+56, yT, w-56-4, 16);
  yT += 16;
  yDigitsLabel->setBounds(xL+4, yT, 56-8, 20);
  yDigitsComboBox->setBounds(xL+56, yT, w-56-4, 20);
  yT += 20;
  yPosLabel->setBounds(xL+4, yT, 56-8, 20);
  yPosComboBox->setBounds(xL+56, yT, w-56-4, 20);
  yT += 20;
  yLogScaleButton->setBounds(xL+w-40-4, yT, 40, 20);

  yT = 0;
  xL += w;

  horizontalGridSetupLabel->setBounds(xL, yT, w, 16);
  yT += 16;
  horizontalCoarseGridButton->setBounds(xL+4, yT, 72-8, 20);
  horizontalCoarseGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);
  yT += 20;
  horizontalFineGridButton->setBounds(xL+4, yT, 72-8, 20);
  horizontalFineGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);
  yT += 24;
  verticalGridSetupLabel->setBounds(xL, yT, w, 16);
  yT += 16;
  verticalCoarseGridButton->setBounds(xL+4, yT, 72-8, 20);
  verticalCoarseGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);
  yT += 20;
  verticalFineGridButton->setBounds(xL+4, yT, 72-8, 20);
  verticalFineGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);

  yT = 0;
  xL += w;

  radialGridSetupLabel->setBounds(xL, yT, w, 16);
  yT += 16;
  radialCoarseGridButton->setBounds(xL+4, yT, 72-8, 20);
  radialCoarseGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);
  yT += 20;
  radialFineGridButton->setBounds(xL+4, yT, 72-8, 20);
  radialFineGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);
  yT += 24;
  angularGridSetupLabel->setBounds(xL, yT, w, 16);
  yT += 16;
  angularCoarseGridButton->setBounds(xL+4, yT, 72-8, 20);
  angularCoarseGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);
  yT += 20;
  angularFineGridButton->setBounds(xL+4, yT, 72-8, 20);
  angularFineGridIntervalLabel->setBounds(xL+72, yT, w-72-4, 20);
}
