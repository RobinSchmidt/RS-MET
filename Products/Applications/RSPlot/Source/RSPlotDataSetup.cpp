#include "RSPlotDataSetup.h"

RSPlotDataSetup::RSPlotDataSetup()
{
  typeOfDataLabel = new Label( String(("typeOfDataLabel")), String(("Type:")) );
  typeOfDataLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( typeOfDataLabel );

  typeOfDataComboBox = new ComboBox(String(("typeOfDataComboBox")));
  typeOfDataComboBox->addItem(String(("Curve Family")),          1);
  typeOfDataComboBox->addItem(String(("Surface Plot")),          2);
  typeOfDataComboBox->setItemEnabled(2, false);
  typeOfDataComboBox->setSelectedId(1);
  addAndMakeVisible( typeOfDataComboBox );

  dataSourceLabel = new Label( String(("dataSourceLabel")), String(("Source:")) );
  dataSourceLabel->setJustificationType(Justification::centred);
  addAndMakeVisible( dataSourceLabel );

  dataSourceComboBox = new ComboBox(String(("dataSourceComboBox")));
  dataSourceComboBox->addItem(String(("Expression")),            1);
  dataSourceComboBox->addItem(String(("Data-File")),             2);
  dataSourceComboBox->setItemEnabled(2, false);
  dataSourceComboBox->setSelectedId(1);
  addAndMakeVisible( dataSourceComboBox );

  //-----------------------------------------------------------------------------------------------
  // widgets for entering the curve- and function-formulas:

  numCurvesLabel = new Label( String(("numCurvesLabel")), String(("Curves:")) );
  numCurvesLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( numCurvesLabel );

  numCurvesSlider = new Slider( String(("numCurvesSlider")) );
  numCurvesSlider->setRange(1.0, 9.0, 1.0);
  numCurvesSlider->setValue(1.0);
  numCurvesSlider->setSliderStyle(Slider::LinearBar);
  addAndMakeVisible( numCurvesSlider );

  xCurveLabel = new Label( String(("xCurveLabel")), String(("x(t)=")) );
  xCurveLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( xCurveLabel );

  xCurveEditLabel = new Label( String(("xCurveEditLabel")), String(("t;")) );
  xCurveEditLabel->setJustificationType(Justification::centredLeft);
  xCurveEditLabel->setColour(Label::outlineColourId, Colours::black);
  xCurveEditLabel->setEditable(true, true);
  addAndMakeVisible( xCurveEditLabel );

  yCurveLabel = new Label( String(("yCurveLabel")), String(("y(x,t)=")) );
  yCurveLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( yCurveLabel );

  yCurveEditLabel = new Label( String(("yCurveEditLabel")), String(("x;")) );
  yCurveEditLabel->setJustificationType(Justification::centredLeft);
  yCurveEditLabel->setColour(Label::outlineColourId, Colours::black);
  yCurveEditLabel->setEditable(true, true);
  addAndMakeVisible( yCurveEditLabel );

  //-----------------------------------------------------------------------------------------------
  // widgets for the parameters a,b,c,d:

  aLabel = new Label( String(("aLabel")), String(("a:")) );
  aLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( aLabel );

  aMinEditLabel = new Label( String(("aMinEditLabel")), String(("0")) );
  aMinEditLabel->setJustificationType(Justification::centredLeft);
  aMinEditLabel->setColour(Label::outlineColourId, Colours::black);
  aMinEditLabel->setEditable(true, true);
  addAndMakeVisible( aMinEditLabel );

  aSlider = new Slider( String(("aSlider")) );
  aSlider->setRange(0.0, 1.0, 0.0);
  aSlider->setValue(0.5);
  aSlider->setSliderStyle(Slider::LinearBar);
  addAndMakeVisible(aSlider);

  aMaxEditLabel = new Label( String(("aMaxEditLabel")), String(("1")) );
  aMaxEditLabel->setJustificationType(Justification::centredLeft);
  aMaxEditLabel->setColour(Label::outlineColourId, Colours::black);
  aMaxEditLabel->setEditable(true, true);
  addAndMakeVisible( aMaxEditLabel );

  bLabel = new Label( String(("bLabel")), String(("b:")) );
  bLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( bLabel );

  bMinEditLabel = new Label( String(("bMinEditLabel")), String(("0")) );
  bMinEditLabel->setJustificationType(Justification::centredLeft);
  bMinEditLabel->setColour(Label::outlineColourId, Colours::black);
  bMinEditLabel->setEditable(true, true);
  addAndMakeVisible( bMinEditLabel );

  bSlider = new Slider( String(("bSlider")) );
  bSlider->setRange(0.0, 1.0, 0.0);
  bSlider->setValue(0.5);
  bSlider->setSliderStyle(Slider::LinearBar);
  addAndMakeVisible(bSlider);

  bMaxEditLabel = new Label( String(("bMaxEditLabel")), String(("1")) );
  bMaxEditLabel->setJustificationType(Justification::centredLeft);
  bMaxEditLabel->setColour(Label::outlineColourId, Colours::black);
  bMaxEditLabel->setEditable(true, true);
  addAndMakeVisible( bMaxEditLabel );

  cLabel = new Label( String(("cLabel")), String(("c:")) );
  cLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( cLabel );

  cMinEditLabel = new Label( String(("cMinEditLabel")), String(("0")) );
  cMinEditLabel->setJustificationType(Justification::centredLeft);
  cMinEditLabel->setColour(Label::outlineColourId, Colours::black);
  cMinEditLabel->setEditable(true, true);
  addAndMakeVisible( cMinEditLabel );

  cSlider = new Slider( String(("cSlider")) );
  cSlider->setRange(0.0, 1.0, 0.0);
  cSlider->setValue(0.5);
  cSlider->setSliderStyle(Slider::LinearBar);
  addAndMakeVisible(cSlider);

  cMaxEditLabel = new Label( String(("cMaxEditLabel")), String(("1")) );
  cMaxEditLabel->setJustificationType(Justification::centredLeft);
  cMaxEditLabel->setColour(Label::outlineColourId, Colours::black);
  cMaxEditLabel->setEditable(true, true);
  addAndMakeVisible( cMaxEditLabel );

  dLabel = new Label( String(("dLabel")), String(("d:")) );
  dLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( dLabel );

  dMinEditLabel = new Label( String(("dMinEditLabel")), String(("0")) );
  dMinEditLabel->setJustificationType(Justification::centredLeft);
  dMinEditLabel->setColour(Label::outlineColourId, Colours::black);
  dMinEditLabel->setEditable(true, true);
  addAndMakeVisible( dMinEditLabel );

  dSlider = new Slider( String(("dSlider")) );
  dSlider->setRange(0.0, 1.0, 0.0);
  dSlider->setValue(0.5);
  dSlider->setSliderStyle(Slider::LinearBar);
  addAndMakeVisible(dSlider);

  dMaxEditLabel = new Label( String(("dMaxEditLabel")), String(("1")) );
  dMaxEditLabel->setJustificationType(Justification::centredLeft);
  dMaxEditLabel->setColour(Label::outlineColourId, Colours::black);
  dMaxEditLabel->setEditable(true, true);
  addAndMakeVisible( dMaxEditLabel );

  //-----------------------------------------------------------------------------------------------
  // widgets for defining the sampling of the t-parameter for the curve or function:

  samplingLabel = new Label( String(("samplingLabel")), String(("Sampling:")) );
  samplingLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( samplingLabel );

  tMinLabel = new Label( String(("tMinLabel")), String(("Min:")) );
  tMinLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( tMinLabel );

  tMinEditLabel = new Label( String(("tMinEditLabel")), String(("0")) );
  tMinEditLabel->setJustificationType(Justification::centredLeft);
  tMinEditLabel->setColour(Label::outlineColourId, Colours::black);
  tMinEditLabel->setEditable(true, true);
  addAndMakeVisible( tMinEditLabel );

  tMaxLabel = new Label( String(("tMaxLabel")), String(("Max:")) );
  tMaxLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( tMaxLabel );

  tMaxEditLabel = new Label( String(("tMaxEditLabel1")), String(("10")) );
  tMaxEditLabel->setJustificationType(Justification::centredLeft);
  tMaxEditLabel->setColour(Label::outlineColourId, Colours::black);
  tMaxEditLabel->setEditable(true, true);
  addAndMakeVisible( tMaxEditLabel );

  numSamplesLabel = new Label( String(("NumSamplesLabel")), String(("Samples:")) );
  numSamplesLabel->setJustificationType(Justification::centredLeft);
  addAndMakeVisible( numSamplesLabel );

  numSamplesSlider = new Slider( String(("numSamplesSlider")) );
  numSamplesSlider->setRange(16.0, 1024.0, 1.0);
  numSamplesSlider->setValue(100.0);
  numSamplesSlider->setSliderStyle(Slider::LinearBar);
  addAndMakeVisible(numSamplesSlider);

  tExponentialSpacingButton = new RButton(String("exp-spacing"));
  tExponentialSpacingButton->setClickingTogglesState(true);
  tExponentialSpacingButton->setToggleState(false, true);
  addAndMakeVisible(tExponentialSpacingButton);

  //-----------------------------------------------------------------------------------------------
  // widgets for the quick-calculator:

  quickCalcLabel = new Label( String(("quickCalcLabel")), String(("Calc: t=")) );
  quickCalcLabel->setJustificationType(Justification::centredLeft);
  quickCalcLabel->setEditable(false, false);
  addAndMakeVisible( quickCalcLabel );

  quickCalcInputLabel = new Label( String(("quickCalcInputLabel")), String(("1")) );
  quickCalcInputLabel->setJustificationType(Justification::centredLeft);
  quickCalcInputLabel->setColour(Label::outlineColourId, Colours::black);
  quickCalcInputLabel->setEditable(true, true);
  addAndMakeVisible( quickCalcInputLabel );

  quickCalcOutputLabelX = new Label( String(("quickCalcOutputLabel")), String(("x=1")) );
  quickCalcOutputLabelX->setJustificationType(Justification::centredLeft);
  quickCalcOutputLabelX->setColour(Label::outlineColourId, Colours::black);
  quickCalcOutputLabelX->setEditable(true, true);
  addAndMakeVisible( quickCalcOutputLabelX );

  quickCalcOutputLabelY = new Label( String(("quickCalcOutputLabel")), String(("y=1")) );
  quickCalcOutputLabelY->setJustificationType(Justification::centredLeft);
  quickCalcOutputLabelY->setColour(Label::outlineColourId, Colours::black);
  quickCalcOutputLabelY->setEditable(true, true);
  addAndMakeVisible( quickCalcOutputLabelY );

  //-----------------------------------------------------------------------------------------------
  // widgets for the data import:

  fileLabel = new Label( String(("fileLabel")), String(("Data-File:")) );
  fileLabel->setJustificationType(Justification::centredLeft);
  addChildComponent(fileLabel);

  fileNameLabel = new Label( String(("fileLabel")), String(("none loaded")) );
  fileNameLabel->setJustificationType(Justification::centredLeft);
  fileNameLabel->setColour(Label::outlineColourId, Colours::black);
  addChildComponent(fileNameLabel);

  loadButton = new RClickButton(String(("Load")));
  addChildComponent(loadButton);
}

RSPlotDataSetup::~RSPlotDataSetup()
{
  deleteAllChildren();
}

void RSPlotDataSetup::setMinMaxA(double newMin, double newMax)
{
  if( newMin < newMax )
  {
    aMinEditLabel->setText(String(newMin), NotificationType::dontSendNotification);
    aMaxEditLabel->setText(String(newMax), NotificationType::dontSendNotification);
    aSlider->setRange(newMin, newMax);
  }
}

void RSPlotDataSetup::setMinMaxB(double newMin, double newMax)
{
  if( newMin < newMax )
  {
    bMinEditLabel->setText(String(newMin), NotificationType::dontSendNotification);
    bMaxEditLabel->setText(String(newMax), NotificationType::dontSendNotification);
    bSlider->setRange(newMin, newMax);
  }
}

void RSPlotDataSetup::setMinMaxC(double newMin, double newMax)
{
  if( newMin < newMax )
  {
    cMinEditLabel->setText(String(newMin), NotificationType::dontSendNotification);
    cMaxEditLabel->setText(String(newMax), NotificationType::dontSendNotification);
    cSlider->setRange(newMin, newMax);
  }
}

void RSPlotDataSetup::setMinMaxD(double newMin, double newMax)
{
  if( newMin < newMax )
  {
    dMinEditLabel->setText(String(newMin), NotificationType::dontSendNotification);
    dMaxEditLabel->setText(String(newMax), NotificationType::dontSendNotification);
    dSlider->setRange(newMin, newMax);
  }
}

void RSPlotDataSetup::setMinMaxT(double newMin, double newMax)
{
  if( newMin < newMax )
  {
    tMinEditLabel->setText(String(newMin), NotificationType::dontSendNotification);
    tMaxEditLabel->setText(String(newMax), NotificationType::dontSendNotification);
  }
}

//-------------------------------------------------------------------------------------------------
// appearance:

void RSPlotDataSetup::paint(Graphics& g)
{
  g.fillAll(Colours::black);
}

void RSPlotDataSetup::resized()
{
  Component::resized();

  int xL, xR, yT, yB; // left, right, top and bottom for the child-component
  int w, h;           // width and height for the child component

  xL = 0;
  xR = getWidth();
  w  = (xR-xL)/3;
  yT = 4;
  yB = getHeight();
  h  = yB-yT;

  dataSourceLabel->setBounds(xL+4, yT, 56, 20);
  dataSourceComboBox->setBounds(xL+56+4, yT, w-56-8, 20);
  xL += w;
  typeOfDataLabel->setBounds(xL+4, yT, 56, 20);
  typeOfDataComboBox->setBounds(xL+56+4, yT, w-56-8, 20);
  xL += w;
  numCurvesLabel->setBounds(xL+4, yT, 56, 20);
  numCurvesSlider->setBounds(xL+56+4, yT, w-56-8, 20);

  xL  = 0;
  w   = getWidth();
  yT += 24;

  xCurveLabel->setBounds(xL+4, yT, 56, 20);
  xCurveEditLabel->setBounds(xL+56+4, yT, w-56-8, 20);
  yT += 24;
  yCurveLabel->setBounds(xL+4, yT, 56, 20);
  yCurveEditLabel->setBounds(xL+56+4, yT, w-56-8, 20);


  w   = getWidth()/4;
  yT += 24;

  aLabel->setBounds(xL+4, yT, 20, 20);
  aMinEditLabel->setBounds(xL+20+4, yT, 32, 20);
  aMaxEditLabel->setBounds(xL+w-32-8, yT, 32, 20);
  aSlider->setBounds(aMinEditLabel->getRight(), yT, 
                     aMaxEditLabel->getX()-aMinEditLabel->getRight(), 20);
  yT += 24;
  bLabel->setBounds(xL+4, yT, 20, 20);
  bMinEditLabel->setBounds(xL+20+4, yT, 32, 20);
  bMaxEditLabel->setBounds(xL+w-32-8, yT, 32, 20);
  bSlider->setBounds(bMinEditLabel->getRight(), yT, 
                     bMaxEditLabel->getX()-bMinEditLabel->getRight(), 20);

  yT  = aLabel->getY();
  xL += w;
  cLabel->setBounds(xL+4, yT, 20, 20);
  cMinEditLabel->setBounds(xL+20+4, yT, 32, 20);
  cMaxEditLabel->setBounds(xL+w-32-8, yT, 32, 20);
  cSlider->setBounds(cMinEditLabel->getRight(), yT, 
                     cMaxEditLabel->getX()-cMinEditLabel->getRight(), 20);
  yT += 24;
  dLabel->setBounds(xL+4, yT, 20, 20);
  dMinEditLabel->setBounds(xL+20+4, yT, 32, 20);
  dMaxEditLabel->setBounds(xL+w-32-8, yT, 32, 20);
  dSlider->setBounds(dMinEditLabel->getRight(), yT, 
                     dMaxEditLabel->getX()-dMinEditLabel->getRight(), 20);

  yT  = aLabel->getY();
  xL += w;

  samplingLabel->setBounds(xL, yT, w, 20);
  tExponentialSpacingButton->setBounds(xL+w-88, yT, 88-8, 20);
  yT += 24;

  tMinEditLabel->setBounds(xL+4, yT, 32, 20);
  tMaxEditLabel->setBounds(xL+w-32-8, yT, 32, 20);
  numSamplesSlider->setBounds(tMinEditLabel->getRight(), yT, 
                              tMaxEditLabel->getX()-tMinEditLabel->getRight(), 20);

  yT  = aLabel->getY();
  xL += w;
  quickCalcLabel->setBounds(xL+4, yT, 56, 20);
  quickCalcInputLabel->setBounds(quickCalcLabel->getRight(), yT, w-56-8, 20);
  yT += 24;
  quickCalcOutputLabelX->setBounds(xL+4, yT, w/2-4, 20);
  xL += w/2-4;
  quickCalcOutputLabelY->setBounds(xL+4, yT, w/2-4, 20);

  yT = numCurvesLabel->getBottom();
  yT += 32;
  fileLabel->setBounds(xL, yT, 80, 20);
  loadButton->setBounds(xR-60, yT-2, 56, 24);
  yT += 24;
  fileNameLabel->setBounds(xL+4, yT, w-4, 20);

}

