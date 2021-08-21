#include "RSPlotContentComponent.h"

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RSPlotContentComponent::RSPlotContentComponent(
  const String &newEditorName) : Component(newEditorName)
{
  setSize(800, 600-24);
  //dataFileName    = String::empty;  // old
  dataFileName = String();  // new

  // initialize the data-arrays:
  xFamilyPointer = new double*[maxNumCurves];
  yFamilyPointer = new double*[maxNumCurves];
  for(int c=0; c<maxNumCurves; c++)
  {
    for(int i=0; i<maxNumValues; i++)
    {
      xValues[c][i] = 0.0;
      yValues[c][i] = 0.0;
    }
    xFamilyPointer[c] = &(xValues[c][0]);
    yFamilyPointer[c] = &(yValues[c][0]);
  }

  // the components for the actual data-visualization:
  /*
  surfacePlot = new SurfacePlot();
  surfacePlot->setRange(-5.0, 5.0, -0.5, 1.5, -5.0, 5.0);
  surfacePlot->useStandardProjection(CoordinateSystem3D::ISOMETRIC);
  addChildComponent(surfacePlot);

  zoomer3D = new CoordinateSystem3DZoomer();
  zoomer3D->setZoomerSize(16);
  zoomer3D->setCoordinateSystem3D(surfacePlot);
  addChildComponent(zoomer3D);
  */

  curveFamilyPlot = new rsDataPlot();
  curveFamilyPlot->setCurrentRange(-5.0, 5.0, -5.0, 5.0);
  curveFamilyPlot->setCurveFamilyValues(100, 1, xFamilyPointer, yFamilyPointer);
  addAndMakeVisible(curveFamilyPlot);

  zoomer2D = new rsPlotZoomer();
  zoomer2D->setZoomerSize(16);
  zoomer2D->setRelativeMargins(0.0, 0.0, 0.0, 0.0);
  zoomer2D->setCoordinateSystem(curveFamilyPlot);
  zoomer2D->setVerticalMouseWheelMode(rsPlotZoomer::zoomViaVerticalMouseWheel);
  addAndMakeVisible(zoomer2D);

  // select the surface-plot component as current plotter by default:
  currentPlotterComponent = CURVE_FAMILY_PLOT;

  // the components for the various tabs:
  axesSetup = new RSPlotAxesSetup();
  dataSetup = new RSPlotDataSetup();

  // the tab-selector:
  setupTabber = new TabbedComponent(TabbedButtonBar::TabsAtLeft);
  setupTabber->setOutline(1);
  setupTabber->addTab(String(("Data")), Colours::white, dataSetup, true);
  setupTabber->addTab(String(("Axes")), Colours::white, axesSetup, true);
  //setupTabber->addTab(String(("File")), Colours::white, fileSetup, true);
  addAndMakeVisible(setupTabber);

  // the component for the file management:
  fileSetup = new RSPlotFileSetup();
  addAndMakeVisible(fileSetup);

  // register ourselves as listener to the widgets which reside in the various
  // tab-components (the widgets are public there):
  axesSetup->xAnnotationEditLabel->addListener(this);
  axesSetup->xMinEditLabel->addListener(this);
  axesSetup->xMaxEditLabel->addListener(this);
  axesSetup->xLogScaleButton->addRButtonListener(this);
  axesSetup->xPosComboBox->addListener(this);
  axesSetup->yAnnotationEditLabel->addListener(this);
  axesSetup->yMinEditLabel->addListener(this);
  axesSetup->yMaxEditLabel->addListener(this);
  axesSetup->yLogScaleButton->addRButtonListener(this);
  axesSetup->yPosComboBox->addListener(this);
  /*
  axesSetup->zAnnotationEditLabel->addListener(this);
  axesSetup->minZEditLabel->addListener(this);
  axesSetup->maxZEditLabel->addListener(this);
  axesSetup->logScaleZButton->addButtonListener(this);
  axesSetup->zAxisComboBox->addListener(this);
  */
  axesSetup->horizontalCoarseGridButton->addRButtonListener(this);
  axesSetup->horizontalCoarseGridIntervalLabel->addListener(this);
  axesSetup->horizontalFineGridButton->addRButtonListener(this);
  axesSetup->horizontalFineGridIntervalLabel->addListener(this);
  axesSetup->verticalCoarseGridButton->addRButtonListener(this);
  axesSetup->verticalCoarseGridIntervalLabel->addListener(this);
  axesSetup->verticalFineGridButton->addRButtonListener(this);
  axesSetup->verticalFineGridIntervalLabel->addListener(this);
  axesSetup->radialCoarseGridButton->addRButtonListener(this);
  axesSetup->radialCoarseGridIntervalLabel->addListener(this);
  axesSetup->radialFineGridButton->addRButtonListener(this);
  axesSetup->radialFineGridIntervalLabel->addListener(this);
  axesSetup->angularCoarseGridButton->addRButtonListener(this);
  axesSetup->angularCoarseGridIntervalLabel->addListener(this);
  axesSetup->angularFineGridButton->addRButtonListener(this);
  axesSetup->angularFineGridIntervalLabel->addListener(this);
  /*
  axesSetup->rotationSliderX->addListener(this);
  axesSetup->rotationSliderY->addListener(this);
  axesSetup->rotationSliderZ->addListener(this);
  */

  dataSetup->dataSourceComboBox->addListener(this);
  dataSetup->typeOfDataComboBox->addListener(this);
  dataSetup->numCurvesSlider->addListener(this);
  dataSetup->xCurveEditLabel->registerTextEntryFieldObserver(this);
  dataSetup->yCurveEditLabel->registerTextEntryFieldObserver(this);
  dataSetup->aMinEditLabel->addListener(this);
  dataSetup->aSlider->addListener(this);
  dataSetup->aMaxEditLabel->addListener(this);
  dataSetup->bMinEditLabel->addListener(this);
  dataSetup->bSlider->addListener(this);
  dataSetup->bMaxEditLabel->addListener(this);
  dataSetup->cMinEditLabel->addListener(this);
  dataSetup->cSlider->addListener(this);
  dataSetup->cMaxEditLabel->addListener(this);
  dataSetup->dMinEditLabel->addListener(this);
  dataSetup->dSlider->addListener(this);
  dataSetup->dMaxEditLabel->addListener(this);
  dataSetup->tMinEditLabel->registerTextEntryFieldObserver(this);
  dataSetup->tMaxEditLabel->registerTextEntryFieldObserver(this);
  dataSetup->numSamplesSlider->addListener(this);
  dataSetup->tExponentialSpacingButton->addRButtonListener(this);
  dataSetup->loadButton->addRButtonListener(this);
  dataSetup->quickCalcInputLabel->addListener(this);

  fileSetup->loadButton->addRButtonListener(this);
  fileSetup->saveButton->addRButtonListener(this);
  fileSetup->plusButton->addRButtonListener(this);
  fileSetup->minusButton->addRButtonListener(this);
  fileSetup->pngExportButton->addRButtonListener(this);
  fileSetup->svgExportButton->addRButtonListener(this);

  // set up the range:
  axesSetup->xMinEditLabel->setText(String(("-2.5")), NotificationType::sendNotification);
  axesSetup->xMaxEditLabel->setText(String(("+2.5")), NotificationType::sendNotification);
  axesSetup->yMinEditLabel->setText(String(("-2.5")), NotificationType::sendNotification);
  axesSetup->yMaxEditLabel->setText(String(("+2.5")), NotificationType::sendNotification);

}

RSPlotContentComponent::~RSPlotContentComponent()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

//void RSPlotContentComponent::buttonClicked(Button *buttonThatWasClicked)
//{
//
//}

void RSPlotContentComponent::rButtonClicked(RButton *buttonThatWasClicked)
{
  int w, h;

  if( buttonThatWasClicked == fileSetup->loadButton )
  {
    StateFileManager::openLoadingDialog();
    fileSetup->fileNameLabel->setText(StateFileManager::getActiveFile().getFileName());
  }
  else if( buttonThatWasClicked == fileSetup->saveButton )
  {
    StateFileManager::openSavingDialog();
    fileSetup->fileNameLabel->setText(StateFileManager::getActiveFile().getFileName());
  }
  else if(buttonThatWasClicked == fileSetup->plusButton )
  {
    StateFileManager::loadNextFile();
    //StateFileManager::loadStateFromXmlFile(StateFileManager::getActiveFile());
    fileSetup->fileNameLabel->setText(StateFileManager::getActiveFile().getFileName());
  }
  else if(buttonThatWasClicked == fileSetup->minusButton )
  {
    StateFileManager::loadPreviousFile();
    //StateFileManager::loadStateFromXmlFile(StateFileManager::getActiveFile());
    fileSetup->fileNameLabel->setText(StateFileManager::getActiveFile().getFileName());
  }
  else if(buttonThatWasClicked == fileSetup->pngExportButton )
  {
    w = fileSetup->imageWidthEditLabel->getText().getIntValue();
    h = fileSetup->imageHeightEditLabel->getText().getIntValue();
    if( currentPlotterComponent == CURVE_FAMILY_PLOT )
      curveFamilyPlot->openExportDialog(w, h, String(("png")), File()); // new
      //curveFamilyPlot->openExportDialog(w, h, String(("png")), File::nonexistent); // old
  }
  else if(buttonThatWasClicked == fileSetup->svgExportButton )
  {
    w = fileSetup->imageWidthEditLabel->getText().getIntValue();
    h = fileSetup->imageHeightEditLabel->getText().getIntValue();
    if( currentPlotterComponent == CURVE_FAMILY_PLOT )
      curveFamilyPlot->openExportDialog(w, h, String(("svg")), File());  // new
      //curveFamilyPlot->openExportDialog(w, h, String(("svg")), File::nonexistent); // old
  }

  else if( buttonThatWasClicked == dataSetup->loadButton )
    openDataLoadingDialog();

  switch( currentPlotterComponent )
  {
  case SURFACE_PLOT:
    {


    } // end of case SURFACE_PLOT:
    break;

  case CURVE_FAMILY_PLOT:
    {
      if( buttonThatWasClicked == axesSetup->xLogScaleButton )
        curveFamilyPlot->useLogarithmicScaleX(axesSetup->xLogScaleButton->getToggleState());
      else if( buttonThatWasClicked == axesSetup->yLogScaleButton )
        curveFamilyPlot->useLogarithmicScaleY(axesSetup->yLogScaleButton->getToggleState());
      else if( buttonThatWasClicked == axesSetup->horizontalCoarseGridButton )
        curveFamilyPlot->setHorizontalCoarseGridVisible(
          axesSetup->horizontalCoarseGridButton->getToggleState() );
      else if( buttonThatWasClicked == axesSetup->horizontalFineGridButton )
        curveFamilyPlot->setHorizontalFineGridVisible(
          axesSetup->horizontalFineGridButton->getToggleState() );
      else if( buttonThatWasClicked == axesSetup->verticalCoarseGridButton )
        curveFamilyPlot->setVerticalCoarseGridVisible(
          axesSetup->verticalCoarseGridButton->getToggleState() );
      else if( buttonThatWasClicked == axesSetup->verticalFineGridButton )
        curveFamilyPlot->setVerticalFineGridVisible(
          axesSetup->verticalFineGridButton->getToggleState() );
      else if( buttonThatWasClicked == axesSetup->radialCoarseGridButton )
        curveFamilyPlot->setRadialCoarseGridVisible(
          axesSetup->radialCoarseGridButton->getToggleState() );
      else if( buttonThatWasClicked == axesSetup->radialFineGridButton )
        curveFamilyPlot->setRadialFineGridVisible(
          axesSetup->radialFineGridButton->getToggleState() );
      else if( buttonThatWasClicked == axesSetup->angularCoarseGridButton )
        curveFamilyPlot->setAngularCoarseGridVisible(
          axesSetup->angularCoarseGridButton->getToggleState() );
      else if( buttonThatWasClicked == axesSetup->angularFineGridButton )
        curveFamilyPlot->setAngularFineGridVisible(
          axesSetup->angularFineGridButton->getToggleState() );

      else if( buttonThatWasClicked == dataSetup->tExponentialSpacingButton )
        calculateData();


    } // end of case CURVE_FAMILY_PLOT:
    break;

  } // end of switch( currentPlotterComponent )
}

void RSPlotContentComponent::comboBoxChanged(ComboBox *comboBoxThatHasChanged)
{
  int id;

  if( comboBoxThatHasChanged == dataSetup->typeOfDataComboBox )
  {
    if( dataSetup->typeOfDataComboBox->getText() == String(("Curve Family")) ||
        dataSetup->typeOfDataComboBox->getText() == String(("Scatter Plot"))    )
    {
      // show the 2D coordinate-system:
      currentPlotterComponent = CURVE_FAMILY_PLOT;
      curveFamilyPlot->setVisible(true);
      zoomer2D->setVisible(true);
      //surfacePlot->setVisible(false);
      //zoomer3D->setVisible(false);
    }
    else if( dataSetup->typeOfDataComboBox->getText() == String(("Surface Plot")) )
    {
      // show the 3d coordinate-system:
      currentPlotterComponent = SURFACE_PLOT;
      curveFamilyPlot->setVisible(false);
      zoomer2D->setVisible(false);
      //surfacePlot->setVisible(true);
      //zoomer3D->setVisible(true);
    }
  }

  switch( currentPlotterComponent )
  {
  case SURFACE_PLOT:
    {


    } // end of case SURFACE_PLOT:
    break;

  case CURVE_FAMILY_PLOT:
    {
      if( comboBoxThatHasChanged == axesSetup->xDigitsComboBox )
      {
        id = axesSetup->xPosComboBox->getSelectedId();
        if( id == 1 )
          curveFamilyPlot->setStringConversionForAxisX(&valueToString0);
        else if( id == 2 )
          curveFamilyPlot->setStringConversionForAxisX(&valueToString1);
        else if( id == 3 )
          curveFamilyPlot->setStringConversionForAxisX(&valueToString2);
        curveFamilyPlot->updatePlotImage(true);
      }
      else if( comboBoxThatHasChanged == axesSetup->xPosComboBox )
      {
        id = axesSetup->xPosComboBox->getSelectedId();
        if( id == 1 )
          curveFamilyPlot->setAxisPositionX(rsPlotSettings::INVISIBLE);
        else if( id == 2 )
          curveFamilyPlot->setAxisPositionX(rsPlotSettings::ZERO);
        else if( id == 3 )
          curveFamilyPlot->setAxisPositionX(rsPlotSettings::TOP);
        else if( id == 4 )
          curveFamilyPlot->setAxisPositionX(rsPlotSettings::BOTTOM);
      }



      else if( comboBoxThatHasChanged == axesSetup->yPosComboBox )
      {
        id = axesSetup->yPosComboBox->getSelectedId();
        if( id == 1 )
          curveFamilyPlot->setAxisPositionY(rsPlotSettings::INVISIBLE);
        else if( id == 2 )
          curveFamilyPlot->setAxisPositionY(rsPlotSettings::ZERO);
        else if( id == 3 )
          curveFamilyPlot->setAxisPositionY(rsPlotSettings::LEFT);
        else if( id == 4 )
          curveFamilyPlot->setAxisPositionY(rsPlotSettings::RIGHT);
      }

      else if( comboBoxThatHasChanged == dataSetup->dataSourceComboBox )
      {
        if( dataSetup->dataSourceComboBox->getText() == String(("File")) )
        {
          dataSetup->xCurveLabel->setVisible(false);
          dataSetup->xCurveEditLabel->setVisible(false);
          dataSetup->yCurveLabel->setVisible(false);
          dataSetup->yCurveEditLabel->setVisible(false);
          dataSetup->aLabel->setVisible(false);
          dataSetup->aMinEditLabel->setVisible(false);
          dataSetup->aSlider->setVisible(false);
          dataSetup->aMaxEditLabel->setVisible(false);
          dataSetup->bLabel->setVisible(false);
          dataSetup->bMinEditLabel->setVisible(false);
          dataSetup->bSlider->setVisible(false);
          dataSetup->bMaxEditLabel->setVisible(false);
          dataSetup->cLabel->setVisible(false);
          dataSetup->cMinEditLabel->setVisible(false);
          dataSetup->cSlider->setVisible(false);
          dataSetup->cMaxEditLabel->setVisible(false);
          dataSetup->dLabel->setVisible(false);
          dataSetup->dMinEditLabel->setVisible(false);
          dataSetup->dSlider->setVisible(false);
          dataSetup->dMaxEditLabel->setVisible(false);
          dataSetup->samplingLabel->setVisible(false);
          dataSetup->tMinLabel->setVisible(false);
          dataSetup->tMinEditLabel->setVisible(false);
          dataSetup->tMaxLabel->setVisible(false);
          dataSetup->tMaxEditLabel->setVisible(false);
          dataSetup->numSamplesLabel->setVisible(false);
          dataSetup->numSamplesSlider->setVisible(false);
          dataSetup->tExponentialSpacingButton->setVisible(false);

          dataSetup->fileLabel->setVisible(true);
          dataSetup->fileNameLabel->setVisible(true);
          dataSetup->loadButton->setVisible(true);
        }
        else
        {
          dataSetup->xCurveLabel->setVisible(true);
          dataSetup->xCurveEditLabel->setVisible(true);
          dataSetup->yCurveLabel->setVisible(true);
          dataSetup->yCurveEditLabel->setVisible(true);
          dataSetup->aLabel->setVisible(true);
          dataSetup->aMinEditLabel->setVisible(true);
          dataSetup->aSlider->setVisible(true);
          dataSetup->aMaxEditLabel->setVisible(true);
          dataSetup->bLabel->setVisible(true);
          dataSetup->bMinEditLabel->setVisible(true);
          dataSetup->bSlider->setVisible(true);
          dataSetup->bMaxEditLabel->setVisible(true);
          dataSetup->cLabel->setVisible(true);
          dataSetup->cMinEditLabel->setVisible(true);
          dataSetup->cSlider->setVisible(true);
          dataSetup->cMaxEditLabel->setVisible(true);
          dataSetup->dLabel->setVisible(true);
          dataSetup->dMinEditLabel->setVisible(true);
          dataSetup->dSlider->setVisible(true);
          dataSetup->dMaxEditLabel->setVisible(true);
          dataSetup->samplingLabel->setVisible(true);
          dataSetup->tMinLabel->setVisible(true);
          dataSetup->tMinEditLabel->setVisible(true);
          dataSetup->tMaxLabel->setVisible(true);
          dataSetup->tMaxEditLabel->setVisible(true);
          dataSetup->numSamplesLabel->setVisible(true);
          dataSetup->numSamplesSlider->setVisible(true);
          dataSetup->tExponentialSpacingButton->setVisible(true);

          dataSetup->fileLabel->setVisible(false);
          dataSetup->fileNameLabel->setVisible(false);
          dataSetup->loadButton->setVisible(false);
        }

      }

    } // end of case CURVE_FAMILY_PLOT:
    break;

  } // end of switch( currentPlotterComponent )
}

void RSPlotContentComponent::textChanged(RTextEntryField *field)
{
  if(field == dataSetup->xCurveEditLabel || field == dataSetup->yCurveEditLabel ||
    field  == dataSetup->tMinEditLabel   || field == dataSetup->tMaxEditLabel)
  {
    calculateData();
  }
}

void RSPlotContentComponent::labelTextChanged(Label *labelThatHasChanged)
{
  double minValue, maxValue;
  double margin = 0.0001;
  if( labelThatHasChanged == dataSetup->aMinEditLabel )
  {
    minValue = dataSetup->aMinEditLabel->getText().getDoubleValue();
    maxValue = dataSetup->aSlider->getMaximum();
    if( minValue < maxValue-margin )
      dataSetup->aSlider->setRange(minValue, maxValue);
    dataSetup->aMinEditLabel->setText(String(dataSetup->aSlider->getMinimum()), NotificationType::dontSendNotification);
    dataSetup->aSlider->repaint();
  }
  else if( labelThatHasChanged == dataSetup->aMaxEditLabel )
  {
    minValue = dataSetup->aSlider->getMinimum();
    maxValue = dataSetup->aMaxEditLabel->getText().getDoubleValue();
    if( minValue < maxValue-margin )
      dataSetup->aSlider->setRange(minValue, maxValue);
    dataSetup->aMaxEditLabel->setText(String(dataSetup->aSlider->getMaximum()), NotificationType::dontSendNotification);
    dataSetup->aSlider->repaint();
  }
  else if( labelThatHasChanged == dataSetup->bMinEditLabel )
  {
    minValue = dataSetup->bMinEditLabel->getText().getDoubleValue();
    maxValue = dataSetup->bSlider->getMaximum();
    if( minValue < maxValue-margin )
      dataSetup->bSlider->setRange(minValue, maxValue);
    dataSetup->bMinEditLabel->setText(String(dataSetup->bSlider->getMinimum()), NotificationType::dontSendNotification);
    dataSetup->bSlider->repaint();

  }
  else if( labelThatHasChanged == dataSetup->bMaxEditLabel )
  {
    minValue = dataSetup->bSlider->getMinimum();
    maxValue = dataSetup->bMaxEditLabel->getText().getDoubleValue();
    if( minValue < maxValue-margin )
      dataSetup->bSlider->setRange(minValue, maxValue);
    dataSetup->bMaxEditLabel->setText(String(dataSetup->bSlider->getMaximum()), NotificationType::dontSendNotification);
    dataSetup->bSlider->repaint();

  }
  else if( labelThatHasChanged == dataSetup->cMinEditLabel )
  {
    minValue = dataSetup->cMinEditLabel->getText().getDoubleValue();
    maxValue = dataSetup->cSlider->getMaximum();
    if( minValue < maxValue-margin )
      dataSetup->cSlider->setRange(minValue, maxValue);
    dataSetup->cMinEditLabel->setText(String(dataSetup->cSlider->getMinimum()), NotificationType::dontSendNotification);
    dataSetup->cSlider->repaint();
  }
  else if( labelThatHasChanged == dataSetup->cMaxEditLabel )
  {
    minValue = dataSetup->cSlider->getMinimum();
    maxValue = dataSetup->cMaxEditLabel->getText().getDoubleValue();
    if( minValue < maxValue-margin )
      dataSetup->cSlider->setRange(minValue, maxValue);
    dataSetup->cMaxEditLabel->setText(String(dataSetup->cSlider->getMaximum()), NotificationType::dontSendNotification);
    dataSetup->cSlider->repaint();
  }
  else if( labelThatHasChanged == dataSetup->dMinEditLabel )
  {
    minValue = dataSetup->dMinEditLabel->getText().getDoubleValue();
    maxValue = dataSetup->dSlider->getMaximum();
    if( minValue < maxValue-margin )
      dataSetup->dSlider->setRange(minValue, maxValue);
    dataSetup->dMinEditLabel->setText(String(dataSetup->dSlider->getMinimum()), NotificationType::dontSendNotification);
    dataSetup->dSlider->repaint();
  }
  else if( labelThatHasChanged == dataSetup->dMaxEditLabel )
  {
    minValue = dataSetup->dSlider->getMinimum();
    maxValue = dataSetup->dMaxEditLabel->getText().getDoubleValue();
    if( minValue < maxValue-margin )
      dataSetup->dSlider->setRange(minValue, maxValue);
    dataSetup->dMaxEditLabel->setText(String(dataSetup->dSlider->getMaximum()), NotificationType::dontSendNotification);
    dataSetup->dSlider->repaint();
  }
  else if( labelThatHasChanged == dataSetup->quickCalcInputLabel )
  {
    updateCalculatorResult();
  }

  switch( currentPlotterComponent )
  {
  case SURFACE_PLOT:
    {
      /*
      if( labelThatHasChanged == axesSetup->xAnnotationEditLabel )
      {
        surfacePlot->setAxisLabel(axesSetup->xAnnotationEditLabel->getText(),
          SurfacePlot::X_AXIS);
      }
      else if( labelThatHasChanged == axesSetup->xMinEditLabel )
      {
        surfacePlot->setMinimum(
          axesSetup->xMinEditLabel->getText().getDoubleValue(), 
          SurfacePlot::X_AXIS);
      }
      else if( labelThatHasChanged == axesSetup->xMaxEditLabel )
      {
        surfacePlot->setMaximum(
          axesSetup->xMaxEditLabel->getText().getDoubleValue(), 
          SurfacePlot::X_AXIS);
      }

      if( labelThatHasChanged == axesSetup->yAnnotationEditLabel )
      {
        surfacePlot->setAxisLabel(axesSetup->yAnnotationEditLabel->getText(),
          SurfacePlot::Y_AXIS);
      }
      else if( labelThatHasChanged == axesSetup->yMinEditLabel )
      {
        surfacePlot->setMinimum(
          axesSetup->yMinEditLabel->getText().getDoubleValue(), 
          SurfacePlot::Y_AXIS);
      }
      else if( labelThatHasChanged == axesSetup->yMaxEditLabel )
      {
        surfacePlot->setMaximum(
          axesSetup->yMaxEditLabel->getText().getDoubleValue(), 
          SurfacePlot::Y_AXIS);
      }
      */

      /*
      else if( labelThatHasChanged == axesSetup->zAnnotationEditLabel )
      {
        surfacePlot->setAxisLabel(axesSetup->zAnnotationEditLabel->getText(),
          SurfacePlot::Z_AXIS);
      }
      else if( labelThatHasChanged == axesSetup->minZEditLabel )
      {
        surfacePlot->setMinimum(
          axesSetup->minZEditLabel->getText().getDoubleValue(), 
          SurfacePlot::Z_AXIS);
      }
      else if( labelThatHasChanged == axesSetup->maxZEditLabel )
      {
        surfacePlot->setMaximum(
          axesSetup->maxZEditLabel->getText().getDoubleValue(), 
          SurfacePlot::Z_AXIS);
      }
      */

    } // end of case SURFACE_PLOT:
    break;

  case CURVE_FAMILY_PLOT:
    {
      if( labelThatHasChanged == axesSetup->xAnnotationEditLabel )
        curveFamilyPlot->setAxisLabelX(axesSetup->xAnnotationEditLabel->getText());
      else if( labelThatHasChanged == axesSetup->xMinEditLabel )
      {
        curveFamilyPlot->setMaximumRangeMinX(axesSetup->xMinEditLabel->getText().getDoubleValue());
        curveFamilyPlot->setCurrentRangeMinX(axesSetup->xMinEditLabel->getText().getDoubleValue());
        axesSetup->xMinEditLabel->setText(String(curveFamilyPlot->getCurrentRangeMinX()), NotificationType::dontSendNotification);
        curveFamilyPlot->updatePlotImage(true);
      }
      else if( labelThatHasChanged == axesSetup->xMaxEditLabel )
      {
        curveFamilyPlot->setMaximumRangeMaxX(axesSetup->xMaxEditLabel->getText().getDoubleValue());
        curveFamilyPlot->setCurrentRangeMaxX(axesSetup->xMaxEditLabel->getText().getDoubleValue());
        axesSetup->xMaxEditLabel->setText(String(curveFamilyPlot->getCurrentRangeMaxX()), NotificationType::dontSendNotification);
        curveFamilyPlot->updatePlotImage(true);
      }
      else if( labelThatHasChanged == axesSetup->yAnnotationEditLabel )
        curveFamilyPlot->setAxisLabelY(axesSetup->yAnnotationEditLabel->getText());
      else if( labelThatHasChanged == axesSetup->yMinEditLabel )
      {
        curveFamilyPlot->setMaximumRangeMinY(axesSetup->yMinEditLabel->getText().getDoubleValue());
        curveFamilyPlot->setCurrentRangeMinY(axesSetup->yMinEditLabel->getText().getDoubleValue());
        axesSetup->yMinEditLabel->setText(String(curveFamilyPlot->getCurrentRangeMinY()), NotificationType::dontSendNotification);
      }
      else if( labelThatHasChanged == axesSetup->yMaxEditLabel )
      {
        curveFamilyPlot->setMaximumRangeMaxY(axesSetup->yMaxEditLabel->getText().getDoubleValue());
        curveFamilyPlot->setCurrentRangeMaxY(axesSetup->yMaxEditLabel->getText().getDoubleValue());
        axesSetup->yMaxEditLabel->setText(String(curveFamilyPlot->getCurrentRangeMaxY()), NotificationType::dontSendNotification);
      }
      else if( labelThatHasChanged == axesSetup->horizontalCoarseGridIntervalLabel )
        curveFamilyPlot->setHorizontalCoarseGridInterval(
        axesSetup->horizontalCoarseGridIntervalLabel->getText().getDoubleValue());
      else if( labelThatHasChanged == axesSetup->horizontalFineGridIntervalLabel )
        curveFamilyPlot->setHorizontalFineGridInterval(
        axesSetup->horizontalFineGridIntervalLabel->getText().getDoubleValue());
      else if( labelThatHasChanged == axesSetup->verticalCoarseGridIntervalLabel )
        curveFamilyPlot->setVerticalCoarseGridInterval(
        axesSetup->verticalCoarseGridIntervalLabel->getText().getDoubleValue());
      else if( labelThatHasChanged == axesSetup->verticalFineGridIntervalLabel )
        curveFamilyPlot->setVerticalFineGridInterval(
        axesSetup->verticalFineGridIntervalLabel->getText().getDoubleValue());
      else if( labelThatHasChanged == axesSetup->radialCoarseGridIntervalLabel )
        curveFamilyPlot->setRadialCoarseGridInterval(
        axesSetup->radialCoarseGridIntervalLabel->getText().getDoubleValue());
      else if( labelThatHasChanged == axesSetup->radialFineGridIntervalLabel )
        curveFamilyPlot->setRadialFineGridInterval(
        axesSetup->radialFineGridIntervalLabel->getText().getDoubleValue());
      else if( labelThatHasChanged == axesSetup->angularCoarseGridIntervalLabel )
        curveFamilyPlot->setAngularCoarseGridInterval(
        axesSetup->angularCoarseGridIntervalLabel->getText().getDoubleValue());
      else if( labelThatHasChanged == axesSetup->angularFineGridIntervalLabel )
        curveFamilyPlot->setAngularFineGridInterval(
        axesSetup->angularFineGridIntervalLabel->getText().getDoubleValue());

    } // end of case CURVE_FAMILY_PLOT:
    break;

  } // end of switch( currentPlotterComponent )

}

void RSPlotContentComponent::sliderValueChanged(Slider *sliderThatHasChanged)
{

  if( sliderThatHasChanged == dataSetup->aSlider           || 
      sliderThatHasChanged == dataSetup->bSlider           || 
      sliderThatHasChanged == dataSetup->cSlider           || 
      sliderThatHasChanged == dataSetup->dSlider           ||
      sliderThatHasChanged == dataSetup->numCurvesSlider   ||
      sliderThatHasChanged == dataSetup->numSamplesSlider      )
  {
    calculateData();
  }

  /*
  switch( currentPlotterComponent )
  {
  case SURFACE_PLOT:
    {


    } // end of case SURFACE_PLOT:
    break;

  case CURVE_FAMILY_PLOT:
    {
      if( sliderThatHasChanged == dataSetup->numSamplesSlider )
        calculateData();
    } // end of case CURVE_FAMILY_PLOT:
    break;

  } // end of switch( currentPlotterComponent )
  */
}

//-------------------------------------------------------------------------------------------------
// state management:

XmlElement* RSPlotContentComponent::getStateAsXml() const
{
  XmlElement* xmlState = new XmlElement(String(("RSPlotState"))); 
    // the XmlElement which stores all the relevant state-information

  xmlState->setAttribute(String(("ExportWidth")),  
    fileSetup->imageWidthEditLabel->getText().getIntValue() );
  xmlState->setAttribute(String(("ExportHeight")), 
    fileSetup->imageHeightEditLabel->getText().getIntValue() );

  xmlState->setAttribute(String(("xAnnotation")), axesSetup->xAnnotationEditLabel->getText());
  xmlState->setAttribute(String(("yAnnotation")), axesSetup->yAnnotationEditLabel->getText());
  //xmlState->setAttribute(String(("zAnnotation")), axesSetup->zAnnotationEditLabel->getText());
  xmlState->setAttribute(String(("xMin")),        
    axesSetup->xMinEditLabel->getText().getDoubleValue());
  xmlState->setAttribute(String(("xMax")),        
    axesSetup->xMaxEditLabel->getText().getDoubleValue());
  xmlState->setAttribute(String(("xLogScale")), axesSetup->xLogScaleButton->getToggleState());
  xmlState->setAttribute(String(("xAxisPosition")), axesSetup->xPosComboBox->getText());
  xmlState->setAttribute(String(("yMin")),        
    axesSetup->yMinEditLabel->getText().getDoubleValue());
  xmlState->setAttribute(String(("yMax")),        
    axesSetup->yMaxEditLabel->getText().getDoubleValue());
  xmlState->setAttribute(String(("yLogScale")), axesSetup->yLogScaleButton->getToggleState());
  xmlState->setAttribute(String(("yPosPosition")), axesSetup->yPosComboBox->getText());
  /*
  xmlState->setAttribute(String(("zMin")),        
    axesSetup->minZEditLabel->getText().getDoubleValue());
  xmlState->setAttribute(String(("zMax")),        
    axesSetup->maxZEditLabel->getText().getDoubleValue());
  xmlState->setAttribute(String(("zLogScale")), axesSetup->logScaleZButton->getToggleState());
  */
  xmlState->setAttribute(String(("HorizontalCoarseGridIsVisible")),
    curveFamilyPlot->isHorizontalCoarseGridVisible());
  xmlState->setAttribute(String(("HorizontalCoarseGridInterval")),
    curveFamilyPlot->getHorizontalCoarseGridInterval());
  xmlState->setAttribute(String(("HorizontalFineGridIsVisible")),      
    curveFamilyPlot->isHorizontalFineGridVisible());
  xmlState->setAttribute(String(("HorizontalFineGridInterval")),  
    curveFamilyPlot->getHorizontalFineGridInterval()); 
  xmlState->setAttribute(String(("VerticalCoarseGridIsVisible")),
    curveFamilyPlot->isVerticalCoarseGridVisible());
  xmlState->setAttribute(String(("VerticalCoarseGridInterval")),
    curveFamilyPlot->getVerticalCoarseGridInterval());
  xmlState->setAttribute(String(("VerticalFineGridIsVisible")),      
    curveFamilyPlot->isVerticalFineGridVisible());
  xmlState->setAttribute(String(("VerticalFineGridInterval")),  
    curveFamilyPlot->getVerticalFineGridInterval()); 
  xmlState->setAttribute(String(("RadialCoarseGridIsVisible")),
    curveFamilyPlot->isRadialCoarseGridVisible());
  xmlState->setAttribute(String(("RadialCoarseGridInterval")),
    curveFamilyPlot->getRadialCoarseGridInterval());
  xmlState->setAttribute(String(("RadialFineGridIsVisible")),      
    curveFamilyPlot->isRadialFineGridVisible());
  xmlState->setAttribute(String(("RadialFineGridInterval")),  
    curveFamilyPlot->getRadialFineGridInterval()); 
  xmlState->setAttribute(String(("AngularCoarseGridIsVisible")),
    curveFamilyPlot->isAngularCoarseGridVisible());
  xmlState->setAttribute(String(("AngularCoarseGridInterval")),
    curveFamilyPlot->getAngularCoarseGridInterval());
  xmlState->setAttribute(String(("AngularFineGridIsVisible")),      
    curveFamilyPlot->isAngularFineGridVisible());
  xmlState->setAttribute(String(("AngularFineGridInterval")),  
    curveFamilyPlot->getAngularFineGridInterval()); 

  // 3D stuff to come .....

  xmlState->setAttribute(String(("DataSource")),  dataSetup->dataSourceComboBox->getText());
  xmlState->setAttribute(String(("DataType")),    dataSetup->typeOfDataComboBox->getText());
  xmlState->setAttribute(String(("NumCurves")),   (int) dataSetup->numCurvesSlider->getValue());
  if( dataSetup->dataSourceComboBox->getText() == String(("Expression")) )
  {
    xmlState->setAttribute(String(("xExpression")),  dataSetup->xCurveEditLabel->getText());
    xmlState->setAttribute(String(("yExpression")),  dataSetup->yCurveEditLabel->getText());
    xmlState->setAttribute(String(("aMin")),         dataSetup->aMinEditLabel->getText());
    xmlState->setAttribute(String(("a")),            dataSetup->aSlider->getValue());
    xmlState->setAttribute(String(("aMax")),         dataSetup->aMaxEditLabel->getText());
    xmlState->setAttribute(String(("bMin")),         dataSetup->bMinEditLabel->getText());
    xmlState->setAttribute(String(("b")),            dataSetup->bSlider->getValue());
    xmlState->setAttribute(String(("bMax")),         dataSetup->bMaxEditLabel->getText());
    xmlState->setAttribute(String(("cMin")),         dataSetup->cMinEditLabel->getText());
    xmlState->setAttribute(String(("c")),            dataSetup->cSlider->getValue());
    xmlState->setAttribute(String(("cMax")),         dataSetup->cMaxEditLabel->getText());
    xmlState->setAttribute(String(("dMin")),         dataSetup->dMinEditLabel->getText());
    xmlState->setAttribute(String(("d")),            dataSetup->dSlider->getValue());
    xmlState->setAttribute(String(("dMax")),         dataSetup->dMaxEditLabel->getText());
    xmlState->setAttribute(String(("tMin")),         dataSetup->tMinEditLabel->getText());
    xmlState->setAttribute(String(("tMax")),         dataSetup->tMaxEditLabel->getText());
    xmlState->setAttribute(String(("tNumSamples")),  (int) dataSetup->numSamplesSlider->getValue());
    xmlState->setAttribute(String(("tExpSpacing")),   
      dataSetup->tExponentialSpacingButton->getToggleState()); 
  }
  else if( dataSetup->dataSourceComboBox->getText() == String(("File")) )
  {
    xmlState->setAttribute(String(("DataFile")), dataFileName);
  }

  return xmlState;
}

bool RSPlotContentComponent::setStateFromXml(const XmlElement &xmlState)
{
 bool success = true; // should report about success, not used yet

 double min, max, coarse, fine;
 bool   logScale;
 String axisPosition;
 int    axisPositionIndex;

 fileSetup->imageWidthEditLabel->setText(
   String(xmlState.getIntAttribute(("ExportWidth"), 320)));
 fileSetup->imageHeightEditLabel->setText(
   String(xmlState.getIntAttribute(("ExportHeight"), 320)));

 // setup the widgets according to the xmlState and thereby we also send the 
 // change-messages which will be picked up by the callback-functions and do something to the
 // CoordinateSystem - therefore we suppress the re-rendering temporarily because we don't want
 // to re-render several tens of times:
 curveFamilyPlot->setAutoReRendering(false);

 // x-axis and vertical grid:
 axesSetup->xAnnotationEditLabel->setText(
   xmlState.getStringAttribute(String(("xAnnotation")), String(("x"))), NotificationType::sendNotification);
 min          = xmlState.getDoubleAttribute(String(("xMin")),   -1.1);
 max          = xmlState.getDoubleAttribute(String(("xMax")),    1.1);
 coarse       = xmlState.getDoubleAttribute(String(("VerticalCoarseGridInterval")), 1.0);
 fine         = xmlState.getDoubleAttribute(String(("VerticalFineGridInterval")),   0.1);
 logScale     = xmlState.getBoolAttribute(String(("xLogScale")), false);
 axisPosition = xmlState.getStringAttribute(String(("xAxisPosition")), String(("Zero")));
 if( axisPosition == String(("Invisible")) )
   axisPositionIndex = rsPlotSettings::INVISIBLE;
 else if( axisPosition == String(("Zero")) )
   axisPositionIndex = rsPlotSettings::ZERO;
 else if( axisPosition == String(("Top")) )
   axisPositionIndex = rsPlotSettings::TOP;
 else if( axisPosition == String(("Bottom")) )
   axisPositionIndex = rsPlotSettings::BOTTOM;
 axesSetup->xMinEditLabel->setText(String(min), NotificationType::dontSendNotification);
 axesSetup->xMaxEditLabel->setText(String(max), NotificationType::dontSendNotification);
 axesSetup->verticalCoarseGridIntervalLabel->setText(String(coarse), NotificationType::dontSendNotification);
 axesSetup->verticalFineGridIntervalLabel->setText(String(fine), NotificationType::dontSendNotification);
 axesSetup->xLogScaleButton->setToggleState(logScale, false);
 axesSetup->xPosComboBox->setText(axisPosition);
 curveFamilyPlot->setupAxisX(min, max, logScale, axisPositionIndex, coarse, fine);
 axesSetup->verticalCoarseGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("VerticalCoarseGridIsVisible")), false), true);
 axesSetup->verticalFineGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("VerticalFineGridIsVisible")), false), true);

 // y-axis and horizontal grid:
 axesSetup->yAnnotationEditLabel->setText(
   xmlState.getStringAttribute(String(("yAnnotation")), String(("y"))), NotificationType::sendNotification);
 min          = xmlState.getDoubleAttribute(String(("yMin")),   -1.1);
 max          = xmlState.getDoubleAttribute(String(("yMax")),    1.1);
 coarse       = xmlState.getDoubleAttribute(String(("HorizontalCoarseGridInterval")), 1.0);
 fine         = xmlState.getDoubleAttribute(String(("HorizontalFineGridInterval")),   0.1);
 logScale     = xmlState.getBoolAttribute(String(("yLogScale")), false);
 axisPosition = xmlState.getStringAttribute(String(("yPosPosition")), String(("Zero")));
 if( axisPosition == String(("Invisible")) )
   axisPositionIndex = rsPlotSettings::INVISIBLE;
 else if( axisPosition == String(("Zero")) )
   axisPositionIndex = rsPlotSettings::ZERO;
 else if( axisPosition == String(("Left")) )
   axisPositionIndex = rsPlotSettings::LEFT;
 else if( axisPosition == String(("Right")) )
   axisPositionIndex = rsPlotSettings::RIGHT;
 axesSetup->yMinEditLabel->setText(String(min), NotificationType::dontSendNotification);
 axesSetup->yMaxEditLabel->setText(String(max), NotificationType::dontSendNotification);
 axesSetup->horizontalCoarseGridIntervalLabel->setText(String(coarse), NotificationType::dontSendNotification);
 axesSetup->horizontalFineGridIntervalLabel->setText(String(fine), NotificationType::dontSendNotification);
 axesSetup->yLogScaleButton->setToggleState(logScale, false);
 axesSetup->yPosComboBox->setText(axisPosition);
 curveFamilyPlot->setupAxisY(min, max, logScale, axisPositionIndex, coarse, fine);
 axesSetup->horizontalCoarseGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("HorizontalCoarseGridIsVisible")), false), true);
 axesSetup->horizontalFineGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("HorizontalFineGridIsVisible")), false), true);

 // radial and angular grids:
 axesSetup->radialCoarseGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("RadialCoarseGridIsVisible")), false), true);
 axesSetup->radialCoarseGridIntervalLabel->setText(
   String(xmlState.getDoubleAttribute(String(("RadialCoarseGridInterval")),  1.0)), NotificationType::sendNotification);
 axesSetup->radialFineGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("RadialFineGridIsVisible")), false), true);
 axesSetup->radialFineGridIntervalLabel->setText(
   String(xmlState.getDoubleAttribute(String(("RadialFineGridInterval")),  0.1)), NotificationType::sendNotification);

 axesSetup->angularCoarseGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("AngularCoarseGridIsVisible")), false), true);
 axesSetup->angularCoarseGridIntervalLabel->setText(
   String(xmlState.getDoubleAttribute(String(("AngularCoarseGridInterval")),  15.0)), NotificationType::sendNotification);
 axesSetup->angularFineGridButton->setToggleState(
   xmlState.getBoolAttribute(String(("AngularFineGridIsVisible")), false), true);
 axesSetup->angularFineGridIntervalLabel->setText(
   String(xmlState.getDoubleAttribute(String(("AngularFineGridInterval")),  5.0)), NotificationType::sendNotification);


 // 3D stuff to come.....


 dataSetup->dataSourceComboBox->setText(
   xmlState.getStringAttribute(String(("DataSource")), String(("Expression"))), NotificationType::dontSendNotification); 
 dataSetup->typeOfDataComboBox->setText(
   xmlState.getStringAttribute(String(("DataType")), String(("Function"))), NotificationType::dontSendNotification); 
 dataSetup->numCurvesSlider->setValue(
   xmlState.getIntAttribute(String(("NumCurves")), 1), NotificationType::sendNotification);

 if( dataSetup->dataSourceComboBox->getText() == String(("Expression")) )
 {
   //dataSetup->xCurveEditLabel->setText(
   //  xmlState.getStringAttribute(String(("xExpression")), String(("t;"))), NotificationType::sendNotification);

   dataSetup->xCurveEditLabel->setText(xmlState.getStringAttribute("xExpression", "t;"));

   dataSetup->yCurveEditLabel->setText(xmlState.getStringAttribute("yExpression", "x;"));

   min = xmlState.getDoubleAttribute(String(("aMin")),    0.0);
   max = xmlState.getDoubleAttribute(String(("aMax")),    1.0);
   dataSetup->setMinMaxA(min, max);
   dataSetup->aSlider->setValue(xmlState.getDoubleAttribute(String(("a")), 0.5), NotificationType::dontSendNotification);

   min = xmlState.getDoubleAttribute(String(("bMin")),    0.0);
   max = xmlState.getDoubleAttribute(String(("bMax")),    1.0);
   dataSetup->setMinMaxB(min, max);
   dataSetup->bSlider->setValue(xmlState.getDoubleAttribute(String(("b")), 0.5), NotificationType::dontSendNotification);

   min = xmlState.getDoubleAttribute(String(("cMin")),    0.0);
   max = xmlState.getDoubleAttribute(String(("cMax")),    1.0);
   dataSetup->setMinMaxC(min, max);
   dataSetup->cSlider->setValue(xmlState.getDoubleAttribute(String(("c")), 0.5), NotificationType::dontSendNotification);

   min = xmlState.getDoubleAttribute(String(("dMin")),    0.0);
   max = xmlState.getDoubleAttribute(String(("dMax")),    1.0);
   dataSetup->setMinMaxD(min, max);
   dataSetup->dSlider->setValue(xmlState.getDoubleAttribute(String(("d")), 0.5), NotificationType::dontSendNotification);

   min = xmlState.getDoubleAttribute(String(("tMin")),     0.0);
   max = xmlState.getDoubleAttribute(String(("tMax")),    10.0);
   dataSetup->setMinMaxT(min, max);
   dataSetup->numSamplesSlider->setValue(
     xmlState.getIntAttribute(String(("tNumSamples")), 100), NotificationType::sendNotification);
   dataSetup->tExponentialSpacingButton->setToggleState(
     xmlState.getBoolAttribute(String(("tExpSpacing")), false), true);
 }
 else if( dataSetup->dataSourceComboBox->getText() == String(("File")) )
 {
   //dataFileName = xmlState.getStringAttribute(String(("DataFile")), String::empty); // old
   dataFileName = xmlState.getStringAttribute(String(("DataFile")), String()); // new
 }

 // set automatic re-rendering back to true and calculate the data:
 curveFamilyPlot->setAutoReRendering(true);
 if( dataSetup->dataSourceComboBox->getText() == String(("Expression")) )
 {
   calculateData();
   curveFamilyPlot->updatePlotImage(true);
 }
 else if( dataSetup->dataSourceComboBox->getText() == String(("File")) )
 {
   File fileToLoadFrom = File(File::getSpecialLocation(
     File::currentExecutableFile).getParentDirectory().getFullPathName() + String(("/")) 
     + dataFileName);
   loadDataFromFile(fileToLoadFrom);
 }

 // adjust the zoomers to maximally zoomed out:
 zoomer2D->zoomToAllXY();
 //zoomer3D->zoomToAllXY();

 return success;
}


//-------------------------------------------------------------------------------------------------
// appearance stuff:

void RSPlotContentComponent::setBounds(int x, int y, 
                                               int width, int height)
{
  Component::setBounds(x, y, width, height);

  int menuHeight = 128;
  int fileWidth  = 128;
  int zoomerSize = zoomer2D->getZoomerSize();

  curveFamilyPlot->setBounds(0, 0, getWidth()-zoomerSize, getHeight()-menuHeight-zoomerSize);
  zoomer2D->alignWidgetsToCoordinateSystem();

  //surfacePlot->setBounds(0, 0, getWidth()-zoomerSize, getHeight()-menuHeight-zoomerSize);
  //zoomer3D->alignWidgetsToCoordinateSystem3D();

  //setupTabber->setBounds(0, 0, getWidth()/4, getHeight());
  setupTabber->setBounds(0, getHeight()-menuHeight, getWidth()-fileWidth, menuHeight);

  fileSetup->setBounds(getWidth()-fileWidth, getHeight()-menuHeight, fileWidth, menuHeight);
}

//-------------------------------------------------------------------------------------------------
// data generation or import:

void RSPlotContentComponent::calculateData()
{
  // calculate the vaules for the quick-calc field:
  updateCalculatorResult();

  // retrieve the number of values and number of curves to be drawn from the widgets:
  int numValues = (int) dataSetup->numSamplesSlider->getValue();
  numValues     = jlimit(1, maxNumValues, numValues);
  int numCurves = (int) dataSetup->numCurvesSlider->getValue();
  numCurves     = jlimit(1, maxNumCurves, numCurves);

  // read out the parameter sliders and assign the parameters accordingly:
  evaluator.assignVariable("a", dataSetup->aSlider->getValue());
  evaluator.assignVariable("b", dataSetup->bSlider->getValue());
  evaluator.assignVariable("c", dataSetup->cSlider->getValue());
  evaluator.assignVariable("d", dataSetup->dSlider->getValue());

  // determine whether we need to draw a bunch of (possibly) unrelated curves or a family of curves
  // indexed by n:
  String xString = dataSetup->xCurveEditLabel->getText();
  String yString = dataSetup->yCurveEditLabel->getText();
  xString = xString.removeCharacters((" "));
  yString = yString.removeCharacters((" "));
  bool curvesAreFamily = true;
  if( xString.contains((";ßx2=")) || yString.contains((";ßy2=")) )
    curvesAreFamily = false;

  // when the curves do form a family, we use the corresponding function and are done:
  if( curvesAreFamily == true )
  {
    fillArraysWithFamilyData();
    curveFamilyPlot->setCurveFamilyValues(numValues, numCurves, xFamilyPointer, yFamilyPointer);
    return;
  }

  // the curves do not form a family ...
  // ... we need to do some string-splitting work, declare a temporary work string:
  String tmpString;
  String startString;  
  String endString;
  int    numExpressions;
  int    i;

  numExpressions = 1;
  if( xString.contains((";ßx2=")) )
  {
    // determine the number of different expressions for x:
    for(i=2; i<=maxNumCurves; i++)
    {
      tmpString = String((";ßx")) + String(i) + String(("="));
      if( xString.contains(tmpString) )
        numExpressions++;
      else
        break; // jump out of the loop for the first missing ;xi=
    }
    numExpressions = jlimit(1, maxNumCurves, numExpressions);

    // extract the successive substrings for the individual expressions for the different x-arrays 
    // and evaluate them separately:
    startString = String(("x1="));
    endString   = String((";ßx2="));
    for(i=1; i<=numExpressions-1; i++)
    {
      if( xString.contains(startString) )
      {
        tmpString = xString.fromFirstOccurrenceOf(startString, false, false);
        tmpString = tmpString.upToFirstOccurrenceOf(endString, false, false);
        tmpString = tmpString + String((";"));
        fillArrayWithDataX(tmpString, i-1);
      }
      startString = String((";ßx")) + String(i+1) + String(("="));
      endString   = String((";ßx")) + String(i+2) + String(("="));
    }
    // the last substring end only with ";" and must be treated separately therfore:
    tmpString = xString.fromFirstOccurrenceOf(startString, false, false);
    tmpString = tmpString.upToLastOccurrenceOf(String((";")), true, false);
    fillArrayWithDataX(tmpString, i-1);
  }
  else
  {
    i = 1;
    fillArrayWithDataX(xString, i-1);
  }
  // fill the remaining x-arrays with the values of the last filled array:
  int j;
  while(i<maxNumCurves)
  {
    for(j=0; j<maxNumValues; j++)
      xValues[i][j] = xValues[i-1][j];
    i++;
  }

  // we have now filled the x-arrays, now we must do the same procedure for the y-arrays:
  numExpressions = 1;
  if( yString.contains((";ßy2=")) )
  {
    for(i=2; i<=maxNumCurves; i++)
    {
      tmpString = String((";ßy")) + String(i) + String(("="));
      if( yString.contains(tmpString) )
        numExpressions++;
      else
        break;
    }
    numExpressions = jlimit(1, maxNumCurves, numExpressions);
    startString    = String(("y1="));
    endString      = String((";ßy2="));
    for(i=1; i<=numExpressions-1; i++)
    {
      if( yString.contains(startString) )
      {
        tmpString = yString.fromFirstOccurrenceOf(startString, false, false);
        tmpString = tmpString.upToFirstOccurrenceOf(endString, false, false);
        tmpString = tmpString + String((";"));
        fillArrayWithDataY(tmpString, i-1);
      }
      startString = String((";ßy")) + String(i+1) + String(("="));
      endString   = String((";ßy")) + String(i+2) + String(("="));
    }
    tmpString = yString.fromFirstOccurrenceOf(startString, false, false);
    tmpString = tmpString.upToLastOccurrenceOf(String((";")), true, false);
    fillArrayWithDataY(tmpString, i-1);

  }
  else
  {
    i = 1;
    fillArrayWithDataY(yString, i-1);
  }
  while(i<maxNumCurves)
  {
    for(j=0; j<maxNumValues; j++)
      yValues[i][j] = yValues[i-1][j];
    i++;
  }

  // pass the new data over to the CurveFamilyPlot object:
  curveFamilyPlot->setCurveFamilyValues(numValues, numCurves, xFamilyPointer, yFamilyPointer);

}

void RSPlotContentComponent::updateCalculatorResult()
{
  double tVal = dataSetup->quickCalcInputLabel->getText().getDoubleValue();
  double xVal, yVal;

  // get the string for x(t) from the editor-field:
  String xString = dataSetup->xCurveEditLabel->getText();
  //xString        = xString.upToFirstOccurrenceOf(String(("ß")), false, false); // triggers jassert
  xString = xString.upToFirstOccurrenceOf(CharPointer_UTF8("ß"), false, false);

  // convert the juce-string to c-string:
  long  length = xString.length();
  char* xStringC = new char[length+1];
  xString.copyToUTF8(xStringC, length+1);

  // create a boolean variable for checking the error state:
  bool evalSuccess = false;

  // parse the expression:
  evalSuccess = evaluator.setExpressionString(xStringC);
  if( !evalSuccess )
  {
    dataSetup->xCurveEditLabel->setColour(Label::backgroundColourId, 
      Colour((uint8)255, (uint8)225, (uint8)225, (uint8)255));
    if(xStringC)
      delete[] xStringC; // delete the char-array before leaving the function
    return;
  }
  else
    dataSetup->xCurveEditLabel->setColour(Label::backgroundColourId, Colours::white);

  // assign the name "t" to the new value:
  evaluator.assignVariable("t", tVal);
  evaluator.assignVariable("n", 1.0);

  // evaluate the expression for this particular t and write the result
  // into xVal (no error checking is done here (!) ):
  xVal = evaluator.evaluateExpression();

  String yString = dataSetup->yCurveEditLabel->getText();
  //yString = yString.upToFirstOccurrenceOf(String(("ß")), false, false);
  yString = yString.upToFirstOccurrenceOf(CharPointer_UTF8("ß"), false, false);
  length = yString.length();
  char* yStringC = new char[length+1];
  yString.copyToUTF8(yStringC, length+1);
  evalSuccess = false;
  evalSuccess = evaluator.setExpressionString(yStringC);
  if( !evalSuccess )
  {
    dataSetup->yCurveEditLabel->setColour(Label::backgroundColourId, 
      Colour((uint8)255, (uint8)225, (uint8)225, (uint8)255));
    if(xStringC)
      delete[] xStringC;
    if(yStringC)
      delete[] yStringC; // delete the char-array before leaving the function
    return;
  }
  else
    dataSetup->yCurveEditLabel->setColour(Label::backgroundColourId, Colours::white);

  // assign the variable-names:
  evaluator.assignVariable("t", tVal);
  evaluator.assignVariable("n", 1.0);
  evaluator.assignVariable("x", xVal);

  // evaluate the expression for this particular t and write the result
  // into xVal (no error checking is done here (!) ):
  yVal = evaluator.evaluateExpression();

  // update the reslut fields:
  dataSetup->quickCalcOutputLabelX->setText(String(("x=")) + String(xVal), NotificationType::dontSendNotification);
  dataSetup->quickCalcOutputLabelY->setText(String(("y=")) + String(yVal), NotificationType::dontSendNotification);

  // delete the created c-strings:
  if(xStringC)
    delete[] xStringC;
  if(yStringC)
    delete[] yStringC;
}

void RSPlotContentComponent::fillArraysWithFamilyData()
{
  // retrieve the number of values and number of curves to be drawn from the widgets:
  int numValues = (int) dataSetup->numSamplesSlider->getValue();
  numValues     = jlimit(1, maxNumValues, numValues);

  int numCurves = (int) dataSetup->numCurvesSlider->getValue();
  numCurves     = jlimit(1, maxNumCurves, numCurves);
 
  // get the string for x(t) from the editor-field:
  String xString = dataSetup->xCurveEditLabel->getText();

  // convert the juce-string to c-string:
  long  length = xString.length();
  char* xStringC = new char[length+1];
  xString.copyToUTF8(xStringC, length+1);

  // create a boolean variable for checking the error state:
  bool evalSuccess = false;

  // parse the expression:
  evalSuccess = evaluator.setExpressionString(xStringC);
  if( !evalSuccess )
  {
    dataSetup->xCurveEditLabel->setColour(Label::backgroundColourId, 
      Colour((uint8)255, (uint8)225, (uint8)225, (uint8)255));
    if(xStringC)
      delete[] xStringC; // delete the char-array before leaving the function
    return;
  }
  else
    dataSetup->xCurveEditLabel->setColour(Label::backgroundColourId, Colours::white);

  // declare/init/calculate the variables for the sampling of the curves:
  double tMin       = dataSetup->tMinEditLabel->getText().getDoubleValue();
  double tMax       = dataSetup->tMaxEditLabel->getText().getDoubleValue();
  bool   expSpacing = dataSetup->tExponentialSpacingButton->getToggleState();
  double tInterval  = (tMax - tMin) / (numValues-1);
  double tFactor    = pow(tMax/tMin, 1.0/(numValues-1));
  double tVal, xVal, yVal;
  int    n, i;

  // loop over the 'numCurves' curves in the curve family:
  for(n=0; n<numCurves; n++)
  {
    // assign n and initialize t:
    evaluator.assignVariable("n", (double) (n+1));
    evaluator.assignVariable("t", tMin);

    // some additional initialization only relavant for exponential spacing:
    tVal = tMin/tFactor; 

    // evaluate the expression (calculate x) for different values of t:
    for(i=0; i<numValues; i++)
    {
      // calculate t-value:
      if( !expSpacing )
        tVal = tMin + i*tInterval; 
      else
        tVal *= tFactor;  

      // reassign the name "t" to the new value:
      evaluator.assignVariable("t", tVal);

      // evaluate the expression for this particular t and write the result
      // into xVal (no error checking is done here (!) ):
      xVal = evaluator.evaluateExpression();

      // write the result into the buffer for the x-values (thereby 
      // typecasting it to float:
      xValues[n][i] = xVal;
    }
  }

  // now we have calculated all the x values, we now do the same procedure for the y-values (this
  // time without the commenting):
  String yString = dataSetup->yCurveEditLabel->getText();
  length = yString.length();
  char* yStringC = new char[length+1];
  yString.copyToUTF8(yStringC, length+1);
  evalSuccess = false;
  evalSuccess = evaluator.setExpressionString(yStringC);
  if( !evalSuccess )
  {
    dataSetup->yCurveEditLabel->setColour(Label::backgroundColourId, 
      Colour((uint8)255, (uint8)225, (uint8)225, (uint8)255));
    if(xStringC)
      delete[] xStringC;
    if(yStringC)
      delete[] yStringC; // delete the char-array before leaving the function
    return;
  }
  else
    dataSetup->yCurveEditLabel->setColour(Label::backgroundColourId, Colours::white);
  for(n=0; n<numCurves; n++)
  {
    evaluator.assignVariable("n", (double) (n+1));
    evaluator.assignVariable("t", tMin);
    tVal = tMin/tFactor; 
    for(i=0; i<numValues; i++)
    {
      if( !expSpacing )
        tVal = tMin + i*tInterval; 
      else
        tVal *= tFactor;  
      evaluator.assignVariable("t", tVal);
      evaluator.assignVariable("x", xValues[n][i]); // (re-)assign the name "x" to the x-value
      yVal = evaluator.evaluateExpression();
      yValues[n][i] = yVal;
    }
  }

  // delete the created c-strings:
  if(xStringC)
    delete[] xStringC;
  if(yStringC)
    delete[] yStringC;
}

void RSPlotContentComponent::fillArrayWithDataX(const String &xString, int curveIndex)
{
  // convert the juce-string to c-string:
  long  length   = xString.length();
  char* xStringC = new char[length+1];
  char  xN[3];   
  xString.copyToUTF8(xStringC, length+1);

  // parse the expression:
  bool evalSuccess = evaluator.setExpressionString(xStringC);
  if( !evalSuccess )
  {
    dataSetup->xCurveEditLabel->setColour(Label::backgroundColourId, 
      Colour((uint8)255, (uint8)225, (uint8)225, (uint8)255));
    if(xStringC)
      delete[] xStringC; // delete the char-array before leaving the function
    return;
  }
  else
    dataSetup->xCurveEditLabel->setColour(Label::backgroundColourId, Colours::white);

  // declare/init/calculate the variables for the sampling of the curves:
  // retrieve the number of values and number of curves to be drawn from the widgets:
  int numValues     = (int) dataSetup->numSamplesSlider->getValue();
  numValues         = jlimit(1, maxNumValues, numValues);
  double tMin       = dataSetup->tMinEditLabel->getText().getDoubleValue();
  double tMax       = dataSetup->tMaxEditLabel->getText().getDoubleValue();
  bool   expSpacing = dataSetup->tExponentialSpacingButton->getToggleState();
  double tInterval  = (tMax - tMin) / (numValues-1);
  double tFactor    = pow(tMax/tMin, 1.0/(numValues-1));
  double tVal, xVal;

  // assign n and initialize t:
  evaluator.assignVariable("n", (double) curveIndex);
  evaluator.assignVariable("t", tMin);

  // some additional initialization only relavant for exponential spacing:
  tVal = tMin/tFactor; 

  // evaluate the expression (calculate x) for different values of t:
  for(int i=0; i<numValues; i++)
  {
    // calculate t-value:
    if( !expSpacing )
      tVal = tMin + i*tInterval; 
    else
      tVal *= tFactor;  

    // reassign the name "t" to the new value:
    evaluator.assignVariable("t", tVal);

    // assign x1, x2, x3, etc. up the xN just below the one, we are currently generating:
    for(int j=0; j<curveIndex; j++)
    {
      (String(("x")) + String(j+1)).copyToUTF8(xN, 3);
      evaluator.assignVariable(xN, xValues[j][i]);
    }

    // evaluate the expression for this particular t and write the result
    // into xVal (no error checking is done here (!) ):
    xVal = evaluator.evaluateExpression();

    // write the result into the buffer for the x-values (thereby 
    // typecasting it to float:
    xValues[curveIndex][i] = xVal;
  }

  // delete the created c-string:
  if(xStringC)
    delete[] xStringC;
}

void RSPlotContentComponent::fillArrayWithDataY(const String &yString, int curveIndex)
{
  // convert the juce-string to c-string:
  long  length   = yString.length();
  char* yStringC = new char[length+1];
  char  yN[3];       
  yString.copyToUTF8(yStringC, length+1);

  // parse the expression:
  bool evalSuccess = evaluator.setExpressionString(yStringC);
  if( !evalSuccess )
  {
    dataSetup->yCurveEditLabel->setColour(Label::backgroundColourId, 
      Colour((uint8)255, (uint8)225, (uint8)225, (uint8)255));
    if(yStringC)
      delete[] yStringC; // delete the char-array before leaving the function
    return;
  }
  else
    dataSetup->yCurveEditLabel->setColour(Label::backgroundColourId, Colours::white);

  // declare/init/calculate the variables for the sampling of the curves:
  // retrieve the number of values and number of curves to be drawn from the widgets:
  int numValues     = (int) dataSetup->numSamplesSlider->getValue();
  numValues         = jlimit(1, maxNumValues, numValues);
  double tMin       = dataSetup->tMinEditLabel->getText().getDoubleValue();
  double tMax       = dataSetup->tMaxEditLabel->getText().getDoubleValue();
  bool   expSpacing = dataSetup->tExponentialSpacingButton->getToggleState();
  double tInterval  = (tMax - tMin) / (numValues-1);
  double tFactor    = pow(tMax/tMin, 1.0/(numValues-1));
  double tVal, yVal;

  // assign n and initialize t:
  evaluator.assignVariable("n", (double) curveIndex);
  evaluator.assignVariable("t", tMin);
  evaluator.assignVariable("x", xValues[curveIndex][0]);

  // some additional initialization only relavant for exponential spacing:
  tVal = tMin/tFactor; 

  // evaluate the expression (calculate x) for different values of t:
  for(int i=0; i<numValues; i++)
  {
    // calculate t-value:
    if( !expSpacing )
      tVal = tMin + i*tInterval; 
    else
      tVal *= tFactor;  

    // reassign the names "t" and "x" to the new values:
    evaluator.assignVariable("t", tVal);
    evaluator.assignVariable("x", xValues[curveIndex][i]); // (re-)assign the name "x" to the x-value

    // assign y1, y2, y3, etc. up the yN just below the one, we are currently generating:
    for(int j=0; j<curveIndex; j++)
    {
      (String(("y")) + String(j+1)).copyToUTF8(yN, 3);
      evaluator.assignVariable(yN, yValues[j][i]);
    }

    // evaluate the expression for this particular t and x and write the result
    // into yVal (no error checking is done here (!) ):  
    yVal = evaluator.evaluateExpression();

    // write the result into the buffer for the y-values:
    yValues[curveIndex][i] = yVal;
  }

  // delete the created c-string:
  if(yStringC)
    delete[] yStringC;
}

void RSPlotContentComponent::openDataLoadingDialog()
{
  FileChooser chooser(("Load Data"), File(currentImageDirectory), String(("*.dat")), true);
  if(chooser.browseForFileToOpen())
    loadDataFromFile(chooser.getResult());
}

void RSPlotContentComponent::loadDataFromFile(const File &dataFileToLoad)
{
  if( !dataFileToLoad.existsAsFile() )
    return;

  String dataString = dataFileToLoad.loadFileAsString();
  dataFileName = dataFileToLoad.getRelativePathFrom(File::getSpecialLocation(
    File::currentExecutableFile).getParentDirectory());
  parseDataString(dataString);
  dataSetup->fileNameLabel->setText(dataFileName, NotificationType::dontSendNotification);
}

void RSPlotContentComponent::parseDataString(const String &dataString)
{
  StringArray lines;
  StringArray tokens1, tokens2;

  double test1, test2;

  lines.addLines(dataString);

  // matlab puts a whitespace in front of the first number and separates values with two 
  // whitespaces - we fix that here by replacing a double whitespace with a single whitespace
  // and also replace tab-stops with whitespaces:
  for(int lineIndex=0; lineIndex<lines.size(); lineIndex++)
  {
    lines.set(lineIndex, lines[lineIndex].trim());
    lines.set(lineIndex, lines[lineIndex].replace(String(("  ")), String((" "))));
    lines.set(lineIndex, lines[lineIndex].replace(String(("\t")), String((" "))));
  }

  int typeOfData = dataSetup->typeOfDataComboBox->getSelectedId();
  switch( typeOfData )
  {
  case FUNCTION_ON_PLANE:
    {
      for(int lineIndex=0; lineIndex < lines.size(); lineIndex++)
      {
        tokens1.clear();
        tokens1.addTokens(lines[lineIndex], true);

        for(int tokenIndex=0; tokenIndex < tokens1.size(); tokenIndex++)
        {
          if( lineIndex == 0 )
            xValues[0][tokenIndex] = tokens1[tokenIndex].getDoubleValue();
          else if( lineIndex == 1 )
            yValues[0][tokenIndex] = tokens1[tokenIndex].getDoubleValue();
        }
        for(int i = tokens1.size(); i < 1024; i++)
        {
          xValues[0][i] = 0.0;
          yValues[0][i] = 0.0;
        }
      }

      // pass the new x- and y-buffers to the funcPlotter and trigger a repaint:
      curveFamilyPlot->setCurveFamilyValues(tokens1.size(), 1, xFamilyPointer, yFamilyPointer);
    }
    break;

  case FUNCTION_FAMILY_ON_PLANE:
    {
      int numFunctions = lines.size()-1; //

      tokens1.clear();
      tokens1.addTokens(lines[0],   true);

      for(int functionIndex=0; functionIndex < numFunctions; functionIndex++)
      {
        tokens2.clear();
        tokens2.addTokens(lines[functionIndex+1], true);

        for(int tokenIndex=0; tokenIndex < tokens2.size(); tokenIndex++)
        {
          test1 = xValues[functionIndex][tokenIndex] = tokens1[tokenIndex].getDoubleValue();
          test2 = yValues[functionIndex][tokenIndex] = tokens2[tokenIndex].getDoubleValue();

          int dummy = 0;
        }
        for(int i = tokens2.size(); i < 1024; i++)
        {
          xValues[functionIndex][i] = 0.0;
          yValues[functionIndex][i] = 0.0;
        }
      }

      // pass the new x- and y-buffers to the funcPlotter and trigger a repaint:
      curveFamilyPlot->setCurveFamilyValues(tokens1.size(), numFunctions, 
        xFamilyPointer, yFamilyPointer);
    }
    break;


  case CURVE_FAMILY_ON_PLANE:
    {
      int numCurves = lines.size()/2; // the number of lines is assumed to be even 
      //int offset =

      for(int curveIndex=0; curveIndex < numCurves; curveIndex++)
      {
        tokens1.clear();
        tokens1.addTokens(lines[curveIndex],   true);
        tokens2.clear();
        tokens2.addTokens(lines[curveIndex+numCurves], true);

        for(int tokenIndex=0; tokenIndex < tokens1.size(); tokenIndex++)
        {
          test1 = xValues[curveIndex][tokenIndex] = tokens1[tokenIndex].getDoubleValue();
          test2 = yValues[curveIndex][tokenIndex] = tokens2[tokenIndex].getDoubleValue();

          int dummy = 0;
        }
        for(int i = tokens1.size(); i < 1024; i++)
        {
          xValues[curveIndex][i] = 0.0;
          yValues[curveIndex][i] = 0.0;
        }
      }
			
      // pass the new x- and y-buffers to the funcPlotter and trigger a repaint:
      curveFamilyPlot->setCurveFamilyValues(tokens1.size(), numCurves, 
        xFamilyPointer, yFamilyPointer);
    }
    break;
		
  } // end of switch     switch( typeOfData )
	
}