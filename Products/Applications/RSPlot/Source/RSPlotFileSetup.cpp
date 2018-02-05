#include "RSPlotFileSetup.h"

RSPlotFileSetup::RSPlotFileSetup()
{
  loadButton = new RButton(String(("Load")) );
  loadButton->setClickingTogglesState(false);
  addAndMakeVisible(loadButton);

  saveButton = new RButton(String(("Save")) );
  saveButton->setClickingTogglesState(false);
  addAndMakeVisible(saveButton);

  plusButton = new RButton(RButton::PLUS);
  plusButton->setClickingTogglesState(false);
  addAndMakeVisible(plusButton);

  minusButton = new RButton(RButton::MINUS);
  minusButton->setClickingTogglesState(false);
  addAndMakeVisible(minusButton);

  imageExportLabel = new RTextField(String(("Export:")));
  imageExportLabel->setJustification(Justification::centredLeft);
  addAndMakeVisible( imageExportLabel );

  imageWidthLabel = new RTextField(String(("W:")));
  imageWidthLabel->setJustification(Justification::centredLeft);
  addAndMakeVisible( imageWidthLabel );

  imageHeightLabel = new RTextField(String(("H:")));
  imageHeightLabel->setJustification(Justification::centredLeft);
  addAndMakeVisible( imageHeightLabel );

  addAndMakeVisible( fileNameLabel = new RTextEntryField() );
  fileNameLabel->setJustification(Justification::centredLeft);
  //fileNameLabel->setColour(Label::backgroundColourId, Colours::white);
  //fileNameLabel->setColour(Label::outlineColourId, Colours::black);
  //addAndMakeVisible( fileNameLabel );

  addAndMakeVisible( imageWidthEditLabel = new RTextEntryField(String(("320"))) );
  //imageWidthEditLabel->setReadOnly(false);
  imageWidthEditLabel->setJustification(Justification::centredLeft);
  //imageWidthEditLabel->setColour(Label::outlineColourId, Colours::black);
  //imageWidthEditLabel->setColour(Label::backgroundColourId, Colours::white);

  addAndMakeVisible( imageHeightEditLabel = new RTextEntryField(String(("320"))) );
  //imageHeightEditLabel->setReadOnly(false);
  imageHeightEditLabel->setJustification(Justification::centredLeft);
  //imageHeightEditLabel->setColour(Label::outlineColourId, Colours::black);
  //imageHeightEditLabel->setColour(Label::backgroundColourId, Colours::white);

  pngExportButton = new RButton(String(("PNG")));
  pngExportButton->setClickingTogglesState(false);
  addAndMakeVisible(pngExportButton);

  svgExportButton = new RButton(String(("SVG")));
  svgExportButton->setClickingTogglesState(false);
  addAndMakeVisible(svgExportButton);
}

RSPlotFileSetup::~RSPlotFileSetup()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// appearance:

void RSPlotFileSetup::paint(Graphics &g)
{
  g.fillAll(Colours::white);
}

void RSPlotFileSetup::resized()
{
  Component::resized();

  int x = 0;
  int w = getWidth();
  int y = 0;
  int h = getHeight();

  fileNameLabel->setBounds(x+4, y+4, w-8, 20);
  y += 24;

  saveButton->setBounds(x+4, y+4, 40, 16); 
  x += 44;
  loadButton->setBounds(x+4, y+4, 40, 16);
  x += 44;
  minusButton->setBounds(x+4, y+4, 16, 16);
  x += 16;
  plusButton->setBounds(x+4, y+4, 16, 16);

  x  = 0;
  y += 28;
  imageExportLabel->setBounds(x+4, y+4, 52, 16);
  x += 52;
  pngExportButton->setBounds(x+4,  y+4, 32, 16);
  x += 32+4;
  svgExportButton->setBounds(x+4,  y+4, 32, 16);

  x  = 0;
  y += 24;
  imageWidthLabel->setBounds(x, y, 24, 16);
  imageWidthEditLabel->setBounds(x+24, y, w/2-24-4, 16);
  x += w/2;
  imageHeightLabel->setBounds(x, y, 24, 16);
  imageHeightEditLabel->setBounds(x+24, y, w/2-24-4, 16);
}
