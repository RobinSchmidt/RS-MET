//#include "rojue_ImageSavingDialog.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ImageSavingDialog::ImageSavingDialog(rsPlot *owner, int defaultWidth, int defaultHeight, 
                                     const String &defaultFormat, const File &defaultTargetFile)
{
  setSize(384, 52);
  ownerSystem = owner;

  int x, y, w, h;
  x = 0;
  y = 0;
  w = getWidth();
  h = getHeight();

  targetFileLabel = new Label(String(("Target File")), String(("File:")) );
  targetFileLabel->setBounds(4, 4, 32, 20);
  addAndMakeVisible(targetFileLabel);

  browseButton = new TextButton(String(("Browse")));
  browseButton->setBounds(w-48-4, 4, 48, 20);
  browseButton->addListener(this);
  addAndMakeVisible(browseButton);

  x = targetFileLabel->getRight();
  w = browseButton->getX() - targetFileLabel->getRight();
  //targetFileEditLabel = new Label(String(("Target File")), String() );
  targetFileEditLabel = new Label(String(("Target File")), defaultTargetFile.getFullPathName() );
  targetFileEditLabel->setColour(Label::backgroundColourId, Colours::white);
  targetFileEditLabel->setColour(Label::outlineColourId, Colours::black);
  targetFileEditLabel->setEditable(true);
  targetFileEditLabel->setBounds(x+4, 4, w-8, 20);
  addAndMakeVisible(targetFileEditLabel);

  x = 0;
  y = targetFileEditLabel->getBottom();
  pixelWidthLabel = new Label(String(("Pixel width")), String(("Width:")) );
  pixelWidthLabel->setBounds(x+4, y+4, 48, 20);
  addAndMakeVisible(pixelWidthLabel);

  x = pixelWidthLabel->getRight();
  //pixelWidthEditLabel = new Label(String(("Pixel width")), String(owner->getWidth()) );
  pixelWidthEditLabel = new Label(String(("Pixel width")), String(defaultWidth) );
  pixelWidthEditLabel->setBounds(x+4, y+4, 48, 20);
  pixelWidthEditLabel->setColour(Label::backgroundColourId, Colours::white);
  pixelWidthEditLabel->setColour(Label::outlineColourId, Colours::black);
  pixelWidthEditLabel->setEditable(true);
  pixelWidthEditLabel->addListener(this);
  addAndMakeVisible(pixelWidthEditLabel);

  x = pixelWidthEditLabel->getRight();
  pixelHeightLabel = new Label(String(("Pixel height")), String(("Height:")) );
  pixelHeightLabel->setBounds(x+4, y+4, 52, 20);
  addAndMakeVisible(pixelHeightLabel);

  x = pixelHeightLabel->getRight();
  //pixelHeightEditLabel = new Label(String(("Pixel height")), String(owner->getHeight()) );
  pixelHeightEditLabel = new Label(String(("Pixel height")), String(defaultHeight) );
  pixelHeightEditLabel->setBounds(x+4, y+4, 48, 20);
  pixelHeightEditLabel->setColour(Label::backgroundColourId, Colours::white);
  pixelHeightEditLabel->setColour(Label::outlineColourId, Colours::black);
  pixelHeightEditLabel->setEditable(true);
  pixelHeightEditLabel->addListener(this);
  addAndMakeVisible(pixelHeightEditLabel);

  x = pixelHeightEditLabel->getRight();
  formatLabel = new Label(String(("Format")), String(("Format:")) );
  formatLabel->setBounds(x+4, y+4, 52, 20);
  addAndMakeVisible(formatLabel);

  x = formatLabel->getRight();
  w = browseButton->getX() - x;
  formatComboBox = new ComboBox(String(("formatComboBox")));
  formatComboBox->setBounds(x+4, y+4, w-8, 20);
  formatComboBox->addItem(String(("svg")),        1);
  formatComboBox->addItem(String(("png")),        2);
  //formatComboBox->setSelectedId(2, true);
  if( defaultFormat == String(("svg")) )
    formatComboBox->setSelectedId(1);
  else if( defaultFormat == String(("png")) )
    formatComboBox->setSelectedId(2);
  else
    formatComboBox->setSelectedId(2);
  addAndMakeVisible(formatComboBox);

  saveButton = new TextButton(String(("Save")));
  saveButton->setBounds(browseButton->getX(), browseButton->getY()+24, 48, 20);
  saveButton->addListener(this);
  addAndMakeVisible(saveButton);

  // initialize the target file field:
  File   thisExeAsFile         = File::getSpecialLocation(File::currentExecutableFile);
  File   thisDirectoryAsFile   = thisExeAsFile.getParentDirectory();
  String thisDirectoryAsString = thisDirectoryAsFile.getFullPathName();
  File   targetFile            = File(thisDirectoryAsString + String(("/ExportedImage")) );
  if( defaultTargetFile != File() )
    targetFile = defaultTargetFile;
  targetFileEditLabel->setText(targetFile.getFullPathName(), juce::NotificationType::dontSendNotification);
}

ImageSavingDialog::~ImageSavingDialog()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int ImageSavingDialog::getSelectedPixelWidth()
{
  int w = pixelWidthEditLabel->getText().getIntValue();
  if( w > 1 )
    return w;
  else
    return 1;
}

int ImageSavingDialog::getSelectedPixelHeight()
{
  int h = pixelHeightEditLabel->getText().getIntValue();
  if( h > 1 )
    return h;
  else
    return 1;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void ImageSavingDialog::buttonClicked(juce::Button *buttonThatWasClicked)
{
  if( buttonThatWasClicked == browseButton )
  {
    FileChooser chooser(String(("Select target file")), 
      File(targetFileEditLabel->getText()), String(), false);
    if(chooser.browseForFileToSave(false))
    {
      File newTargetFile = chooser.getResult();
      targetFileEditLabel->setText(newTargetFile.getFullPathName(), 
        NotificationType::dontSendNotification);
      saveNow();
    }
  }
  else if( buttonThatWasClicked == saveButton )
  {
    saveNow();
  }
}

void ImageSavingDialog::comboBoxChanged(juce::ComboBox *comboBoxThatHasChanged)
{

}

void ImageSavingDialog::labelTextChanged(juce::Label *labelThatHasChanged)
{

}

//-------------------------------------------------------------------------------------------------
// saving:

void ImageSavingDialog::saveNow()
{
  int w = getSelectedPixelWidth();
  int h = getSelectedPixelHeight();
  File fileToSaveTo = File(targetFileEditLabel->getText());

  if( formatComboBox->getText() == String(("png")) )
  {
    // append the proper extension, if not already there:
    if( !fileToSaveTo.hasFileExtension(String(("png"))) )
      fileToSaveTo = fileToSaveTo.withFileExtension(String(("png")));

    // ask user for overwriting when the file already exists:
    if( fileToSaveTo.existsAsFile() )
    {
      bool overwrite = AlertWindow::showOkCancelBox(
        AlertWindow::WarningIcon, 
        String(("Overwrite Warning")),                        
        String(("File already exists. Overwrite?")),
        String(("yes")),
        String(("no"))  );

      if( overwrite )
        fileToSaveTo.deleteFile();
      else
        return;
    }

    Image* im = ownerSystem->getPlotAsImage(w, h);

    // create a PNGImagefileFormat object:
    PNGImageFormat pngFormat;

    // create the file output stream:
    FileOutputStream fileStream(fileToSaveTo);

    bool success = false;
    success = pngFormat.writeImageToStream(*im, fileStream);

    delete im;
  }
  else if( formatComboBox->getText() == String(("svg")) )
  {
    // append the proper extension, if not already there:
    if( !fileToSaveTo.hasFileExtension(String(("svg"))) )
      fileToSaveTo = fileToSaveTo.withFileExtension(String(("svg")));

    // ask user for overwriting when the file already exists:
    if( fileToSaveTo.existsAsFile() )
    {
      bool overwrite = AlertWindow::showOkCancelBox(
        AlertWindow::WarningIcon, 
        String(("Overwrite Warning")),                        
        String(("File already exists. Overwrite?")),
        String(("yes")),
        String(("no"))  );

      if( overwrite )
        fileToSaveTo.deleteFile();
      else
        return;
    }

    // create a dummy-image - this has the side-effect to make the svg-drawing inside the 
    // CoordinatSystem to have the right dimensions - elegantize this someday...
    //const Image* dummyImage = ownerSystem->getPlotAsImage(w, h);
    //delete dummyImage;

    XmlElement* theSVG = ownerSystem->getPlotAsSVG(w, h);

    /*String myXmlDoc = theSVG->createDocument(String());*/
    String myXmlDoc = theSVG->toString();
    fileToSaveTo.create();
    fileToSaveTo.appendText(myXmlDoc);

    delete theSVG;
  }
}




