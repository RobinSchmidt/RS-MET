
//-------------------------------------------------------------------------------------------------
// construction/destruction:

ImageFileManager::ImageFileManager() 
{
  imageIsUnsaved               = false;
  currentImageFileName         = String::empty;
  currentImageFileNameWithStar = String::empty;
}

ImageFileManager::~ImageFileManager()
{

}

//-------------------------------------------------------------------------------------------------
// others:

/*
Image* ImageFileManager::loadImageFromFile(const File& fileToLoadFrom)
{

XmlDocument myDocument(fileToLoadFrom);
XmlElement* xmlState = myDocument.getDocumentElement();

if( xmlState != NULL )
{
setStateFromXml(*xmlState);
delete xmlState;
}

imageIsUnsaved = false;

return new Image; // something to do here....
}
*/

void ImageFileManager::saveImageToFile(const File& fileToSaveTo, const Image* imageToSave)
{
  // the user has already confirmed to overwrite the file, so we do it:
  if( fileToSaveTo.existsAsFile() )
    fileToSaveTo.deleteFile();

  // create a PNGImagefileFormat object:
  PNGImageFormat pngFormat;

  // create the file output stream:
  FileOutputStream fileStream(fileToSaveTo);

  bool success = false;
  success = pngFormat.writeImageToStream(*imageToSave, fileStream);

  imageIsUnsaved = false;
}

/*
void ImageFileManager::openImageLoadingDialog()
{

FileChooser chooser(T("Load Image"), File(currentImageDirectory), String(T("*.png")), true);

if(chooser.browseForFileToOpen())
{
File fileToLoad = chooser.getResult();
loadImageFromFile(fileToLoad);

// remember where (in what directory) we are:
currentImageDirectory = fileToLoad.getParentDirectory().getFullPathName(); 

// update the image-field on the GUI
currentImageFileName         = fileToLoad.getFileNameWithoutExtension();
currentImageFileNameWithStar = currentImageFileName + String(T("*"));

// update the array with the images which are in the current 
// image-directory (which is needed for image increment/decrement):
imagesInCurrentImageDirectory.clear();
File(currentImageDirectory).findChildFiles(imagesInCurrentImageDirectory, 
File::findFiles, 
false, 
JUCE_T("*.png"));
numImagesInCurrentImageDirectory = imagesInCurrentImageDirectory.size();
imagesInCurrentImageDirectory.sort(fileComparator);

// find the index of the currently loaded file inside the current directory:
currentImageIndex = 0;
for(int i=0; i<numImagesInCurrentImageDirectory; i++)
{
if(imagesInCurrentImageDirectory[i]->getFileNameWithoutExtension() 
== currentImageFileName)
{
currentImageIndex = i;
break;
}
}
}
}
*/

void ImageFileManager::openImageSavingDialog(const Image* imageToSave)
{
  FileChooser chooser("Save Image", File(currentImageDirectory), String("*.png"), true);

  if(chooser.browseForFileToSave(true))
  {
    File fileToSaveTo = chooser.getResult();

    // auto append the the extension .png if the user did not type it in
    // manually:
    if ( !fileToSaveTo.hasFileExtension(String(".png")) )
      fileToSaveTo = fileToSaveTo.withFileExtension(String(".png")) ;

    saveImageToFile(fileToSaveTo, imageToSave);

    // remember where (in what directory) we are:
    currentImageDirectory = fileToSaveTo.getParentDirectory().getFullPathName();

    // update the image-field on the GUI
    currentImageFileName         = fileToSaveTo.getFileNameWithoutExtension();
    currentImageFileNameWithStar = currentImageFileName + String("*");

    // update the array with the images which are in the current 
    // image-directory (which is needed for image increment/decrement):
    /*
    imagesInCurrentImageDirectory.clear();
    File(currentImageDirectory).findChildFiles(imagesInCurrentImageDirectory, 
      File::findFiles, 
      false, 
      JUCE_T("*.png"));
    numImagesInCurrentImageDirectory = imagesInCurrentImageDirectory.size();
    imagesInCurrentImageDirectory.sort(fileComparator);
    */

    // find the index of the currently loaded file inside the current directory:
    currentImageIndex = 0;
    for(int i=0; i<numImagesInCurrentImageDirectory; i++)
    {
      if(imagesInCurrentImageDirectory[i].getFileNameWithoutExtension() == currentImageFileName)
      {
        currentImageIndex = i;
        break;
      }
    }
  }
}
