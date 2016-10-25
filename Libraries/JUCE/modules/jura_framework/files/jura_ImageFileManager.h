#ifndef jura_ImageFileManager_h
#define jura_ImageFileManager_h

/** This class manages image-files.  .... a lot of nonsense is still going on here... */

class JUCE_API ImageFileManager : virtual public FileManager

{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ImageFileManager();

  /** Destructor. */
  virtual ~ImageFileManager();

protected:

  // open a dialog-box to load/save the parameters from/to a xml file:
  //virtual Image* openImageLoadingDialog();
  virtual void   openImageSavingDialog(const Image* imageToSave);

  // the actual load/save functions:
  //virtual Image* loadImageFromFile(const File& fileToLoadFrom);
  virtual void   saveImageToFile(const juce::File& fileToSaveTo, const Image* imageToSave);

  // for remebering the directories in which was most recntly browsed for 
  // images:
  bool imageIsUnsaved;
  juce::String currentImageFileNameWithStar;
  juce::String currentImageDirectory;
  juce::String currentImageFileName;
  juce::Array<juce::File> imagesInCurrentImageDirectory;  // holds all image files in the current directory
  int currentImageIndex;                                  // index of the current image in the array
  int numImagesInCurrentImageDirectory;
  FileComparator fileComparator;                          // this object is needed to sort the image-array

  juce_UseDebuggingNewOperator;
};

#endif 