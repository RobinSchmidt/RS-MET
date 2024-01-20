#ifndef jura_FileManager_h
#define jura_FileManager_h

// ToDo: move them into a FileTools.h/cpp file, Maybe wrap them as static functions into a class.

/** Returns the directory of the current application (or .dll) as String. */
JUCE_API juce::String getApplicationDirectory();

/** Returns the user's application data directory. */
JUCE_API juce::String getUserAppDataDirectory();

/** Returns the folder used for documents that are shared among all users of the machine. */
JUCE_API juce::String getCommonDocumentsDirectory();

/** Returns the directory, where the support files such as presets, samples, etc. are supposed to 
be found. */
JUCE_API juce::String getSupportDirectory();
 // maybe move these 4 functions as static functions into FileManager

/** Given a path that may contain backslashes, this function returns a version of that path where 
the backslashes have been replaced by forward slashes. */
JUCE_API juce::File convertBackslashToSlash(const juce::File& path);
// not yet tested

/** Returns a file object if the file with the path given by 'path' exists and is a valid
audio file, otherwise it returns File(). */
JUCE_API juce::File getAudioFileFromPath(const juce::String& path);

/** Checks whether a directory or any of its subdirectories (if any) has files in it and returns
true, if so and false otherwise. */
JUCE_API bool hasDirectoryFiles(const juce::File& directoryToCheck);

// maybe factor out into a file FileTools - put also the XmlTools into this file (they are also
// related to load/save xml files - maybe that should be als a class and these functions should
// be static member functions


//=================================================================================================

/** Compares the two files according to their filenames lexicographically. ... */

class JUCE_API FileComparator
{

public:

  static int compareElements(const juce::File first, const juce::File second) throw()
  {
    if(first.getFileName() < second.getFileName())
      return -1;
    else if(first.getFileName() > second.getFileName())
      return 1;
    else
      return 0;
  }
  // todo: Check, if arguments can be const references. If so, use them!

  juce_UseDebuggingNewOperator;
};

class FileManager;

//=================================================================================================

/** This class can be used to keep track of the currently active file in some FileManager object.
The FileManager will invoke the virtual method activeFileChanged in any of its attached listeners
whenever the currently active file was changed. */

class JUCE_API FileManagerListener
{

public:

  virtual ~FileManagerListener() {}

  /** The callback that is called when the currently active file was changed in some FileManager
  to which we have registered as listener. */
  virtual void activeFileChanged(FileManager *fileManagerThatHasChanged) = 0;

  // virtual void activeFileBecameDirty() {}  

  juce_UseDebuggingNewOperator;

};
// todo: rename to FileManagerObserver, absorb into FileManager as nested class

//=================================================================================================

/** This class can be used to manage a bunch of files in a directory to allow for loading, saving,
skipping to next/previous, etc. */

class JUCE_API FileManager
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  FileManager();

  /** Destructor. */
  virtual ~FileManager();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the permissible patterns which files must adhere to in order to be considered
  as permissible. The patterns should be delimited by a semicolon - for example, to support .wav
  and .flac files, you would call it like:
  setPermissibleWildcardPatters( String(T("*.wav;*.flac")) );  */
  virtual void setPermissibleWildcardPatterns(const juce::String& newPatterns)
  {
    wildcardPatterns = newPatterns; updateFileList();
  }

  /** Selects, whether or not subdirectories should be scanned recursively when updating the
  fileList. */
  virtual void setRecurseSubDirectories(bool shouldRecurse)
  {
    recurseSubDirectories = shouldRecurse; updateFileList();
  }

  //-----------------------------------------------------------------------------------------------
  // file-management:

  /** Sets some nominal root-directory. Relative paths returned by getRelativePath() are with
  respect to this directory. Optionally we may also move into the root directory (making it the
  active directory. */
  virtual bool setRootDirectory(const juce::File& newRootDirectory, bool gotoRootDirectory);

  /** Sets the active directory and updates the filelist so as to contain the contents of the
  new directory. */
  virtual bool setActiveDirectory(const juce::File& newActiveDirectory);

  /** Marks a a new file as active and informs whether or not this was operation was
  successful. */
  virtual bool setActiveFile(const juce::File& newActiveFile);

  /** Sets the currently active file, but only if the file is present in the fileList. */
  virtual bool setActiveFileIfInList(const juce::File& newActiveFile);

  /** Marks the current file as 'saved', indicating that the file on the disk is in sync with
  the currently visible state. When the user performs some action that makes the current state
  different from the contents of the file, pass false here such that we can show a star next to
  the filename. */
  virtual void markFileAsClean(bool shouldBeMarkedAsClean);

  /** Opposite of markFileAsClean. */
  virtual void markFileAsDirty() { markFileAsClean(false); }
  // maybe remove markFileAsClean and keep only markfileAsDirty

  /** Override this function to load a new file. */
  virtual bool loadFile(const juce::File& fileToLoad) = 0;

  /** Override this function to save into a file. */
  virtual bool saveToFile(const juce::File& fileToSaveTo) = 0;

  /** Loads the next file in the current directory. */
  virtual bool loadNextFile();

  /** Loads the previous file in the current directory. */
  virtual bool loadPreviousFile();

  /** Opens a dialog-box to load from a file and returns true when a file was loaded and false
  when the dialog was closed without loading anything. */
  virtual bool openLoadingDialog(const juce::String &dialogTitle = juce::String("Load File"));

  /** Opens a dialog-box to save to a fileand returns true when a file was loaded and false
  when the dialog was closed without saving anything. */
  virtual bool openSavingDialog(const juce::String &dialogTitle = juce::String("Save File"));

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the currently marked file (in order to load it or whatever else. */
  virtual juce::File getActiveFile();

  /** Returns the relative path of the active file with respect to the nominal root directory
  (which is the application directory by default and can be set via setRootDirectory). */
  virtual juce::String getActiveFileRelativePath();

  /** Returns the directory in which the active file resides. */
  virtual juce::File getActiveDirectory();

  /** Returns the directory in which the active file resides as String. */
  virtual juce::String getActiveDirectoryAsString();

  /** Returns the list of files as juce::Array. */
  virtual juce::Array<juce::File> getFileList();

  /** Returns the list of files as a String. */
  virtual juce::String getFileListAsString(
    const juce::String &fileSeparatorString = juce::String("\n"));

  /** Returns the number of files in the list. If you use this in a loop in conjunction with
  getFileByIndex(int), you may want to wrap the loop into calls to lockFileList/unlockFileList
  to make it thread-safe. */
  virtual int getNumFilesInList();

  /** Returns a file form the list at the given index - if the index is out of range, it will
  return File(). */
  virtual juce::File getFileByIndex(int index);

  /** Returns the index of the given file inside our list - when the file is not in the list,
  it returns -1. */
  virtual int getFileIndexInList(const juce::File& fileToLookFor) const;

  /** Returns true when the given file is in our list, false otherwise. */
  virtual bool isFileInList(const juce::File& fileToLookFor) const
  {
    return getFileIndexInList(fileToLookFor) != -1;
  }

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Locks the member array fileList. */
  virtual void lockFileList();

  /** Unlocks the member array fileList. */
  virtual void unlockFileList();

  /** Adds a listener that will be informed about changes in the currently active file via a
  callback to FileManagerListener::activeFileChanged. */
  virtual void addFileManagerListener(FileManagerListener *listenerToAdd);

  /** Removes a FileManagerListener. */
  virtual void removeFileManagerListener(FileManagerListener *listenerToRemove);

  /** Updates the array which contains all the relevant files in order to skip through them
  via plus/minus buttons. */
  virtual void updateFileList();

  /** A static string member used to locate the company-specific support file folder. By default,
  it's set to "RS-MET", but since it's public, you can change that from any convenient place.  */
  static juce::String companyName;
  // Maybe rename to something more appropriate like supportFileFolder. Why do we need this at all
  // when we also have the rootDirectory member?

protected:

  /** Sends out a message to our listeners to inform them that the currently active file was
  changed. */
  virtual void notifyListeners();
  // rename to sendActiveFileChangeNotification and add sendActiveFileBecameDirtyNotification


  bool         fileIsDirty;
  bool         recurseSubDirectories;
  juce::File   rootDirectory;                   // a root directory for relative paths
  juce::File   activeDirectory;                 // directory where we currently are
  juce::String defaultExtension;
  juce::String wildcardPatterns;
  juce::Array<juce::File, CriticalSection> fileList;  // holds all relevant files in a flat array
  int activeFileIndex;                                // index of the currently active file in the array
  FileComparator fileComparator;                      // this object is needed to sort the file-array
  juce::Array<FileManagerListener*, CriticalSection> listeners;  // our listeners

  // Try to switch from juce::Array to std::vector for the fileList and listeners

  juce_UseDebuggingNewOperator;
};




#endif
