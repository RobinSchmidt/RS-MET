String getApplicationDirectory()
{
  File thisExeAsFile       = File::getSpecialLocation(File::currentApplicationFile);
  File thisDirectoryAsFile = thisExeAsFile.getParentDirectory();

  // debug code due to problem on mac:
  //bool isDir = thisDirectoryAsFile.isDirectory();
  if( !thisDirectoryAsFile.isDirectory() )
  {
    jassertfalse;
  }

  String thisDirectoryAsString = thisDirectoryAsFile.getFullPathName();
  return thisDirectoryAsString;
}

juce::String getUserAppDataDirectory()
{
  return File::getSpecialLocation(juce::File::userApplicationDataDirectory).getFullPathName();
}

juce::String getCommonDocumentsDirectory()
{
  return File::getSpecialLocation(juce::File::commonDocumentsDirectory).getFullPathName();
}

juce::String getSupportDirectory()
{
  //return getCommonDocumentsDirectory() + File::getSeparatorString() + FileManager::companyName;
  return getUserAppDataDirectory() + File::getSeparatorString() + FileManager::companyName;
}

juce::File convertBackslashToSlash(const juce::File& path)
{
  String str = path.getFullPathName();
  str = str.replaceCharacter('\\', '/');
  return File(str);
}

File getAudioFileFromPath(const String& path)
{
  File file = File(path);
  if( !file.existsAsFile() )
    return File();

  // the file exists - now check whether it is a valid audio file:
  AudioFileInfo info(file);
  if( !info.isValidAudioFile )
    return File();

  return file;
}

bool hasDirectoryFiles(const File& directoryToCheck)
{
  // when the passed file is not a directory, we consider it as having itself as file:
  if( !directoryToCheck.isDirectory() )
    return true;

  //OwnedArray<File> foundFiles;
  //directoryToCheck.findChildFiles(foundFiles, File::findFiles, true);

  juce::Array<File> foundFiles;
  directoryToCheck.findChildFiles(foundFiles, File::findFiles, true);
  if( foundFiles.size() > 0 )
    return true;
  else
    return false;
  // Note: we can't use File::getNumberOfChildFiles here because it doesn't scan recursively
}

//=================================================================================================

juce::String FileManager::companyName = "RS-MET";

// construction/destruction:

FileManager::FileManager()
{
  fileIsDirty           = false;
  recurseSubDirectories = false;
  rootDirectory         = getApplicationDirectory();
  activeDirectory       = rootDirectory;
  defaultExtension      = String();
  wildcardPatterns      = String("*");
  activeFileIndex       = -1;
}

FileManager::~FileManager()
{

}

//-------------------------------------------------------------------------------------------------
// file-management:

bool FileManager::setRootDirectory(const File& newRootDirectory, bool gotoRootDirectory)
{
  if( newRootDirectory.exists() && newRootDirectory.isDirectory() )
  {
    rootDirectory = newRootDirectory;
    if( gotoRootDirectory == true )
      return setActiveDirectory(rootDirectory);
    else
      return false;
  }
  else
    return false;
}

bool FileManager::setActiveDirectory(const File &newActiveDirectory)
{
  if( newActiveDirectory.exists() && newActiveDirectory.isDirectory() )
  {
    if( newActiveDirectory != activeDirectory )
    {
      activeDirectory = newActiveDirectory;
      updateFileList();
    }
    return true;
  }
  return false;
}

bool FileManager::setActiveFile(const File& newActiveFile)
{
  fileList.getLock().enter();
  if( getActiveDirectory() != newActiveFile.getParentDirectory() )
  {
    int oldIndex = activeFileIndex;
    if( setActiveDirectory(newActiveFile.getParentDirectory()) == false )
    {
      activeFileIndex = oldIndex;
      fileList.getLock().exit();
      return false;
    }
  }
  for(int i=0; i<fileList.size(); i++)
  {
    if(fileList[i] == newActiveFile)
    {
      activeFileIndex = i;
      fileList.getLock().exit();
      notifyListeners();
      return true;
    }
  }
  fileList.getLock().exit();
  return false;
}

bool FileManager::setActiveFileIfInList(const File& newActiveFile)
{
  if( fileList.contains(newActiveFile) )
    return setActiveFile(newActiveFile);
  else
    return false;
}

void FileManager::markFileAsClean(bool shouldBeMarkedAsClean)
{
  fileIsDirty = !shouldBeMarkedAsClean;
}

/*
bool FileManager::loadFile(const File& fileToLoad, bool sendActiveFileChangeMessage)
{
  bool result = setActiveFile(fileToLoad);
  if( result == true )
    fileIsUnsaved = false;
  return result;
}

bool FileManager::saveToFile(const File& fileToSaveTo, bool sendActiveFileChangeMessage)
{
  bool result = setActiveFile(fileToSaveTo);
  if( result == true )
    fileIsUnsaved = false;
  return result;
}
*/

bool FileManager::loadNextFile()
{
  fileList.getLock().enter();
  if(fileList.size() <= 0)
  {
    fileList.getLock().exit();
    return false;
  }
  activeFileIndex += 1;
  while( activeFileIndex > fileList.size()-1 )
    activeFileIndex -= fileList.size();
  bool result = loadFile(fileList[activeFileIndex]);
  fileList.getLock().exit();
  return result;
}

bool FileManager::loadPreviousFile()
{
  fileList.getLock().enter();
  if(fileList.size() <= 0)
  {
    fileList.getLock().exit();
    return false;
  }
  activeFileIndex -= 1;
  while( activeFileIndex < 0 )
    activeFileIndex += fileList.size();
  bool result = loadFile(fileList[activeFileIndex]);
  fileList.getLock().exit();
  return result;
}

bool FileManager::openLoadingDialog(const String &dialogTitle)
{
  FileChooser chooser(dialogTitle, getActiveDirectory(), wildcardPatterns, true);
  if( chooser.browseForFileToOpen() )
  {
    File fileToLoad = chooser.getResult();
    bool result     = loadFile(fileToLoad);
    if( result == true )
    {
      //setActiveFile(fileToLoad); 
      // BUG (possibly): I think, this call may be superfluous. The loadFile call did already call
      // that. Calling it here again may lead to redundant activeFileChanged callbacks.
      // ...OK...commenting it out seems to have fixed the issue. More tests are needed. If all
      // works well, we may delete it finally (and then simplify the remaining code). Try loading
      // xml presets in various modules, waveforms ins Straightliner's oscs, colormaps in scope

      return true;
    }
    else
      return false;
  }
  else
    return false;
}

bool FileManager::openSavingDialog(const String &dialogTitle)
{
  FileChooser chooser(dialogTitle, getActiveDirectory(), wildcardPatterns, true);
  if(chooser.browseForFileToSave(true))
  {
    File fileToSaveTo = chooser.getResult();
    if ( !fileToSaveTo.hasFileExtension(defaultExtension) && defaultExtension != String() )
      fileToSaveTo = fileToSaveTo.withFileExtension( defaultExtension ) ;
    bool result = saveToFile(fileToSaveTo);
    if( result == true )
    {
      setActiveFile(fileToSaveTo);
      return true;
    }
    else
      return false;
  }
  else
    return false;
}

void FileManager::lockFileList()
{
  fileList.getLock().enter();
}

void FileManager::unlockFileList()
{
  fileList.getLock().exit();
}

void FileManager::addFileManagerListener(FileManagerListener *listenerToAdd)
{
  listeners.getLock().enter();
  listeners.addIfNotAlreadyThere(listenerToAdd);
  listeners.getLock().exit();
}

void FileManager::removeFileManagerListener(FileManagerListener *listenerToRemove)
{
  listeners.getLock().enter();
  listeners.removeFirstMatchingValue(listenerToRemove);
  listeners.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

File FileManager::getActiveFile()
{
  if( activeFileIndex == -1 )
    return File();
  else
    return fileList[activeFileIndex];
}

String FileManager::getActiveFileRelativePath()
{
  File activeFile = getActiveFile();
  return activeFile.getRelativePathFrom(rootDirectory);
}

File FileManager::getActiveDirectory()
{
  return activeDirectory;
}

String FileManager::getActiveDirectoryAsString()
{
  return getActiveDirectory().getFullPathName();
}

juce::Array<File> FileManager::getFileList()
{
  juce::Array<File> result;
  fileList.getLock().enter();
  for(int i=0; i<fileList.size(); i++)
    result.add(fileList[i]);
  fileList.getLock().exit();
  return result;
}

String FileManager::getFileListAsString(const String &fileSeparatorString)
{
  String fileListString;
  fileList.getLock().enter();
  for(int f=0; f<fileList.size(); f++)
    fileListString += fileList[f].getFullPathName() + fileSeparatorString;
  fileList.getLock().exit();
  return fileListString;
}

int FileManager::getNumFilesInList()
{
  fileList.getLock().enter();
  int result = fileList.size();
  fileList.getLock().exit();
  return result;
}

File FileManager::getFileByIndex(int index)
{
  fileList.getLock().enter();
  File result = File();
  if( index >= 0 && index < fileList.size() )
    result = fileList[index];
  fileList.getLock().exit();
  return result;
}

int FileManager::getFileIndexInList(const File& fileToLookFor) const
{
  fileList.getLock().enter();
  for(int i=0; i<fileList.size(); i++)
  {
    if(fileList[i] == fileToLookFor)
    {
      fileList.getLock().exit();
      return i;
    }
  }
  fileList.getLock().exit();
  return -1;
}

//-------------------------------------------------------------------------------------------------
// others:

void FileManager::updateFileList()
{
  fileList.getLock().enter();
  if( getActiveDirectory().exists() && getActiveDirectory().isDirectory() )
  {
    String oldFileName = getActiveFile().getFileName();
    fileList.clear();
    juce::Array<File> tmpFileList;
    File(getActiveDirectory()).findChildFiles(tmpFileList, File::findFiles, recurseSubDirectories);
    WildcardFileFilter fileFilter(wildcardPatterns, String(), String("Wildcard filter") );
    int i;
    for(i=0; i<tmpFileList.size(); i++)
    {
      if( fileFilter.isFileSuitable(tmpFileList[i]) )
        fileList.add(tmpFileList[i]);
    }
    fileList.sort(fileComparator);
    activeFileIndex = -1;
    for(i=0; i<fileList.size(); i++)
    {
      if(fileList[i].getFileName() == oldFileName)
      {
        activeFileIndex = i;
        break;
      }
    }
  }
  fileList.getLock().exit();
}

void FileManager::notifyListeners()
{
  listeners.getLock().enter();
  for(int i=0; i<listeners.size(); i++)
    listeners[i]->activeFileChanged(this);
  listeners.getLock().exit();
}
