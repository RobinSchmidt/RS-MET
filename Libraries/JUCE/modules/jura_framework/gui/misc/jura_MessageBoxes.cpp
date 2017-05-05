

// maybe move to header file - at least a declaration
void JUCE_CALLTYPE showWarningBox(const String& title, const String& message) // rename to showWarningBox
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, title, message, "OK");
}


void showMemoryAllocationErrorBox(String sourceFunctionName)
{
  showWarningBox("Memory Allocation Error", 
    "Memory could not be allocated in: " + sourceFunctionName);
}

void showAudioFileInvalidErrorBox(String fileName)
{
  showWarningBox("Audiofile invalid", 
    "The file: " + fileName + " is not a supported audiofile.");
}

void showAudioFileTooLargeErrorBox()
{
  showWarningBox("Audiofile too large",  
    "Audiofiles with more than " + String(INT_MAX) + " samples are not supported");
}

void showAudioFileWriteErrorBox(String fileName)
{
  showWarningBox("Audiofile write error", 
    "The audio-file: " + fileName + " could not be written");
}

void showDirectoryCouldNotBeCreatedBox(File directoryThatCouldNotBeCreated)
{
  showWarningBox("Directory could not be created", 
    "The directory: " +  directoryThatCouldNotBeCreated.getFullPathName()
    + " could not be created");
}

void showEnterNameErrorBox()
{
  showWarningBox("Please enter a name", 
    "Please enter a name for the new project in the textfield.");
}

void showFileCouldNotBeCopiedBox(const String& fileName, const String& attemptedDirectoryName)
{
  showWarningBox("File could not be copied", 
    "The file: " + fileName + " could not be copied into the directory" + attemptedDirectoryName);
}

void showFileCouldNotBeDeletedBox(const String& fileName)
{
  showWarningBox("File could not be deleted", 
    "The file: " + fileName + " could not be deleted");
}

void showFileCouldNotBeMovedBox(const String& fileName, const String& attemptedDirectoryName)
{
  showWarningBox("File could not be moved", 
    "The file: " + fileName + " could not be moved into the directory" + attemptedDirectoryName);
}

void showFileIsNoDirectoryErrorBox(const File& file)
{
  showWarningBox("Not a directory", "The chosen file: " + file.getFileName() + " is not a directory");
}

void showFileNotFoundOrInvalidAudioFileBox(const String& fileName)
{
  showWarningBox("File not found or invalid", 
    "The file: " + fileName
    + " was not found or is not a valid audiofile - clips using this file will be omitted");
}

void showInvalidPathNameErrorBox(const String& theInvalidName)
{
  showWarningBox("Invalid path name", "The chosen project directory name: " + theInvalidName 
    + " is not a valid name for a directory. Please enter another name." );
}

void showPleaseSelectProjectFileBox()
{
  showWarningBox("Please select a project file", 
    "Please select a valid Mixsonic project file to load.");
}

void showProjectCouldNotBeCreatedBox(const String& projectName)
{
  showWarningBox("Project could not be created", 
    "The project: " + projectName + " could not be created");
}

void showProjectDirectoryAlreadyExistsErrorBox(File existingProjectDirectory)
{
  showWarningBox("Directory already exists", 
    "The chosen project directory: " + existingProjectDirectory.getFullPathName() 
    + " already exists. Please enter another name.");
}

void showProjectsParentDirectoryInvalidBox(const String& pathString)
{
  showWarningBox("Invalid path for projects", 
    "The parent directory for the projects: " + pathString + " does not exist or is invalid.");
}

void showSampleContentPathInvalidBox(const String& pathString)
{
  showWarningBox("Invalid sample content path", 
    "The path: " + pathString + " for the sample content does not exist or is invalid.");
}

void showSettingsFileIsMissingBox()
{
  showWarningBox("Settings file missing", 
    "The global application settings file: MixsonicSettings.xml was not found in the \
    application directory - Mixsonic will use fallback settings, which are likely to \
be meaningless. Please re-install to resolve this issue." );
}

void showSettingsFileIsInvalidBox()
{
  showWarningBox("Settings file missing", 
    "The global application settings file: MixsonicSettings.xml is corrupted - Mixsonic will\
 use fallback settings, which are likely to be meaningless. Please re-install to resolve\
 this issue.");
}

// some of them are app-specific - move them into appropriate place

// OK/Cancel boxes:

bool JUCE_CALLTYPE showOkCancelWarningBox(const String& title, const String& message) // rename to showWarningBox
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, title, message, "OK", "Cancel");
}

bool showDirectoryDeleteWarningBox(String directoryName)
{
  return showOkCancelWarningBox("Delete directory?", 
    "The directory: " + directoryName +
    " will be deleted from the project directory and all clips using a sample from this\
 directory will be removed - proceed?");
}

bool showOverwriteAudioFileWarningBox(String fileName)
{
  return showOkCancelWarningBox("Overwrite audio file?", 
    "The audio-file: " + fileName + " will be overwritten - proceed?");
}

bool showOverwriteSongFileWarningBox(String fileName)
{
  return showOkCancelWarningBox("Overwrite song file?", 
    "The song-file: " + fileName + " will be overwritten - proceed?");
}

bool showSampleDeleteWarningBox(String fileName)
{
  return showOkCancelWarningBox("Delete sample?", 
    "The sample file: " + fileName 
    + " will be deleted from the project directory and all clips using that sample will be removed - proceed?");
}

bool showSongUnsavedWarningBox()
{
  return showOkCancelWarningBox("Song not saved!", 
    "The current song is not saved and will be discarded - proceed?");
}

